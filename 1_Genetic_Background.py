from Bio import Entrez
from Bio import SeqIO
import time
import os
import re
from io import StringIO
from typing import Dict, Optional, Tuple
import validators
from http.client import IncompleteRead
import socket
from urllib.error import HTTPError, URLError
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('1_Genetic_Background.log', encoding='utf-8')
    ]
)

class ValidationError(Exception):
    pass

def validate_email(email: str) -> bool:
    return bool(validators.email(email))

def validate_api_key(api_key: str) -> bool:
    return bool(re.match(r'^[a-zA-Z0-9]{36}$', api_key))

def validate_taxid(taxid: str) -> bool:
    return bool(re.match(r'^\d+$', taxid))

def safe_write_to_file(filepath: str, content: str) -> None:
    try:
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except OSError as e:
        raise IOError(f"Failed to write to {filepath}: {e}")

def fetch_batch(db: str, webenv: str, query_key: str, start: int, 
                batch_size: int, rettype: str, retmode: str, retries: int = 3) -> Optional[str]:
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(
                db=db, webenv=webenv, query_key=query_key,
                retstart=start, retmax=batch_size, 
                rettype=rettype, retmode=retmode
            )
            return handle.read()
            
        except HTTPError as e:
            if hasattr(e, 'code') and e.code == 400:
                logging.error(f"Bad Request Error: {e}")
                return None
            logging.error(f"HTTP Error on attempt {attempt + 1}: {e}")
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            
        except (URLError, IncompleteRead, socket.timeout) as e:
            logging.error(f"Network error on attempt {attempt + 1}: {e}")
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
                
        except Exception as e:
            logging.error(f"Unexpected error on attempt {attempt + 1}: {e}")
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
    
    return None

def validate_fasta_file(filepath: str) -> bool:
    if not os.path.exists(filepath):
        return True
    try:
        with open(filepath, 'r', encoding='utf-8') as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            return len(records) > 0
    except Exception:
        return False

def normalize_gene_name(gene_name: str) -> str:
    return gene_name.upper().replace(' ', '_').replace('-', '_').replace('/', '_')

def perform_search(term: str, retries: int = 3) -> Optional[Tuple[int, str, str]]:
    for attempt in range(retries):
        try:
            search_handle = Entrez.esearch(db="nucleotide", term=term, usehistory="y", retmax=0)
            results = Entrez.read(search_handle)
            count = int(results["Count"])
            return count, results["WebEnv"], results["QueryKey"]
        except Exception as e:
            if attempt == retries - 1:
                logging.error(f"Search failed after {retries} attempts: {e}")
                return None
            time.sleep(2 * (attempt + 1))
    return None

def process_taxid(taxid: str, output_dir: str) -> None:
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        organism_name = records[0]["ScientificName"]
        logging.info(f"\n{'='*50}\nProcessing {organism_name} (TaxID: {taxid})\n{'='*50}")
    except HTTPError as e:
        if hasattr(e, 'code') and e.code == 400:
            logging.error(f"Invalid TaxID {taxid}: Bad Request")
            return
        logging.error(f"Error fetching taxonomy data for {taxid}: {e}")
        return
    except Exception as e:
        logging.error(f"Error fetching taxonomy data for {taxid}: {e}")
        return

    search_term = f"txid{taxid}[Organism]"
    search_result = perform_search(search_term)
    if not search_result:
        return
    
    count, webenv, query_key = search_result
    logging.info(f"Found {count} total sequences for {organism_name}")
    if count == 0:
        logging.info(f"No sequences found for taxid {taxid}")
        return

    gene_counts: Dict[str, int] = {}
    feature_types = set()
    batch_size = 100
    total_records_processed = 0
    total_gene_features_found = 0
    
    logging.info(f"\nProcessing sequences in batches of {batch_size}:")
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        logging.info(f"\nFetching and analyzing records {start+1} to {end} ({end-start} records)")
        
        data = fetch_batch(
            db="nucleotide", webenv=webenv, query_key=query_key,
            start=start, batch_size=batch_size, rettype="gb", retmode="text"
        )
        
        if not data:
            logging.warning(f"Failed to fetch batch {start+1} to {end}, skipping...")
            continue

        try:
            records = SeqIO.parse(StringIO(data), "gb")
            batch_records_processed = 0
            batch_gene_features_found = 0
            
            for record in records:
                batch_records_processed += 1
                for feature in record.features:
                    feature_types.add(feature.type)
                    
                    gene_name = None
                    
                    for qualifier in ['gene', 'locus_tag', 'product']:
                        if qualifier in feature.qualifiers:
                            gene_name = feature.qualifiers[qualifier][0]
                            break
                    
                    if gene_name:
                        batch_gene_features_found += 1
                        gene_name = normalize_gene_name(gene_name)
                        gene_counts[gene_name] = gene_counts.get(gene_name, 0) + 1
            
            total_records_processed += batch_records_processed
            total_gene_features_found += batch_gene_features_found
            logging.info(f"Processed {batch_records_processed} records in this batch")
            logging.info(f"Found {batch_gene_features_found} gene features in this batch")
            
        except Exception as e:
            logging.error(f"Error processing batch: {e}")
            continue

        time.sleep(0.1)

    logging.info(f"\nProcessing complete for {organism_name}:")
    logging.info(f"Total records processed: {total_records_processed}")
    logging.info(f"Total gene features found: {total_gene_features_found}")
    logging.info(f"Feature types found: {', '.join(sorted(feature_types))}")
    logging.info(f"Total unique genes found: {len(gene_counts)}")

    filtered_genes = {gene: count for gene, count in gene_counts.items() if count >= 8}
    sorted_genes = sorted(filtered_genes.items(), key=lambda x: x[1], reverse=True)

    if not sorted_genes:
        logging.info(f"No genes with sequence count above 2 for TaxID {taxid}")
        return

    stats_file = os.path.join(output_dir, f"{taxid}_gene_stats.txt")
    try:
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write(f"Organism: {organism_name}\n")
            f.write(f"TaxID: {taxid}\n")
            f.write(f"Total sequences analyzed: {count}\n")
            f.write(f"Total unique genes found: {len(gene_counts)}\n")
            f.write(f"Feature types found: {', '.join(sorted(feature_types))}\n\n")
            f.write("Gene counts (≥2 sequences):\n")
            for gene, count in sorted_genes:
                f.write(f"{gene}\t{count}\n")
        logging.info(f"\nGene statistics saved to: {stats_file}")
    except IOError as e:
        logging.error(f"Failed to save statistics: {e}")
        return

    logging.info("\nFetching individual sequences for each gene...")
    for gene, count in sorted_genes:
        logging.info(f"\nProcessing gene '{gene}' ({count} sequences in total)")
        fasta_file = os.path.join(output_dir, f"{taxid}_{gene}.fasta")
        
        if os.path.exists(fasta_file):
            if not validate_fasta_file(fasta_file):
                logging.warning(f"Existing FASTA file {fasta_file} is corrupted, will recreate")
                try:
                    os.remove(fasta_file)
                except OSError as e:
                    logging.error(f"Failed to remove corrupted file: {e}")
                    continue
            else:
                logging.info(f"File {fasta_file} already exists, will append new sequences")
        else:
            logging.info(f"Creating new file: {fasta_file}")

        gene_query = f"txid{taxid}[Organism] AND {gene}[All Fields]"
        gene_search_result = perform_search(gene_query)
        
        if not gene_search_result:
            continue
            
        gene_count, gene_webenv, gene_query_key = gene_search_result
        if gene_count == 0:
            logging.warning(f"No sequences found for gene {gene}")
            continue

        logging.info(f"Found {gene_count} sequences for gene {gene}")
        try:
            with open(fasta_file, 'a' if os.path.exists(fasta_file) else 'w', encoding='utf-8') as f:
                sequences_saved = 0
                for start in range(0, gene_count, batch_size):
                    end = min(gene_count, start + batch_size)
                    logging.info(f"Fetching sequences {start+1}-{end} for gene '{gene}'")
                    data = fetch_batch(
                        db="nucleotide", webenv=gene_webenv,
                        query_key=gene_query_key, start=start,
                        batch_size=batch_size, rettype="fasta",
                        retmode="text"
                    )
                    
                    if data:
                        f.write(data)
                        sequences_saved += end - start
                    
                    time.sleep(0.1)
            logging.info(f"Successfully saved {sequences_saved} sequences to {fasta_file}")
                    
        except Exception as e:
            logging.error(f"Error saving sequences for gene {gene}: {e}")

def main():
    try:
        email = input("Enter your email address: ").strip()
        if not validate_email(email):
            raise ValidationError("Invalid email format")

        api_key = input("Enter your NCBI API key: ").strip()
        if not validate_api_key(api_key):
            raise ValidationError("Invalid API key format")

        user_data_dir = os.path.join(os.getcwd(), 'user_data')
        safe_write_to_file(os.path.join(user_data_dir, 'email.txt'), email)
        safe_write_to_file(os.path.join(user_data_dir, 'api_key.txt'), api_key)

        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'Genetic_Background'

        output_dir = os.path.join(os.getcwd(), '1_')
        os.makedirs(output_dir, exist_ok=True)

        taxids_input = input("Enter TaxID number(s) separated by commas: ").strip()
        taxid_list = [taxid.strip() for taxid in taxids_input.split(',')]
        
        invalid_taxids = [taxid for taxid in taxid_list if not validate_taxid(taxid)]
        if invalid_taxids:
            raise ValidationError(f"Invalid TaxID format: {', '.join(invalid_taxids)}")

        for taxid in taxid_list:
            process_taxid(taxid, output_dir)

    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

if __name__ == '__main__':
    main()
