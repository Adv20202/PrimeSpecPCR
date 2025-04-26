#!/usr/bin/env python3
"""
PrimeSpecPCR - Genetic Background Module (Optimized)
Improved batch processing and file-based input handling for GUI compatibility
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time
import os
import re
import sys
from io import StringIO
from typing import Dict, Optional, Tuple, List, Any
import validators
from http.client import IncompleteRead
import socket
from urllib.error import HTTPError, URLError
import logging
import sys

# Configure logging with explicit flush after each message
class FlushingStreamHandler(logging.StreamHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()

# Custom formatter to remove "- INFO -" from messages
class CustomFormatter(logging.Formatter):
    def format(self, record):
        # Skip timestamp, just return the message
        return record.getMessage()

# Set up logging with immediate flushing and modified time format
formatter = CustomFormatter('%(message)s')
handler1 = FlushingStreamHandler(sys.stdout)
handler1.setFormatter(formatter)
handler2 = logging.FileHandler('1_Genetic_Background.log', encoding='utf-8')
handler2.setFormatter(formatter)

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(handler1)
logger.addHandler(handler2)

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

# Function to get accurate gene count for a specific gene
def get_gene_sequence_count(taxid: str, gene_name: str) -> int:
    """Get the actual sequence count for a specific gene using precise search"""
    gene_query = f"txid{taxid}[Organism] AND {gene_name}[All Fields]"
    gene_search_result = perform_search(gene_query)
    
    if not gene_search_result:
        return 0
    
    count, _, _ = gene_search_result
    return count

# Function for screening genes for a given TaxID
def screen_taxid(taxid: str, output_dir: str) -> Tuple[Optional[str], Optional[Dict[str, int]]]:
    try:
        logging.info(f"Fetching taxonomy information for TaxID: {taxid}")
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        organism_name = records[0]["ScientificName"]
        logging.info(f"\n{'='*50}\nScreening {organism_name} (TaxID: {taxid})\n{'='*50}")
    except HTTPError as e:
        if hasattr(e, 'code') and e.code == 400:
            logging.error(f"Invalid TaxID {taxid}: Bad Request")
            return None, None
        logging.error(f"Error fetching taxonomy data for {taxid}: {e}")
        return None, None
    except Exception as e:
        logging.error(f"Error fetching taxonomy data for {taxid}: {e}")
        return None, None

    search_term = f"txid{taxid}[Organism]"
    logging.info(f"Searching for sequences with term: {search_term}")
    search_result = perform_search(search_term)
    if not search_result:
        return None, None
    
    count, webenv, query_key = search_result
    logging.info(f"Found {count} total sequences for {organism_name}, please wait..")
    if count == 0:
        logging.info(f"No sequences found for taxid {taxid}")
        return organism_name, {}

    gene_counts: Dict[str, int] = {}
    feature_types = set()
    batch_size = 100
    total_records_processed = 0
    total_gene_features_found = 0
    
    # Process more sequences - instead of just 10 batches/1000 sequences, process up to 5000
    # but still have a reasonable limit to avoid excessive API calls for very large datasets
    max_sequences = min(50000, count)
    max_batches = (max_sequences + batch_size - 1) // batch_size
    
    # No intermediate progress updates - just process all batches silently
    for batch_idx in range(max_batches):
        start = batch_idx * batch_size
        end = min(max_sequences, start + batch_size)
        
        data = fetch_batch(
            db="nucleotide", webenv=webenv, query_key=query_key,
            start=start, batch_size=batch_size, rettype="gb", retmode="text"
        )
        
        if not data:
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
            
        except Exception as e:
            logging.error(f"Error processing batch: {e}")
            continue

        # More careful with rate limits for longer runs
        time.sleep(0.1)
    
    logging.info(f"\nScreening complete for {organism_name}:")
    logging.info(f"Total records processed: {total_records_processed}")
    # Usunięto informację o cechach genowych (gene features)
    logging.info(f"Total unique genes found: {len(gene_counts)}")

    # Verify gene counts with precise search
    verified_gene_counts = {}
    
    # Sort genes by occurrence count (descending)
    sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
    
    # For efficiency, only get precise counts for top genes (max 50)
    top_genes = sorted_genes[:50]
    
    for gene, approx_count in top_genes:
        # Get accurate sequence count with precise query
        precise_count = get_gene_sequence_count(taxid, gene)
        if precise_count > 0:
            verified_gene_counts[gene] = precise_count
    
    # Sort the verified counts
    sorted_verified_genes = sorted(verified_gene_counts.items(), key=lambda x: x[1], reverse=True)

    stats_file = os.path.join(output_dir, f"{taxid}_gene_stats.txt")
    try:
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write(f"Organism: {organism_name}\n")
            f.write(f"TaxID: {taxid}\n")
            f.write(f"Total sequences analyzed: {count}\n")
            f.write(f"Total unique genes found: {len(verified_gene_counts)}\n")
            f.write(f"Feature types found: {', '.join(sorted(feature_types))}\n\n")
            f.write("Gene counts:\n")
            
            # Use verified gene counts for the stats file with continuous numbering
            for i, (gene, count) in enumerate(sorted_verified_genes, 1):
                # Start numbering from 301 if not among first 300 genes
                gene_num = i if i <= 300 else i
                f.write(f"{gene_num}. {gene}\t{count}\n")
                
        logging.info(f"\nGene statistics saved to: {stats_file}")
    except IOError as e:
        logging.error(f"Failed to save statistics: {e}")
        return organism_name, None

    return organism_name, verified_gene_counts

# Function to fetch reference sequence for BLAST
def fetch_reference_sequence(taxid: str, gene_name: str, sequence_number: int) -> Optional[Tuple[str, str]]:
    """Fetch a specific sequence to be used as a reference for BLAST."""
    try:
        gene_query = f"txid{taxid}[Organism] AND {gene_name}[All Fields]"
        logging.info(f"Searching for reference sequences with query: {gene_query}")
        
        search_result = perform_search(gene_query)
        if not search_result:
            logging.error(f"Failed to search for reference sequence: {gene_name}")
            return None
            
        count, webenv, query_key = search_result
        if count == 0:
            logging.error(f"No sequences found for {gene_name}")
            return None
            
        if sequence_number > count:
            logging.error(f"Requested sequence number {sequence_number} exceeds available sequences ({count})")
            return None
        
        # Fetch the specific sequence
        data = fetch_batch(
            db="nucleotide", webenv=webenv, query_key=query_key,
            start=sequence_number-1, batch_size=1, rettype="gb",
            retmode="text"
        )
        
        if not data:
            logging.error(f"Failed to fetch reference sequence data")
            return None
            
        # Parse the GenBank record
        record = next(SeqIO.parse(StringIO(data), "gb"))
        sequence = str(record.seq)
        accession = record.id
        
        # Return the sequence and its accession number
        return sequence, accession
        
    except Exception as e:
        logging.error(f"Error fetching reference sequence: {e}")
        return None

# Function to run BLAST search and download homologous regions
def blast_and_download_homologs(taxid: str, gene_name: str, ref_sequence: str, ref_accession: str, output_dir: str, 
                              max_hits: int = 1000, evalue: float = 0.1) -> None:
    """
    Use BLAST to identify homologous sequences and download only the homologous regions.
    """
    blast_file = os.path.join(output_dir, f"{taxid}_{gene_name}_blast_results.xml")
    fasta_file = os.path.join(output_dir, f"{taxid}_{gene_name}.fasta")
    
    logging.info(f"Running BLAST search for gene {gene_name} using reference sequence {ref_accession}...")
    
    # Check if BLAST results already exist
    if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
        logging.info(f"Using existing BLAST results from {blast_file}")
        with open(blast_file, 'r') as f:
            blast_results = f.read()
    else:
        try:
            # Run BLAST search
            entrez_query = f"txid{taxid}[ORGN]"  # Limit to the specified taxon
            
            # Run BLAST with retries
            for attempt in range(3):
                try:
                    result_handle = NCBIWWW.qblast(
                        program="blastn", 
                        database="nt", 
                        sequence=ref_sequence,
                        entrez_query=entrez_query,
                        expect=evalue,
                        hitlist_size=max_hits
                    )
                    blast_results = result_handle.read()
                    break
                except Exception as e:
                    logging.error(f"BLAST search attempt {attempt+1} failed: {e}")
                    if attempt < 2:  # On 3rd attempt, just give up
                        time.sleep(10 * (attempt + 1))
                    else:
                        raise
            
            # Save BLAST results for future use
            with open(blast_file, 'w') as f:
                f.write(blast_results)
                
            logging.info(f"BLAST search completed and results saved to {blast_file}")
                
        except Exception as e:
            logging.error(f"Error performing BLAST search: {e}")
            return
    
    # Create or append to the FASTA output file
    mode = 'a' if os.path.exists(fasta_file) else 'w'
    
    try:
        # Parse BLAST results
        blast_records = NCBIXML.parse(StringIO(blast_results))
        record = next(blast_records)
        
        total_hits = len(record.alignments)
        if total_hits == 0:
            logging.warning(f"No BLAST hits found for {gene_name}")
            return
            
        logging.info(f"Found {total_hits} BLAST hits for {gene_name}, downloading homologous regions...")
        
        # Open output file for writing FASTA sequences
        with open(fasta_file, mode, encoding='utf-8') as f:
            # Add reference sequence first
            f.write(f">{ref_accession} [REFERENCE] {gene_name}\n{ref_sequence}\n")
            sequences_saved = 1
            
            # Process each alignment
            for alignment in record.alignments:
                # Skip the reference sequence itself
                if alignment.accession == ref_accession:
                    continue
                    
                # Get hit details
                for hsp in alignment.hsps:
                    # Get the subject sequence
                    subject_seq = hsp.sbjct
                    
                    # Create FASTA header
                    header = f">{alignment.accession} {alignment.hit_def} [E={hsp.expect:.2e}] [identity={hsp.identities}/{hsp.align_length}] {gene_name}"
                    
                    # Write to FASTA file
                    f.write(f"{header}\n{subject_seq}\n")
                    sequences_saved += 1
                    
                    # Only take the best HSP for each alignment
                    #break
                    
                # Add a small delay to be gentle on API
                time.sleep(0.1)
        
        logging.info(f"Successfully saved {sequences_saved} sequences to {fasta_file}")
        
    except Exception as e:
        logging.error(f"Error processing BLAST results: {e}")

# Function for downloading FASTA sequences for selected genes using BLAST
def download_gene_sequences(taxid: str, organism_name: str, selected_genes: List[str], output_dir: str) -> None:
    logging.info(f"\nDownloading sequences for {len(selected_genes)} selected genes...")
    
    for gene_idx, gene in enumerate(selected_genes):
        # Add progress reporting
        logging.info(f"Processing gene {gene_idx+1} of {len(selected_genes)} ({gene})")
        
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
                logging.info(f"File {fasta_file} already exists, will use existing sequences")
                continue
        
        gene_query = f"txid{taxid}[Organism] AND {gene}[All Fields]"
        gene_search_result = perform_search(gene_query)
        
        if not gene_search_result:
            continue
            
        gene_count, gene_webenv, gene_query_key = gene_search_result
        if gene_count == 0:
            logging.warning(f"No sequences found for gene {gene}")
            continue

        # Get sequences for selection
        logging.info(f"Found {gene_count} potential reference sequences for gene {gene}")
        
        # Fetch all sequences for display
        # Let's handle chunking to avoid overloading the API
        display_count = gene_count
        batch_size = 100
        total_batches = (display_count + batch_size - 1) // batch_size

        logging.info(f"Retrieving {display_count} sequences for gene {gene}...")

        all_records = []
        for batch_idx in range(total_batches):
            start = batch_idx * batch_size
            current_batch_size = min(batch_size, display_count - start)
            
            # Fetch batch silently (removed the logging message)
            data = fetch_batch(
                db="nucleotide", webenv=gene_webenv, query_key=gene_query_key,
                start=start, batch_size=current_batch_size, rettype="gb", retmode="text"
            )
            
            if not data:
                logging.error(f"Failed to fetch sequence data for gene {gene} batch {batch_idx+1}")
                continue
                
            batch_records = list(SeqIO.parse(StringIO(data), "gb"))
            all_records.extend(batch_records)
            
            # Add a small delay to avoid rate limits
            time.sleep(0.1)
        
        if not all_records:
            logging.error(f"Failed to fetch any sequence data for gene {gene}")
            continue
        
        # Parse and display reference sequence options
        try:
            logging.info(f"\nSelect a reference sequence for {gene}:")
            
            for i, record in enumerate(all_records, 1):
                # Create a more informative description
                description = record.description
                if hasattr(record, 'features') and record.features:
                    for feature in record.features:
                        if feature.type == 'source' and 'organism' in feature.qualifiers:
                            description = f"{feature.qualifiers['organism'][0]} - {description}"
                            break
                
                sequence_length = len(record.seq)
                # Modified format
                logging.info(f"{i}. {record.id}: {description}: {sequence_length} bp")
            
            # Ask user to select reference sequence
            selected_seq = int(get_console_input(f"Select reference sequence number (1-{len(all_records)} or 0 if you skip this): "))
            
            if selected_seq < 1 or selected_seq > len(all_records):
                logging.error(f"Invalid selection: {selected_seq}")
                continue
            
            # Use the directly selected record from all_records instead of fetching again
            record = all_records[selected_seq-1]
            ref_sequence = str(record.seq)
            ref_accession = record.id
            
            # Use BLAST to find homologs and download them
            blast_and_download_homologs(taxid, gene, ref_sequence, ref_accession, output_dir)
                
        except Exception as e:
            logging.error(f"Error processing sequences for gene {gene}: {e}")

def get_user_data():
    """Get email and API key from user data directory or use defaults."""
    user_data_dir = os.path.join(os.getcwd(), 'user_data')
    email = None
    api_key = None
    
    if os.path.exists(user_data_dir):
        email_file = os.path.join(user_data_dir, 'email.txt')
        api_key_file = os.path.join(user_data_dir, 'api_key.txt')
        
        if os.path.exists(email_file):
            with open(email_file, 'r') as f:
                email = f.read().strip()
                
        if os.path.exists(api_key_file):
            with open(api_key_file, 'r') as f:
                api_key = f.read().strip()
    
    return email, api_key

def get_taxids():
    """Get TaxID from user_config.ini or use default."""
    taxids = []
    config_file = os.path.join(os.getcwd(), 'user_config.ini')
    
    if os.path.exists(config_file):
        import configparser
        config = configparser.ConfigParser()
        config.read(config_file)
        
        if 'User' in config and 'taxids' in config['User']:
            taxids_str = config['User']['taxids']
            taxids = [taxid.strip() for taxid in taxids_str.split(',') if taxid.strip()]
    
    return taxids

def auto_select_genes(gene_counts: Dict[str, int], min_count: int = 5) -> List[str]:
    """Automatically select genes based on minimum sequence count."""
    selected_genes = []
    
    # Sort genes by sequence count (descending)
    sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Select genes with at least min_count sequences
    for gene, count in sorted_genes:
        if count >= min_count:
            selected_genes.append(gene)
    
    # Limit to top 5 genes if too many are selected
    if len(selected_genes) > 5:
        selected_genes = selected_genes[:5]
        
    return selected_genes

# File-based input handling method
def get_console_input(prompt, default=None):
    """
    File-based input function for better GUI compatibility
    """
    # Create an input request file to signal GUI
    input_request_file = os.path.join(os.getcwd(), 'input_request.txt')
    input_response_file = os.path.join(os.getcwd(), 'input_response.txt')
    
    # Delete response file if it exists
    if os.path.exists(input_response_file):
        try:
            os.remove(input_response_file)
        except:
            pass
    
    # Only log the waiting message, don't print to avoid duplication
    logging.info(f"Waiting for input: {prompt}")
    
    # Write the prompt to the request file
    with open(input_request_file, 'w', encoding='utf-8') as f:
        f.write(prompt)
    
    # Wait for response file to appear
    timeout = 60000
    start_time = time.time()
    while not os.path.exists(input_response_file):
        time.sleep(0.5)
        if time.time() - start_time > timeout:
            logging.warning("Input timeout. Using default value.")
            return default if default is not None else ""
    
    # Wait a moment to ensure file is completely written
    time.sleep(0.1)
    
    # Read response
    try:
        with open(input_response_file, 'r', encoding='utf-8') as f:
            user_input = f.read().strip()
        
        # Delete the response file after reading
        try:
            os.remove(input_response_file)
        except:
            pass
            
        # Also delete request file
        try:
            os.remove(input_request_file)
        except:
            pass
            
        return user_input
    except Exception as e:
        logging.error(f"Error reading input: {e}")
        return default if default is not None else ""

def cleanup_previous_run():
    """Remove log file and directories/files created by previous run of module 1"""
    logging.info("Cleaning up files from previous run...")
    
    # Remove log file if it exists
    log_file = '1_Genetic_Background.log'
    if os.path.exists(log_file):
        try:
            os.remove(log_file)
            logging.info(f"Removed log file: {log_file}")
        except OSError as e:
            logging.error(f"Failed to remove log file: {e}")
    
    # Remove output directory and its contents
    output_dir = os.path.join(os.getcwd(), '1_')
    if os.path.exists(output_dir):
        try:
            for root, dirs, files in os.walk(output_dir, topdown=False):
                for file in files:
                    file_path = os.path.join(root, file)
                    os.remove(file_path)
                    logging.info(f"Removed file: {file_path}")
                for dir in dirs:
                    dir_path = os.path.join(root, dir)
                    os.rmdir(dir_path)
                    logging.info(f"Removed directory: {dir_path}")
            os.rmdir(output_dir)
            logging.info(f"Removed output directory: {output_dir}")
        except OSError as e:
            logging.error(f"Failed to remove output directory or its contents: {e}")
    
    # Remove input request file if it exists
    input_request_file = os.path.join(os.getcwd(), 'input_request.txt')
    if os.path.exists(input_request_file):
        try:
            os.remove(input_request_file)
            logging.info(f"Removed input request file: {input_request_file}")
        except OSError as e:
            logging.error(f"Failed to remove input request file: {e}")
    
    # Remove input response file if it exists
    input_response_file = os.path.join(os.getcwd(), 'input_response.txt')
    if os.path.exists(input_response_file):
        try:
            os.remove(input_response_file)
            logging.info(f"Removed input response file: {input_response_file}")
        except OSError as e:
            logging.error(f"Failed to remove input response file: {e}")
    
    logging.info("Cleanup completed")

def main():
    try:
        # Clean up files from previous run
        cleanup_previous_run()
        
        # Check command line arguments
        force_non_interactive = "--non-interactive" in sys.argv
        interactive_mode = "--interactive" in sys.argv or not force_non_interactive
        
        # Get email and API key from user_data directory
        email, api_key = get_user_data()
        
        # When run from GUI, we should always use saved credentials
        # Only ask for input in interactive standalone mode
        if interactive_mode and (not email or not api_key):
            # Interactive mode for email/API key
            email = get_console_input("Enter your email address: ")
            if not validate_email(email):
                raise ValidationError("Invalid email format")

            api_key = get_console_input("Enter your NCBI API key: ")
            if not validate_api_key(api_key):
                raise ValidationError("Invalid API key format")

            user_data_dir = os.path.join(os.getcwd(), 'user_data')
            os.makedirs(user_data_dir, exist_ok=True)
            safe_write_to_file(os.path.join(user_data_dir, 'email.txt'), email)
            safe_write_to_file(os.path.join(user_data_dir, 'api_key.txt'), api_key)
        else:
            if not email or not api_key:
                logging.error("Email or API key not found in user data. Please provide them in the GUI.")
                return
                
            logging.info(f"Using saved email: {email}")
            logging.info(f"Using saved API key: {api_key[:5]}...{api_key[-5:]}")

        # Set up Entrez
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'PrimeSpecPCR'

        output_dir = os.path.join(os.getcwd(), '1_')
        os.makedirs(output_dir, exist_ok=True)

        # Get TaxIDs
        taxids = get_taxids()
        
        # Only ask for TaxID in interactive standalone mode, not when run from GUI
        if interactive_mode and not taxids:
            taxids_input = get_console_input("Enter TaxID number(s) separated by commas: ")
            taxid_list = [taxid.strip() for taxid in taxids_input.split(',')]
        else:
            if not taxids:
                logging.error("No TaxIDs found in configuration. Please provide them in the GUI.")
                return
                
            taxid_list = taxids
            logging.info(f"Using saved TaxIDs: {', '.join(taxid_list)}")
        
        invalid_taxids = [taxid for taxid in taxid_list if not validate_taxid(taxid)]
        if invalid_taxids:
            raise ValidationError(f"Invalid TaxID format: {', '.join(invalid_taxids)}")

        # Stage 1: Gene screening for each taxon
        taxid_results = {}
        for taxid_idx, taxid in enumerate(taxid_list):
            logging.info(f"Processing TaxID {taxid_idx+1} of {len(taxid_list)}: {taxid}")
            organism_name, gene_counts = screen_taxid(taxid, output_dir)
            if organism_name and gene_counts:
                taxid_results[taxid] = (organism_name, gene_counts)

        # Stage 2: Gene selection and download
        for taxid_idx, (taxid, (organism_name, gene_counts)) in enumerate(taxid_results.items()):
            logging.info(f"\n{'='*50}")
            logging.info(f"Results for {organism_name} (TaxID: {taxid}) [{taxid_idx+1}/{len(taxid_results)}]")
            
            sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
            
            # Display genes with sequence counts (for logging)
            logging.info("\nTop genes by sequence count:")
            for i, (gene, count) in enumerate(sorted_genes[:300], 1):
                logging.info(f"{i}. {gene}: {count} sequences")

            if len(sorted_genes) > 300:
                logging.info(f"...and {len(sorted_genes) - 300} more genes (see {taxid}_gene_stats.txt for complete list)")
            
            # Choose selection mode based on arguments
            if interactive_mode:
                # Interactive mode - automatically use option 1 (select specific genes by number)
                logging.info("\n")
                
                selected_genes = []
                # Select specific genes
                selections = get_console_input("Enter gene numbers separated by commas (ranges like '5-7' are allowed): ")
                try:
                    parts = selections.split(',')
                    for part in parts:
                        if '-' in part:
                            start, end = map(int, part.split('-'))
                            for i in range(start, end + 1):
                                if 1 <= i <= len(sorted_genes):
                                    selected_genes.append(sorted_genes[i-1][0])
                        else:
                            i = int(part)
                            if 1 <= i <= len(sorted_genes):
                                selected_genes.append(sorted_genes[i-1][0])
                except ValueError:
                    logging.error("Invalid input format. No genes selected.")            
            else:
                # Non-interactive mode (GUI) - auto select genes
                selected_genes = auto_select_genes(gene_counts)
                logging.info(f"Automatically selected {len(selected_genes)} genes with at least 8 sequences")
                
                # Print selected genes for the log
                if selected_genes:
                    logging.info("Selected genes:")
                    for i, gene in enumerate(selected_genes, 1):
                        count = gene_counts.get(gene, 0)
                        logging.info(f"  {i}. {gene} ({count} sequences)")
                        
            # Download sequences for selected genes
            if selected_genes:
                download_gene_sequences(taxid, organism_name, selected_genes, output_dir)
            else:
                logging.info("No genes selected for download.")

        logging.info("\nAll tasks completed!")

    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        import traceback
        logging.error(traceback.format_exc())

if __name__ == '__main__':
    main()
