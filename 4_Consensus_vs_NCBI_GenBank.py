import os 
import glob
import csv
import time
import threading
import logging
import sys
from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez

class ValidationError(Exception):
    pass

@dataclass
class BlastParameters:
    program: str = "blastn"
    database: str = "nt"
    megablast: bool = True
    hitlist_size: int = 5000
    expect: float = 10
    input_dir: str = "3_"
    output_dir: str = "4_"
    format_type: str = "XML"
    max_wait_time: int = 600
    word_size: int = 11

class BlastThread(threading.Thread):
    def __init__(self, params: BlastParameters, sequence: str, entrez_query: Optional[str] = None):
        super().__init__()
        self.params = params
        self.sequence = sequence
        self.entrez_query = entrez_query
        self.result_handle: Optional[Any] = None
        self.error: Optional[Exception] = None

    def run(self):
        try:
            self.result_handle = NCBIWWW.qblast(
                self.params.program,
                self.params.database,
                self.sequence,
                megablast=self.params.megablast,
                hitlist_size=self.params.hitlist_size,
                entrez_query=self.entrez_query,
                format_type=self.params.format_type,
                expect=self.params.expect,
                word_size=self.params.word_size
            )
        except Exception as e:
            self.error = e

def setup_logging(script_name: str) -> None:
    log_name = f"{os.path.splitext(script_name)[0]}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_name, encoding='utf-8')
        ]
    )

def update_progress(prefix: str, elapsed: int) -> None:
    sys.stdout.write(f"\r{prefix} Time elapsed: {elapsed} seconds")
    sys.stdout.flush()

def find_consensus_files(directory: str) -> List[str]:
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist")
    files = glob.glob(os.path.join(directory, "*consensus*.txt"))
    return sorted(files)

def read_sequence(file_path: str) -> Optional[str]:
    try:
        with open(file_path, 'r', encoding='utf-8') as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if not records:
                logging.warning(f"No sequences found in file {file_path}")
                return None
            return str(records[0].seq)
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        return None

def verify_taxid(accession: str, target_taxid: str) -> bool:
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
        record = Entrez.read(handle)[0]
        handle.close()
        
        features = record.get('GBSeq_feature-table', [])
        for feature in features:
            if feature.get('GBFeature_key') == 'source':
                quals = feature.get('GBFeature_quals', [])
                for qual in quals:
                    if qual.get('GBQualifier_name') == 'db_xref':
                        value = qual.get('GBQualifier_value', '')
                        if value.startswith('taxon:'):
                            organism_taxid = value.split(':')[1]
                            return organism_taxid == target_taxid
        return False
    except Exception as e:
        logging.error(f"Error verifying taxid for accession {accession}: {e}")
        return False

def perform_blast_with_timer(params: BlastParameters, sequence: str, 
                           entrez_query: Optional[str] = None) -> Optional[Any]:
    max_retries = 3
    current_try = 0
    
    while current_try < max_retries:
        blast_thread = BlastThread(params, sequence, entrez_query)
        blast_thread.start()
        
        prefix = f"Performing BLAST search {f'with query: {entrez_query}' if entrez_query else 'without restrictions'}..."
        start_time = time.time()
        
        while blast_thread.is_alive():
            elapsed = int(time.time() - start_time)
            update_progress(prefix, elapsed)
                
            if elapsed >= params.max_wait_time:
                print(f"\nBLAST search timed out (attempt {current_try + 1}/{max_retries})")
                blast_thread.join(1)
                break
                
            time.sleep(1)
        
        print()
        
        if blast_thread.error:
            logging.error(f"Error performing BLAST (attempt {current_try + 1}/{max_retries}): {blast_thread.error}")
            current_try += 1
            if current_try < max_retries:
                logging.info("Retrying BLAST search...")
                time.sleep(5)
            continue
            
        if blast_thread.result_handle:
            return blast_thread.result_handle
            
        current_try += 1
        
    logging.error("All BLAST attempts failed")
    return None

def parse_blast_results(result_handle: Any, taxid: Optional[str] = None, used_accessions: Optional[set] = None) -> List[Dict[str, str]]:
    records = []
    if used_accessions is None:
        used_accessions = set()
        
    try:
        blast_record = NCBIXML.read(result_handle)
        for alignment in blast_record.alignments:
            if alignment.accession in used_accessions:
                continue
                
            for hsp in alignment.hsps:
                species = alignment.hit_def.split('[')[-1].strip(']') if '[' in alignment.hit_def else 'Unknown'
                
                percent_identity = (hsp.identities / hsp.align_length) * 100
                
                if taxid:
                    try:
                        if not verify_taxid(alignment.accession, taxid):
                            logging.debug(f"Skipping {alignment.accession} - different taxid")
                            continue
                    except Exception as e:
                        logging.error(f"Error verifying taxid for {alignment.accession}: {e}")
                        continue

                record = {
                    'Species': species,
                    'Record ID': alignment.accession,
                    'Description': alignment.hit_def,
                    'Max Score': str(hsp.score),
                    'Total Score': str(hsp.bits),
                    'Query Cover': f"{(hsp.align_length / len(blast_record.query)) * 100:.2f}%",
                    'E value': str(hsp.expect),
                    'Per. Ident': f"{percent_identity:.2f}%",
                    'Alignment Length': hsp.align_length
                }
                
                if hsp.expect < 1 and percent_identity > 65.0:
                    records.append(record)
                    used_accessions.add(alignment.accession)
                    break
                    
        return records
    except Exception as e:
        logging.error(f"Error parsing BLAST results: {e}")
        return []

def write_csv(output_dir: str, base_file_name: str, records: List[Dict[str, str]], suffix: str = "") -> None:
    base_name = os.path.splitext(os.path.basename(base_file_name))[0]
    csv_file = f"{base_name}_{suffix}.csv" if suffix else f"{base_name}.csv"
    csv_file_path = os.path.join(output_dir, csv_file)
    
    fieldnames = ['Species', 'Record ID', 'Description', 'Max Score', 'Total Score', 
                 'Query Cover', 'E value', 'Per. Ident']
    
    try:
        with open(csv_file_path, 'w', newline='', encoding='utf-8') as csvhandle:
            writer = csv.DictWriter(csvhandle, fieldnames=fieldnames)
            writer.writeheader()
            for record in records:
                record_copy = record.copy()
                record_copy.pop('Alignment Length', None)
                writer.writerow(record_copy)
        logging.info(f"Results saved to {csv_file_path}")
    except Exception as e:
        logging.error(f"Error saving CSV file {csv_file_path}: {e}")

def summarize_records(taxid_records: List[Dict[str, str]], 
                     no_taxid_records: List[Dict[str, str]]) -> Dict[str, int]:
    summary = {
        'num_taxid': len(taxid_records),
        'num_no_taxid': len(no_taxid_records),
        'difference': len(no_taxid_records) - len(taxid_records),
        'identical_per_ident': sum(1 for record in no_taxid_records 
                                 if record['Per. Ident'] in {r['Per. Ident'] for r in taxid_records})
    }
    return summary

def read_user_data() -> Tuple[str, str]:
    user_data = {'email': '', 'api_key': ''}
    user_data_dir = os.path.join(os.getcwd(), 'user_data')
    
    try:
        with open(os.path.join(user_data_dir, 'email.txt'), 'r', encoding='utf-8') as email_file:
            user_data['email'] = email_file.read().strip()
        with open(os.path.join(user_data_dir, 'api_key.txt'), 'r', encoding='utf-8') as api_file:
            user_data['api_key'] = api_file.read().strip()
        
        if not user_data['email'] or not user_data['api_key']:
            raise ValidationError("Email or API key is empty")
            
        return user_data['email'], user_data['api_key']
    except Exception as e:
        raise ValidationError(f"Error reading user data files: {e}")

def process_consensus_file(file: str, params: BlastParameters, taxid: Optional[str]) -> None:
    logging.info(f"\nProcessing file: {file}")
    
    sequence = read_sequence(file)
    if not sequence:
        return

    used_accessions = set()

    if taxid:
        result_handle_taxid = perform_blast_with_timer(params, sequence, f"txid{taxid}[Organism:exp]")
        
        if result_handle_taxid:
            logging.info("Retrieving BLAST results with TaxID...")
            taxid_records = parse_blast_results(result_handle_taxid, taxid, used_accessions)
            write_csv(params.output_dir, file, taxid_records, f"taxid_{taxid}")
            result_handle_taxid.close()
            
            used_accessions.update(record['Record ID'] for record in taxid_records)
        else:
            logging.warning("No BLAST results with TaxID")
            taxid_records = []

    query = f"NOT txid{taxid}[Organism:exp]" if taxid else None
    result_handle_other = perform_blast_with_timer(params, sequence, query)
    
    if result_handle_other:
        logging.info(f"Retrieving BLAST results {'excluding specified TaxID' if taxid else 'without restrictions'}...")
        other_records = parse_blast_results(result_handle_other, None, used_accessions)
        write_csv(params.output_dir, file, other_records, "no_taxid")
        result_handle_other.close()
        
        if taxid:
            summary = summarize_records(taxid_records, other_records)
            logging.info(f"\n--- Summary for file: {file} ---")
            logging.info(f"Number of records with TaxID {taxid}: {summary['num_taxid']}")
            logging.info(f"Number of records without TaxID: {summary['num_no_taxid']}")
            logging.info(f"Difference: {summary['difference']}")
            logging.info(f"Number of sequences with identical 'Per. Ident': {summary['identical_per_ident']}")
    else:
        logging.warning(f"No BLAST results {'excluding TaxID' if taxid else 'without restrictions'}")
        logging.info(f"\n--- Summary for file: {file} ---")
        logging.info("No records found in second search")

def main():
    try:
        setup_logging(__file__)
        
        email: str
        api_key: str
        email, api_key = read_user_data()
        
        Entrez.email = email
        Entrez.api_key = api_key

        params = BlastParameters()
        os.makedirs(params.output_dir, exist_ok=True)

        taxid_input = input("Enter TaxID to limit BLAST search (or leave blank to skip): ").strip()
        if taxid_input and not taxid_input.isdigit():
            raise ValidationError("Invalid TaxID format. Must be a number.")
        taxid = taxid_input if taxid_input else None

        consensus_files = find_consensus_files(params.input_dir)
        if not consensus_files:
            raise ValidationError("No consensus files found in the input directory")

        logging.info(f"Found {len(consensus_files)} files to process")

        for file in consensus_files:
            process_consensus_file(file, params, taxid)

        logging.info("\nProcessing completed")

    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

if __name__ == '__main__':
    main()
