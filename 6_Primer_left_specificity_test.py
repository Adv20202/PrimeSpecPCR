import os
import re
import csv
import time
import threading
import logging
import sys
from typing import Dict, List, Optional, Set, Any, Tuple
from dataclasses import dataclass
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
from io import StringIO

class ValidationError(Exception):
    pass

@dataclass
class BlastParameters:
    program: str = "blastn"
    database: str = "nt"
    megablast: bool = False
    expect: float = 100
    word_size: int = 7
    hitlist_size: int = 500
    alignments: int = 500
    descriptions: int = 500
    format_type: str = "XML"
    gapcosts: str = "5 2"
    max_records: int = 500
    max_wait_time: int = 600
    input_dir: str = "5_"
    output_dir: str = "6_"
    min_identity: float = 85.0
    min_coverage: float = 85.0

class BlastThread(threading.Thread):
    def __init__(self, sequence: str, params: BlastParameters):
        super().__init__()
        self.sequence = sequence
        self.params = params
        self.result_handle: Optional[Any] = None
        self.error: Optional[Exception] = None

    def run(self) -> None:
        try:
            self.result_handle = NCBIWWW.qblast(
                self.params.program,
                self.params.database,
                self.sequence,
                megablast=self.params.megablast,
                expect=self.params.expect,
                word_size=self.params.word_size,
                hitlist_size=self.params.hitlist_size,
                alignments=self.params.alignments,
                descriptions=self.params.descriptions,
                format_type=self.params.format_type,
                gapcosts=self.params.gapcosts
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

def find_primer_files(directory: str) -> List[str]:
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist")
    pattern = re.compile(r'^MAFFT_\d{1,2}_consensus_primers\.csv$')
    files = [f for f in os.listdir(directory) if pattern.match(f)]
    return sorted(files)

def get_completed_sequences(output_dir: str) -> Set[str]:
    """Get set of sequences that have already been processed based on existing files."""
    completed = set()
    if os.path.exists(output_dir):
        for filename in os.listdir(output_dir):
            if filename.endswith('.fasta'):
                sequence = filename[:-6].replace('_', '')
                completed.add(sequence)
    return completed

def sanitize_filename(sequence: str) -> str:
    return re.sub(r'[^\w\-_.]', '_', sequence)

def update_progress(prefix: str, elapsed: int) -> None:
    sys.stdout.write(f"\r{prefix} Time elapsed: {elapsed} seconds")
    sys.stdout.flush()

def process_sequence(idx: int, total: int) -> None:
    sys.stdout.write(f"\rProcessing sequence {idx}/{total}")
    sys.stdout.flush()

def fetch_sequence_from_entrez(accession: str, sequence_cache: Dict[str, str]) -> Tuple[str, str]:
    attempts = 0
    while attempts < 3:
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text"
            )
            fasta_data = handle.read()
            handle.close()
            
            record = SeqIO.read(StringIO(fasta_data), "fasta")
            full_sequence = str(record.seq)
            full_seq_id = str(record.id)
            
            sequence_cache[accession] = full_sequence
            sequence_cache[f"{accession}_id"] = full_seq_id
            
            time.sleep(0.1)
            return full_sequence, full_seq_id
            
        except HTTPError as e:
            attempts += 1
            logging.warning(f"HTTPError while fetching {accession}: {e}. Attempt {attempts}/3")
            time.sleep(5)
        except Exception as e:
            logging.error(f"Error fetching sequence {accession}: {e}")
            raise
    
    raise ValidationError(f"Failed to fetch sequence after 3 attempts: {accession}")

def process_blast_results(result_handle: Any, fasta_filename: str, sequence_cache: Dict[str, str], params: BlastParameters) -> None:
    try:
        blast_record = NCBIXML.read(result_handle)
        accession_start_dict: Dict[str, List[int]] = {}
        
        for alignment in blast_record.alignments:
            accession = alignment.accession
            for hsp in alignment.hsps:
                query_coverage = (hsp.align_length / len(blast_record.query)) * 100
                identity = (hsp.identities / hsp.align_length) * 100
                
                if (identity >= params.min_identity and 
                    query_coverage >= params.min_coverage and 
                    hsp.expect <= params.expect):
                    if accession in accession_start_dict:
                        accession_start_dict[accession].append(hsp.sbjct_start)
                    else:
                        accession_start_dict[accession] = [hsp.sbjct_start]

        total_records = len(accession_start_dict)
        logging.info(f"Records found after filtering: {total_records}")

        if total_records == 0:
            return

        if total_records > params.max_records:
            logging.info(f"Limiting to top {params.max_records} records")
            accession_start_dict = dict(list(accession_start_dict.items())[:params.max_records])
            total_records = params.max_records

        with open(fasta_filename, 'w', encoding='utf-8') as fasta_file:
            for idx, (accession, start_positions) in enumerate(accession_start_dict.items(), 1):
                process_sequence(idx, total_records)
                
                if accession in sequence_cache:
                    full_sequence = sequence_cache[accession]
                    full_seq_id = sequence_cache[f"{accession}_id"]
                else:
                    full_sequence, full_seq_id = fetch_sequence_from_entrez(accession, sequence_cache)

                for subject_start in start_positions:
                    seq_start = max(1, subject_start)
                    seq_end = seq_start + 999
                    trimmed_sequence = full_sequence[seq_start - 1:seq_end]
                    fasta_file.write(f">{full_seq_id}_{seq_start}_{seq_end}\n{trimmed_sequence}\n")
                
        print()
        
    except Exception as e:
        logging.error(f"Error processing BLAST results: {e}")
        raise

def perform_blast_with_timer(sequence: str, params: BlastParameters) -> Optional[Any]:
    max_retries = 3
    current_try = 0
    
    while current_try < max_retries:
        logging.info(f"Starting BLAST for sequence: {sequence}")
        
        blast_thread = BlastThread(sequence, params)
        blast_thread.start()
        
        start_time = time.time()
        while blast_thread.is_alive():
            elapsed = int(time.time() - start_time)
            if elapsed >= params.max_wait_time:
                print(f"\nBLAST search timed out (attempt {current_try + 1}/{max_retries})")
                blast_thread.join(1)
                break
            update_progress("BLAST search", elapsed)
            time.sleep(1)
            
        print()
        
        if blast_thread.error:
            logging.error(f"Error in BLAST search (attempt {current_try + 1}/{max_retries}): {blast_thread.error}")
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

def process_primer_file(file_path: str, processed_sequences: Set[str], sequence_cache: Dict[str, str], params: BlastParameters) -> None:
    logging.info(f"Processing file: {file_path}")
    
    try:
        with open(file_path, 'r', newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            if 'PRIMER_LEFT_SEQUENCE' not in reader.fieldnames:
                raise ValidationError(f"Column PRIMER_LEFT_SEQUENCE not found in {file_path}")

            for row in reader:
                sequence = row.get('PRIMER_LEFT_SEQUENCE', '').strip().upper()
                if not sequence:
                    logging.warning("Empty primer sequence found, skipping")
                    continue

                if sequence in processed_sequences:
                    logging.info(f"Sequence {sequence} already processed, skipping")
                    continue

                fasta_filename = os.path.join(params.output_dir, f"{sanitize_filename(sequence)}.fasta")

                if os.path.exists(fasta_filename):
                    logging.info(f"Output file already exists for sequence {sequence}, skipping")
                    processed_sequences.add(sequence)
                    continue

                result_handle = perform_blast_with_timer(sequence, params)
                if result_handle:
                    xml_filename = fasta_filename.replace('.fasta', '_blast_result.xml')
                    with open(xml_filename, 'w', encoding='utf-8') as xml_file:
                        xml_file.write(result_handle.read())
                    result_handle.seek(0)
                    
                    process_blast_results(result_handle, fasta_filename, sequence_cache, params)
                    result_handle.close()
                    processed_sequences.add(sequence)

    except Exception as e:
        logging.error(f"Error processing file {file_path}: {e}")
        raise

def read_user_data() -> Tuple[str, str]:
    try:
        with open(os.path.join('user_data', 'email.txt'), 'r', encoding='utf-8') as email_file:
            email = email_file.read().strip()
        with open(os.path.join('user_data', 'api_key.txt'), 'r', encoding='utf-8') as api_file:
            api_key = api_file.read().strip()
            
        if not email or not api_key:
            raise ValidationError("Email or API key is empty")
            
        return email, api_key
    except Exception as e:
        raise ValidationError(f"Error reading user data: {e}")

def main() -> None:
    try:
        setup_logging(__file__)
        
        email, api_key = read_user_data()
        Entrez.email = email
        Entrez.api_key = api_key

        params = BlastParameters()
        os.makedirs(params.output_dir, exist_ok=True)

        primer_files = find_primer_files(params.input_dir)
        if not primer_files:
            raise ValidationError("No primer files found matching pattern MAFFT_<index>_consensus_primers.csv")

        logging.info(f"Found {len(primer_files)} files to process")

        processed_sequences = get_completed_sequences(params.output_dir)
        if processed_sequences:
            logging.info(f"Found {len(processed_sequences)} already processed sequences")
            
        sequence_cache: Dict[str, str] = {}

        for file in primer_files:
            file_path = os.path.join(params.input_dir, file)
            process_primer_file(file_path, processed_sequences, sequence_cache, params)

        logging.info("Processing completed successfully")

    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

if __name__ == '__main__':
    main()
