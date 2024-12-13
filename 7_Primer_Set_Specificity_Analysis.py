import os
import re
import csv
import sys
import logging
from dataclasses import dataclass
from typing import Dict, Optional, List
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('primer_specificity.log', encoding='utf-8')
    ]
)

class ValidationError(Exception):
    pass

@dataclass
class PrimerSet:
    set_number: int
    left_sequence: str
    right_sequence: str
    probe_sequence: str

@dataclass
class AnalysisResult:
    left_annotation: str
    probe_annotation: str
    right_annotation: str

def safe_write_to_file(filepath: str, content: str) -> None:
    """Safely write content to a file, creating directories if needed."""
    try:
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except OSError as e:
        raise IOError(f"Failed to write to {filepath}: {e}")

def find_csv_files(directory: str) -> List[str]:
    """Find all CSV files matching the MAFFT pattern in the directory."""
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist")
    
    pattern = re.compile(r'^MAFFT_\d+.*\.csv$')
    files = [f for f in os.listdir(directory) if pattern.match(f)]
    if not files:
        return []
    files.sort(key=lambda x: int(re.search(r'MAFFT_(\d+)', x).group(1) if re.search(r'MAFFT_(\d+)', x) else 0))
    return files

def sanitize_filename(filename: str) -> str:
    """Sanitize filename to be safe for filesystem."""
    return re.sub(r'[^\w\-_\. ]', '_', filename)

def update_progress(current: int, total: int) -> None:
    """Display progress bar in the console."""
    sys.stdout.write(f'\rProcessing primer set {current}/{total}')
    sys.stdout.flush()

def read_primer_sets(csv_file: str) -> Dict[int, PrimerSet]:
    """Read primer sets from CSV file."""
    primer_sets = {}
    try:
        with open(csv_file, 'r', newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                set_number = int(row.get('PRIMER_SET_NUMBER', '0'))
                primer_sets[set_number] = PrimerSet(
                    set_number=set_number,
                    left_sequence=row.get('PRIMER_LEFT_SEQUENCE', '').strip().upper(),
                    right_sequence=row.get('PRIMER_RIGHT_SEQUENCE', '').strip().upper(),
                    probe_sequence=row.get('PROBE_SEQUENCE', '').strip().upper()
                )
    except (ValueError, KeyError) as e:
        raise ValidationError(f"Invalid CSV format in {csv_file}: {e}")
    return primer_sets

def count_mismatches(sequence: str, fragment: str) -> int:
    """Count minimum number of mismatches between sequence and fragment."""
    if not fragment:
        return 0
    min_mismatches = len(fragment)
    for i in range(len(sequence) - len(fragment) + 1):
        mismatches = sum(1 for a, b in zip(sequence[i:i+len(fragment)], fragment) if a != b)
        if mismatches < min_mismatches:
            min_mismatches = mismatches
        if min_mismatches == 0:
            break
    return min_mismatches

def analyze_primer(sequence: str, primer_seq: str, is_left: bool = True) -> str:
    """Analyze primer sequence against target sequence."""
    primer_seq = primer_seq.upper()
    sequence = sequence.upper()

    if primer_seq in sequence:
        return 'A'

    if is_left:
        primer_fragment = primer_seq[-14:]
        fragment_last4 = primer_fragment[-4:]
        fragment_next10 = primer_fragment[:-4]
        remaining_seq = primer_seq[:-14]
    else:
        primer_fragment = primer_seq[:14]
        fragment_last4 = primer_fragment[:4]
        fragment_next10 = primer_fragment[4:]
        remaining_seq = primer_seq[14:]

    annotation = ''

    annotation += 'B' if fragment_last4 in sequence else 'C'

    mismatches = count_mismatches(sequence, fragment_next10)
    if mismatches == 0:
        annotation += 'A'
    elif mismatches == 1:
        annotation += 'B'
    else:
        annotation += 'C'

    if remaining_seq:
        mismatches = count_mismatches(sequence, remaining_seq)
        if mismatches == 0:
            annotation += 'A'
        elif mismatches <= 2:
            annotation += 'B'
        else:
            annotation += 'C'

    return annotation

def analyze_probe(sequence: str, probe_seq: str) -> str:
    """Analyze probe sequence against target sequence."""
    probe_seq = probe_seq.upper()
    sequence = sequence.upper()

    if probe_seq in sequence:
        return 'A'

    first7 = probe_seq[:7]
    rest = probe_seq[7:]

    annotation = 'B' if rest in sequence else 'C'
    
    mismatches = count_mismatches(sequence, first7)
    if mismatches == 0:
        annotation += 'A'
    elif mismatches == 1:
        annotation += 'B'
    else:
        annotation += 'C'

    return annotation

def analyze_sequence(sequence: str, primer_set: PrimerSet) -> AnalysisResult:
    """Analyze a sequence against a primer set."""
    right_primer_rc = str(Seq(primer_set.right_sequence).reverse_complement())
    
    return AnalysisResult(
        left_annotation=analyze_primer(sequence, primer_set.left_sequence, is_left=True),
        probe_annotation=analyze_probe(sequence, primer_set.probe_sequence),
        right_annotation=analyze_primer(sequence, right_primer_rc, is_left=False)
    )

def process_primer_set(primer_set: PrimerSet, fasta_dir: str, output_dir: str) -> Optional[bool]:
    """Process a single primer set and create output FASTA file."""
    sanitized_left_seq = sanitize_filename(primer_set.left_sequence)
    fasta_filename = os.path.join(fasta_dir, f"{sanitized_left_seq}.fasta")
    
    if not os.path.exists(fasta_filename):
        return None
        
    output_fasta = os.path.join(output_dir, f"primer_set_{primer_set.set_number}.fasta")
    
    try:
        sequences = list(SeqIO.parse(fasta_filename, 'fasta'))
        for seq_record in sequences:
            analysis_result = analyze_sequence(str(seq_record.seq), primer_set)
            seq_record.id = f"{seq_record.id}__{analysis_result.left_annotation}_{analysis_result.probe_annotation}_{analysis_result.right_annotation}"
            seq_record.description = ''
        
        SeqIO.write(sequences, output_fasta, 'fasta')
        
        with open(output_fasta, 'r', encoding='utf-8') as f:
            headers = [line for line in f if line.startswith('>')]
        with open(output_fasta, 'w', encoding='utf-8') as f:
            f.writelines(headers)
            
        return True
    except Exception as e:
        logging.error(f"Error processing primer set {primer_set.set_number}: {e}")
        return False

def main():
    try:
        input_dir = Path(os.getcwd()) / '5_'
        fasta_dir = Path(os.getcwd()) / '6_'
        output_dir = Path(os.getcwd()) / '7_'
        output_dir.mkdir(exist_ok=True)

        csv_files = find_csv_files(str(input_dir))
        if not csv_files:
            raise ValidationError('No CSV files starting with "MAFFT_" found in directory "5_"')

        total_missing_files = 0
        for csv_index, csv_file in enumerate(csv_files, 1):
            logging.info(f"\nProcessing CSV file {csv_index}/{len(csv_files)}: {csv_file}")
            
            try:
                primer_sets = read_primer_sets(str(input_dir / csv_file))
                missing_files = 0
                
                for set_index, (set_number, primer_set) in enumerate(sorted(primer_sets.items()), 1):
                    update_progress(set_index, len(primer_sets))
                    
                    if not all([primer_set.left_sequence, primer_set.right_sequence, primer_set.probe_sequence]):
                        continue
                        
                    result = process_primer_set(primer_set, str(fasta_dir), str(output_dir))
                    if result is None:
                        missing_files += 1
                
                total_missing_files += missing_files
                if missing_files > 0:
                    logging.info(f"\nMissing FASTA files in {csv_file}: {missing_files}")
                
            except Exception as e:
                logging.error(f"Error processing CSV file {csv_file}: {e}")
                continue

        if total_missing_files > 0:
            logging.info(f"\nTotal missing FASTA files: {total_missing_files}")
        else:
            logging.info("\nAll FASTA files were processed successfully.")

    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

if __name__ == '__main__':
    main()
