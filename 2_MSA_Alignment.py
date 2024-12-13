import os
import glob
import subprocess
import logging
from typing import List, Tuple, Optional
from dataclasses import dataclass
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('2_MSA_Alignment.log', encoding='utf-8')
    ]
)

class ValidationError(Exception):
    pass

@dataclass
class SequenceGroup:
    indices: List[int]
    sequences: List[Tuple[str, str]]
    reference_idx: Optional[int] = None

def safe_write_to_file(filepath: str, content: str) -> None:
    try:
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except OSError as e:
        raise IOError(f"Failed to write to {filepath}: {e}")

def find_fasta_files(directory: str) -> List[str]:
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist")
    
    extensions = ['*.fasta', '*.fa', '*.fna', '*.ffn', '*.faa', '*.frn']
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(directory, ext)))
    return sorted(files)

def validate_group_indices(groups: List[str], num_files: int) -> List[List[int]]:
    group_indices: List[List[int]] = []
    
    for group in groups:
        indices_str = group.strip().split(',')
        try:
            indices = [int(i.strip()) for i in indices_str]
            if any(i < 0 or i >= num_files for i in indices):
                raise ValidationError(f"Index out of range in group: {group}")
            if len(set(indices)) != len(indices):
                raise ValidationError(f"Duplicate indices in group: {group}")
            group_indices.append(indices)
        except ValueError:
            raise ValidationError(f"Invalid number format in group: {group}")
    
    return group_indices

def read_fasta(file: str) -> List[Tuple[str, str]]:
    try:
        sequences = []
        for record in SeqIO.parse(file, "fasta"):
            sequences.append((str(record.description), str(record.seq)))
        return sequences
    except Exception as e:
        raise IOError(f"Error reading FASTA file {file}: {e}")

def validate_reference_index(reference_idx: int, num_sequences: int) -> bool:
    return 0 <= reference_idx < num_sequences

def prepare_temporary_fasta(sequences: List[Tuple[str, str]], reference_idx: int, filename: str) -> None:
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            for idx, (title, seq) in enumerate(sequences):
                header = f"{title}_ref" if idx == reference_idx else title
                f.write(f">{header}\n{seq}\n")
    except Exception as e:
        raise IOError(f"Error writing temporary FASTA file: {e}")

def run_mafft(input_file: str, output_file: str, retries: int = 3) -> Optional[str]:
    for attempt in range(retries):
        try:
            result = subprocess.run(
                ['mafft', '--quiet', input_file],
                capture_output=True,
                text=True,
                check=True
            )
            safe_write_to_file(output_file, result.stdout)
            return result.stdout
        except subprocess.CalledProcessError as e:
            logging.error(f"MAFFT error on attempt {attempt + 1}: {e.stderr}")
            if attempt == retries - 1:
                return None
        except Exception as e:
            logging.error(f"Unexpected error running MAFFT on attempt {attempt + 1}: {e}")
            if attempt == retries - 1:
                return None
    return None

def process_sequence_group(group: SequenceGroup, output_dir: str) -> Optional[str]:
    try:
        logging.info(f"\nProcessing group with indices: {group.indices}")
        
        filtered_sequences = [(title, seq) for title, seq in group.sequences if len(seq) < 3000]
        if not filtered_sequences:
            logging.info("No sequences shorter than 3000 nucleotides in this group. Skipping MSA analysis.")
            return None
        
        logging.info(f"Using {len(filtered_sequences)} sequences for MSA analysis")
        
        temp_fasta = os.path.join(output_dir, "temp_group.fasta")
        prepare_temporary_fasta(filtered_sequences, group.reference_idx or 0, temp_fasta)
        
        output_file = os.path.join(output_dir, f"MAFFT_{'_'.join(map(str, group.indices))}.txt")
        result = run_mafft(temp_fasta, output_file)
        
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
            
        if result:
            logging.info(f"MSA analysis saved to: {output_file}")
            return result
        else:
            logging.error("MSA analysis failed")
            return None
            
    except Exception as e:
        logging.error(f"Error processing sequence group: {e}")
        return None

def main():
    try:
        input_dir = os.path.join(os.getcwd(), '1_')
        output_dir = os.path.join(os.getcwd(), '2_')
        os.makedirs(output_dir, exist_ok=True)

        fasta_files = find_fasta_files(input_dir)
        if not fasta_files:
            raise ValidationError("No FASTA files found in the input directory")

        logging.info("\nAvailable FASTA files:")
        for idx, file in enumerate(fasta_files):
            logging.info(f"{idx}: {os.path.basename(file)}")

        groups_input = input("\nEnter groups of file indices for MSA analysis, separated by semicolons (e.g., 0,1,3;2,4,5): ")
        groups = validate_group_indices(groups_input.strip().split(';'), len(fasta_files))

        msa_results = []
        for group_indices in groups:
            group_sequences = []
            for idx in group_indices:
                sequences = read_fasta(fasta_files[idx])
                group_sequences.extend(sequences)

            logging.info("\nAvailable sequences (sorted by length):")
            sorted_sequences = sorted(group_sequences, key=lambda x: len(x[1]))
            for idx, (title, seq) in enumerate(sorted_sequences):
                logging.info(f"{idx}: {title} (length: {len(seq)})")

            while True:
                try:
                    ref_idx = int(input(f"\nChoose reference sequence (0-{len(sorted_sequences)-1}): "))
                    if validate_reference_index(ref_idx, len(sorted_sequences)):
                        break
                    logging.error("Invalid reference index")
                except ValueError:
                    logging.error("Please enter a valid number")

            group = SequenceGroup(group_indices, sorted_sequences, ref_idx)
            result = process_sequence_group(group, output_dir)
            if result:
                msa_results.append((group_indices, result))

        if msa_results:
            logging.info("\nMSA analysis completed successfully")
            for indices, result in msa_results:
                logging.info(f"\nResults for group {indices}:")
                logging.info(result[:500] + "..." if len(result) > 500 else result)
        else:
            logging.warning("\nNo successful MSA analyses to report")

    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

if __name__ == '__main__':
    main()
