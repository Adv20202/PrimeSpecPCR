import os
import glob
import logging
from typing import Dict, List, Optional
from dataclasses import dataclass

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('3_Consensus.log', encoding='utf-8')
    ]
)

class ValidationError(Exception):
    pass

@dataclass
class ConsensusParameters:
    threshold: float
    line_length: int = 60

def safe_write_to_file(filepath: str, content: str) -> None:
    try:
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except OSError as e:
        raise IOError(f"Failed to write to {filepath}: {e}")

def find_mafft_files(directory: str) -> List[str]:
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist")
    
    files = glob.glob(os.path.join(directory, "MAFFT*.txt"))
    return sorted(files)

def read_fasta_from_mafft(file_path: str) -> Optional[List[str]]:
    sequences = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            title = ''
            seq = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if title and seq:
                        sequences.append(seq)
                    title = line
                    seq = ''
                else:
                    seq += line.upper()
            if title and seq:
                sequences.append(seq)
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        return None
    return sequences

def create_consensus(sequences: List[str], params: ConsensusParameters) -> Optional[str]:
    if not sequences:
        return None

    num_sequences = len(sequences)
    length = len(sequences[0])

    # Validate sequence lengths
    if any(len(seq) != length for seq in sequences):
        logging.error("Not all sequences have the same length")
        return None

    consensus = []
    for pos in range(length):
        nucleotides: Dict[str, int] = {}
        for seq in sequences:
            nuc = seq[pos]
            if nuc in ['A', 'C', 'G', 'T']:
                nucleotides[nuc] = nucleotides.get(nuc, 0) + 1
        
        if not nucleotides:
            consensus.append('N')
            continue
            
        nucleotide_items = nucleotides.items()
        most_common = max(nucleotide_items, key=lambda item: item[1])
        if most_common[1] / num_sequences >= params.threshold:
            consensus.append(most_common[0])
        else:
            consensus.append('N')
            
    return ''.join(consensus)

def format_sequence(sequence: str, line_length: int) -> str:
    return '\n'.join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))

def save_consensus(input_file: str, consensus: str, output_dir: str, params: ConsensusParameters) -> None:
    output_filename = f"{os.path.splitext(os.path.basename(input_file))[0]}_consensus.txt"
    output_path = os.path.join(output_dir, output_filename)
    
    try:
        formatted_sequence = format_sequence(consensus, params.line_length)
        content = f">Consensus_sequence\n{formatted_sequence}\n"
        safe_write_to_file(output_path, content)
        logging.info(f"Consensus sequence saved to: {output_path}")
    except Exception as e:
        logging.error(f"Error saving consensus to {output_path}: {e}")
        raise

def process_mafft_file(input_file: str, output_dir: str, params: ConsensusParameters) -> None:
    logging.info(f"\nProcessing file: {input_file}")
    
    sequences = read_fasta_from_mafft(input_file)
    if not sequences:
        logging.warning(f"No valid sequences in {input_file}")
        return
        
    logging.info(f"Found {len(sequences)} sequences")
    
    consensus = create_consensus(sequences, params)
    if not consensus:
        logging.error(f"Failed to create consensus for {input_file}")
        return
        
    logging.info(f"Created consensus sequence (length: {len(consensus)})")
    save_consensus(input_file, consensus, output_dir, params)

def main():
    try:
        input_dir = os.path.join(os.getcwd(), '2_')
        output_dir = os.path.join(os.getcwd(), '3_')
        os.makedirs(output_dir, exist_ok=True)

        params = ConsensusParameters(threshold=0.8)

        mafft_files = find_mafft_files(input_dir)
        if not mafft_files:
            raise ValidationError("No MAFFT*.txt files found in the input directory")

        logging.info(f"Found {len(mafft_files)} MAFFT files for analysis")

        for file in mafft_files:
            process_mafft_file(file, output_dir, params)

        logging.info("\nConsensus analysis completed successfully")

    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

if __name__ == '__main__':
    main()
