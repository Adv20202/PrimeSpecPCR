import os
import glob
import primer3
import csv
import time
import re
import logging
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from tqdm import tqdm

@dataclass
class PrimerParameters:
    min_size: int
    max_size: int
    opt_size: int = 20
    min_tm: float = 58.0
    opt_tm: float = 60.0
    max_tm: float = 62.0
    min_gc: float = 30.0
    max_gc: float = 70.0
    max_poly_x: int = 5
    salt_monovalent: float = 50.0
    salt_divalent: float = 1.5
    dna_conc: float = 50.0
    max_ns_accepted: int = 0
    input_dir: str = "3_"
    output_dir: str = "5_"
    max_consecutive_failures: int = 10

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

def find_consensus_files(directory: str) -> List[str]:
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist")
    files = glob.glob(os.path.join(directory, "*consensus*.txt"))
    return sorted(files)

def read_sequence(file_path: str) -> Optional[str]:
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            sequence = ""
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                sequence += line.upper()
            return sequence if sequence else None
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        return None

def get_primer_parameters() -> Optional[PrimerParameters]:
    try:
        user_input = input("Enter amplicon length range (e.g., 100-250): ").strip()
        match = re.match(r'^(\d+)-(\d+)$', user_input)
        if not match:
            raise ValidationError("Invalid format. Use 'min-max', e.g., '100-250'.")
        
        min_size, max_size = map(int, match.groups())
        if min_size >= max_size:
            raise ValidationError("Minimum size must be less than maximum size.")
            
        return PrimerParameters(min_size=min_size, max_size=max_size)
    except ValidationError as e:
        logging.error(f"Validation error: {e}")
        return None
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        return None

def get_primer3_settings(params: PrimerParameters) -> Dict:
    return {
        'PRIMER_OPT_SIZE': params.opt_size,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SIZE': 25,
        'PRIMER_INTERNAL_MIN_SIZE': 18,
        'PRIMER_INTERNAL_OPT_SIZE': 20,
        'PRIMER_NUM_RETURN': 2000,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 24,
        'PRIMER_MIN_TM': params.min_tm,
        'PRIMER_OPT_TM': params.opt_tm,
        'PRIMER_MAX_TM': params.max_tm,
        'PRIMER_MIN_GC': params.min_gc,
        'PRIMER_MAX_GC': params.max_gc,
        'PRIMER_MAX_POLY_X': params.max_poly_x,
        'PRIMER_INTERNAL_MAX_POLY_X': 4,
        'PRIMER_SALT_MONOVALENT': params.salt_monovalent,
        'PRIMER_SALT_DIVALENT': params.salt_divalent,
        'PRIMER_DNA_CONC': params.dna_conc,
        'PRIMER_MAX_NS_ACCEPTED': params.max_ns_accepted,
        'PRIMER_PRODUCT_SIZE_RANGE': [[params.min_size, params.max_size]],
        'PRIMER_EXPLAIN_FLAG': 1,
        'PRIMER_MAX_SELF_ANY_TH': 35.0,
        'PRIMER_MAX_SELF_END_TH': 35.0,
        'PRIMER_MAX_HAIRPIN_TH': 24.0,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 35.0,
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
    }

def design_primers(sequence: str, settings: Dict, excluded_regions: List[List[int]]) -> Optional[Dict]:
    try:
        sequence_args = {
            'SEQUENCE_ID': 'consensus',
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': [0, len(sequence)],
        }
        
        if excluded_regions:
            sequence_args['SEQUENCE_EXCLUDED_REGION'] = excluded_regions
        
        start_time = time.time()
        result = primer3.bindings.designPrimers(sequence_args, settings)
        design_time = time.time() - start_time
        
        if design_time > 1:
            logging.info(f"Primer design took {design_time:.2f} seconds.")
            
        return result
    except Exception as e:
        logging.debug(f"Primer design iteration ended: {e}")
        return None

def process_primer_set(result: Dict, i: int) -> Optional[Dict]:
    primer_left_seq = result.get(f'PRIMER_LEFT_{i}_SEQUENCE', '')
    primer_right_seq = result.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')
    probe_seq = result.get(f'PRIMER_INTERNAL_{i}_SEQUENCE', '')
    
    if not all([primer_left_seq, primer_right_seq, probe_seq]):
        return None
        
    if len({primer_left_seq, primer_right_seq, probe_seq}) != 3:
        return None
        
    return {
        'SEQUENCE_ID': 'Consensus_sequence',
        'PRIMER_SET_NUMBER': i + 1,
        'PRIMER_LEFT_SEQUENCE': primer_left_seq,
        'PRIMER_LEFT_TM': result.get(f'PRIMER_LEFT_{i}_TM', ''),
        'PRIMER_LEFT_GC%': result.get(f'PRIMER_LEFT_{i}_GC_PERCENT', ''),
        'PRIMER_LEFT_SELF_ANY_TH': result.get(f'PRIMER_LEFT_{i}_SELF_ANY_TH', ''),
        'PRIMER_LEFT_SELF_END_TH': result.get(f'PRIMER_LEFT_{i}_SELF_END_TH', ''),
        'PRIMER_LEFT_HAIRPIN_TH': result.get(f'PRIMER_LEFT_{i}_HAIRPIN_TH', ''),
        'PRIMER_RIGHT_SEQUENCE': primer_right_seq,
        'PRIMER_RIGHT_TM': result.get(f'PRIMER_RIGHT_{i}_TM', ''),
        'PRIMER_RIGHT_GC%': result.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', ''),
        'PRIMER_RIGHT_SELF_ANY_TH': result.get(f'PRIMER_RIGHT_{i}_SELF_ANY_TH', ''),
        'PRIMER_RIGHT_SELF_END_TH': result.get(f'PRIMER_RIGHT_{i}_SELF_END_TH', ''),
        'PRIMER_RIGHT_HAIRPIN_TH': result.get(f'PRIMER_RIGHT_{i}_HAIRPIN_TH', ''),
        'PROBE_SEQUENCE': probe_seq,
        'PROBE_TM': result.get(f'PRIMER_INTERNAL_{i}_TM', ''),
        'PROBE_GC%': result.get(f'PRIMER_INTERNAL_{i}_GC_PERCENT', ''),
        'PROBE_SELF_ANY_TH': result.get(f'PRIMER_INTERNAL_{i}_SELF_ANY_TH', ''),
        'PROBE_SELF_END_TH': result.get(f'PRIMER_INTERNAL_{i}_SELF_END_TH', ''),
        'PROBE_HAIRPIN_TH': result.get(f'PRIMER_INTERNAL_{i}_HAIRPIN_TH', ''),
        'PRODUCT_SIZE': result.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', '')
    }

def save_primers_to_csv(output_path: str, primer_data: List[Dict]) -> None:
    if not primer_data:
        logging.warning("No primer data to save")
        return
        
    headers = [
        'SEQUENCE_ID', 'PRIMER_SET_NUMBER',
        'PRIMER_LEFT_SEQUENCE', 'PRIMER_LEFT_TM', 'PRIMER_LEFT_GC%',
        'PRIMER_LEFT_SELF_ANY_TH', 'PRIMER_LEFT_SELF_END_TH', 'PRIMER_LEFT_HAIRPIN_TH',
        'PRIMER_RIGHT_SEQUENCE', 'PRIMER_RIGHT_TM', 'PRIMER_RIGHT_GC%',
        'PRIMER_RIGHT_SELF_ANY_TH', 'PRIMER_RIGHT_SELF_END_TH', 'PRIMER_RIGHT_HAIRPIN_TH',
        'PROBE_SEQUENCE', 'PROBE_TM', 'PROBE_GC%',
        'PROBE_SELF_ANY_TH', 'PROBE_SELF_END_TH', 'PROBE_HAIRPIN_TH',
        'PRODUCT_SIZE'
    ]
    
    try:
        with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()
            writer.writerows(primer_data)
        logging.info(f"Results saved to: {output_path}")
        logging.info(f"Saved {len(primer_data)} primer sets")
    except Exception as e:
        logging.error(f"Error saving CSV file {output_path}: {e}")

def process_consensus_file(file_path: str, params: PrimerParameters, primer3_settings: Dict) -> None:
    logging.info(f"\nProcessing file: {file_path}")
    
    sequence = read_sequence(file_path)
    if not sequence:
        logging.warning(f"No valid sequence in file: {file_path}")
        return
        
    primer_data: List[Dict] = []
    designed_sets = set()
    excluded_regions: List[List[int]] = []
    consecutive_failures = 0
    
    output_path = os.path.join(
        params.output_dir,
        f"{os.path.splitext(os.path.basename(file_path))[0]}_primers.csv"
    )
    
    with tqdm(total=params.max_consecutive_failures) as pbar:
        while consecutive_failures < params.max_consecutive_failures:
            try:
                result = design_primers(sequence, primer3_settings, excluded_regions)
                if not result:
                    logging.info("No more primer sets can be designed")
                    break
                    
                left_primers = [k for k in result.keys() if re.match(r'PRIMER_LEFT_\d+_SEQUENCE', k)]
                if not left_primers:
                    consecutive_failures += 1
                    pbar.update(1)
                    continue
                    
                new_set_found = False
                for i in range(len(left_primers)):
                    primer_set = process_primer_set(result, i)
                    if not primer_set:
                        continue
                        
                    set_tuple = (
                        primer_set['PRIMER_LEFT_SEQUENCE'],
                        primer_set['PRIMER_RIGHT_SEQUENCE'],
                        primer_set['PROBE_SEQUENCE']
                    )
                    
                    if set_tuple in designed_sets:
                        continue
                    
                    primer_data.append(primer_set)
                    designed_sets.add(set_tuple)
                    new_set_found = True
                    
                    for pos_info in [
                        result.get(f'PRIMER_LEFT_{i}'),
                        result.get(f'PRIMER_RIGHT_{i}'),
                        result.get(f'PRIMER_INTERNAL_{i}')
                    ]:
                        if pos_info:
                            excluded_regions.append([pos_info[0], pos_info[1]])
                            
                if new_set_found:
                    consecutive_failures = 0
                    save_primers_to_csv(output_path, primer_data)
                else:
                    consecutive_failures += 1
                    logging.info(f"No new sets found. Attempt {consecutive_failures} of {params.max_consecutive_failures}")
                
                pbar.update(1)
                
            except Exception as e:
                logging.debug(f"Iteration error: {e}")
                break
    
    if primer_data:
        save_primers_to_csv(output_path, primer_data)
    else:
        logging.warning(f"Failed to design primer sets for file: {file_path}")

def main():
    try:
        setup_logging(__file__)
        
        params = get_primer_parameters()
        if not params:
            raise ValidationError("Failed to get valid primer parameters")
            
        os.makedirs(params.output_dir, exist_ok=True)
        
        primer3_settings = get_primer3_settings(params)
        
        consensus_files = find_consensus_files(params.input_dir)
        if not consensus_files:
            raise ValidationError("No consensus files found in the input directory")
            
        logging.info(f"Found {len(consensus_files)} files to process")
        
        for file in consensus_files:
            process_consensus_file(file, params, primer3_settings)
            
        logging.info("\nPrimer design completed successfully")
            
    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

if __name__ == '__main__':
    main()
