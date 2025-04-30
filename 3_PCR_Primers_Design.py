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
import sys
import traceback

class ValidationError(Exception):
    pass

@dataclass
class PrimerParameters:
    min_size: int
    max_size: int
    max_primer_sets: int = 10
    opt_size: int = 20
    min_tm: float = 57.0
    opt_tm: float = 59.0
    max_tm: float = 62.0
    min_gc: float = 30.0
    max_gc: float = 70.0
    max_poly_x: int = 4
    salt_monovalent: float = 50.0
    salt_divalent: float = 1.5
    dna_conc: float = 50.0
    max_ns_accepted: int = 0
    input_dir: str = "2_"
    output_dir: str = "3_"
    max_consecutive_failures: int = 10


def setup_logging(script_name: str) -> None:
    log_name = f"{os.path.splitext(script_name)[0]}.log"
    
    # Create a custom formatter without timestamp and level
    class CustomFormatter(logging.Formatter):
        def format(self, record):
            return record.getMessage()
    
    # Set up handlers with the custom formatter
    formatter = CustomFormatter()
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    
    file_handler = logging.FileHandler(log_name, encoding='utf-8')
    file_handler.setFormatter(formatter)
    
    # Configure root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Remove any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Add our custom handlers
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

def find_consensus_files(directory: str) -> List[str]:
    """
    Find CONSENSUS fasta files in the specified directory.
    Raises FileNotFoundError if directory doesn't exist.
    """
    if not os.path.exists(directory):
        logging.error(f"Directory {directory} does not exist")
        raise FileNotFoundError(f"Directory {directory} does not exist")
    
    files = glob.glob(os.path.join(directory, "*CONSENSUS*.fasta"))
    if not files:
        logging.warning(f"No consensus files found in {directory}")
    
    return sorted(files)

def read_sequence(file_path: str) -> Optional[Tuple[str, str]]:
    """
    Read DNA sequence from FASTA file.
    Returns a tuple of (sequence, header) or (None, header) if no valid sequence.
    """
    try:
        # Check if file exists
        if not os.path.isfile(file_path):
            logging.error(f"File does not exist: {file_path}")
            return (None, "")
            
        # First read the header to extract for display
        header = ""
        with open(file_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
        
        # Now read the actual sequence
        with open(file_path, 'r', encoding='utf-8') as f:
            sequence = ""
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                sequence += line.upper()
                
            # Handle non-standard nucleotides by replacing them with N
            sequence = re.sub(r'[^ACGTN]', 'N', sequence)
            
            # Verify sequence is not empty
            if not sequence:
                logging.error(f"No valid sequence found in {file_path}")
                return (None, header)
                
            # Remove any digits, spaces or other non-DNA characters
            sequence = re.sub(r'[^ACGTN]', '', sequence)
            
            # Log sequence length for debugging
            logging.info(f"Read sequence of length {len(sequence)} bp")
            
            # For direct comparison with Primer3 web output, log first 20 and last 20 bases
            if len(sequence) >= 40:
                logging.info(f"Sequence starts with: {sequence[:20]}")
                logging.info(f"Sequence ends with: {sequence[-20:]}")
            else:
                logging.info(f"Full sequence: {sequence}")
            
            return (sequence, header)
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        logging.error(traceback.format_exc())
        return (None, "")

def get_primer_parameters(interactive_mode=False, file_path="") -> Optional[PrimerParameters]:
    """
    Get primer design parameters from the user or defaults.
    In interactive mode, displays file info before asking for input.
    """
    try:
        if interactive_mode:
            # Display file name before asking for input
            logging.info(f"Processing file: {file_path}")
            
            # Read and display the header
            with open(file_path, 'r', encoding='utf-8') as f:
                header = f.readline().strip()
                logging.info(header)
            
            # Get amplicon size range
            user_input = get_console_input("Enter amplicon length range (e.g., 100-250): ").strip()
        else:
            user_input = input("Enter amplicon length range (e.g., 100-250): ").strip()
            
        match = re.match(r'^(\d+)-(\d+)$', user_input)
        if not match:
            logging.error("Invalid format. Use 'min-max', e.g., '100-250'.")
            return None
        
        min_size, max_size = map(int, match.groups())
        if min_size >= max_size:
            logging.error("Minimum size must be less than maximum size.")
            return None
        
        # Ask for maximum number of primer sets
        if interactive_mode:
            max_sets_input = get_console_input("Enter maximum number of primer sets to generate: ").strip()
        else:
            max_sets_input = input("Enter maximum number of primer sets to generate: ").strip()
            
        try:
            max_primer_sets = int(max_sets_input)
            if max_primer_sets <= 0:
                logging.warning("Invalid number, using default of 10 primer sets")
                max_primer_sets = 10
            # Cap the maximum number to prevent GUI freezing
            if max_primer_sets > 2000:
                logging.warning(f"Limiting to maximum of 2000 primer sets to prevent performance issues")
                max_primer_sets = 2000
        except ValueError:
            logging.warning("Invalid input, using default of 10 primer sets")
            max_primer_sets = 10
            
        # Create output directory
        output_dir = "3_"
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
                logging.info(f"Created output directory: {output_dir}")
            except Exception as e:
                logging.error(f"Failed to create output directory: {e}")
                return None
            
        # Return parameters matching Primer3 web defaults
        return PrimerParameters(
            min_size=min_size, 
            max_size=max_size, 
            max_primer_sets=max_primer_sets,
            opt_size=20,
            min_tm=57.0,
            opt_tm=59.0,
            max_tm=62.0,
            min_gc=30.0,
            max_gc=70.0,
            max_poly_x=4,
            output_dir=output_dir
        )
    except ValidationError as e:
        logging.error(f"Validation error: {e}")
        return None
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        logging.error(traceback.format_exc())
        return None

def load_primer_settings(settings_file: str = "PCR_primer_settings.txt") -> Dict:
    """
    Load primer design settings from settings file.
    Returns a dictionary with the settings.
    """
    settings = {}
    
    # Check if the settings file exists
    if not os.path.exists(settings_file):
        logging.warning(f"Settings file {settings_file} not found. Creating with default values.")
        create_default_settings_file(settings_file)
    
    try:
        with open(settings_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                    
                # Parse setting key and value
                if '=' in line:
                    key, value = [part.strip() for part in line.split('=', 1)]
                    
                    # Convert value to appropriate type
                    if value.replace('.', '', 1).isdigit():
                        if '.' in value:
                            settings[key] = float(value)
                        else:
                            settings[key] = int(value)
                    elif value.lower() in ['true', 'false']:
                        settings[key] = value.lower() == 'true'
                    else:
                        settings[key] = value
    
    except Exception as e:
        logging.error(f"Error loading settings from {settings_file}: {e}")
        logging.error(traceback.format_exc())
        # Return empty dict so caller can use defaults
        return {}
    
    return settings

def create_default_settings_file(settings_file: str) -> None:
    """
    Create a default settings file with description and default values.
    """
    default_settings = """# PCR Primer Design Settings
# Modify values as needed, but maintain the format: setting_name = value

# Primer size parameters
PRIMER_OPT_SIZE = 20
PRIMER_MIN_SIZE = 18
PRIMER_MAX_SIZE = 23

# Tm parameters (melting temperature)
PRIMER_MIN_TM = 57.0
PRIMER_OPT_TM = 59.0
PRIMER_MAX_TM = 62.0
PRIMER_MAX_DIFF_TM = 5.0

# GC content parameters
PRIMER_MIN_GC = 30.0
PRIMER_OPT_GC_PERCENT = 50.0
PRIMER_MAX_GC = 70.0

# Probe parameters (internal oligo)
PRIMER_PICK_INTERNAL_OLIGO = 1
PRIMER_INTERNAL_MIN_SIZE = 18
PRIMER_INTERNAL_OPT_SIZE = 20
PRIMER_INTERNAL_MAX_SIZE = 27
PRIMER_INTERNAL_MIN_TM = 57.0
PRIMER_INTERNAL_OPT_TM = 60.0
PRIMER_INTERNAL_MAX_TM = 63.0
PRIMER_INTERNAL_MIN_GC = 20.0
PRIMER_INTERNAL_OPT_GC_PERCENT = 50.0
PRIMER_INTERNAL_MAX_GC = 80.0

# Poly-X parameters
PRIMER_MAX_POLY_X = 4
PRIMER_INTERNAL_MAX_POLY_X = 5

# Salt and DNA concentration parameters
PRIMER_SALT_MONOVALENT = 50.0
PRIMER_SALT_DIVALENT = 1.5
PRIMER_DNA_CONC = 50.0
PRIMER_MAX_NS_ACCEPTED = 0
PRIMER_SALT_CORRECTIONS = 1

# Thermodynamic alignment parameters
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 1
PRIMER_MAX_SELF_ANY_TH = 45.0
PRIMER_MAX_SELF_END_TH = 35.0
PRIMER_PAIR_MAX_COMPL_ANY_TH = 45.0
PRIMER_PAIR_MAX_COMPL_END_TH = 35.0
PRIMER_MAX_HAIRPIN_TH = 24.0

# Probe thermodynamic parameters
PRIMER_INTERNAL_MAX_SELF_ANY_TH = 47.0
PRIMER_INTERNAL_MAX_SELF_END_TH = 47.0
PRIMER_INTERNAL_MAX_HAIRPIN_TH = 47.0

# Template mispriming
PRIMER_MAX_TEMPLATE_MISPRIMING_TH = 40.0
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH = 70.0

# Additional parameters
PRIMER_EXPLAIN_FLAG = 1
PRIMER_LIBERAL_BASE = 1
PRIMER_FIRST_BASE_INDEX = 1
PRIMER_MAX_END_STABILITY = 9.0
PRIMER_MAX_LIBRARY_MISPRIMING = 12.00
PRIMER_PAIR_MAX_LIBRARY_MISPRIMING = 20.00
PRIMER_INSIDE_PENALTY = -1.0
PRIMER_OUTSIDE_PENALTY = 0.0
"""
    
    try:
        with open(settings_file, 'w') as f:
            f.write(default_settings)
        logging.info(f"Created default settings file: {settings_file}")
    except Exception as e:
        logging.error(f"Failed to create default settings file: {e}")

def display_current_settings(settings: Dict) -> None:
    """
    Display current primer design settings to the user.
    """
    logging.info("\nCurrent Primer Design Settings:")
    logging.info("===============================")
    
    # Display settings in categories
    categories = {
        "Primer Size": ["PRIMER_MIN_SIZE", "PRIMER_OPT_SIZE", "PRIMER_MAX_SIZE"],
        "Melting Temperature": ["PRIMER_MIN_TM", "PRIMER_OPT_TM", "PRIMER_MAX_TM", "PRIMER_MAX_DIFF_TM"],
        "GC Content": ["PRIMER_MIN_GC", "PRIMER_OPT_GC_PERCENT", "PRIMER_MAX_GC"],
        "Probe Parameters": ["PRIMER_INTERNAL_MIN_SIZE", "PRIMER_INTERNAL_OPT_SIZE", "PRIMER_INTERNAL_MAX_SIZE",
                             "PRIMER_INTERNAL_MIN_TM", "PRIMER_INTERNAL_OPT_TM", "PRIMER_INTERNAL_MAX_TM",
                             "PRIMER_INTERNAL_MIN_GC", "PRIMER_INTERNAL_OPT_GC_PERCENT", "PRIMER_INTERNAL_MAX_GC"],
        "Structural Constraints": ["PRIMER_MAX_POLY_X", "PRIMER_MAX_SELF_ANY_TH", "PRIMER_MAX_SELF_END_TH", 
                               "PRIMER_MAX_HAIRPIN_TH", "PRIMER_INTERNAL_MAX_POLY_X",
                               "PRIMER_INTERNAL_MAX_SELF_ANY_TH", "PRIMER_INTERNAL_MAX_SELF_END_TH", 
                               "PRIMER_INTERNAL_MAX_HAIRPIN_TH"],
        "Chemical Parameters": ["PRIMER_SALT_MONOVALENT", "PRIMER_SALT_DIVALENT", "PRIMER_DNA_CONC"]
    }
    
    for category, keys in categories.items():
        logging.info(f"\n{category}:")
        for key in keys:
            if key in settings:
                logging.info(f"  {key} = {settings[key]}")
    
    logging.info("\nFor a complete list of all settings, please see PCR_primer_settings.txt")

def get_primer3_settings(params: PrimerParameters) -> Dict:
    """
    Returns a dictionary with Primer3 settings based on the Primer3 web interface parameters
    and the settings file. File settings take precedence over defaults.
    """
    # Load settings from file
    file_settings = load_primer_settings()
    
    # Display current settings to the user
    display_current_settings(file_settings)
    
    # Base settings dictionary with defaults
    settings = {
        # Primer size parameters
        'PRIMER_OPT_SIZE': params.opt_size,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 23,
        
        # Internal oligo (probe) parameters - enable internal oligo selection
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MIN_SIZE': 18,
        'PRIMER_INTERNAL_OPT_SIZE': 20,
        'PRIMER_INTERNAL_MAX_SIZE': 27,
        
        # Number of primers to return - match Primer3 web exactly
        'PRIMER_NUM_RETURN': params.max_primer_sets,
        
        # Tm parameters (melting temperature)
        'PRIMER_MIN_TM': params.min_tm,
        'PRIMER_OPT_TM': params.opt_tm,
        'PRIMER_MAX_TM': params.max_tm,
        'PRIMER_MAX_DIFF_TM': 5.0,  # Max difference in Tm between primers
        
        # Internal oligo Tm parameters
        'PRIMER_INTERNAL_MIN_TM': 57.0,
        'PRIMER_INTERNAL_OPT_TM': 60.0,
        'PRIMER_INTERNAL_MAX_TM': 63.0,
        
        # GC content parameters
        'PRIMER_MIN_GC': params.min_gc,
        'PRIMER_OPT_GC_PERCENT': 50.0,
        'PRIMER_MAX_GC': params.max_gc,
        
        # Internal oligo GC content parameters
        'PRIMER_INTERNAL_MIN_GC': 20.0,
        'PRIMER_INTERNAL_OPT_GC_PERCENT': 50.0,
        'PRIMER_INTERNAL_MAX_GC': 80.0,
        
        # Max poly-X (e.g., AAAAA or CCCCC)
        'PRIMER_MAX_POLY_X': params.max_poly_x,
        'PRIMER_INTERNAL_MAX_POLY_X': 5,
        
        # Salt and DNA concentration parameters
        'PRIMER_SALT_MONOVALENT': params.salt_monovalent,
        'PRIMER_SALT_DIVALENT': params.salt_divalent,
        'PRIMER_DNA_CONC': params.dna_conc,
        'PRIMER_MAX_NS_ACCEPTED': params.max_ns_accepted,
        'PRIMER_SALT_CORRECTIONS': 1,  # SantaLucia 1998
        
        # Internal oligo salt and DNA concentration
        'PRIMER_INTERNAL_SALT_MONOVALENT': params.salt_monovalent,
        'PRIMER_INTERNAL_SALT_DIVALENT': params.salt_divalent,
        'PRIMER_INTERNAL_DNA_CONC': params.dna_conc,
        
        # Product size range - use user-provided values
        'PRIMER_PRODUCT_SIZE_RANGE': [[params.min_size, params.max_size]],
        
        # Thermodynamic alignment parameters (secondary structures)
        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
        'PRIMER_MAX_SELF_ANY_TH': 45.0,
        'PRIMER_MAX_SELF_END_TH': 35.0,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0, 
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
        'PRIMER_MAX_HAIRPIN_TH': 24.0,
        
        # Internal oligo thermodynamic parameters
        'PRIMER_INTERNAL_MAX_SELF_ANY_TH': 47.0,
        'PRIMER_INTERNAL_MAX_SELF_END_TH': 47.0,
        'PRIMER_INTERNAL_MAX_HAIRPIN_TH': 47.0,
        
        # Template mispriming
        'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 40.0,
        'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': 70.0,
        
        # Additional parameters
        'PRIMER_EXPLAIN_FLAG': 1,
        'PRIMER_LIBERAL_BASE': 1,
        'PRIMER_FIRST_BASE_INDEX': 1,
        'PRIMER_MIN_END_QUALITY': 0,
        'PRIMER_MIN_QUALITY': 0,

        # End stability
        'PRIMER_MAX_END_STABILITY': 9.0,
        
        # Library mispriming
        'PRIMER_MAX_LIBRARY_MISPRIMING': 12.00,
        'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': 20.00,
        
        # Target penalty parameters
        'PRIMER_INSIDE_PENALTY': -1.0,
        'PRIMER_OUTSIDE_PENALTY': 0.0
    }
    
    # Override default settings with values from the settings file
    for key, value in file_settings.items():
        if key in settings:
            settings[key] = value
    
    return settings

def design_primers_simplified(sequence: str, params: PrimerParameters) -> Optional[Dict]:
    """A simplified approach to design primers when the standard method fails."""
    try:
        logging.info("Attempting primer design with simplified settings...")
        
        # Basic settings that should work with most primer3 versions
        basic_settings = {
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 23,
            'PRIMER_PRODUCT_SIZE_RANGE': [[params.min_size, params.max_size]],
            'PRIMER_LIBERAL_BASE': 1,
            'PRIMER_NUM_RETURN': params.max_primer_sets
        }
        
        # Sequence args with simplified parameters
        sequence_args = {
            'SEQUENCE_ID': 'consensus',
            'SEQUENCE_TEMPLATE': sequence
        }
        
        # Redirect stdout to capture primer3 output
        import io
        import sys
        original_stdout = sys.stdout
        output_capture = io.StringIO()
        sys.stdout = output_capture
        
        # Try the primer design using the most basic approach
        result = None
        try:
            logging.info("Trying design_primers in primer3 module directly...")
            if hasattr(primer3, 'design_primers'):
                result = primer3.design_primers(sequence_args, basic_settings)
            else:
                logging.info("Using primer3.bindings.designPrimers...")
                result = primer3.bindings.designPrimers(sequence_args, basic_settings)
        except Exception as e:
            logging.error(f"Simplified primer design approach failed: {e}")
            logging.error(traceback.format_exc())
            result = None
            
        # Restore stdout
        sys.stdout = original_stdout
        
        # Check if we got any primers
        if result and any(key.startswith('PRIMER_LEFT_0_') for key in result.keys()):
            return result
        else:
            logging.warning("Simplified approach failed to design primers")
            return None
    except Exception as e:
        logging.error(f"Error in simplified primer design: {e}")
        logging.error(traceback.format_exc())
        return None

def design_primers(sequence: str, settings: Dict, excluded_regions: List[List[int]] = None) -> Optional[Dict]:
    """
    Uses primer3 to design primers based on the provided settings.
    Returns the primer3 result dictionary or None if design fails.
    """
    try:
        # Enhanced logging for debugging
        logging.info(f"Starting primer design on sequence of length {len(sequence)} bp")
        
        # Validate inputs before proceeding
        if not sequence or len(sequence) == 0:
            logging.error("Empty sequence provided to primer design function")
            return None
        
        if excluded_regions is None:
            excluded_regions = []
            
        # Check for valid excluded regions
        valid_excluded_regions = []
        for region in excluded_regions:
            if len(region) == 2 and region[0] >= 0 and region[1] > 0 and region[0] + region[1] <= len(sequence):
                valid_excluded_regions.append(region)
            else:
                logging.warning(f"Skipping invalid excluded region: {region}")
                
        # Use only valid excluded regions
        excluded_regions = valid_excluded_regions
        
        sequence_args = {
            'SEQUENCE_ID': 'consensus',
            'SEQUENCE_TEMPLATE': sequence
        }
        
        # The SEQUENCE_INCLUDED_REGION parameter seems to be causing issues,
        # so let's not include it in our sequence_args
        
        # Ensure sequence length is valid
        if len(sequence) <= 0:
            logging.error("Invalid sequence length")
            return None
        
        if excluded_regions:
            sequence_args['SEQUENCE_EXCLUDED_REGION'] = excluded_regions
            logging.info(f"Using {len(excluded_regions)} excluded regions")
        
        # Log the actual parameters being sent to primer3
        logging.info(f"Primer3 args: SEQUENCE_INCLUDED_REGION=[0, {len(sequence)}], " +
                    f"excluded regions count: {len(excluded_regions)}")
        
        # Redirect stdout to capture primer3 output
        import io
        import sys
        original_stdout = sys.stdout
        output_capture = io.StringIO()
        sys.stdout = output_capture
        
        start_time = time.time()
        
        # Try different primer3 API approaches
        result = None
        try:
            # Try first with the newer design_primers function if available
            logging.info("Trying primer3.design_primers...")
            result = primer3.design_primers(sequence_args, settings)
        except (AttributeError, TypeError, OSError) as e:
            # Fall back to the original bindings if the newer API fails
            logging.info(f"First method failed: {e}")
            logging.info("Falling back to primer3.bindings.designPrimers API")
            try:
                result = primer3.bindings.designPrimers(sequence_args, settings)
            except Exception as inner_e:
                logging.error(f"designPrimers failed: {inner_e}")
                result = None
            
        design_time = time.time() - start_time
        
        # Restore standard output and get any output from primer3
        sys.stdout = original_stdout
        primer3_output = output_capture.getvalue()
        if primer3_output.strip():
            logging.info(f"Primer3 output: {primer3_output[:200]}...")
        
        if result is None:
            logging.error("Primer design failed with all methods")
            # Try simplified approach as a fallback
            return design_primers_simplified(sequence, PrimerParameters(
                min_size=settings['PRIMER_PRODUCT_SIZE_RANGE'][0][0],
                max_size=settings['PRIMER_PRODUCT_SIZE_RANGE'][0][1],
                max_primer_sets=settings['PRIMER_NUM_RETURN']
            ))
            
        # Log the number of results obtained
        left_primers = [k for k in result.keys() if k.startswith('PRIMER_LEFT_') and k.endswith('_SEQUENCE')]
        logging.info(f"Primer design completed in {design_time:.2f} seconds, found {len(left_primers)} primer pairs")
            
        # Check if we got any results
        has_primers = any(key.startswith('PRIMER_LEFT_0_') for key in result.keys())
        if not has_primers:
            logging.info("No suitable primers found with current settings")
            # Try simplified approach as a fallback
            return design_primers_simplified(sequence, PrimerParameters(
                min_size=settings['PRIMER_PRODUCT_SIZE_RANGE'][0][0],
                max_size=settings['PRIMER_PRODUCT_SIZE_RANGE'][0][1],
                max_primer_sets=settings['PRIMER_NUM_RETURN']
            ))
            
        return result
    except Exception as e:
        # Enhanced error logging with full traceback
        logging.error(f"Primer design failed: {str(e)}")
        logging.error(traceback.format_exc())
        # Try simplified approach as a fallback
        return design_primers_simplified(sequence, PrimerParameters(
            min_size=settings['PRIMER_PRODUCT_SIZE_RANGE'][0][0],
            max_size=settings['PRIMER_PRODUCT_SIZE_RANGE'][0][1],
            max_primer_sets=settings['PRIMER_NUM_RETURN']
        ))

def process_primer_set(result: Dict, i: int) -> Optional[Dict]:
    """Process a single primer set from the Primer3 results."""
    try:
        primer_left_seq = result.get(f'PRIMER_LEFT_{i}_SEQUENCE', '')
        primer_right_seq = result.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')
        probe_seq = result.get(f'PRIMER_INTERNAL_{i}_SEQUENCE', '')
        
        # Skip if any primer or probe is missing
        if not all([primer_left_seq, primer_right_seq, probe_seq]):
            return None
            
        # Skip if all are identical (indicates a problem)
        if len({primer_left_seq, primer_right_seq, probe_seq}) != 3:
            return None
            
        # Also store the position information which is important for exact matching
        left_pos = result.get(f'PRIMER_LEFT_{i}', [0, 0])
        right_pos = result.get(f'PRIMER_RIGHT_{i}', [0, 0])
        internal_pos = result.get(f'PRIMER_INTERNAL_{i}', [0, 0])
            
        return {
            'SEQUENCE_ID': 'Consensus_sequence',
            'PRIMER_SET_NUMBER': i + 1,
            'PRIMER_LEFT_SEQUENCE': primer_left_seq,
            'PRIMER_LEFT_TM': result.get(f'PRIMER_LEFT_{i}_TM', ''),
            'PRIMER_LEFT_GC%': result.get(f'PRIMER_LEFT_{i}_GC_PERCENT', ''),
            'PRIMER_LEFT_SELF_ANY_TH': result.get(f'PRIMER_LEFT_{i}_SELF_ANY_TH', ''),
            'PRIMER_LEFT_SELF_END_TH': result.get(f'PRIMER_LEFT_{i}_SELF_END_TH', ''),
            'PRIMER_LEFT_HAIRPIN_TH': result.get(f'PRIMER_LEFT_{i}_HAIRPIN_TH', ''),
            'PRIMER_LEFT_POS': left_pos,
            'PRIMER_RIGHT_SEQUENCE': primer_right_seq,
            'PRIMER_RIGHT_TM': result.get(f'PRIMER_RIGHT_{i}_TM', ''),
            'PRIMER_RIGHT_GC%': result.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', ''),
            'PRIMER_RIGHT_SELF_ANY_TH': result.get(f'PRIMER_RIGHT_{i}_SELF_ANY_TH', ''),
            'PRIMER_RIGHT_SELF_END_TH': result.get(f'PRIMER_RIGHT_{i}_SELF_END_TH', ''),
            'PRIMER_RIGHT_HAIRPIN_TH': result.get(f'PRIMER_RIGHT_{i}_HAIRPIN_TH', ''),
            'PRIMER_RIGHT_POS': right_pos,
            'PROBE_SEQUENCE': probe_seq,
            'PROBE_TM': result.get(f'PRIMER_INTERNAL_{i}_TM', ''),
            'PROBE_GC%': result.get(f'PRIMER_INTERNAL_{i}_GC_PERCENT', ''),
            'PROBE_SELF_ANY_TH': result.get(f'PRIMER_INTERNAL_{i}_SELF_ANY_TH', ''),
            'PROBE_SELF_END_TH': result.get(f'PRIMER_INTERNAL_{i}_SELF_END_TH', ''),
            'PROBE_HAIRPIN_TH': result.get(f'PRIMER_INTERNAL_{i}_HAIRPIN_TH', ''),
            'PROBE_POS': internal_pos,
            'PRODUCT_SIZE': result.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', ''),
            'PAIR_PENALTY': result.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0)  # Store penalty for sorting
        }
    except Exception as e:
        logging.error(f"Error processing primer set {i}: {e}")
        return None

def save_primers_to_csv(output_path: str, primer_data: List[Dict]) -> bool:
    """
    Save the designed primer sets to a CSV file.
    Returns True if successful, False otherwise.
    """
    if not primer_data:
        logging.warning("No primer data to save")
        return False
        
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
    
    # Filter out position and penalty data before saving to CSV
    filtered_data = []
    for item in primer_data:
        filtered_item = {k: v for k, v in item.items() if k in headers}
        filtered_data.append(filtered_item)
    
    try:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()
            writer.writerows(filtered_data)
        
        # Verify the file was actually created
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            logging.info(f"Results saved to: {output_path}")
            logging.info(f"Saved {len(primer_data)} primer sets")
            return True
        else:
            logging.error(f"File was not created or is empty: {output_path}")
            return False
    except Exception as e:
        logging.error(f"Error saving CSV file {output_path}: {e}")
        logging.error(traceback.format_exc())
        return False

def process_consensus_file(file_path: str, params: PrimerParameters, primer3_settings: Dict, interactive_mode: bool = False) -> bool:
    """
    Process a single consensus file to design primers and probes.
    Returns True if successful, False otherwise.
    """
    # Check if file exists
    if not os.path.isfile(file_path):
        logging.error(f"File does not exist: {file_path}")
        return False
    
    # Extract gene name for logging
    gene_name = extract_gene_name_from_file(file_path)
    
    # Read the sequence and header
    sequence_info = read_sequence(file_path)
    if not sequence_info or not sequence_info[0]:
        logging.warning(f"No valid sequence in file: {file_path}")
        return False
        
    sequence, header = sequence_info
    
    # Log hash of the sequence to ensure we're using the exact same input
    import hashlib
    seq_hash = hashlib.md5(sequence.encode()).hexdigest()
    logging.info(f"Sequence hash: {seq_hash}")
    
    # Validate sequence length
    if len(sequence) < params.min_size:
        logging.error(f"Sequence too short ({len(sequence)} bp) for minimum amplicon size ({params.min_size} bp)")
        return False
    
    # Count Ns in sequence to give warning if too many
    n_count = sequence.count('N')
    n_percent = (n_count / len(sequence)) * 100 if sequence else 0
    
    if n_percent > 10:
        logging.warning(f"Sequence contains {n_percent:.1f}% ambiguous bases (N), which may limit primer design.")
    
    # Define the output path
    output_path = os.path.join(
        params.output_dir,
        f"{os.path.splitext(os.path.basename(file_path))[0]}_primers.csv"
    )
    
    primer_data: List[Dict] = []
    
    # Initialize progress bar
    pbar = None
    
    try:
        # Only show tqdm progress bar in terminal mode, not in GUI
        if not interactive_mode:
            pbar = tqdm(total=params.max_primer_sets, desc=f"Designing primers for {gene_name}")
        
        logging.info(f"Designing up to {params.max_primer_sets} primer sets for {gene_name}...")
        
        # Try direct primer design
        result = design_primers(sequence, primer3_settings)
        
        # Process results if we got any
        if result:
            # Extract primer information
            for i in range(params.max_primer_sets):
                if f'PRIMER_LEFT_{i}_SEQUENCE' not in result:
                    break
                    
                primer_set = process_primer_set(result, i)
                if primer_set:
                    primer_data.append(primer_set)
                    
                    if pbar:
                        pbar.update(1)
            
            # Save results only if we found some primers
            if primer_data:
                success = save_primers_to_csv(output_path, primer_data)
                if success:
                    # Double-check that the file exists and has the right content
                    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                        logging.info(f"Successfully designed {len(primer_data)} primer sets for {gene_name}")
                        logging.info(f"Saved {len(primer_data)} primer sets to {output_path}")
                        return True
                    else:
                        logging.error(f"Expected output file {output_path} doesn't exist or is empty")
                        return False
                else:
                    logging.error(f"Failed to save primer sets to {output_path}")
                    return False
            else:
                logging.warning(f"No valid primer sets found for {gene_name}")
                return False
        else:
            logging.warning(f"Primer3 did not return any results for {gene_name}")
            return False
            
    finally:
        if pbar:
            pbar.close()
    
    # If we reach here without returning, it means we didn't succeed
    logging.warning(f"Failed to design primer sets for file: {file_path}")
    if n_percent > 20:
        logging.warning(f"This may be due to the high percentage of ambiguous bases ({n_percent:.1f}%).")
    return False

def extract_gene_name_from_file(file_path: str) -> str:
    """Extract gene name from FASTA header or filename."""
    try:
        if not os.path.isfile(file_path):
            return os.path.basename(file_path)
            
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
            if first_line.startswith('>'):
                # Extract gene name from header
                header = first_line[1:]
                # Try to extract gene name
                if '*' in header:
                    # Extract text between asterisks if present
                    gene_parts = header.split('*')
                    if len(gene_parts) >= 2:
                        return gene_parts[1].strip()
                # Try other parsing methods
                header_parts = header.split('_')
                for part in header_parts:
                    if 'gene' in part.lower() or 'rna' in part.lower() or 'tubulin' in part.lower():
                        return part
                
                # Use part after organism name as fallback
                if '_' in header and header.count('_') >= 2:
                    return header.split('_')[1]
    except Exception as e:
        logging.error(f"Error extracting gene name: {e}")
    
    # Fallback to filename parsing
    filename = os.path.basename(file_path)
    if '_' in filename:
        parts = filename.split('_')
        # Try to find 'gene' or 'consensus' in parts
        for i, part in enumerate(parts):
            if part.lower() == 'consensus' and i+1 < len(parts):
                return parts[i+1]
    
    return os.path.splitext(os.path.basename(file_path))[0]

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
        except Exception as e:
            logging.debug(f"Failed to remove response file: {e}")
    
    logging.info(f"Waiting for input: {prompt}")
    
    # Write the prompt to the request file
    try:
        with open(input_request_file, 'w', encoding='utf-8') as f:
            f.write(prompt)
    except Exception as e:
        logging.error(f"Failed to create input request file: {e}")
        # Fall back to default value if we can't create the request file
        return default if default is not None else ""
    
    # Wait for response file to appear
    timeout = 600
    start_time = time.time()
    while not os.path.exists(input_response_file):
        time.sleep(0.5)
        # Check if the input request file still exists
        if not os.path.exists(input_request_file):
            logging.warning("Input request file was removed unexpectedly")
            return default if default is not None else ""
        if time.time() - start_time > timeout:
            logging.warning("Input timeout. Using default value.")
            try:
                os.remove(input_request_file)  # Clean up request file on timeout
            except Exception:
                pass
            return default if default is not None else ""
    
    # Wait a moment to ensure file is completely written
    time.sleep(0.1)
    
    # Read response
    try:
        with open(input_response_file, 'r', encoding='utf-8') as f:
            user_input = f.read().strip()
        
        # Delete the files after reading
        cleanup_files = [input_response_file]
        if os.path.exists(input_request_file):
            cleanup_files.append(input_request_file)
            
        for file_path in cleanup_files:
            try:
                os.remove(file_path)
            except Exception as e:
                logging.debug(f"Failed to remove file {file_path}: {e}")
            
        return user_input
    except Exception as e:
        logging.error(f"Error reading input response: {e}")
        # Try to clean up the request file if it still exists
        if os.path.exists(input_request_file):
            try:
                os.remove(input_request_file)
            except Exception:
                pass
        return default if default is not None else ""

def cleanup_previous_run(output_dir: str) -> None:
    """
    Clean up files and directories from previous run of Module 3.
    Removes all files in the output directory and recreates it.
    """
    if os.path.exists(output_dir):
        try:
            logging.info(f"Cleaning up previous run files in {output_dir}...")
            # Remove all files in the directory
            for filename in os.listdir(output_dir):
                file_path = os.path.join(output_dir, filename)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                        logging.info(f"Deleted file: {filename}")
                    elif os.path.isdir(file_path):
                        import shutil
                        shutil.rmtree(file_path)
                        logging.info(f"Deleted directory: {filename}")
                except Exception as e:
                    logging.error(f"Failed to delete {file_path}: {e}")
            
            # Remove the directory itself
            os.rmdir(output_dir)
            logging.info(f"Removed directory: {output_dir}")
        except Exception as e:
            logging.error(f"Error during cleanup: {e}")

def main():
    try:
        setup_logging(__file__)
        
        # Check command line arguments
        interactive_mode = "--interactive" in sys.argv
        
        # Clean up previous run
        output_dir = "3_"
        cleanup_previous_run(output_dir)
        
        # Create output directory
        try:
            os.makedirs(output_dir, exist_ok=True)
            logging.info(f"Using output directory: {output_dir}")
        except Exception as e:
            logging.error(f"Failed to create output directory {output_dir}: {e}")
            return

        # Create or load the settings file for primer design
        settings_file = "PCR_primer_settings.txt"
        if not os.path.exists(settings_file):
            create_default_settings_file(settings_file)
            logging.info(f"Created default primer settings file: {settings_file}")
        
        # Prompt user about primer settings
        if interactive_mode:
            logging.info(f"\nPrimer design settings are stored in {settings_file}")
            logging.info("You can modify the settings by editing this file.")
            modify_settings = get_console_input("Do you want to modify primer design settings? (yes/no): ").lower()
            
            if modify_settings in ["yes", "y"]:
                # Show current settings before user modifies them
                settings = load_primer_settings(settings_file)
                display_current_settings(settings)
                
                logging.info(f"\nPlease edit {settings_file} with your desired parameter values.")
                logging.info("Save the file when you're done, then return to this program.")
                
                # Wait for user to confirm they've edited the file
                proceed = get_console_input("Have you finished editing the settings file? (yes/no): ").lower()
                
                if proceed not in ["yes", "y"]:
                    logging.info("Primer design canceled by user.")
                    return
                
                # Reload settings to pick up any changes
                settings = load_primer_settings(settings_file)
                logging.info("Updated settings loaded successfully.")
        
        # Find consensus files
        try:
            consensus_files = find_consensus_files("2_")
            if not consensus_files:
                logging.error("No consensus files found in the input directory (2_).")
                return
                
            logging.info(f"Found {len(consensus_files)} consensus files to process")
        except FileNotFoundError as e:
            logging.error(f"Error finding consensus files: {e}")
            return
        
        # Process each consensus file separately
        successful_files = 0
        for file_idx, file_path in enumerate(consensus_files):
            try:
                logging.info(f"\n{'='*60}")
                logging.info(f"Processing file {file_idx+1} of {len(consensus_files)}: {file_path}")
                
                # Get parameters for this specific file
                file_params = get_primer_parameters(interactive_mode, file_path)
                if not file_params:
                    logging.warning(f"Skipping file {file_path} due to invalid parameters")
                    continue
                
                # Configure primer3 settings based on the settings file
                primer3_settings = get_primer3_settings(file_params)
                
                # Design primers for this consensus
                result = process_consensus_file(file_path, file_params, primer3_settings, interactive_mode)
                
                # Verify the output file exists after processing
                output_path = os.path.join(
                    file_params.output_dir,
                    f"{os.path.splitext(os.path.basename(file_path))[0]}_primers.csv"
                )
                
                if result and os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                    successful_files += 1
                else:
                    logging.error(f"Failed to create valid output file for {file_path}")
                
            except Exception as e:
                logging.error(f"Error processing file {file_path}: {e}")
                # Print a brief traceback but not the full one to avoid flooding GUI
                tb = traceback.format_exc().split('\n')
                if len(tb) > 5:
                    tb = tb[:3] + ['...'] + tb[-2:]
                logging.error('\n'.join(tb))
        
        # Summary of results
        if successful_files > 0:
            logging.info(f"\nPrimer design successfully completed for {successful_files} out of {len(consensus_files)} files.")
        else:
            logging.error("\nPrimer design failed for all files. Please check the logs for details.")
        
    except Exception as e:
        logging.error(f"Unexpected error in main function: {e}")
        logging.error(traceback.format_exc())

if __name__ == '__main__':
    main()
