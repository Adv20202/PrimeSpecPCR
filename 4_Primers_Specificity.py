#!/usr/bin/env python3
"""
PrimeSpecPCR - Primer Specificity Test Module
Tests the specificity of designed primer sets by comparing them against GenBank sequences.
"""

import os
import glob
import csv
import re
import time
import logging
import webbrowser
import subprocess
import sys
import json
import hashlib
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

# Import required libraries with automatic installation if missing
try:
    from Bio import Entrez, SeqIO
    from Bio.Align import PairwiseAligner
except ImportError:
    # Suppress stdout during installation
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython>=1.79"], 
                              stdout=devnull, stderr=devnull)
    from Bio import Entrez, SeqIO
    from Bio.Align import PairwiseAligner

try:
    import concurrent.futures
except ImportError:
    # Suppress stdout during installation
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "futures"], 
                              stdout=devnull, stderr=devnull)
    import concurrent.futures

try:
    import requests
except ImportError:
    # Suppress stdout during installation
    with open(os.devnull, 'w') as devnull:
        try:
            # First try regular pip install
            subprocess.check_call([sys.executable, "-m", "pip", "install", "requests"], 
                                  stdout=devnull, stderr=devnull)
        except subprocess.CalledProcessError:
            # If that fails, try with --user flag
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", "requests"], 
                                      stdout=devnull, stderr=devnull)
            except subprocess.CalledProcessError:
                # If that also fails, try with --break-system-packages
                try:
                    subprocess.check_call([sys.executable, "-m", "pip", "install", "requests", "--break-system-packages"], 
                                          stdout=devnull, stderr=devnull)
                except subprocess.CalledProcessError:
                    print("WARNING: Could not install requests library automatically.")
                    print("Please run: pip install --user requests")
                    print("Or create a virtual environment.")
                    sys.exit(1)
    import requests

import traceback
from urllib.error import HTTPError, URLError
from http.client import IncompleteRead
import socket
import copy
import xml.etree.ElementTree as ET
from collections import defaultdict

# Configure logging
class CustomFormatter(logging.Formatter):
    def format(self, record):
        return record.getMessage()

def setup_logging(script_name: str) -> None:
    log_name = f"{os.path.splitext(script_name)[0]}.log"
    
    formatter = CustomFormatter()
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    
    file_handler = logging.FileHandler(log_name, encoding='utf-8')
    file_handler.setFormatter(formatter)
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Remove any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Add our custom handlers
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

# Data structures
@dataclass
class PrimerPair:
    id: str
    left_primer: str
    right_primer: str
    probe: str
    product_size: int
    left_tm: float
    right_tm: float
    probe_tm: float
    left_gc: float = 0.0
    right_gc: float = 0.0
    probe_gc: float = 0.0
    amplicon_sequence: str = ""

@dataclass
class SpecificitySettings:
    max_hits: int = 1000
    min_identity: float = 85.0
    max_e_value: float = 100
    word_size: int = 11
    fetch_limit: int = 1000
    max_total_mismatches: int = 3
    max_3prime_mismatches: int = 2
    bases_from_3prime: int = 6
    max_mismatches_probe: int = 3
    cache_dir: str = "4_cache"
    output_dir: str = "4_"
    report_filename: str = "specificity_report.html"
    max_threads: int = 4  # Maximum number of parallel threads to use
    
    def validate(self) -> bool:
        """Validate settings are within reasonable bounds."""
        if self.max_hits < 1 or self.max_hits > 50000:
            logging.error(f"Invalid max_hits value: {self.max_hits}. Must be between 1 and 50000.")
            return False
        if self.min_identity < 50 or self.min_identity > 100:
            logging.error(f"Invalid min_identity value: {self.min_identity}. Must be between 50 and 100.")
            return False
        if self.max_e_value < 0:
            logging.error(f"Invalid max_e_value value: {self.max_e_value}. Must be positive.")
            return False
        if self.word_size < 4 or self.word_size > 20:
            logging.error(f"Invalid word_size value: {self.word_size}. Must be between 4 and 20.")
            return False
        if self.fetch_limit < 1 or self.fetch_limit > 30000:
            logging.error(f"Invalid fetch_limit value: {self.fetch_limit}. Must be between 1 and 30000.")
            return False
        if self.max_total_mismatches < 0:
            logging.error(f"Invalid max_total_mismatches value: {self.max_total_mismatches}. Must be non-negative.")
            return False
        if self.max_3prime_mismatches < 0:
            logging.error(f"Invalid max_3prime_mismatches value: {self.max_3prime_mismatches}. Must be non-negative.")
            return False
        if self.bases_from_3prime < 1:
            logging.error(f"Invalid bases_from_3prime value: {self.bases_from_3prime}. Must be at least 1.")
            return False
        if self.max_mismatches_probe < 0:
            logging.error(f"Invalid max_mismatches_probe value: {self.max_mismatches_probe}. Must be non-negative.")
            return False
        if self.max_threads < 1:
            logging.error(f"Invalid max_threads value: {self.max_threads}. Must be at least 1.")
            return False
        return True

@dataclass
class SpecificityMatch:
    accession: str
    organism: str
    description: str
    left_match: bool = False
    right_match: bool = False
    probe_match: bool = False
    amplicon_size: Optional[int] = None
    left_position: Optional[int] = None
    right_position: Optional[int] = None
    probe_position: Optional[int] = None
    left_alignment: str = ""
    right_alignment: str = ""
    probe_alignment: str = ""
    left_mismatches: int = 0
    right_mismatches: int = 0
    probe_mismatches: int = 0
    left_3prime_mismatches: int = 0
    right_3prime_mismatches: int = 0
    taxid: str = ""
    
    @property
    def all_match(self) -> bool:
        return self.left_match and self.right_match and self.probe_match
    
    @property
    def match_score(self) -> int:
        """Return a score based on matching components (0-3)."""
        return sum([self.left_match, self.right_match, self.probe_match])

def get_user_data() -> Tuple[str, str]:
    """Get email and API key from user data directory."""
    user_data_dir = os.path.join(os.getcwd(), 'user_data')
    email = ""
    api_key = ""
    
    if os.path.exists(user_data_dir):
        email_file = os.path.join(user_data_dir, 'email.txt')
        api_key_file = os.path.join(user_data_dir, 'api_key.txt')
        
        if os.path.exists(email_file):
            with open(email_file, 'r', encoding='utf-8') as f:
                email = f.read().strip()
                
        if os.path.exists(api_key_file):
            with open(api_key_file, 'r', encoding='utf-8') as f:
                api_key = f.read().strip()
    
    return email, api_key

def get_taxid_for_accession(accession: str, email: str, api_key: str) -> str:
    """Get taxonomy ID for a GenBank accession number using NCBI E-utilities."""
    if not accession:
        return ""
        
    # Configure Entrez
    Entrez.email = email
    Entrez.api_key = api_key
    Entrez.tool = "PrimeSpecPCR"
    
    try:
        # Use elink to get taxonomy ID from nucleotide accession
        handle = Entrez.elink(dbfrom="nuccore", db="taxonomy", id=accession)
        record = Entrez.read(handle)
        handle.close()
        
        # Extract taxid from the response
        if record and len(record[0]["LinkSetDb"]) > 0:
            links = record[0]["LinkSetDb"][0]["Link"]
            if links:
                return links[0]["Id"]  # Return the first taxid
    except Exception as e:
        logging.debug(f"Error getting taxid for {accession}: {e}")
    
    return ""

def get_taxonomy_lineage(taxid: str, email: str, api_key: str) -> Dict[str, str]:
    """Get full taxonomy lineage for a given taxonomy ID."""
    if not taxid:
        return {}
        
    # Configure Entrez
    Entrez.email = email
    Entrez.api_key = api_key
    Entrez.tool = "PrimeSpecPCR"
    
    try:
        # Fetch taxonomy information
        handle = Entrez.efetch(db="taxonomy", id=taxid)
        records = Entrez.read(handle)
        handle.close()
        
        if records and len(records) > 0:
            record = records[0]
            lineage_info = {
                "taxid": taxid,
                "scientific_name": record.get("ScientificName", ""),
                "lineage": record.get("Lineage", ""),
                "rank": record.get("Rank", "")
            }
            
            # Add lineage taxids
            lineage_ex = record.get("LineageEx", [])
            for entry in lineage_ex:
                rank = entry.get("Rank", "").lower()
                if rank and "taxid" in entry and "ScientificName" in entry:
                    lineage_info[f"{rank}_taxid"] = entry["TaxId"]
                    lineage_info[f"{rank}_name"] = entry["ScientificName"]
            
            return lineage_info
    except Exception as e:
        logging.debug(f"Error getting taxonomy lineage for {taxid}: {e}")
    
    return {}

def get_user_target_taxids() -> List[str]:
    """Ask user for target taxonomy IDs to evaluate specificity against."""
    logging.info("\nSpecificity Evaluation Against Target Taxonomy")
    logging.info("============================================")
    
    target_taxids = []
    taxid_input = get_console_input(
        "Enter NCBI Taxonomy ID(s) to evaluate primer specificity against (comma-separated): ",
        ""
    )
    
    if taxid_input:
        # Parse and validate taxids
        for taxid in taxid_input.split(","):
            taxid = taxid.strip()
            if taxid.isdigit():
                target_taxids.append(taxid)
            else:
                logging.warning(f"Invalid taxonomy ID format: {taxid} - skipping")
    
    if target_taxids:
        logging.info(f"Using target taxonomy IDs: {', '.join(target_taxids)}")
    else:
        logging.info("No target taxonomy IDs provided. Proceeding without taxonomy specificity evaluation.")
    
    return target_taxids

def get_console_input(prompt: str, default: Optional[str] = None) -> str:
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
    
    logging.info(f"Waiting for input: {prompt}")
    
    # Write the prompt to the request file
    with open(input_request_file, 'w', encoding='utf-8') as f:
        f.write(prompt)
    
    # Wait for response file to appear
    timeout = 600  # 10 minutes timeout
    start_time = time.time()
    while not os.path.exists(input_response_file):
        time.sleep(0.1)
        if time.time() - start_time > timeout:
            logging.warning("Input timeout. Using default value.")
            return default if default is not None else ""
    
    # Wait a moment to ensure file is completely written
    time.sleep(0.1)
    
    # Read response
    try:
        with open(input_response_file, 'r', encoding='utf-8') as f:
            user_input = f.read().strip()
        
        # Delete the files after reading
        try:
            os.remove(input_response_file)
            os.remove(input_request_file)
        except:
            pass
            
        return user_input
    except Exception as e:
        logging.error(f"Error reading input: {e}")
        return default if default is not None else ""

def load_specificity_settings(interactive_mode: bool = False) -> SpecificitySettings:
    """Load specificity test settings from a file or create default settings."""
    settings_file = "primer_specificity_settings.txt"
    settings = SpecificitySettings()
    
    # Check if settings file exists
    if os.path.exists(settings_file):
        try:
            with open(settings_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    if '=' in line:
                        key, value = [part.strip() for part in line.split('=', 1)]
                        if key == 'MAX_HITS':
                            settings.max_hits = int(value)
                        elif key == 'MIN_IDENTITY':
                            settings.min_identity = float(value)
                        elif key == 'MAX_E_VALUE':
                            settings.max_e_value = float(value)
                        elif key == 'WORD_SIZE':
                            settings.word_size = int(value)
                        elif key == 'FETCH_LIMIT':
                            settings.fetch_limit = int(value)
                        elif key == 'MAX_TOTAL_MISMATCHES':
                            settings.max_total_mismatches = int(value)
                        elif key == 'MAX_3PRIME_MISMATCHES':
                            settings.max_3prime_mismatches = int(value)
                        elif key == 'BASES_FROM_3PRIME':
                            settings.bases_from_3prime = int(value)
                        elif key == 'MAX_MISMATCHES_PROBE':
                            settings.max_mismatches_probe = int(value)
                        elif key == 'CACHE_DIR':
                            settings.cache_dir = value
                        elif key == 'OUTPUT_DIR':
                            settings.output_dir = value
            
            logging.info(f"Loaded settings from {settings_file}")
        except Exception as e:
            logging.error(f"Error loading settings from {settings_file}: {e}")
            # Create default settings file
            create_default_settings_file(settings_file)
    else:
        logging.info(f"Settings file {settings_file} not found, creating with default values")
        create_default_settings_file(settings_file)
    
    # In interactive mode, ask if user wants to modify settings
    if interactive_mode:
        display_current_settings(settings)
        modify = get_console_input("Do you want to modify specificity test settings? (yes/no): ", "no").lower()
        
        if modify in ['yes', 'y']:
            logging.info(f"\nPlease edit {settings_file} with your desired parameter values.")
            logging.info("Save the file when you're done, then return to this program.")
            
            proceed = get_console_input("Have you finished editing the settings file? (yes/no): ", "yes").lower()
            
            if proceed in ["yes", "y"]:
                # Reload settings
                new_settings = load_specificity_settings(False)  # Avoid recursion
                settings = new_settings
                logging.info("Updated settings loaded successfully.")
                display_current_settings(settings)
    
    # Create directories
    os.makedirs(settings.cache_dir, exist_ok=True)
    os.makedirs(settings.output_dir, exist_ok=True)
    
    return settings

def create_default_settings_file(settings_file: str) -> None:
    """Create a default settings file with description and default values."""
    default_settings = """# Primer Specificity Test Settings
# Modify values as needed, but maintain the format: SETTING_NAME = value
# Lines starting with # are comments and will be ignored by the program

# =====================================================================
# Search parameters
# =====================================================================

# MAX_HITS: Maximum number of BLAST hits to retrieve per primer
MAX_HITS = 1000

# MIN_IDENTITY: Minimum percent identity for a valid match (%)
MIN_IDENTITY = 85.0

# MAX_E_VALUE: Maximum E-value for BLAST searches
MAX_E_VALUE = 100

# WORD_SIZE: Word size for BLAST searches
WORD_SIZE = 11

# FETCH_LIMIT: Maximum number of sequences to fetch from GenBank
FETCH_LIMIT = 1000

# =====================================================================
# Primer specificity stringency parameters
# =====================================================================

# MAX_TOTAL_MISMATCHES: Maximum number of total mismatches for a primer to be specific
MAX_TOTAL_MISMATCHES = 3

# MAX_3PRIME_MISMATCHES: Maximum number of mismatches within the last X base pairs at the 3' end
MAX_3PRIME_MISMATCHES = 2

# BASES_FROM_3PRIME: Number of base pairs from the 3' end to check for mismatches
BASES_FROM_3PRIME = 6

# MAX_MISMATCHES_PROBE: Maximum number of mismatches within the probe sequence
MAX_MISMATCHES_PROBE = 3

# =====================================================================
# Directories
# =====================================================================

# CACHE_DIR: Directory to store cached GenBank sequences
CACHE_DIR = 4_cache

# OUTPUT_DIR: Directory to store output files
OUTPUT_DIR = 4_
"""
    
    try:
        with open(settings_file, 'w', encoding='utf-8') as f:
            f.write(default_settings)
        logging.info(f"Created default settings file: {settings_file}")
    except Exception as e:
        logging.error(f"Failed to create default settings file: {e}")

def display_current_settings(settings: SpecificitySettings) -> None:
    """Display current specificity test settings to the user."""
    logging.info("\nCurrent Primer Specificity Test Settings:")
    logging.info("=======================================")
    logging.info(f"Maximum BLAST hits: {settings.max_hits}")
    logging.info(f"Minimum identity (%): {settings.min_identity}")
    logging.info(f"Maximum E-value: {settings.max_e_value}")
    logging.info(f"Word size: {settings.word_size}")
    logging.info(f"Fetch limit: {settings.fetch_limit}")
    logging.info("\nPrimer Specificity Stringency:")
    logging.info(f"Maximum total mismatches for specificity: {settings.max_total_mismatches}")
    logging.info(f"Maximum 3' end mismatches: {settings.max_3prime_mismatches}")
    logging.info(f"3' end region size (bp): {settings.bases_from_3prime}")
    logging.info(f"Maximum probe mismatches: {settings.max_mismatches_probe}")
    logging.info("\nDirectories:")
    logging.info(f"Cache directory: {settings.cache_dir}")
    logging.info(f"Output directory: {settings.output_dir}")
    logging.info("=======================================")

def find_primer_files(directory: str = "3_") -> List[str]:
    """Find CSV files containing primer data in the specified directory."""
    if not os.path.exists(directory):
        logging.error(f"Directory {directory} does not exist")
        return []
    
    files = glob.glob(os.path.join(directory, "*primers.csv"))
    if not files:
        logging.warning(f"No primer files found in {directory}")
    
    return sorted(files)

def read_primer_pairs(csv_file: str) -> List[PrimerPair]:
    """Read primer pairs from a CSV file."""
    primer_pairs = []
    
    try:
        with open(csv_file, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            
            # Validate that the CSV has the required headers
            required_headers = [
                'PRIMER_SET_NUMBER', 
                'PRIMER_LEFT_SEQUENCE', 
                'PRIMER_RIGHT_SEQUENCE',
                'PROBE_SEQUENCE',
                'PRODUCT_SIZE'
            ]
            
            headers = reader.fieldnames if reader.fieldnames else []
            missing_headers = [h for h in required_headers if h not in headers]
            
            if missing_headers:
                logging.error(f"CSV file {csv_file} is missing required headers: {', '.join(missing_headers)}")
                return []
            
            # Read each row
            for row in reader:
                try:
                    primer_id = str(row.get('PRIMER_SET_NUMBER', ''))
                    
                    # Basic validation
                    left_primer = row.get('PRIMER_LEFT_SEQUENCE', '').upper()
                    right_primer = row.get('PRIMER_RIGHT_SEQUENCE', '').upper()
                    probe = row.get('PROBE_SEQUENCE', '').upper()
                    
                    # Skip if any primer is missing
                    if not (left_primer and right_primer and probe):
                        logging.warning(f"Skipping incomplete primer set {primer_id} in {csv_file}")
                        continue
                    
                    # Get additional data if available
                    try:
                        product_size = int(row.get('PRODUCT_SIZE', '0'))
                    except ValueError:
                        product_size = 0
                        
                    try:
                        left_tm = float(row.get('PRIMER_LEFT_TM', '0'))
                    except ValueError:
                        left_tm = 0.0
                        
                    try:
                        right_tm = float(row.get('PRIMER_RIGHT_TM', '0'))
                    except ValueError:
                        right_tm = 0.0
                        
                    try:
                        probe_tm = float(row.get('PROBE_TM', '0'))
                    except ValueError:
                        probe_tm = 0.0
                        
                    try:
                        left_gc = float(row.get('PRIMER_LEFT_GC%', '0'))
                    except ValueError:
                        left_gc = 0.0
                        
                    try:
                        right_gc = float(row.get('PRIMER_RIGHT_GC%', '0'))
                    except ValueError:
                        right_gc = 0.0
                        
                    try:
                        probe_gc = float(row.get('PROBE_GC%', '0'))
                    except ValueError:
                        probe_gc = 0.0
                    
                    primer_pair = PrimerPair(
                        id=primer_id,
                        left_primer=left_primer,
                        right_primer=right_primer,  # Keep the original reverse primer sequence
                        probe=probe,
                        product_size=product_size,
                        left_tm=left_tm,
                        right_tm=right_tm,
                        probe_tm=probe_tm,
                        left_gc=left_gc,
                        right_gc=right_gc,
                        probe_gc=probe_gc
                    )
                    
                    primer_pairs.append(primer_pair)
                except Exception as e:
                    logging.error(f"Error processing primer set {row.get('PRIMER_SET_NUMBER', 'unknown')}: {e}")
                    continue
    except Exception as e:
        logging.error(f"Error reading CSV file {csv_file}: {e}")
        return []
    
    logging.info(f"Read {len(primer_pairs)} primer pairs from {csv_file}")
    return primer_pairs

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                       'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                       'N': 'N', 'n': 'n', '-': '-', '.': '.'}
    
    return ''.join(complement_dict.get(base, 'N') for base in reversed(sequence))

def find_cache_file(primer: str, cache_dir: str) -> Optional[str]:
    """Find a cached file for the given primer sequence."""
    # Create a hash of the primer sequence to use as filename
    primer_hash = hashlib.md5(primer.encode()).hexdigest()
    cache_file = os.path.join(cache_dir, f"{primer_hash}.fasta")
    
    # Only use cache if file exists and has content (non-empty)
    if os.path.exists(cache_file) and os.path.getsize(cache_file) > 10:  # Require at least some meaningful content
        # Check if the file is too old (more than 7 days)
        file_age = time.time() - os.path.getmtime(cache_file)
        if file_age > (7 * 24 * 60 * 60):  # 7 days in seconds
            logging.info(f"Cache file is older than 7 days, will perform new search")
            return None
        return cache_file
    
    return None

def direct_blast_search(primer: str, settings: SpecificitySettings, email: str, api_key: str) -> List[Dict]:
    """
    Perform a BLAST search directly using the NCBI BLAST REST API.
    Returns a list of hit dictionaries with hit positions.
    """
    logging.info(f"Starting BLAST search for primer: {primer}")
    
    # Set up parameters for BLASTN search
    blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    # Put request - initiates the BLAST search
    put_params = {
        'CMD': 'Put',
        'PROGRAM': 'blastn',
        'DATABASE': 'nt',
        'QUERY': primer,
        'HITLIST_SIZE': settings.max_hits,
        'EXPECT': settings.max_e_value,      # Much more permissive E-value to find more distant matches
        'WORD_SIZE': settings.word_size,       # Reduced word size to increase sensitivity
        'FILTER': 'F',
        'FORMAT_TYPE': 'XML',
        'EMAIL': email,
        'TOOL': 'PrimeSpecPCR',
        'MEGABLAST': 'off',
        'THRESHOLD': 100,     # More permissive threshold
        'GAPCOSTS': '5 2',    # Adjust gap costs to allow more alignment flexibility
        'MATCH_SCORES': '1,-1' # More permissive match/mismatch scoring
    }
    
    if api_key:
        put_params['API_KEY'] = api_key
    
    try:
        # Submit the search
        response = requests.post(blast_url, data=put_params)
        response.raise_for_status()
        
        logging.info(f"BLAST search response status: {response.status_code}")
        logging.info(f"Response text excerpt: {response.text[:200]}...")
        
        # Extract RID (Request ID) from the response
        match = re.search(r'(?:RID|QBlastInfoBegin)\s*=\s*([A-Za-z0-9-]+)', response.text)
        if not match:
            # Try alternative pattern used in newer BLAST API responses
            match = re.search(r'<p>Request ID: ([A-Za-z0-9-]+)</p>', response.text)
            
        if not match:
            logging.error("Failed to extract RID from BLAST response")
            logging.debug(f"Full response: {response.text}")
            return []
            
        rid = match.group(1).strip()
        logging.info(f"Extracted RID: {rid}")
        
        # Poll for results
        max_attempts = 30
        attempts = 0
        while attempts < max_attempts:
            attempts += 1
            
            # Wait between polling attempts
            time.sleep(10)
            
            # Get request - checks the status and retrieves results when ready
            get_params = {
                'CMD': 'Get',
                'FORMAT_TYPE': 'XML',
                'RID': rid
            }
            
            logging.info(f"Polling BLAST results, attempt {attempts}/{max_attempts}...")
            response = requests.get(blast_url, params=get_params)
            response.raise_for_status()
            
            # Check if search is still running
            if "Status=WAITING" in response.text:
                logging.info("BLAST search still in progress, waiting...")
                continue
                
            # Check if search is done
            if "Status=READY" in response.text or "<BlastOutput>" in response.text:
                # Parse XML results using ElementTree
                hits = []
                
                try:
                    # Check for valid XML with BlastOutput
                    if "<BlastOutput>" not in response.text:
                        if "Status=WAITING" in response.text or "Status=UNKNOWN" in response.text:
                            logging.warning("BLAST search still in progress or unknown status")
                            continue
                        
                        logging.warning("Unexpected response format from BLAST")
                        logging.debug(f"Response content: {response.text[:500]}...")
                        return []
                    
                    # Parse the XML response
                    try:
                        # Clean response text if needed
                        xml_text = response.text
                        if "&" in xml_text and not "&amp;" in xml_text:
                            xml_text = xml_text.replace('&', '&amp;')
                        
                        root = ET.fromstring(xml_text)
                    except ET.ParseError as e:
                        logging.error(f"XML parsing error: {e}")
                        # Try to fix common XML issues
                        xml_text = re.sub(r'[^\x09\x0A\x0D\x20-\uD7FF\uE000-\uFFFD\U00010000-\U0010FFFF]', '', xml_text)
                        try:
                            root = ET.fromstring(xml_text)
                        except ET.ParseError as e2:
                            logging.error(f"Failed to parse XML even after cleanup: {e2}")
                            return []
                    
                    # Process each hit from the XML
                    hit_elements = root.findall(".//Hit")
                    if not hit_elements:
                        logging.info("No hits found in BLAST results")
                        return []
                        
                    for hit in hit_elements:
                        try:
                            hit_id_elem = hit.find(".//Hit_accession")
                            hit_def_elem = hit.find(".//Hit_def")
                            
                            if hit_id_elem is None or hit_def_elem is None:
                                continue
                                
                            hit_id = hit_id_elem.text
                            hit_def = hit_def_elem.text
                            organism = "Unknown"
                            
                            # Extract organism name
                            org_match = re.search(r'\[(.*?)\]', hit_def)
                            if org_match:
                                organism = org_match.group(1)
                            
                            # Process the HSPs (High-scoring Segment Pairs)
                            hsps = hit.findall(".//Hsp")
                            if not hsps:
                                continue
                                
                            for hsp in hsps:
                                identity_elem = hsp.find(".//Hsp_identity")
                                align_len_elem = hsp.find(".//Hsp_align-len")
                                hit_from_elem = hsp.find(".//Hsp_hit-from")
                                
                                if identity_elem is None or align_len_elem is None or hit_from_elem is None:
                                    continue
                                    
                                try:
                                    identity = int(identity_elem.text)
                                    align_len = int(align_len_elem.text)
                                    hit_position = int(hit_from_elem.text)
                                    identity_percent = (identity / align_len) * 100
                                    
                                    if identity_percent >= settings.min_identity:
                                        hits.append({
                                            'accession': hit_id,
                                            'description': hit_def,
                                            'organism': organism,
                                            'identity_percent': identity_percent,
                                            'hit_position': hit_position
                                        })
                                        break  # Only add each hit once
                                except (ValueError, TypeError) as e:
                                    logging.debug(f"Error processing HSP data: {e}")
                                    continue
                        except Exception as e:
                            logging.debug(f"Error processing hit element: {e}")
                            continue
                        
                        # Limit the number of hits
                        if len(hits) >= settings.fetch_limit:
                            break
                            
                    return hits
                    
                except Exception as e:
                    logging.error(f"Error parsing BLAST XML: {e}")
                    return []
            
            # If we've reached maximum attempts without getting results
            logging.warning(f"BLAST search timed out after {max_attempts} attempts")
            
            # As a last resort, try using the Entrez.esearch directly to find sequences
            try:
                logging.info("Attempting direct Entrez search as fallback...")
                # Create search term based on primer sequence with wildcards to allow mismatches
                search_term = primer
                # Add organism qualifier if it would help (example with "bacteria" as fallback)
                handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=settings.fetch_limit)
                record = Entrez.read(handle)
                handle.close()
                
                if record["Count"] != "0":
                    ids = record["IdList"]
                    logging.info(f"Found {len(ids)} sequences via direct Entrez search")
                    
                    # Format results like BLAST hits
                    fallback_hits = []
                    for sequence_id in ids[:settings.fetch_limit]:
                        fallback_hits.append({
                            'accession': sequence_id,
                            'description': "Found via Entrez direct search",
                            'organism': "Unknown",
                            'identity_percent': 85.0,  # Default value
                            'hit_position': 1  # Default position
                        })
                    return fallback_hits
            except Exception as e:
                logging.error(f"Fallback Entrez search failed: {e}")
            
            return []
        
    except Exception as e:
        logging.error(f"Error performing BLAST search: {e}")
        return []

def fetch_matching_sequences(primer: str, settings: SpecificitySettings, email: str, api_key: str) -> Optional[str]:
    """
    Fetch sequences from GenBank that match the given primer.
    Returns the path to the fasta file containing the sequences.
    Limits each sequence to 1000 nucleotides starting from the primer match position.
    """
    # Check cache first
    cache_file = find_cache_file(primer, settings.cache_dir)
    if cache_file:
        return cache_file
    
    # Set up Entrez email and API key
    Entrez.email = email
    Entrez.api_key = api_key
    Entrez.tool = "PrimeSpecPCR"
    
    # Create a hash of the primer sequence to use as filename
    primer_hash = hashlib.md5(primer.encode()).hexdigest()
    cache_file = os.path.join(settings.cache_dir, f"{primer_hash}.fasta")
    
    try:
        # Temporarily reduce logging level for BLAST operations
        original_log_level = logging.getLogger().level
        logging.getLogger().setLevel(logging.WARNING)
        
        # Get BLAST hits with initial settings
        hits = direct_blast_search(primer, settings, email, api_key)
        
        if not hits:
            logging.info("No BLAST hits found, attempting Entrez.esearch fallback for approximate matches")
            # Ustawienia Entrez
            Entrez.email = email
            Entrez.api_key = api_key
            # Szukaj za pomocą esearch
            try:
                handle = Entrez.esearch(
                    db="nucleotide",
                    term=primer,
                    retmax=settings.fetch_limit
                )
                record = Entrez.read(handle)
                handle.close()
                ids = record.get("IdList", [])
                if ids:
                    logging.info(f"Found {len(ids)} sequences via Entrez.esearch fallback")
                    hits = [{
                        'accession': seq_id,
                        'description': 'Found via Entrez.esearch fallback',
                        'organism': 'Unknown',
                        'identity_percent': settings.min_identity,
                        'hit_position': 1
                    } for seq_id in ids]
                else:
                    logging.info("Entrez.esearch fallback returned no matches.")
            except Exception as e:
                logging.warning(f"Entrez.esearch fallback failed: {e}")


        # If no hits or few hits (< 5), retry with more permissive settings
        if not hits or len(hits) < 5:
            logging.info(f"Only {len(hits) if hits else 0} hits found with initial parameters, trying more permissive settings...")
            modified_settings = copy.deepcopy(settings)
            modified_settings.word_size = 4  # Most permissive word size
            modified_settings.min_identity = 60.0  # Much lower identity threshold
            modified_settings.max_e_value = 100000  # Very high E-value
            
            # Try again with modified settings
            additional_hits = direct_blast_search(primer, modified_settings, email, api_key)
            
            # Merge unique hits from both searches
            if additional_hits:
                existing_accessions = {hit['accession'] for hit in hits} if hits else set()
                for hit in additional_hits:
                    if hit['accession'] not in existing_accessions:
                        if hits:
                            hits.append(hit)
                        else:
                            hits = [hit]
                        existing_accessions.add(hit['accession'])
            
            # Log the total number of hits after retrying
            if hits:
                logging.info(f"Found {len(hits)} total hits after trying permissive settings")
            else:
                logging.warning(f"No hits found for primer even with permissive settings: {primer}")
                logging.getLogger().setLevel(original_log_level)  # Restore log level
                return None
        
        # Create dictionary of hit positions by accession
        hit_positions = {hit['accession']: hit.get('hit_position', 1) for hit in hits}
        
        # Extract accession IDs from hits
        hit_ids = [hit['accession'] for hit in hits]
        
        # Fetch the sequences from GenBank
        sequences = []
        batch_size = 100  # Process in batches to avoid E-utility limits
        max_sequence_length = 1000  # Limit sequence length to 1000 nucleotides from primer match
        
        for i in range(0, len(hit_ids), batch_size):
            batch_ids = hit_ids[i:i+batch_size]
            
            # Try multiple times in case of network errors
            success = False
            for attempt in range(3):
                try:
                    # Use efetch to get sequences directly rather than using qblast
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=",".join(batch_ids),
                        rettype="fasta",
                        retmode="text"
                    )
                    
                    batch_data = fetch_handle.read()
                    fetch_handle.close()
                    
                    # Verify we have valid FASTA data
                    if batch_data and (">" in batch_data):
                        # Process and validate FASTA entries
                        valid_entries = []
                        entries = batch_data.split("\n>")
                        
                        # Process first entry separately
                        first_entry = entries[0]
                        if first_entry.startswith(">"):
                            # Parse header and sequence
                            lines = first_entry.split("\n")
                            header = lines[0]
                            seq = "".join(lines[1:])
                            
                            # Extract accession from header
                            accession = header.split()[0].strip(">")
                            
                            # Get hit position for this accession
                            hit_pos = hit_positions.get(accession, 1)
                            hit_pos = max(1, hit_pos - 1)  # Convert to 0-based and ensure >= 1
                            
                            # Trim sequence to max_sequence_length starting from hit position
                            if hit_pos < len(seq):
                                trimmed_seq = seq[hit_pos:hit_pos + max_sequence_length]
                            else:
                                trimmed_seq = seq[:max_sequence_length]  # Fallback if position is invalid
                            
                            # Recreate FASTA entry with trimmed sequence
                            trimmed_entry = header + "\n" + trimmed_seq
                            valid_entries.append(trimmed_entry)
                        
                        # Process remaining entries
                        for entry in entries[1:]:
                            # Parse header and sequence
                            lines = entry.split("\n")
                            header = ">" + lines[0]
                            seq = "".join(lines[1:])
                            
                            # Extract accession from header
                            accession = header.split()[0].strip(">")
                            
                            # Get hit position for this accession
                            hit_pos = hit_positions.get(accession, 1)
                            hit_pos = max(1, hit_pos - 1)  # Convert to 0-based and ensure >= 1
                            
                            # Trim sequence to max_sequence_length starting from hit position
                            if hit_pos < len(seq):
                                trimmed_seq = seq[hit_pos:hit_pos + max_sequence_length]
                            else:
                                trimmed_seq = seq[:max_sequence_length]  # Fallback if position is invalid
                            
                            # Recreate FASTA entry with trimmed sequence
                            trimmed_entry = header + "\n" + trimmed_seq
                            valid_entries.append(trimmed_entry)
                        
                        # Add valid entries to sequences list
                        if valid_entries:
                            sequences.append("\n".join(valid_entries))
                            success = True
                            break
                    
                except (HTTPError, URLError, IncompleteRead, socket.timeout) as e:
                    time.sleep(3 * (attempt + 1))  # Exponential backoff
                
                # Rate limiting - wait between requests
                time.sleep(0.1)
                
        if not hits and not sequences:
            hits = fallback_blast_search(primer, settings, email, api_key)
    
        # Write sequences to cache file
        if sequences:
            with open(cache_file, 'w', encoding='utf-8') as f:
                f.write('\n'.join(sequences))
            
            logging.getLogger().setLevel(original_log_level)  # Restore log level
            return cache_file
        else:
            # Create an empty cache file to avoid repeated failed searches
            with open(cache_file, 'w', encoding='utf-8') as f:
                pass
            logging.getLogger().setLevel(original_log_level)  # Restore log level
            return None
        
    except Exception as e:
        logging.error(f"Error fetching sequences for primer {primer}: {e}")
        return None

def create_alignment_string(primer: str, target: str) -> str:
    """
    Create a Primer-BLAST style alignment string with dots for matches
    and letters for mismatches.
    """
    alignment = ""
    min_len = min(len(primer), len(target))
    
    for i in range(min_len):
        if primer[i] == target[i]:  
            alignment += "."  # Match symbol (dot for Primer-BLAST style)
        else:
            alignment += target[i]  # Show the target base at mismatch positions
    
    return alignment

def count_mismatches(primer: str, target: str, settings: SpecificitySettings) -> Tuple[int, int, str]:
    """
    Count total mismatches and 3' end mismatches between primer and target.
    Uses Bio.Align.PairwiseAligner instead of deprecated pairwise2.
    Returns (total mismatches, 3' end mismatches, alignment string)
    """
    # Simple direct comparison for short primers - much faster than alignment
    if len(primer) <= 30 and len(target) <= 30 and abs(len(primer) - len(target)) <= 3:
        # Direct comparison
        min_len = min(len(primer), len(target))
        total_mismatches = sum(1 for i in range(min_len) if primer[i] != target[i])
        
        # Count 3' end mismatches
        three_prime_mismatches = 0
        three_prime_region = min(settings.bases_from_3prime, len(primer))
        
        for i in range(1, three_prime_region + 1):
            p_idx = len(primer) - i
            t_idx = len(target) - i
            if p_idx >= 0 and t_idx >= 0:  # Check indices are valid
                if primer[p_idx] != target[t_idx]:
                    three_prime_mismatches += 1
        
        alignment = create_alignment_string(primer, target)
        return total_mismatches, three_prime_mismatches, alignment
    
    # Use the new PairwiseAligner for more complex cases
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    try:
        # Limit to one alignment for speed
        alignments = aligner.align(primer, target)
        alignment_obj = next(iter(alignments))
        
        # Get the aligned sequences
        alignment_str = str(alignment_obj)
        alignment_lines = alignment_str.split('\n')
        
        if len(alignment_lines) >= 3:
            primer_aligned = alignment_lines[0]
            target_aligned = alignment_lines[2]
            
            # Count total mismatches (excluding gaps)
            total_mismatches = 0
            for p, t in zip(primer_aligned, target_aligned):
                if p != t and p != '-' and t != '-':
                    total_mismatches += 1
            
            # Count 3' end mismatches
            three_prime_mismatches = 0
            three_prime_region = settings.bases_from_3prime
            
            # Find the actual primer sequence (no gaps) in the alignment
            primer_aligned_no_gaps = ''.join(p for p in primer_aligned if p != '-')
            
            # Find the last non-gap position of the primer
            last_pos = -1
            for i in range(len(primer_aligned)-1, -1, -1):
                if primer_aligned[i] != '-':
                    last_pos = i
                    break
            
            # Count mismatches in the 3' region
            if last_pos >= 0:
                count = 0
                for i in range(last_pos, -1, -1):
                    if primer_aligned[i] != '-':
                        if count < three_prime_region:
                            if (i < len(target_aligned) and 
                                primer_aligned[i] != target_aligned[i] and 
                                target_aligned[i] != '-'):
                                three_prime_mismatches += 1
                            count += 1
                        else:
                            break
            
            # Create Primer-BLAST style alignment string
            alignment = create_alignment_string(primer, target)
            return total_mismatches, three_prime_mismatches, alignment
    
    except Exception as e:
        # Fallback to simple comparison if alignment fails
        logging.debug(f"Alignment failed: {e}. Using simple comparison.")
        min_len = min(len(primer), len(target))
        total_mismatches = sum(1 for i in range(min_len) if primer[i] != target[i])
        
        # Count 3' end mismatches
        three_prime_mismatches = 0
        three_prime_region = min(settings.bases_from_3prime, len(primer))
        for i in range(1, three_prime_region + 1):
            if len(primer) - i < len(primer) and len(target) - i < len(target):
                if primer[len(primer) - i] != target[len(target) - i]:
                    three_prime_mismatches += 1
        
        alignment = create_alignment_string(primer, target)
        return total_mismatches, three_prime_mismatches, alignment

def fallback_blast_search(primer, settings, email, api_key):
    """Fallback BLAST search using Biopython's NCBIWWW module"""
    from Bio.Blast import NCBIWWW, NCBIXML
    
    logging.info(f"Using fallback BLAST method for primer: {primer}")
    
    try:
        # Set NCBI options
        Entrez.email = email
        Entrez.api_key = api_key
        
        # Perform BLAST search
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=primer,
            expect=1000,
            word_size=7,
            megablast=True,
            hitlist_size=settings.max_hits
        )
        
        # Parse results
        blast_records = NCBIXML.parse(result_handle)
        record = next(blast_records)
        
        # Extract hits
        hits = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                identity_percent = (hsp.identities / hsp.align_length) * 100
                if identity_percent >= settings.min_identity:
                    # Extract organism name if possible
                    organism = "Unknown"
                    org_match = re.search(r'\[(.*?)\]', alignment.title)
                    if org_match:
                        organism = org_match.group(1)
                    
                    hits.append({
                        'accession': alignment.accession,
                        'description': alignment.title,
                        'organism': organism,
                        'identity_percent': identity_percent,
                        'hit_position': hsp.sbjct_start  # Add position information
                    })
                    break  # Only count each hit once
            
            if len(hits) >= settings.fetch_limit:
                break
                
        return hits
        
    except Exception as e:
        logging.error(f"Fallback BLAST search failed: {e}")
        return []

def process_sequence(seq_record, primer_pair, settings):
    """Process a single sequence for primer matches - runs in a separate thread"""
    # Extract the organism name from the description if possible
    description = seq_record.description
    organism = "Unknown"
    
    # Try to extract organism name using common patterns
    org_match = re.search(r'\[(.*?)\]', description)
    if org_match:
        organism = org_match.group(1)
        
    # Convert sequence to uppercase string for searching
    seq_str = str(seq_record.seq).upper()
    
    # Look for matches to the left primer
    left_primer = primer_pair.left_primer
    right_rc = reverse_complement(primer_pair.right_primer)  # Already in forward orientation
    probe = primer_pair.probe
    
    # Create a match object
    match = SpecificityMatch(
        accession=seq_record.id,
        organism=organism,
        description=description.replace(seq_record.id, '').strip()
    )
    
    # Search for approximate matches to left primer using sliding window
    left_best_pos = -1
    left_best_mismatches = len(left_primer)  # Initialize with worst case
    left_best_3prime_mismatches = settings.bases_from_3prime  # Initialize with worst case
    left_best_alignment = ""
    
    # Skip if sequence is too short
    if len(seq_str) < len(left_primer):
        return None
    
    # Use a more efficient approach - sample positions rather than checking every position
    # This dramatically speeds up processing for long sequences
    seq_len = len(seq_str)
    step_size = 1  # Use 1 for short sequences
    
    if seq_len > 100000:      # sequences longer than 100 kb
        step_size = 10         # sample every 10th base
    else:
        step_size = 1          # full scan for sequences <=100 kb
    
    # Try positions with sampling to speed up
    for i in range(0, len(seq_str) - len(left_primer) + 1, step_size):
        target_region = seq_str[i:i+len(left_primer)]
        total_mm, three_prime_mm, alignment = count_mismatches(left_primer, target_region, settings)
        
        # Update if this is a better match
        if (total_mm < left_best_mismatches or 
            (total_mm == left_best_mismatches and three_prime_mm < left_best_3prime_mismatches)):
            left_best_mismatches = total_mm
            left_best_3prime_mismatches = three_prime_mm
            left_best_pos = i
            left_best_alignment = alignment
    
    # Check if we found a sufficiently good match for left primer
    if left_best_pos >= 0:
        match.left_position = left_best_pos
        match.left_mismatches = left_best_mismatches
        match.left_3prime_mismatches = left_best_3prime_mismatches
        match.left_alignment = left_best_alignment
        
        # Apply the specificity criteria
        if (left_best_mismatches <= settings.max_total_mismatches + 1 or
            (left_best_3prime_mismatches <= settings.max_3prime_mismatches + 1 and 
             left_best_mismatches <= settings.max_total_mismatches + 2)):
            match.left_match = True
        
        # Now check for right primer in expected amplicon region
        expected_amplicon_max = 2000  # Maximum reasonable amplicon size to search for
        search_end = min(left_best_pos + expected_amplicon_max, len(seq_str))
        expected_amplicon_region = seq_str[left_best_pos:search_end]
        
        # Search for approximate matches to right primer
        right_best_pos = -1
        right_best_mismatches = len(right_rc)
        right_best_3prime_mismatches = settings.bases_from_3prime
        right_best_alignment = ""
        
        if len(expected_amplicon_region) >= len(right_rc):
            for i in range(0, len(expected_amplicon_region) - len(right_rc) + 1, step_size):
                target_region = expected_amplicon_region[i:i+len(right_rc)]
                total_mm, three_prime_mm, alignment = count_mismatches(right_rc, target_region, settings)
                
                # Update if this is a better match
                if (total_mm < right_best_mismatches or 
                    (total_mm == right_best_mismatches and three_prime_mm < right_best_3prime_mismatches)):
                    right_best_mismatches = total_mm
                    right_best_3prime_mismatches = three_prime_mm
                    right_best_pos = i
                    right_best_alignment = alignment
        
        # Check if we found a match for right primer
        if right_best_pos >= 0:
            right_absolute_pos = left_best_pos + right_best_pos
            match.right_position = right_absolute_pos
            match.right_mismatches = right_best_mismatches
            match.right_3prime_mismatches = right_best_3prime_mismatches
            match.right_alignment = right_best_alignment
            
            # Apply the specificity criteria for right primer
            if (right_best_mismatches <= settings.max_total_mismatches and
                right_best_3prime_mismatches <= settings.max_3prime_mismatches):
                match.right_match = True
            
            # Calculate amplicon size
            amplicon_size = right_absolute_pos + len(right_rc) - left_best_pos
            match.amplicon_size = amplicon_size
        
        # Search for approximate matches to probe
        probe_best_pos = -1
        probe_best_mismatches = len(probe)
        probe_best_alignment = ""
        
        if len(expected_amplicon_region) >= len(probe):
            for i in range(0, len(expected_amplicon_region) - len(probe) + 1, step_size):
                target_region = expected_amplicon_region[i:i+len(probe)]
                
                # For probe, we only care about total mismatches, not 3' mismatches
                # Using a simplified version of count_mismatches
                total_mm = sum(1 for j in range(len(probe)) if j < len(target_region) and probe[j] != target_region[j])
                alignment = create_alignment_string(probe, target_region)
                
                # Update if this is a better match
                if total_mm < probe_best_mismatches:
                    probe_best_mismatches = total_mm
                    probe_best_pos = i
                    probe_best_alignment = alignment
        
        # Check if we found a match for probe
        if probe_best_pos >= 0:
            probe_absolute_pos = left_best_pos + probe_best_pos
            match.probe_position = probe_absolute_pos
            match.probe_mismatches = probe_best_mismatches
            match.probe_alignment = probe_best_alignment
            
            # Apply the specificity criteria for probe
            if probe_best_mismatches <= settings.max_mismatches_probe:
                match.probe_match = True
        
        # Return the match if we found at least one component match
        if match.left_match or match.right_match or match.probe_match:
            return match
    
    return None

def align_sequences(primer_pair: PrimerPair, sequence_file: str, settings: Optional[SpecificitySettings] = None) -> List[SpecificityMatch]:
    """
    Align primer pair components to sequences and identify matches based on specificity criteria.
    Uses parallel processing for efficiency.
    """
    matches = []
    
    # Use default settings if none provided
    if settings is None:
        settings = SpecificitySettings()
    
    try:
        if not os.path.exists(sequence_file):
            return []
            
        # Parse sequences from the FASTA file
        sequences = []
        with open(sequence_file, 'r', encoding='utf-8') as f:
            sequences = list(SeqIO.parse(f, 'fasta'))
            
        if not sequences:
            return []
        
        # Process sequences in parallel using a thread pool
        with concurrent.futures.ThreadPoolExecutor(max_workers=settings.max_threads) as executor:
            # Submit all sequences for processing
            future_to_seq = {executor.submit(process_sequence, seq, primer_pair, settings): seq for seq in sequences}
            
            # Process as they complete
            for future in concurrent.futures.as_completed(future_to_seq):
                try:
                    # Get the result and add to matches if valid
                    match = future.result()
                    if match is not None:
                        matches.append(match)
                except Exception:
                    pass
        
    except Exception:
        pass
    
    return matches

def create_html_report(all_results: Dict[str, List[SpecificityMatch]], primer_pairs: Dict[str, PrimerPair], settings: SpecificitySettings, target_taxids: List[str] = None, input_file: str = None) -> str:
    """Create an HTML report with JavaScript for interactive features in a Primer-BLAST style."""
    # Create output directory
    os.makedirs(settings.output_dir, exist_ok=True)
    
    # Prepare data for the report
    report_data = []
    
    # Format the current date and time - only once when the report is generated
    report_generation_time = time.strftime("%d.%m.%Y, %H:%M:%S")
    
    # Process each primer set
    for primer_id, matches in all_results.items():
        if primer_id not in primer_pairs:
            continue
            
        pair = primer_pairs[primer_id]
        
        # Count matching sequences based on new specificity criteria
        valid_matches = [
            m for m in matches if (
                m.left_match and 
                m.right_match and 
                m.probe_match
            )
        ]
        
        total_matches = len(valid_matches)
        
        # Count matches to target taxids if provided
        target_taxid_matches = 0
        if target_taxids and valid_matches:
            target_taxid_matches = sum(1 for m in valid_matches if m.taxid in target_taxids)
        
        # Sort matches by match score (all three components first)
        sorted_matches = sorted(matches, key=lambda m: (-m.match_score, m.organism))
        
        # Add to report data
        report_data.append({
            'primer_id': primer_id,
            'left_primer': pair.left_primer,
            'right_primer': pair.right_primer,  # Keep in forward orientation
            'probe': pair.probe,
            'product_size': pair.product_size,
            'left_tm': pair.left_tm,
            'right_tm': pair.right_tm,
            'probe_tm': pair.probe_tm,
            'left_gc': pair.left_gc,
            'right_gc': pair.right_gc,
            'probe_gc': pair.probe_gc,
            'total_matches': total_matches,
            'target_taxid_matches': target_taxid_matches,
            'target_taxids': target_taxids if target_taxids else [],
            'detailed_matches': [
                {
                    'accession': m.accession,
                    'organism': m.organism,
                    'description': m.description,
                    'left_match': m.left_match,
                    'right_match': m.right_match,
                    'probe_match': m.probe_match,
                    'amplicon_size': m.amplicon_size if m.amplicon_size else 0,
                    'left_position': m.left_position,
                    'right_position': m.right_position,
                    'probe_position': m.probe_position,
                    'left_mismatches': m.left_mismatches,
                    'right_mismatches': m.right_mismatches,
                    'probe_mismatches': m.probe_mismatches,
                    'left_3prime_mismatches': m.left_3prime_mismatches,
                    'right_3prime_mismatches': m.right_3prime_mismatches,
                    'left_alignment': m.left_alignment,
                    'right_alignment': m.right_alignment,
                    'probe_alignment': m.probe_alignment,
                    'taxid': m.taxid,
                    'is_target_taxid': m.taxid in target_taxids if target_taxids else False
                } for m in sorted_matches
            ]
        })
    
    # Sort report data by primer ID
    report_data.sort(key=lambda x: int(x['primer_id']) if x['primer_id'].isdigit() else x['primer_id'])
    
    # Generate filename based on input file if provided
    report_filename = "specificity_report.html"
    if input_file:
        base_name = os.path.basename(input_file)
        if base_name.endswith("_primers.csv"):
            report_name_part = base_name.replace("_primers.csv", "")
            report_filename = f"specificity_report_{report_name_part}_primers.html"
    
    report_path = os.path.join(settings.output_dir, report_filename)
    
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Primer Specificity Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            color: #333;
        }
        h1, h2, h3, h4 {
            color: #2c3e50;
        }
        .summary {
            background-color: #f8f9fa;
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 20px;
            border-left: 5px solid #4caf50;
        }
        .stats-container {
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            margin-bottom: 30px;
        }
        .stat-card {
            flex: 1;
            min-width: 250px;
            background-color: #fff;
            border-radius: 5px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            padding: 15px;
        }
        .stat-value {
            font-size: 24px;
            font-weight: bold;
            color: #3498db;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 20px;
        }
        th, td {
            padding: 8px 12px;
            border-bottom: 1px solid #ddd;
            text-align: left;
        }
        th {
            background-color: #f2f2f2;
            font-weight: bold;
            position: sticky;
            top: 0;
        }
        .sortable {
            cursor: pointer;
        }
        .sortable:after {
            content: ' ⇵';
            font-size: 12px;
            color: #999;
        }
        .sequence-alignment {
            font-family: 'Courier New', Courier, monospace;
            white-space: pre;
            margin: 10px 0;
            font-size: 14px;
            letter-spacing: 0;
            line-height: 1.5;
            tab-size: 4;
        }
        tr:hover {
            background-color: #f5f5f5;
        }
        .primer-details {
            display: none;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
            margin-top: 10px;
            border: 1px solid #ddd;
        }
        .toggle-details {
            cursor: pointer;
            color: #007bff;
            text-decoration: underline;
        }
        .filter-container {
            margin-bottom: 20px;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
        }
        .filter-container input, .filter-container select {
            padding: 8px;
            margin-right: 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
        }
        .filter-container button {
            padding: 8px 15px;
            background-color: #4caf50;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }
        .filter-container button:hover {
            background-color: #45a049;
        }
        .sequence-alignment {
            font-family: 'Courier New', Courier, monospace;
            white-space: pre;
            margin: 10px 0;
            font-size: 14px;
            letter-spacing: 0;
        }
        .product-details {
            margin-top: 5px;
            margin-bottom: 15px;
            padding-left: 20px;
        }
        .product-template {
            margin-bottom: 30px;
            border-left: 3px solid #4caf50;
            padding-left: 15px;
        }
        .no-matches {
            color: #856404;
            background-color: #fff3cd;
            padding: 10px;
            border-radius: 5px;
            margin-top: 10px;
        }
        .primer-table {
            width: 100%;
            margin-bottom: 20px;
        }
        .primer-table th, .primer-table td {
            padding: 8px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        .primer-table th {
            background-color: #f2f2f2;
        }
        .settings-info {
            background-color: #e3f2fd;
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 20px;
        }
        .settings-info h3 {
            margin-top: 0;
            color: #0d47a1;
        }
        .settings-table {
            width: 100%;
            max-width: 600px;
        }
        .settings-table th {
            width: 60%;
        }
    </style>
</head>
<body>
    <h1>Primer Specificity Analysis Report</h1>
    
    <div class="summary">
        <h2>Summary</h2>
        <p>Analysis of <span id="primerCount">0</span> primer pairs for cross-reactivity with sequences in GenBank.</p>
        <p>Generated on <span id="reportDate">REPORT_DATE_PLACEHOLDER</span></p>
    </div>
    
    <div class="settings-info">
        <h3>Specificity Settings</h3>
        <table class="settings-table">
            <tr>
                <th>Maximum total mismatches allowed</th>
                <td>MAX_TOTAL_MISMATCHES_PLACEHOLDER</td>
            </tr>
            <tr>
                <th>Maximum 3' end mismatches allowed</th>
                <td>MAX_3PRIME_MISMATCHES_PLACEHOLDER</td>
            </tr>
            <tr>
                <th>3' end region size (bases)</th>
                <td>BASES_FROM_3PRIME_PLACEHOLDER</td>
            </tr>
            <tr>
                <th>Maximum probe mismatches allowed</th>
                <td>MAX_MISMATCHES_PROBE_PLACEHOLDER</td>
            </tr>
            <tr>
                <th>Target NCBI Taxonomy IDs</th>
                <td>TARGET_TAXIDS_PLACEHOLDER</td>
            </tr>
        </table>
    </div>
    
    <div class="stat-card">
        <h3>Primers Analyzed</h3>
        <div class="stat-value" id="totalPrimers">0</div>
    </div>
    
    <div class="filter-container">
        <h3>Filter Primers</h3>
        <div style="display: flex; flex-wrap: wrap; gap: 10px; margin-bottom: 10px;">
            <div>
                <label for="searchInput">Search by ID or sequence:</label>
                <input type="text" id="searchInput" placeholder="Enter primer ID or sequence">
            </div>
            
    <div>
        <label for="minAmpliconSize">Amplicon size range:</label>
        <input type="number" id="minAmpliconSize" placeholder="Min" min="0" style="width: 70px;">
        <label for="maxAmpliconSize">to</label>
        <input type="number" id="maxAmpliconSize" placeholder="Max" min="0" style="width: 70px;">
    </div>

    </div>
    <div style="margin-top: 10px;">
        <button onclick="applyFilters()">Apply Filters</button>
        <button onclick="resetFilters()">Reset</button>
    </div>
    </div>  <!-- Closing div.filter-container -->
        
    <h2>Primer Sets</h2>
    <table id="primersTable">
        <thead>
            <tr>
                <th>ID</th>
                <th>Forward Primer</th>
                <th>Reverse Primer</th>
                <th>Probe</th>
                <th>Product Size</th>
                <th class="sortable" onclick="sortTable(5)">Matched NCBI GenBank Records</th>
                <th class="sortable" onclick="sortTable(6)">Matches to Target Taxa</th>
                <th>Details</th>
            </tr>
        </thead>
        <tbody id="primersTableBody">
            <!-- Table rows will be inserted here by JavaScript -->
        </tbody>
    </table>

    <script>
        // Report data
        const reportData = REPORT_DATA_PLACEHOLDER;
        
        // Update summary statistics
        document.getElementById('primerCount').textContent = reportData.length;
        document.getElementById('totalPrimers').textContent = reportData.length;
        document.getElementById('reportDate').textContent = "REPORT_DATE_PLACEHOLDER";
        
        // Populate the primers table
        function populateTable(data) {
            const tableBody = document.getElementById('primersTableBody');
            tableBody.innerHTML = '';
            
            data.forEach((primer, index) => {
                const row = document.createElement('tr');
                
                // Handle case where no target taxids were specified
                const targetTaxaCell = primer.target_taxids && primer.target_taxids.length > 0 
                    ? `<td>${primer.target_taxid_matches} / ${primer.total_matches}</td>`
                    : `<td>N/A</td>`;
                
                row.innerHTML = `
                    <td>${primer.primer_id}</td>
                    <td>${primer.left_primer}</td>
                    <td>${primer.right_primer}</td>
                    <td>${primer.probe}</td>
                    <td>${primer.product_size}</td>
                    <td>${primer.total_matches}</td>
                    ${targetTaxaCell}
                    <td><span class="toggle-details" onclick="toggleDetails(${index})">Show Details</span></td>
                `;
                
                tableBody.appendChild(row);
                
                // Create details row
                const detailsRow = document.createElement('tr');
                detailsRow.id = `details-${index}`;
                
                const detailsCell = document.createElement('td');
                detailsCell.colSpan = 8;
                detailsCell.innerHTML = `
                    <div class="primer-details">
                        <h3>Primer pair ${primer.primer_id}</h3>
                        
                        <table class="primer-table">
                            <thead>
                                <tr>
                                    <th>Sequence (5'->3')</th>
                                    <th>Length</th>
                                    <th>Tm</th>
                                    <th>GC%</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td>Forward primer: ${primer.left_primer}</td>
                                    <td>${primer.left_primer.length}</td>
                                    <td>${primer.left_tm.toFixed(2)}</td>
                                    <td>${primer.left_gc.toFixed(2)}</td>
                                </tr>
                                <tr>
                                    <td>Reverse primer: ${primer.right_primer}</td>
                                    <td>${primer.right_primer.length}</td>
                                    <td>${primer.right_tm.toFixed(2)}</td>
                                    <td>${primer.right_gc.toFixed(2)}</td>
                                </tr>
                                <tr>
                                    <td>Probe: ${primer.probe}</td>
                                    <td>${primer.probe.length}</td>
                                    <td>${primer.probe_tm.toFixed(2)}</td>
                                    <td>${primer.probe_gc.toFixed(2)}</td>
                                </tr>
                            </tbody>
                        </table>
                        
                        <h4>Products on target templates</h4>
                        ${primer.detailed_matches.length > 0 
                            ? formatDetailedMatches(primer) 
                            : '<div class="no-matches">No matching sequences found.</div>'}
                    </div>
                `;
                
                detailsCell.querySelector('.primer-details').style.display = 'none';
                detailsRow.appendChild(detailsCell);
                tableBody.appendChild(detailsRow);
            });
        }
        
        // Format detailed matches in Primer-BLAST style
        function formatDetailedMatches(primer) {
            let matchesHtml = '';
            
            // Count valid matches (those where all three components match)
            const validMatches = primer.detailed_matches.filter(match => 
                match.left_match && match.right_match && match.probe_match
            );
            
            if (validMatches.length === 0) {
                return '<div class="no-matches">No sequences meet all specificity criteria.</div>';
            }
            
            // Add CSS for target taxids
            matchesHtml += `
                <style>
                    .target-taxid {
                        background-color: #e8f5e9;
                        border-left: 3px solid #2e7d32;
                    }
                </style>
            `;
            
            // DNA complement lookup table
            const complementMap = {
                'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                '.': '.'  // Keep dots as dots
            };
            
            // Function to complement a DNA sequence
            const complement = (seq) => {
                return seq.split('').map(base => complementMap[base] || base).join('');
            };
            
            validMatches.forEach(match => {
                // Calculate the product
                const productLength = match.amplicon_size;
                
                // Format positions with proper padding
                const leftPosDisplay = match.left_position !== null ? match.left_position + 1 : '?';
                const rightPosDisplay = match.right_position !== null ? match.right_position + 1 : '?';
                const probePosDisplay = match.probe_position !== null ? match.probe_position + 1 : '?';
                
                // Calculate end positions
                const leftEndPos = match.left_position !== null ? match.left_position + primer.left_primer.length : '?';
                const rightEndPos = match.right_position !== null ? match.right_position + primer.right_primer.length : '?';
                const probeEndPos = match.probe_position !== null ? match.probe_position + primer.probe.length : '?';
                
                // Process the right alignment to make complementary
                const rightAlignment = match.right_alignment || '.'.repeat(primer.right_primer.length);
                const complementaryRightAlignment = complement(rightAlignment).split('').reverse().join('');
                
                // Add class for target taxids
                const targetClass = match.is_target_taxid ? 'target-taxid' : '';
                const taxidInfo = match.taxid ? ` [TaxID: ${match.taxid}]` : '';
                
                
                // Calculate padding for position numbers to handle very large GenBank positions
                const maxPosLength = Math.max(
                    leftPosDisplay.toString().length,
                    leftEndPos.toString().length,
                    rightPosDisplay.toString().length,
                    rightEndPos.toString().length,
                    probePosDisplay.toString().length,
                    probeEndPos.toString().length
                );

                // Ensure we have proper monospace alignment by using pre-formatted text with consistent spacing
                matchesHtml += `
                    <div class="product-template ${targetClass}">
                        <div>&gt;${match.accession} ${match.description}${taxidInfo}</div>
                        <div>product length = ${productLength}</div>

                        <pre class="sequence-alignment">Forward primer  1${' '.repeat(maxPosLength-1)}    ${primer.left_primer}  ${primer.left_primer.length}
Template        ${leftPosDisplay.toString().padStart(maxPosLength, ' ')}    ${match.left_alignment || '.'.repeat(primer.left_primer.length)}  ${leftEndPos}</pre>

            <pre class="sequence-alignment">Reverse primer  1${' '.repeat(maxPosLength-1)}    ${primer.right_primer}  ${primer.right_primer.length}
Template        ${rightEndPos.toString().padStart(maxPosLength, ' ')}    ${complementaryRightAlignment}  ${rightPosDisplay}</pre>

            <pre class="sequence-alignment">Probe           1${' '.repeat(maxPosLength-1)}    ${primer.probe}  ${primer.probe.length}
Template        ${probePosDisplay.toString().padStart(maxPosLength, ' ')}    ${match.probe_alignment || '.'.repeat(primer.probe.length)}  ${probeEndPos}</pre>
                    </div>
                `;
            });
            
            return matchesHtml;
        }
        
        // Toggle details visibility
        function toggleDetails(index) {
            const detailsRow = document.getElementById(`details-${index}`);
            const detailsDiv = detailsRow.querySelector('.primer-details');
            
            if (detailsDiv.style.display === 'none') {
                detailsDiv.style.display = 'block';
            } else {
                detailsDiv.style.display = 'none';
            }
        }

        // Filter functions
        function applyFilters() {
            const searchText = document.getElementById('searchInput').value.toLowerCase();
            const minAmpliconSize = parseInt(document.getElementById('minAmpliconSize').value) || 0;
            const maxAmpliconSize = parseInt(document.getElementById('maxAmpliconSize').value) || Number.MAX_SAFE_INTEGER;
            
            let filteredData = [...reportData];
            
            // Apply text search if specified
            if (searchText) {
                filteredData = filteredData.filter(primer => 
                    primer.primer_id.toString().toLowerCase().includes(searchText) ||
                    primer.left_primer.toLowerCase().includes(searchText) ||
                    primer.right_primer.toLowerCase().includes(searchText) ||
                    primer.probe.toLowerCase().includes(searchText)
                );
            }
            
            // Apply amplicon size filter
            if (minAmpliconSize > 0 || maxAmpliconSize < Number.MAX_SAFE_INTEGER) {
                filteredData = filteredData.filter(primer => 
                    primer.product_size >= minAmpliconSize && 
                    primer.product_size <= maxAmpliconSize
                );
            }
            
            populateTable(filteredData);
        }
        
        function resetFilters() {
            document.getElementById('searchInput').value = '';
            document.getElementById('minAmpliconSize').value = '';
            document.getElementById('maxAmpliconSize').value = '';
            populateTable(reportData);
        }
        
        // Table sorting functionality
        let sortAscending = true;
        
        function sortTable(columnIndex) {
            sortAscending = !sortAscending;
            
            // Sort the actual data
            const sortedData = [...reportData].sort((a, b) => {
                let aValue, bValue;
                
                if (columnIndex === 5) {  // Matched NCBI GenBank Records
                    aValue = parseInt(a.total_matches) || 0;
                    bValue = parseInt(b.total_matches) || 0;
                } 
                else if (columnIndex === 6) {  // Matches to Target Taxa
                    aValue = parseInt(a.target_taxid_matches) || 0;
                    bValue = parseInt(b.target_taxid_matches) || 0;
                }
                
                return sortAscending ? aValue - bValue : bValue - aValue;
            });
            
            // Repopulate the table with sorted data
            populateTable(sortedData);
        }
        
        // Initialize the page
        populateTable(reportData);
    </script>
</body>
</html>"""

# Replace placeholders with actual data
    html_content = html_template.replace('REPORT_DATA_PLACEHOLDER', json.dumps(report_data))
    html_content = html_content.replace('REPORT_DATE_PLACEHOLDER', report_generation_time)
    html_content = html_content.replace('MAX_TOTAL_MISMATCHES_PLACEHOLDER', str(settings.max_total_mismatches))
    html_content = html_content.replace('MAX_3PRIME_MISMATCHES_PLACEHOLDER', str(settings.max_3prime_mismatches))
    html_content = html_content.replace('BASES_FROM_3PRIME_PLACEHOLDER', str(settings.bases_from_3prime))
    html_content = html_content.replace('MAX_MISMATCHES_PROBE_PLACEHOLDER', str(settings.max_mismatches_probe))
    
    # Format target taxids for display
    target_taxids_display = "None specified" if not target_taxids or len(target_taxids) == 0 else ", ".join(target_taxids)
    html_content = html_content.replace('TARGET_TAXIDS_PLACEHOLDER', target_taxids_display)
    
    # Write the HTML file
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logging.info(f"Created HTML report at {report_path}")
    return report_path

def open_in_browser(file_path: str) -> None:
    """Open an HTML file in the default web browser."""
    try:
        # Ensure the file exists
        if not os.path.exists(file_path):
            logging.error(f"File {file_path} does not exist")
            return
            
        # Convert to absolute path (required for some systems)
        file_path = os.path.abspath(file_path)
        
        # Use webbrowser module to open file
        logging.info(f"Opening report in default browser: {file_path}")
        webbrowser.open(f'file://{file_path}')
        
    except Exception as e:
        logging.error(f"Error opening browser: {e}")
        # Try alternative methods if webbrowser module fails
        try:
            if sys.platform == 'darwin':  # macOS
                subprocess.run(['open', file_path], check=True)
            elif sys.platform == 'win32':  # Windows
                os.startfile(file_path)
            elif sys.platform.startswith('linux'):  # Linux
                subprocess.run(['xdg-open', file_path], check=True)
        except Exception as e2:
            logging.error(f"Alternative method to open browser failed: {e2}")
            logging.info(f"Please open the report manually at: {file_path}")

def main():
    try:
        setup_logging(__file__)
        
        logging.info("PrimeSpecPCR - Primer Specificity Test Module")
        logging.info("============================================")
        
        # Check if running in interactive mode
        interactive_mode = "--interactive" in sys.argv
        
        # Get email and API key
        email, api_key = get_user_data()
        if not email or not api_key:
            logging.error("Email or API key not found. Please ensure they are set in the GUI.")
            return
            
        logging.info(f"Using email: {email}")
        logging.info(f"Using API key: {api_key[:5]}...{api_key[-5:]}")
        
        # Ask for target taxonomy IDs
        target_taxids = get_user_target_taxids()
        
        # Load specificity test settings
        settings = load_specificity_settings(interactive_mode)
        if not settings.validate():
            logging.error("Invalid settings. Please correct the settings file and try again.")
            return

        # Ask if user wants to clear cache
        if interactive_mode:
            clear_cache = get_console_input("Do you want to clear the cache and force new BLAST searches? (yes/no): ", "no").lower()
            if clear_cache in ['yes', 'y']:
                cache_dir = settings.cache_dir
                if os.path.exists(cache_dir):
                    logging.info(f"Clearing cache directory: {cache_dir}")
                    for file in os.listdir(cache_dir):
                        try:
                            os.remove(os.path.join(cache_dir, file))
                        except Exception as e:
                            logging.error(f"Error removing cache file {file}: {e}")
                    logging.info("Cache cleared successfully")
        
        # Find primer files
        primer_files = find_primer_files()
        if not primer_files:
            logging.error("No primer files found in the input directory (3_).")
            return
            
        logging.info(f"Found {len(primer_files)} primer files")
        
        # Ask user to select files to analyze
        if interactive_mode and len(primer_files) > 1:
            logging.info("\nAvailable primer files:")
            for i, file in enumerate(primer_files):
                logging.info(f"{i+1}: {os.path.basename(file)}")
                
            file_selection = get_console_input(
                "Enter file numbers to analyze (comma-separated, or 'all'): ",
                "all"
            )
            
            if file_selection.lower() != "all":
                try:
                    indices = [int(i.strip()) - 1 for i in file_selection.split(",")]
                    primer_files = [primer_files[i] for i in indices if 0 <= i < len(primer_files)]
                except (ValueError, IndexError):
                    logging.warning("Invalid selection. Analyzing all files.")
        
        # Process each primer file
        taxid_cache = {}  # Cache for taxonomy IDs to avoid repeated lookups

        for file_idx, file_path in enumerate(primer_files):
            logging.info(f"\nProcessing file {file_idx + 1} of {len(primer_files)}: {os.path.basename(file_path)}")
            
            # Create new dictionaries for each file to avoid mixing results
            file_primer_pairs = {}  # Dictionary to store primer information by ID for this file
            file_results = {}  # Dictionary to store specificity results by primer ID for this file
            
            # Read primer pairs from file
            primer_pairs = read_primer_pairs(file_path)
            if not primer_pairs:
                logging.warning(f"No valid primer pairs found in {file_path}")
                continue
                
            # Store primer info for this file
            for pair in primer_pairs:
                file_primer_pairs[pair.id] = pair
            
            # Process primer pairs by reusing results for identical primers
            processed_left_primers = {}  # Cache for primers already processed
            
            for pair_idx, pair in enumerate(primer_pairs):
                logging.info(f"Processing primer pair {pair.id} ({pair_idx + 1} of {len(primer_pairs)})")
                
                # Check if we've already processed this left primer
                if pair.left_primer in processed_left_primers:
                    logging.info(f"Using cached results for left primer: {pair.left_primer}")
                    sequence_file = processed_left_primers[pair.left_primer]
                else:
                    # Fetch sequences matching the left primer
                    sequence_file = fetch_matching_sequences(pair.left_primer, settings, email, api_key)
                    processed_left_primers[pair.left_primer] = sequence_file
                
                if not sequence_file:
                    logging.warning(f"No sequences found for left primer: {pair.left_primer}")
                    file_results[pair.id] = []
                    continue

                # Align sequences to the primer pair with specificity settings - but suppress detailed logging
                original_log_level = logging.getLogger().level
                logging.getLogger().setLevel(logging.WARNING)  # Temporarily increase log level
                matches = align_sequences(pair, sequence_file, settings)
                logging.getLogger().setLevel(original_log_level)  # Restore original log level

                # Count matches meeting the new specificity criteria
                valid_matches = [m for m in matches if m.left_match and m.right_match and m.probe_match]
                valid_count = len(valid_matches)
                logging.info(f"Found {valid_count} matching sequences.")

                # If we have target taxids, fetch taxonomy info for only valid matches
                if target_taxids and valid_matches:
                    logging.info(f"Fetching taxonomy information for {valid_count} matches, please wait..")
                    
                    for match in valid_matches:
                        if match.accession in taxid_cache:
                            match.taxid = taxid_cache[match.accession]
                        else:
                            # Fetch taxid for this accession without detailed logging
                            taxid = get_taxid_for_accession(match.accession, email, api_key)
                            match.taxid = taxid
                            taxid_cache[match.accession] = taxid
                            
                            # Rate limiting - prevent overloading NCBI servers
                            time.sleep(0.1)
                    
                    logging.info("Success")
                
                # Store results for this file
                file_results[pair.id] = matches
                
                # Progress update with percentage
                progress = (pair_idx + 1) / len(primer_pairs) * 100
                logging.info(f"Progress: {progress:.1f}%")
            
            # Create report for this file
            if file_results:
                report_path = create_html_report(file_results, file_primer_pairs, settings, target_taxids, file_path)
                
                # Open report in browser
                open_in_browser(report_path)
                
                logging.info(f"\nSpecificity analysis complete! Report saved to {report_path}")
                logging.info(f"You can open an HTML file in the  web browser.")
            else:
                logging.error(f"No specificity results to report for file {os.path.basename(file_path)}.")
    
    except Exception as e:
        logging.error(f"Unexpected error in main function: {e}")
        logging.error(traceback.format_exc())

if __name__ == "__main__":
    main()
