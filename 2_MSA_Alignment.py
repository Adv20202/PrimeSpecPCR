import sys
import os
import glob
import subprocess
import logging
import time
import tempfile
import platform
import shutil
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
from io import StringIO
import re
import collections

# Use Biopython for sequence handling
try:
    from Bio import SeqIO, AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython>=1.79", "--upgrade"])
    from Bio import SeqIO, AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

# Standard logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
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
    """Safely write content to a file, creating directories if needed."""
    try:
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except OSError as e:
        raise IOError(f"Failed to write to {filepath}: {e}")

def find_fasta_files(directory: str) -> List[str]:
    """Find all FASTA files in the specified directory."""
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist")
    
    extensions = ['*.fasta', '*.fa', '*.fna', '*.ffn', '*.faa', '*.frn']
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(directory, ext)))
    return sorted(files)

def validate_group_indices(groups: List[str], num_files: int) -> List[List[int]]:
    """Validate the group indices provided by the user."""
    group_indices: List[List[int]] = []
    
    for group in groups:
        indices_str = group.strip().split(',')
        try:
            # Odejmij 1 od każdego indeksu, aby przekształcić indeksy od 1 do indeksów od 0
            indices = [int(i.strip()) - 1 for i in indices_str]
            if any(i < 0 or i >= num_files for i in indices):
                raise ValidationError(f"Index out of range in group: {group}")
            if len(set(indices)) != len(indices):
                raise ValidationError(f"Duplicate indices in group: {group}")
            group_indices.append(indices)
        except ValueError:
            raise ValidationError(f"Invalid number format in group: {group}")
    
    return group_indices

def read_fasta(file: str) -> List[Tuple[str, str]]:
    """Read sequences from a FASTA file."""
    try:
        sequences = []
        for record in SeqIO.parse(file, "fasta"):
            sequences.append((str(record.description), str(record.seq)))
        return sequences
    except Exception as e:
        raise IOError(f"Error reading FASTA file {file}: {e}")

def validate_reference_index(reference_idx: int, num_sequences: int) -> bool:
    """Validate the reference sequence index."""
    return 0 <= reference_idx < num_sequences

def create_temporary_fasta(sequences: List[Tuple[str, str]], reference_idx: int) -> str:
    """Creates a temporary FASTA file with all sequences and returns the path."""
    try:
        fd, temp_path = tempfile.mkstemp(suffix='.fasta')
        with os.fdopen(fd, 'w') as f:
            # Reference sequence first
            ref_title, ref_seq = sequences[reference_idx]
            f.write(f">Query\n{ref_seq}\n")
            
            # Other sequences
            for i, (title, seq) in enumerate(sequences):
                if i != reference_idx:
                    # Use JX IDs similar to the example
                    f.write(f">JX{533400+i}\n{seq}\n")
        
        return temp_path
    except Exception as e:
        logging.error(f"Error creating temporary FASTA file: {e}")
        raise e

def locate_mafft():
    """
    Locate the MAFFT executable based on OS and installation method.
    According to MAFFT documentation for different operating systems.
    """
    system = platform.system()
    current_dir = os.getcwd()
    
    # First check if mafft is in PATH - common for all systems when installed properly
    mafft_in_path = shutil.which('mafft')
    if mafft_in_path:
        logging.info(f"Found MAFFT in PATH at {mafft_in_path}")
        return mafft_in_path
    
    # Check for local installation in the current directory
    if system == "Windows":
        # On Windows, look for mafft.bat in the project directory
        local_mafft = os.path.join(current_dir, "mafft.bat")
        if os.path.exists(local_mafft):
            logging.info(f"Found MAFFT batch file at {local_mafft}")
            return local_mafft
            
        # Check in mafft directory structure
        mafft_dir = os.path.join(current_dir, "mafft")
        if os.path.exists(mafft_dir):
            mafft_bat = os.path.join(mafft_dir, "mafft-win", "mafft.bat")
            if os.path.exists(mafft_bat):
                logging.info(f"Found MAFFT in project directory at {mafft_bat}")
                return mafft_bat
    
    elif system == "Darwin":  # macOS
        # On macOS, according to the documentation, MAFFT is run using mafft.bat
        mafft_dir = os.path.join(current_dir, "mafft")
        if os.path.exists(mafft_dir):
            # Per documentation, use mafft.bat in mafft-mac folder
            mafft_bat = os.path.join(mafft_dir, "mafft-mac", "mafft.bat")
            if os.path.exists(mafft_bat):
                # Check if the file is executable, if not, try to make it executable
                if not os.access(mafft_bat, os.X_OK):
                    try:
                        os.chmod(mafft_bat, 0o755)  # rwxr-xr-x
                        logging.info(f"Made {mafft_bat} executable with permissions 0o755")
                    except Exception as e:
                        logging.error(f"Failed to make {mafft_bat} executable: {e}")
                        # Continue anyway as we'll test execution later
                logging.info(f"Found MAFFT batch file at {mafft_bat}")
                return mafft_bat
    
    else:  # Linux
        # First look for mafft executable directly
        local_mafft = os.path.join(current_dir, "mafft")
        if os.path.exists(local_mafft) and os.access(local_mafft, os.X_OK):
            logging.info(f"Found local MAFFT executable at {local_mafft}")
            return local_mafft

        # Check in various known Linux locations based on the screenshot
        mafft_dir = os.path.join(current_dir, "mafft")
        if os.path.exists(mafft_dir):
            # Look in mafftdir/bin/mafft
            mafft_bin = os.path.join(mafft_dir, "mafftdir", "bin", "mafft")
            if os.path.exists(mafft_bin):
                if not os.access(mafft_bin, os.X_OK):
                    try:
                        os.chmod(mafft_bin, 0o755)
                        logging.info(f"Made {mafft_bin} executable with permissions 0o755")
                    except Exception as e:
                        logging.error(f"Failed to make MAFFT executable: {e}")
                logging.info(f"Found MAFFT in project directory at {mafft_bin}")
                return mafft_bin
            
            # Also check directly in the bin directory which is shown in your screenshot
            bin_dir = os.path.join(mafft_dir, "bin")
            if os.path.exists(bin_dir):
                mafft_in_bin = os.path.join(bin_dir, "mafft")
                if os.path.exists(mafft_in_bin):
                    if not os.access(mafft_in_bin, os.X_OK):
                        try:
                            os.chmod(mafft_in_bin, 0o755)
                            logging.info(f"Made {mafft_in_bin} executable with permissions 0o755")
                        except Exception as e:
                            logging.error(f"Failed to make MAFFT executable: {e}")
                    logging.info(f"Found MAFFT in bin directory at {mafft_in_bin}")
                    return mafft_in_bin
                    
            # Check if there's a mafft file in the root mafft directory
            root_mafft = os.path.join(mafft_dir, "mafft")
            if os.path.exists(root_mafft):
                if not os.access(root_mafft, os.X_OK):
                    try:
                        os.chmod(root_mafft, 0o755)
                        logging.info(f"Made {root_mafft} executable with permissions 0o755")
                    except Exception as e:
                        logging.error(f"Failed to make MAFFT executable: {e}")
                logging.info(f"Found MAFFT in mafft directory at {root_mafft}")
                return root_mafft
    
    logging.warning("MAFFT executable not found - MSA will not be possible")
    return None

def perform_alignment_mafft(sequences: List[Tuple[str, str]], reference_idx: int) -> Optional[List[Tuple[str, str]]]:
    """Perform sequence alignment using MAFFT, properly handling different OS configurations."""
    try:
        logging.info("Performing MSA using MAFFT")
        
        # Find MAFFT executable
        mafft_path = locate_mafft()
        if not mafft_path:
            logging.error("MAFFT executable not found - cannot perform alignment")
            # Try to run setup_mafft_from_local() from install_dependencies.py if available
            try:
                logging.info("Attempting to setup MAFFT from local installation files...")
                
                # Import the function if the module exists
                setup_mafft_success = False
                try:
                    import importlib.util
                    spec = importlib.util.spec_from_file_location("install_dependencies", 
                                                                 os.path.join(os.getcwd(), "install_dependencies.py"))
                    if spec and spec.loader:
                        install_deps = importlib.util.module_from_spec(spec)
                        spec.loader.exec_module(install_deps)
                        if hasattr(install_deps, "setup_mafft_from_local"):
                            setup_mafft_success = install_deps.setup_mafft_from_local()
                except Exception as setup_err:
                    logging.error(f"Error importing setup_mafft_from_local: {setup_err}")
                
                if setup_mafft_success:
                    logging.info("MAFFT setup completed successfully, trying to locate it again")
                    mafft_path = locate_mafft()
                    if not mafft_path:
                        logging.error("Still could not find MAFFT after setup")
                        return None
                else:
                    logging.error("Failed to setup MAFFT from local files")
                    return None
            except Exception as e:
                logging.error(f"Failed to setup MAFFT: {e}")
                return None
        
        # Create a temporary input FASTA file
        temp_fasta_path = create_temporary_fasta(sequences, reference_idx)
        output_file = tempfile.mktemp(suffix=".afa")
        
        system = platform.system()
        aligned_seqs = None  # Initialize as None for error cases
        
        # SPECIAL HANDLING FOR LINUX SYSTEMS
        if system == "Linux":
            logging.info("Using Linux-specific fallback method for running MAFFT")
            
            # Try to install MAFFT with apt-get if it's not executable
            if not os.access(mafft_path, os.X_OK):
                logging.warning(f"{mafft_path} is not executable")
                try:
                    # Try to install system MAFFT
                    logging.info("Attempting to install system MAFFT...")
                    subprocess.run(["sudo", "apt-get", "update"], check=False)
                    subprocess.run(["sudo", "apt-get", "install", "-y", "mafft"], check=False)
                    system_mafft = shutil.which("mafft")
                    if system_mafft:
                        logging.info(f"Successfully installed system MAFFT at {system_mafft}")
                        mafft_path = system_mafft
                except Exception as install_err:
                    logging.error(f"Failed to install system MAFFT: {install_err}")
            
            # 1. Ensure the MAFFT file is executable if it's our local one
            try:
                os.chmod(mafft_path, 0o755)
                logging.info(f"Set permissions on {mafft_path}")
            except Exception as perm_err:
                logging.warning(f"Error setting permissions (may be normal for system MAFFT): {perm_err}")
            
            # 2. Try direct execution first - often simpler
            logging.info(f"Attempting direct execution: {mafft_path} with {temp_fasta_path}")
            try:
                with open(output_file, 'w') as out_f:
                    proc = subprocess.run(
                        [mafft_path, "--quiet", "--auto", temp_fasta_path],
                        stdout=out_f,
                        stderr=subprocess.PIPE,
                        text=True
                    )
                
                if proc.returncode == 0:
                    logging.info("Direct execution succeeded")
                    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                        # SUCCESS CASE 1: Direct execution worked
                        aligned_seqs = []
                        for record in SeqIO.parse(output_file, "fasta"):
                            aligned_seqs.append((record.id, str(record.seq)))
                        logging.info(f"Successfully parsed {len(aligned_seqs)} aligned sequences")
                else:
                    logging.warning(f"Direct execution failed with code {proc.returncode}")
                    if proc.stderr:
                        logging.warning(f"Error: {proc.stderr}")
            
            except Exception as direct_err:
                logging.warning(f"Direct execution error: {direct_err}")
            
            # If direct execution failed, try with a wrapper script
            if not aligned_seqs:
                script_path = tempfile.mktemp(suffix=".sh")
                try:
                    with open(script_path, 'w') as script_file:
                        script_file.write("#!/bin/bash\n")
                        script_file.write(f"'{mafft_path}' --quiet --auto '{temp_fasta_path}' > '{output_file}'\n")
                    
                    # Make script executable    
                    os.chmod(script_path, 0o755)
                    logging.info(f"Created wrapper script at {script_path}")
                    
                    # Run script with bash
                    logging.info(f"Running: bash {script_path}")
                    proc = subprocess.run(
                        ["bash", script_path],
                        stderr=subprocess.PIPE,
                        text=True
                    )
                    
                    if proc.returncode == 0:
                        logging.info("Script execution succeeded")
                        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                            # SUCCESS CASE 2: Script execution worked
                            aligned_seqs = []
                            for record in SeqIO.parse(output_file, "fasta"):
                                aligned_seqs.append((record.id, str(record.seq)))
                            logging.info(f"Successfully parsed {len(aligned_seqs)} aligned sequences")
                    else:
                        logging.error(f"Script execution failed with code {proc.returncode}")
                        logging.error(f"Error: {proc.stderr}")
                
                except Exception as script_err:
                    logging.error(f"Script execution error: {script_err}")
                finally:
                    # Clean up script file
                    if os.path.exists(script_path):
                        os.remove(script_path)
            
            # As a last resort, try system-installed mafft
            if not aligned_seqs:
                system_mafft = shutil.which("mafft")
                if system_mafft and system_mafft != mafft_path:
                    logging.info(f"Trying system MAFFT at {system_mafft}")
                    try:
                        with open(output_file, 'w') as out_f:
                            proc = subprocess.run(
                                [system_mafft, "--quiet", "--auto", temp_fasta_path],
                                stdout=out_f,
                                stderr=subprocess.PIPE,
                                text=True
                            )
                        
                        if proc.returncode == 0:
                            logging.info("System MAFFT execution succeeded")
                            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                                # SUCCESS CASE 3: System MAFFT worked
                                aligned_seqs = []
                                for record in SeqIO.parse(output_file, "fasta"):
                                    aligned_seqs.append((record.id, str(record.seq)))
                                logging.info(f"Successfully parsed {len(aligned_seqs)} aligned sequences")
                        else:
                            logging.error(f"System MAFFT execution failed with code {proc.returncode}")
                            logging.error(f"Error: {proc.stderr}")
                    
                    except Exception as sys_err:
                        logging.error(f"System MAFFT error: {sys_err}")
                        
        # For macOS specifically use mafft.bat with shell=True
        elif system == "Darwin" and "mafft.bat" in mafft_path:
            logging.info(f"Running macOS-specific MAFFT command with mafft.bat")
            
            # Execute mafft.bat with file redirection, which requires shell=True
            cmd = f"{mafft_path} --quiet --auto {temp_fasta_path} > {output_file}"
            logging.info(f"Running command: {cmd}")
            
            try:
                # Use shell=True for macOS to handle redirection properly
                result = subprocess.run(
                    cmd,
                    shell=True,  # Shell required for redirection on macOS
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
                
                if result.stderr:
                    logging.info(f"MAFFT stderr: {result.stderr[:200]}...")
                
                # Parse the alignment results
                aligned_seqs = []
                for record in SeqIO.parse(output_file, "fasta"):
                    aligned_seqs.append((record.id, str(record.seq)))
                logging.info(f"Successfully parsed {len(aligned_seqs)} aligned sequences")
                    
            except subprocess.CalledProcessError as e:
                logging.error(f"MAFFT command failed with error code {e.returncode}")
                if e.stderr:
                    logging.error(f"MAFFT stderr: {e.stderr[:200]}...")
                aligned_seqs = None
                
        else:
            # For Windows and fallback for other systems
            mafft_cmd = [mafft_path, "--quiet", "--auto", temp_fasta_path]
            logging.info(f"Running standard command: {' '.join(mafft_cmd)}")
            
            try:
                with open(output_file, 'w') as outfile:
                    result = subprocess.run(
                        mafft_cmd,
                        stdout=outfile,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True
                    )
                
                if result.stderr:
                    logging.info(f"MAFFT stderr: {result.stderr[:200]}...")
                
                # Parse the alignment results
                aligned_seqs = []
                for record in SeqIO.parse(output_file, "fasta"):
                    aligned_seqs.append((record.id, str(record.seq)))
                logging.info(f"Successfully parsed {len(aligned_seqs)} aligned sequences")
                    
            except subprocess.CalledProcessError as e:
                logging.error(f"MAFFT command failed with error code {e.returncode}")
                if e.stderr:
                    logging.error(f"MAFFT stderr: {e.stderr[:200]}...")
                aligned_seqs = None
                
        # Final check across all methods
        if aligned_seqs is None or len(aligned_seqs) == 0:
            logging.error("Failed to generate valid alignment with any method")
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                try:
                    with open(output_file, 'r') as f:
                        content = f.read(500)
                        logging.error(f"Output file exists but parsing failed, content: {content}")
                except Exception:
                    pass
            return None
            
        # Clean up temporary files
        try:
            if os.path.exists(temp_fasta_path):
                os.remove(temp_fasta_path)
            if os.path.exists(output_file):
                os.remove(output_file)
        except Exception as cleanup_err:
            logging.warning(f"Error cleaning up temporary files: {cleanup_err}")
            
        return aligned_seqs
            
    except Exception as e:
        logging.error(f"Alignment error: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return None

def format_blast_like_alignment(aligned_sequences: List[Tuple[str, str]], block_size: int = 60) -> str:
    """Format alignment in BLAST-like output with dots for matching residues."""
    if not aligned_sequences:
        return ""
    
    # Get the reference sequence (first in list)
    ref_id, ref_seq = aligned_sequences[0]
    
    # Adjust sequences to the same length if needed
    max_len = max(len(seq) for _, seq in aligned_sequences)
    ref_seq = ref_seq.ljust(max_len, '-')
    
    # Format the alignment in blocks
    result = []
    sequence_positions = {ref_id: 0}
    
    for start_pos in range(0, len(ref_seq), block_size):
        end_pos = min(start_pos + block_size, len(ref_seq))
        
        # Count non-gap characters in the reference segment
        ref_segment = ref_seq[start_pos:end_pos]
        non_gap_count_ref = sum(1 for c in ref_segment if c != '-')
        
        # Calculate query positions
        if ref_id not in sequence_positions:
            sequence_positions[ref_id] = 0
        query_start_pos = sequence_positions[ref_id] + 1
        query_end_pos = query_start_pos + non_gap_count_ref - 1
        sequence_positions[ref_id] += non_gap_count_ref
        
        # Add the reference sequence block
        result.append(f"{ref_id} {query_start_pos:<3} {ref_segment} {query_end_pos}")
        
        # Add each subject sequence with dots for matches
        for seq_id, seq in aligned_sequences[1:]:
            if seq_id not in sequence_positions:
                sequence_positions[seq_id] = 0
                
            # Ensure sequence is long enough
            seq = seq.ljust(max_len, '-')
            
            # Extract segment
            subj_segment = seq[start_pos:end_pos]
            
            # Count non-gap characters for position calculation
            non_gap_count_subj = sum(1 for c in subj_segment if c != '-')
            
            subj_start_pos = sequence_positions[seq_id] + 1
            subj_end_pos = subj_start_pos + non_gap_count_subj - 1
            if non_gap_count_subj == 0:  # Handle case with only gaps
                subj_end_pos = subj_start_pos
            sequence_positions[seq_id] += non_gap_count_subj
            
            # Create formatted segment with appropriate notation
            formatted_segment = ""
            for i in range(len(subj_segment)):
                if i < len(ref_segment):
                    # Only show a dot if the nucleotides match and both are not gaps
                    if ref_segment[i] == subj_segment[i] and ref_segment[i] != '-' and subj_segment[i] != '-':
                        formatted_segment += "."  # Match - use a dot
                    else:
                        formatted_segment += subj_segment[i]  # Mismatch - show the actual nucleotide
                else:
                    formatted_segment += " "  # Padding if subject is shorter
            
            # Add the subject line
            result.append(f"{seq_id:<10} {subj_start_pos:<3} {formatted_segment} {subj_end_pos}")
        
        result.append("")  # Empty line between blocks
    
    return '\n'.join(result)

def extract_gene_name(sequence_titles: List[str], group_indices: List[int] = None) -> str:
    """Extract gene name from sequence titles and file groups for the consensus sequence header."""
    # Look for common gene names in sequence titles
    gene_patterns = [
        r'(\w+(?:\s*\d*)?\s*(?:ribosomal\s*RNA|rRNA))',
        r'((?:elongation\s*factor|tef)\s*\d*\s*(?:alpha)?)',
        r'((?:beta|β)[\-\s]*tubulin\s*\d*)',
        r'((?:cytochrome\s*oxidase|cox)\s*(?:subunit\s*)?[IV\d]+)',
        r'(internal\s*transcribed\s*spacer\s*\d*)'
    ]
    
    potential_genes = []
    for title in sequence_titles:
        for pattern in gene_patterns:
            matches = re.findall(pattern, title, re.IGNORECASE)
            if matches:
                for match in matches:
                    if isinstance(match, tuple):  # Some regex patterns return tuples
                        potential_genes.extend([m for m in match if m])
                    else:
                        potential_genes.append(match)
    
    # Group similar gene names (normalize case, remove spaces)
    normalized_genes = {}
    for gene in potential_genes:
        key = gene.lower().replace(" ", "")
        if key in normalized_genes:
            normalized_genes[key].append(gene)
        else:
            normalized_genes[key] = [gene]
    
    # Get the most common gene name from each group
    most_common_genes = []
    for gene_group in normalized_genes.values():
        counter = collections.Counter(gene_group)
        most_common = counter.most_common(1)[0][0]
        most_common_genes.append(most_common)
    
    if most_common_genes:
        # Join multiple gene names
        return "_".join(most_common_genes)
    
    # Fallback: look for common gene abbreviations
    common_genes = {
        'tef1': 'elongation_factor_1_alpha',
        'TEF1': 'elongation_factor_1_alpha',
        'ITS': 'internal_transcribed_spacer',
        'COX1': 'cytochrome_oxidase_subunit_I',
        'cox1': 'cytochrome_oxidase_subunit_I',
        'rRNA': 'ribosomal_RNA',
        'rDNA': 'ribosomal_DNA',
        'beta-tubulin': 'beta_tubulin'
    }
    
    for title in sequence_titles:
        for abbr, full_name in common_genes.items():
            if abbr in title:
                return full_name
    
    # Last resort - use "consensus" if nothing else works
    return "consensus"

def generate_consensus_sequence(aligned_sequences: List[Tuple[str, str]], threshold_percent: float) -> str:
    """
    Generate a consensus sequence from aligned sequences based on the specified threshold percentage.
    
    Args:
        aligned_sequences: List of tuples containing (sequence_id, aligned_sequence)
        threshold_percent: Percentage threshold (0-100) for including a nucleotide in the consensus
    
    Returns:
        String representation of the consensus sequence
    """
    if not aligned_sequences:
        return ""
    
    # Get the length of the alignment (all sequences should have the same length after alignment)
    alignment_length = max(len(seq) for _, seq in aligned_sequences)
    
    # Initialize consensus sequence with empty characters
    consensus = []
    
    # Count nucleotides at each position
    for position in range(alignment_length):
        # Get nucleotides at this position, excluding gaps
        nucleotides = [seq[position] for _, seq in aligned_sequences if position < len(seq) and seq[position] != '-']
        
        # If no non-gap nucleotides at this position, add a gap to the consensus
        if not nucleotides:
            consensus.append('-')
            continue
            
        # Count occurrences of each nucleotide
        nucleotide_counts = collections.Counter(nucleotides)
        
        # Get the most common nucleotide and its count
        most_common = nucleotide_counts.most_common(1)[0]
        most_common_nucleotide, most_common_count = most_common
        
        # Calculate percentage of the most common nucleotide
        percent = (most_common_count / len(nucleotides)) * 100
        
        # Add the nucleotide to the consensus if it meets the threshold, otherwise add 'N'
        if percent >= threshold_percent:
            consensus.append(most_common_nucleotide)
        else:
            consensus.append('N')
    
    return ''.join(consensus)

def save_consensus_to_fasta(consensus_sequence: str, gene_name: str, group_indices: List[int], output_dir: str) -> str:
    """
    Save the consensus sequence to a FASTA file with an appropriate header.
    """
    # Remove any gaps from the final consensus sequence
    # consensus_sequence = consensus_sequence.replace('-', '')
    
    # Generate filename
    indices_str = '_'.join(map(str, group_indices))
    consensus_file = os.path.join(output_dir, f"CONSENSUS_{indices_str}.fasta")
    
    # Create the FASTA file with more descriptive header
    with open(consensus_file, 'w') as f:
        f.write(f">Puccinia_recondita_{gene_name}_consensus_group_{indices_str}\n")
        
        # Write sequence in blocks of 60 characters
        for i in range(0, len(consensus_sequence), 60):
            f.write(f"{consensus_sequence[i:i+60]}\n")
    
    return consensus_file

def process_sequence_group(group: SequenceGroup, output_dir: str, interactive_mode: bool = False) -> Optional[str]:
    """Process a group of sequences to generate a multiple sequence alignment using MAFFT."""
    try:
        logging.info(f"\nProcessing group with indices: {group.indices}")
        
        filtered_sequences = [(title, seq) for title, seq in group.sequences if len(seq) < 3000]
        if not filtered_sequences:
            logging.info("No sequences shorter than 3000 nucleotides in this group. Skipping MSA analysis.")
            return None
        
        logging.info(f"Using {len(filtered_sequences)} sequences for MSA analysis")
        
        # Get reference sequence
        ref_idx = group.reference_idx or 0
        ref_title, ref_sequence = filtered_sequences[ref_idx]
        logging.info(f"Using sequence '{ref_title}' as reference")
        
        # Perform the alignment using MAFFT - no fallback methods anymore
        logging.info("Performing alignment with MAFFT...")
        aligned_sequences = perform_alignment_mafft(filtered_sequences, ref_idx)
        
        if not aligned_sequences:
            logging.error("MAFFT alignment failed. Cannot continue without MAFFT.")
            return None
        
        # Format results in BLAST-like format
        blast_like_alignment = format_blast_like_alignment(aligned_sequences)
        
        # Save results
        mafft_output_file = os.path.join(output_dir, f"MAFFT_{'_'.join(map(str, group.indices))}.txt")
        blast_output_file = os.path.join(output_dir, f"BLAST_{'_'.join(map(str, group.indices))}.txt")
        
        # Save MSA alignment to file
        with open(mafft_output_file, 'w') as f:
            for seq_id, sequence in aligned_sequences:
                f.write(f">{seq_id}\n{sequence}\n")
        
        # Save BLAST-like alignment to file
        with open(blast_output_file, 'w') as f:
            f.write(blast_like_alignment)
        
        logging.info(f"MSA analysis saved to: {mafft_output_file}")
        logging.info(f"BLAST-like alignment saved to: {blast_output_file}")
        
        # Verify the alignment quality
        ref_seq = aligned_sequences[0][1]
        match_percentages = []
        
        for seq_id, seq in aligned_sequences[1:]:
            matches = sum(1 for i in range(min(len(ref_seq), len(seq))) 
                         if i < len(ref_seq) and i < len(seq) and
                         ref_seq[i] == seq[i] and ref_seq[i] != '-' and seq[i] != '-')
            total = sum(1 for i in range(min(len(ref_seq), len(seq))) 
                       if i < len(ref_seq) and i < len(seq) and
                       ref_seq[i] != '-' and seq[i] != '-')
            
            if total > 0:
                match_percentage = (matches / total) * 100
                match_percentages.append(match_percentage)
                logging.info(f"Sequence {seq_id} matches reference at {match_percentage:.1f}% of positions")
        
        if match_percentages:
            avg_match = sum(match_percentages) / len(match_percentages)
            logging.info(f"Average match percentage: {avg_match:.1f}%")
            
            if avg_match < 70:
                logging.warning("Average sequence similarity is below 70%. Alignment may be suboptimal.")
        
        # Ask user for consensus threshold
        default_threshold = 80
        threshold_prompt = f"\nEnter consensus threshold percentage for group {group.indices} (default {default_threshold}%): "
        
        if interactive_mode:
            threshold_input = get_console_input(threshold_prompt, str(default_threshold))
        else:
            threshold_input = input(threshold_prompt)
            if not threshold_input.strip():
                threshold_input = str(default_threshold)
        
        try:
            threshold = float(threshold_input.strip())
            if threshold < 0 or threshold > 100:
                logging.warning(f"Invalid threshold value: {threshold}. Using default of {default_threshold}%")
                threshold = default_threshold
        except ValueError:
            logging.warning(f"Invalid threshold format: {threshold_input}. Using default of {default_threshold}%")
            threshold = default_threshold
            
        # Generate consensus sequence
        logging.info(f"Generating consensus sequence with {threshold}% threshold...")
        consensus_sequence = generate_consensus_sequence(aligned_sequences, threshold)
        
        if consensus_sequence:
            # Extract gene name from sequence titles
            sequence_titles = [title for title, _ in filtered_sequences]
            gene_name = extract_gene_name(sequence_titles, group.indices)
            
            # Save consensus to FASTA file
            consensus_file = save_consensus_to_fasta(consensus_sequence, gene_name, group.indices, output_dir)
            logging.info(f"Consensus sequence saved to: {consensus_file}")
            
            # Count nucleotides in consensus
            non_gap_count = sum(1 for nt in consensus_sequence if nt != '-')
            n_count = consensus_sequence.count('N')
            resolved_percent = ((non_gap_count - n_count) / non_gap_count) * 100 if non_gap_count > 0 else 0
            
            logging.info(f"Consensus sequence statistics:")
            logging.info(f"  - Total length: {len(consensus_sequence)}")
            logging.info(f"  - Non-gap nucleotides: {non_gap_count}")
            logging.info(f"  - Ambiguous (N) positions: {n_count}")
            logging.info(f"  - Resolved positions: {non_gap_count - n_count} ({resolved_percent:.1f}%)")
        else:
            logging.warning("Failed to generate consensus sequence")
            
        return blast_like_alignment
    except Exception as e:
        logging.error(f"Error processing sequence group: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return None

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
    
    logging.info(f"Waiting for input: {prompt}")
    
    # Write the prompt to the request file
    with open(input_request_file, 'w', encoding='utf-8') as f:
        f.write(prompt)
    
    # Wait for response file to appear
    timeout = 600  # 10 minutes timeout
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

def cleanup_previous_run():
    """Remove log files and other files/directories from previous run of module 2."""
    try:
        # Remove log file if it exists
        log_file = '2_MSA_Alignment.log'
        if os.path.exists(log_file):
            os.remove(log_file)
            print(f"Removed previous log file: {log_file}")
        
        # Remove output directory and its contents if it exists
        output_dir = os.path.join(os.getcwd(), '2_')
        if os.path.exists(output_dir):
            for root, dirs, files in os.walk(output_dir, topdown=False):
                for file in files:
                    file_path = os.path.join(root, file)
                    os.remove(file_path)
                    print(f"Removed file: {file_path}")
                for dir in dirs:
                    dir_path = os.path.join(root, dir)
                    os.rmdir(dir_path)
                    print(f"Removed directory: {dir_path}")
            
            # Remove the main output directory after its contents are removed
            if os.path.exists(output_dir):  # Check again in case it was already deleted
                os.rmdir(output_dir)
                print(f"Removed output directory: {output_dir}")
        
        # Remove input request and response files if they exist
        for file in ['input_request.txt', 'input_response.txt']:
            if os.path.exists(file):
                os.remove(file)
                print(f"Removed file: {file}")
                
    except Exception as e:
        print(f"Error during cleanup: {e}")
        # Continue execution even if cleanup fails

def main():
    """Main function to handle sequence alignment."""
    try:
        # Clean up files from previous run
        cleanup_previous_run()
        
        # Check if running in interactive mode
        interactive_mode = "--interactive" in sys.argv
        
        # Check for MAFFT and report
        mafft_path = locate_mafft()
        if mafft_path:
            # Check if executable is properly accessible based on OS
            if platform.system() == "Windows":
                logging.info(f"Found MAFFT at {mafft_path}")
            elif os.access(mafft_path, os.X_OK):
                logging.info(f"Found MAFFT at {mafft_path} (executable: Yes)")
            else:
                logging.warning(f"Found MAFFT at {mafft_path} but it is not executable")
                try:
                    os.chmod(mafft_path, 0o755)
                    logging.info(f"Made {mafft_path} executable with permissions 0o755")
                except Exception as e:
                    logging.error(f"Failed to make MAFFT executable: {e}")
                    logging.warning("You may need to manually set the executable permissions")
        else:
            logging.error("MAFFT not found. Please make sure it's installed properly.")
            logging.error("You can install it by running install_dependencies.py")
            sys.exit(1)
        
        input_dir = os.path.join(os.getcwd(), '1_')
        output_dir = os.path.join(os.getcwd(), '2_')
        os.makedirs(output_dir, exist_ok=True)

        fasta_files = find_fasta_files(input_dir)
        if not fasta_files:
            raise ValidationError("No FASTA files found in the input directory")

        logging.info("\nAvailable FASTA files:")
        for idx, file in enumerate(fasta_files):
            logging.info(f"{idx+1}: {os.path.basename(file)}")

        # Get group indices from user
        if interactive_mode:
            groups_input = get_console_input(
                "\nEnter groups of file indices for MSA analysis, separated by semicolons (e.g., 0,1,3;2,4,5): "
            )
        else:
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
                logging.info(f"{idx+1}: {title} (length: {len(seq)})")

            # Get reference sequence index
            if interactive_mode:
                ref_input = get_console_input(
                    f"\nChoose reference sequence (1-{len(sorted_sequences)}): "
                )
            else:
                ref_input = input(f"\nChoose reference sequence (0-{len(sorted_sequences)-1}): ")
                
            try:
                ref_idx = int(ref_input) - 1
                if not validate_reference_index(ref_idx, len(sorted_sequences)):
                    logging.error("Invalid reference index")
                    continue
            except ValueError:
                logging.error("Please enter a valid number")
                continue

            group = SequenceGroup(group_indices, sorted_sequences, ref_idx)
            result = process_sequence_group(group, output_dir, interactive_mode)
            
            if result:
                msa_results.append((group_indices, result))

        if msa_results:
            logging.info("\nMSA analysis, BLAST-like alignment, and consensus generation completed successfully")
            for indices, result in msa_results:
                logging.info(f"\nResults for group {indices}:")
                # Show just a preview of the result (first few lines)
                preview = result.split('\n')[:10]
                logging.info('\n'.join(preview) + '\n...')
        else:
            logging.warning("\nNo successful MSA analyses to report")
            
    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        import traceback
        logging.error(traceback.format_exc())

if __name__ == '__main__':
    main()
