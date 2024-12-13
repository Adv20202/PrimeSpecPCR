import os
import re
import csv
import time 
import logging
import shutil
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Any, Tuple
from Bio import Entrez, SeqIO
from tqdm import tqdm
import sys
import gc

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="resource_tracker")

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('8_Species_Specificity_Test.log', encoding='utf-8')
    ]
)

class ValidationError(Exception):
    pass

@dataclass
class PrimerInfo:
    sequence: str
    tm: float
    gc_percent: float
    self_any_th: float
    self_end_th: float
    hairpin_th: float

@dataclass
class PrimerSet:
    set_number: int
    product_size: str
    left_primer: PrimerInfo
    right_primer: PrimerInfo
    probe: PrimerInfo

@dataclass
class SpeciesCounts:
    full_match: Dict[str, int]
    partial_match: Dict[str, int]
    no_match: Dict[str, int]
    species_taxid_mapping: Dict[str, Optional[int]]

@dataclass
class AnalysisResult:
    primer_set: PrimerSet
    ncbi_records_found: int
    species_counts: SpeciesCounts
    csv_file: str

class SpeciesAnalyzer:
    def __init__(self, email: str, api_key: str):
        self.email = email
        self.api_key = api_key
        self.taxid_species_cache: Dict[int, str] = {}
        self.accession_species_cache: Dict[str, Dict[str, Any]] = {}
        Entrez.email = email
        Entrez.api_key = api_key

    def get_species_name(self, accession: str, verbose: bool = False) -> Tuple[str, Optional[int]]:
        max_retries = 3
        retry_count = 0
        
        if accession in self.accession_species_cache:
            cached = self.accession_species_cache[accession]
            return cached['species'], cached['taxid']

        while retry_count < max_retries:
            try:
                if verbose:
                    logging.info(f"Fetching species name for accession {accession}...")
                handle = Entrez.esummary(db="nucleotide", id=accession, retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                
                if not records:
                    retry_count += 1
                    time.sleep(0.1)
                    continue
                    
                taxid = int(records[0]['TaxId'])
                species_name = self.get_species_name_from_taxid(taxid, verbose)
                
                if species_name != 'Unknown':
                    self.accession_species_cache[accession] = {
                        'species': species_name,
                        'taxid': taxid
                    }
                    time.sleep(0.1)
                    return species_name, taxid
                
                retry_count += 1
                time.sleep(0.1)
                
            except Exception as e:
                if verbose:
                    logging.error(f"Error fetching species name for accession {accession}: {e}")
                retry_count += 1
                time.sleep(0.1)
                
        try:
            handle = Entrez.esearch(db="taxonomy", term=accession)
            record = Entrez.read(handle)
            handle.close()
            
            if record["Count"] != "0":
                taxid = int(record["IdList"][0])
                species_name = self.get_species_name_from_taxid(taxid, verbose)
                
                self.accession_species_cache[accession] = {
                    'species': species_name,
                    'taxid': taxid
                }
                return species_name, taxid
                
        except Exception:
            pass
            
        return self.get_fallback_species_name(accession)

    def get_fallback_species_name(self, accession: str) -> Tuple[str, Optional[int]]:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            for feature in record.features:
                if feature.type == "source":
                    if "organism" in feature.qualifiers:
                        species_name = feature.qualifiers["organism"][0]
                        try:
                            handle = Entrez.esearch(db="taxonomy", term=species_name)
                            result = Entrez.read(handle)
                            handle.close()
                            if result["Count"] != "0":
                                taxid = int(result["IdList"][0])
                                return species_name, taxid
                        except:
                            pass
                        return species_name, None
            
        except Exception:
            pass
            
        return "Unidentified species", None

    def get_species_name_from_taxid(self, taxid: int, verbose: bool = False) -> str:
        max_retries = 3
        retry_count = 0
        
        if taxid in self.taxid_species_cache:
            return self.taxid_species_cache[taxid]

        while retry_count < max_retries:
            try:
                if verbose:
                    logging.info(f"Retrieving species name for TaxID {taxid}...")
                handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                
                if records:
                    species_name = records[0]['ScientificName']
                    self.taxid_species_cache[taxid] = species_name
                    time.sleep(0.1)
                    return species_name
                    
                retry_count += 1
                time.sleep(0.1)
                
            except Exception as e:
                if verbose:
                    logging.error(f"Error fetching species name for TaxID {taxid}: {e}")
                retry_count += 1
                time.sleep(0.1)
                
        return "Unidentified species"

class DataProcessor:
    def __init__(self, analyzer: SpeciesAnalyzer):
        self.analyzer = analyzer

    def read_primer_set(self, row: Dict[str, str]) -> PrimerSet:
        def parse_float(value: str) -> float:
            try:
                return round(float(value), 2)
            except (ValueError, TypeError):
                return 0.0

        return PrimerSet(
            set_number=int(row.get('PRIMER_SET_NUMBER', '0')),
            product_size=row.get('PRODUCT_SIZE', ''),
            left_primer=PrimerInfo(
                sequence=row.get('PRIMER_LEFT_SEQUENCE', '').strip().upper(),
                tm=parse_float(row.get('PRIMER_LEFT_TM', '0')),
                gc_percent=parse_float(row.get('PRIMER_LEFT_GC%', '0')),
                self_any_th=parse_float(row.get('PRIMER_LEFT_SELF_ANY_TH', '0')),
                self_end_th=parse_float(row.get('PRIMER_LEFT_SELF_END_TH', '0')),
                hairpin_th=parse_float(row.get('PRIMER_LEFT_HAIRPIN_TH', '0'))
            ),
            right_primer=PrimerInfo(
                sequence=row.get('PRIMER_RIGHT_SEQUENCE', '').strip().upper(),
                tm=parse_float(row.get('PRIMER_RIGHT_TM', '0')),
                gc_percent=parse_float(row.get('PRIMER_RIGHT_GC%', '0')),
                self_any_th=parse_float(row.get('PRIMER_RIGHT_SELF_ANY_TH', '0')),
                self_end_th=parse_float(row.get('PRIMER_RIGHT_SELF_END_TH', '0')),
                hairpin_th=parse_float(row.get('PRIMER_RIGHT_HAIRPIN_TH', '0'))
            ),
            probe=PrimerInfo(
                sequence=row.get('PROBE_SEQUENCE', '').strip().upper(),
                tm=parse_float(row.get('PROBE_TM', '0')),
                gc_percent=parse_float(row.get('PROBE_GC%', '0')),
                self_any_th=parse_float(row.get('PROBE_SELF_ANY_TH', '0')),
                self_end_th=parse_float(row.get('PROBE_SELF_END_TH', '0')),
                hairpin_th=parse_float(row.get('PROBE_HAIRPIN_TH', '0'))
            )
        )

    def categorize_scoring(self, scoring: str) -> str:
        if scoring == 'A_A_A':
            return 'A'
        elif 'C' in scoring:
            return 'C'
        else:
            return 'B'

    def analyze_sequence_record(self, record: Any) -> Tuple[str, str, Optional[int]]:
        header = str(record.description)
        accession_match = re.match(r'^>?([\w\.]+)', header)
        accession = accession_match.group(1).split('.')[0] if accession_match else "Unknown"
        
        scoring_match = re.search(r'__(\S+)$', header)
        scoring = scoring_match.group(1) if scoring_match else ''
            
        species_name, taxid = self.analyzer.get_species_name(accession)
        return species_name, self.categorize_scoring(scoring), taxid

    def process_fasta_file(self, fasta_path: str, taxid_list: List[int]) -> SpeciesCounts:
        species_counts = SpeciesCounts(
            full_match={},
            partial_match={},
            no_match={},
            species_taxid_mapping={}
        )

        if not os.path.exists(fasta_path):
            return species_counts

        try:
            sequences = list(SeqIO.parse(fasta_path, 'fasta'))
            total = len(sequences)
            terminal_width = shutil.get_terminal_size().columns
                
            for idx, record in enumerate(sequences, 1):
                sys.stdout.write('\r' + ' ' * terminal_width + '\r')
                progress_msg = f"Analyzing sequence {idx}/{total}: {record.id}"
                sys.stdout.write(progress_msg)
                sys.stdout.flush()
                
                species_name, category, taxid = self.analyze_sequence_record(record)
                
                if taxid is not None:
                    species_counts.species_taxid_mapping[species_name] = taxid
                
                if category == 'A':
                    species_counts.full_match[species_name] = species_counts.full_match.get(species_name, 0) + 1
                elif category == 'B':
                    species_counts.partial_match[species_name] = species_counts.partial_match.get(species_name, 0) + 1
                else:
                    species_counts.no_match[species_name] = species_counts.no_match.get(species_name, 0) + 1

            sys.stdout.write('\r' + ' ' * terminal_width + '\r')
            sys.stdout.flush()

        except Exception as e:
            logging.error(f"Error processing FASTA file {fasta_path}: {e}")

        return species_counts

class DashboardGenerator:
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir

    def sort_primer_sets(self, results: List[AnalysisResult]) -> List[AnalysisResult]:
        def get_sort_key(result: AnalysisResult) -> Tuple[int, int, int]:
            full_match_max = max(result.species_counts.full_match.values()) if result.species_counts.full_match else 0
            partial_match_max = max(result.species_counts.partial_match.values()) if result.species_counts.partial_match else 0
            no_match_max = max(result.species_counts.no_match.values()) if result.species_counts.no_match else 0
            return (-full_match_max, -partial_match_max, -no_match_max)
            
        return sorted(results, key=get_sort_key)

    def create_dashboard(self, results: List[AnalysisResult], taxid_list: List[int]) -> None:
        results_by_file: Dict[str, List[AnalysisResult]] = {}
        for result in results:
            if result.csv_file not in results_by_file:
                results_by_file[result.csv_file] = []
            results_by_file[result.csv_file].append(result)

        for file_results in results_by_file.values():
            file_results.sort(key=lambda x: (
                -max(x.species_counts.full_match.values()) if x.species_counts.full_match else 0,
                -max(x.species_counts.partial_match.values()) if x.species_counts.partial_match else 0
            ))

        html_content = """<!DOCTYPE html>
<html>
<head>
    <title>Species Specificity Analysis</title>
    <style>
        body { font-family: Arial, sans-serif; padding: 20px; line-height: 1.6; }
        .primer-set { border: 1px solid #ccc; margin: 20px 0; padding: 15px; }
        .file-section { margin-bottom: 40px; }
        .file-header { background: #f0f0f0; padding: 10px; margin-bottom: 20px; }
        .match-category { margin: 10px 0; }
        .target-species { color: green; }
        .species-list { margin-left: 20px; }
        .primer-info { font-family: monospace; background: #f5f5f5; padding: 10px; margin: 5px 0; }
    </style>
</head>
<body>"""

        for csv_file, file_results in results_by_file.items():
            html_content += f"""
    <div class='file-section'>
        <div class='file-header'>
            <h2>Results for {csv_file}</h2>
        </div>"""

            for result in file_results:
                html_content += f"""
        <div class='primer-set'>
            <div class='primer-info'>
                Primer Set {result.primer_set.set_number} | 
                Left: {result.primer_set.left_primer.sequence} | 
                Right: {result.primer_set.right_primer.sequence} | 
                Probe: {result.primer_set.probe.sequence} | 
                Product Size: {result.primer_set.product_size}
            </div>
            
            <div class='match-category'>
                <strong>Full match between primers and NCBI records:</strong><br/>
                <div class='species-list'>"""

                for species, count in sorted(result.species_counts.full_match.items(), 
                                          key=lambda x: (-x[1], x[0])):
                    species_marker = "*" if result.species_counts.species_taxid_mapping.get(species) in taxid_list else ""
                    html_content += f"{species}{species_marker}: {count} records; "

                html_content += """</div>
            </div>
            
            <div class='match-category'>
                <strong>Slight nucleotide difference (high probability of amplification):</strong><br/>
                <div class='species-list'>"""

                for species, count in sorted(result.species_counts.partial_match.items(), 
                                          key=lambda x: (-x[1], x[0])):
                    species_marker = "*" if result.species_counts.species_taxid_mapping.get(species) in taxid_list else ""
                    html_content += f"{species}{species_marker}: {count} records; "

                html_content += """</div>
            </div>
            
            <div class='match-category'>
                <strong>Significant differences (low probability of amplification):</strong><br/>
                <div class='species-list'>"""

                for species, count in sorted(result.species_counts.no_match.items(), 
                                          key=lambda x: (-x[1], x[0])):
                    species_marker = "*" if result.species_counts.species_taxid_mapping.get(species) in taxid_list else ""
                    html_content += f"{species}{species_marker}: {count} records; "

                html_content += """</div>
            </div>
            <br/>
        </div>"""

            html_content += """
    </div>"""

        html_content += """
</body>
</html>"""

        output_path = self.output_dir / 'species_specificity_dashboard.html'
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        import webbrowser
        webbrowser.open(str(output_path))

def main():
    try:
        current_dir = Path.cwd()
        input_dir = current_dir / '5_'
        fasta_dir = current_dir / '7_'
        output_dir = current_dir / '8_'
        output_dir.mkdir(exist_ok=True)

        try:
            email, api_key = read_user_data()
        except Exception as e:
            raise ValidationError(f"Error reading user data: {e}")

        analyzer = SpeciesAnalyzer(email, api_key)
        processor = DataProcessor(analyzer)
        dashboard = DashboardGenerator(output_dir)

        taxid_list = get_taxids_from_user(analyzer.taxid_species_cache)

        csv_files = find_csv_files(str(input_dir))
        if not csv_files:
            raise ValidationError('No CSV files starting with "MAFFT_" found in directory "5_"')

        all_results: List[AnalysisResult] = []
        
        for csv_index, csv_file in enumerate(csv_files, 1):
            logging.info(f"\nProcessing CSV file {csv_index}/{len(csv_files)}: {csv_file}")
            
            try:
                results = process_csv_file(
                    csv_file=csv_file,
                    input_dir=input_dir,
                    fasta_dir=fasta_dir,
                    processor=processor,
                    taxid_list=taxid_list
                )
                all_results.extend(results)
                
            except Exception as e:
                logging.error(f"Error processing CSV file {csv_file}: {e}")
                continue

        if all_results:
            dashboard.create_dashboard(all_results, taxid_list)
        else:
            logging.warning("No results to display in dashboard")

    except ValidationError as e:
        logging.error(f"Validation error: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

def read_user_data() -> Tuple[str, str]:
    user_data_dir = Path('user_data')
    email_file = user_data_dir / 'email.txt'
    api_key_file = user_data_dir / 'api_key.txt'

    if not email_file.exists() or not api_key_file.exists():
        raise ValidationError("Missing user data files")

    email = email_file.read_text().strip()
    api_key = api_key_file.read_text().strip()

    if not email or not api_key:
        raise ValidationError("Empty email or API key")

    return email, api_key

def get_taxids_from_user(taxid_species_cache: Dict[int, str]) -> List[int]:
    print("\nPlease enter the TaxID(s) for which the primers should be specific.")
    taxid_input = input("Enter TaxID(s) separated by commas (e.g., 12345,67890): ").strip()
    
    taxid_list = []
    for taxid_str in taxid_input.split(','):
        taxid_str = taxid_str.strip()
        if taxid_str.isdigit():
            taxid = int(taxid_str)
            taxid_list.append(taxid)
            if taxid not in taxid_species_cache:
                taxid_species_cache[taxid] = 'Unknown'
        else:
            logging.warning(f'Invalid TaxID: {taxid_str}')

    if not taxid_list:
        raise ValidationError('No valid TaxIDs provided')

    return taxid_list

def find_csv_files(directory: str) -> List[str]:
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist")
    
    pattern = re.compile(r'^MAFFT_\d+.*\.csv$')
    files = [f for f in os.listdir(directory) if pattern.match(f)]
    
    def extract_number(filename: str) -> int:
        match = re.search(r'MAFFT_(\d+)', filename)
        return int(match.group(1)) if match else 0
    
    return sorted(files, key=extract_number)

def process_csv_file(
    csv_file: str,
    input_dir: Path,
    fasta_dir: Path,
    processor: DataProcessor,
    taxid_list: List[int]
) -> List[AnalysisResult]:
    results = []
    pbar = None
    csv_path = input_dir / csv_file

    try:
        with open(csv_path, 'r', newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            total_rows = len(rows)
            
            terminal_width = shutil.get_terminal_size().columns - 20
            
            pbar = tqdm(total=total_rows, 
                       desc=f"Processing {csv_file}",
                       unit="sets",
                       position=0,
                       leave=True,
                       ncols=min(100, terminal_width))
            
            for row in rows:
                primer_set = processor.read_primer_set(row)
                if not all([primer_set.left_primer.sequence, 
                          primer_set.right_primer.sequence, 
                          primer_set.probe.sequence]):
                    pbar.update(1)
                    continue

                fasta_filename = f"primer_set_{primer_set.set_number}.fasta"
                fasta_path = fasta_dir / fasta_filename

                if not fasta_path.exists():
                    logging.warning(f"FASTA file not found: {fasta_filename}")
                    pbar.update(1)
                    continue

                try:
                    num_records = sum(1 for _ in SeqIO.parse(fasta_path, 'fasta'))
                    species_counts = processor.process_fasta_file(str(fasta_path), taxid_list)
                    
                    result = AnalysisResult(
                        primer_set=primer_set,
                        ncbi_records_found=num_records,
                        species_counts=species_counts,
                        csv_file=csv_file
                    )
                    
                    results.append(result)
                    
                except Exception as e:
                    logging.error(f"Error processing FASTA file {fasta_filename}: {e}")
                finally:
                    pbar.update(1)

    except Exception as e:
        logging.error(f"Error processing CSV file {csv_file}: {e}")
        raise
    finally:
        if pbar is not None:
            pbar.close()
            pbar.clear()
            del pbar
        gc.collect()
    
    return results

if __name__ == '__main__':
    main()
