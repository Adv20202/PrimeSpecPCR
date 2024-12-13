# PrimeSpecPCR: species-specific primer design toolkit

## Table of Contents

- [Introduction](#introduction)
- [Key Features](#key-features)
- [Installation](#installation)
- [Workflow](#workflow)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [License](#license)
- [Disclaimer](#disclaimer)

## Introduction

This repository contains a suite of Python scripts designed for the analysis and processing of genetic data. The scripts facilitate key tasks such as retrieving genetic sequences, performing multiple sequence alignments, generating consensus sequences, designing primers, and evaluating primer specificity. Each script is modular, enabling focused execution of specific stages in the genetic data workflow.

The repository is tailored for research workflows that require rigorous sequence validation, alignment accuracy, and primer design reliability. The primary aim is to streamline computational tasks associated with genetic sequence analysis and enhance the reproducibility of bioinformatics pipelines.

## Key Features

- **Genetic Sequence Retrieval**: Automates the extraction of nucleotide sequences from NCBI databases using specified TaxIDs and query terms.

- **Multiple Sequence Alignment (MSA)**: Implements MAFFT for aligning sequences with customizable input and output configurations, ensuring consistent results for downstream analyses.

- **Consensus Sequence Generation**: Produces high-confidence consensus sequences from aligned data, with adjustable thresholds for nucleotide agreement.

- **Primer Design**: Utilizes the Primer3 library for efficient primer generation, tailored to user-specified parameters such as melting temperature, GC content, and amplicon size.

- **Primer Specificity Testing**: Evaluates primer sequences against genetic databases to identify potential cross-reactivity and optimize target specificity.

- **Species-Specificity Analysis**: Integrates BLAST searches and annotation tools to assess the taxonomic specificity of designed primers, aiding in the identification of species-level targets.

- **Comprehensive Logging**: Provides detailed logging at every step of the workflow, ensuring traceability and facilitating troubleshooting.

- **Flexible Modular Design**: Each script operates independently, allowing users to customize workflows according to specific research needs.

## Installation

To set up the environment and install the necessary dependencies, follow these steps:

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/Adv20202/PrimeSpecPCR.git
   cd PrimeSpecPCR
   ```

2. **Create a Virtual Environment** (optional but recommended):

   ```bash
   python3 -m venv env
   source env/bin/activate  # On Windows: .\env\Scripts\activate
   ```

3. **Install Dependencies**:

   Install the required Python libraries as specified in the `requirements.txt` file:

   ```bash
   pip install -r requirements.txt
   ```

4. **Set Up MAFFT**:

   Ensure that MAFFT is installed and added to your system's PATH:
   - Download and install MAFFT from [MAFFT Official Site](https://mafft.cbrc.jp/alignment/software/).

5. **Obtain an NCBI API Key**:

   To access NCBI services, you must have an API key. Visit the [NCBI API Key Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25497/) for guidance on obtaining one. The scripts will prompt you for this key when needed.

6. **Run a Test Script**:

   Verify the installation by running one of the provided scripts:

   ```bash
   python 1_Genetic_Background.py
   ```

This setup ensures that all tools and dependencies are correctly configured for the analysis pipeline.

## Workflow

The analysis pipeline consists of the following steps:

1. **Sequence Retrieval**:

   - Input TaxIDs or query terms to download sequences from NCBI.
   - Validate and process retrieved genetic data.

2. **Multiple Sequence Alignment (MSA)**:

   - Align sequences using MAFFT.
   - Export alignment files for further processing.

3. **Consensus Sequence Generation**:

   - Create consensus sequences from the aligned data.
   - Specify nucleotide agreement thresholds as needed.

4. **Consensus vs. NCBI GenBank Analysis**:

   - Compare consensus sequences with NCBI GenBank entries.
   - Identify potential matches and mismatches.

5. **Primer Design**:

   - Generate primers using Primer3.
   - Specify parameters such as GC content, melting temperature, and amplicon size.

6. **Primer Specificity Testing**:

   - Evaluate primers against genetic databases.
   - Ensure high specificity to target sequences.

7. **Primer Set Specificity Analysis**:

   - Analyze specificity of primer sets across sequences.
   - Generate detailed annotations for each primer set.

8. **Species-Specificity Analysis**:

   - Assess primer matches to specific taxa.
   - Generate detailed reports for interpretation.

## Usage

Run the scripts in sequential order to perform a complete analysis. For example:

```bash
python 1_Genetic_Background.py

python 2_MSA_Alignment.py

python 3_Consensus.py

python 4_Consensus_vs_NCBI_GenBank.py

python 5_PCR_primers_design.py

python 6_Primer_left_specificity_test.py

python 7_Primer_Set_Specificity_Analysis.py

python 8_Species_Specificity_Test.py
```

Each script will prompt for required input and provide detailed logs of the operations performed.

## Dependencies

- **Python**: Version 3.8 or higher.
- **BioPython**: For sequence handling and database queries.
- **Primer3-py**: For primer design.
- **MAFFT**: For multiple sequence alignment.
- **NCBI Entrez**: For sequence retrieval.
- **Other Libraries**:
  - `numpy`
  - `pandas`
  - `tqdm`
  - `logging`

Install all Python dependencies using:

```bash
pip install -r requirements.txt
```

## License

This repository is licensed under the MIT License. See the `LICENSE` file for details.

## Disclaimer

This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability arising from the use of this software.
