# PrimeSpecPCR Quick Start Guide

## Overview

PrimeSpecPCR is a toolkit for designing species-specific PCR primers. It automates the process from retrieving genetic sequences through primer design and specificity testing.

## Prerequisites

Before using PrimeSpecPCR, you need:

1. **NCBI Account and API Key**:
   - Register at NCBI: https://www.ncbi.nlm.nih.gov/
   - Get an API key: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

2. **TaxID Number(s)**:
   - Identify the NCBI Taxonomy ID for your target organism(s).
   - Example TaxIDs: 
     - 36050 for *Fusarium poae* (wheat pathogen)
     - 27348 for *Puccinia recondita* (wheat leaf rust fungus)
     - 4787 for *Rhizopus stolonifer* (black bread mold)

## Getting Started

### Option 1: Using the GUI (recommended)

1. **Launch PrimeSpecPCR**:
   ```
   python PrimeSpecPCR.py
   ```

2. **Enter your credentials**:
   - Fill in your NCBI Email
   - Enter your NCBI API Key
   - Enter TaxID number(s)
   - Click "Save Configuration"

3. **Run the modules in sequence**:
   - Start with Module 1 (Genetic Sequence Retrieval)
   - Continue with each module in order (1 through 8)
   - The interface will show progress and log output

### Option 2: Using Individual Scripts

For advanced users who prefer command-line operation:

1. First, set up your configuration:
   ```
   python PrimeSpecPCR.py
   ```
   (Enter your credentials and save)

2. Then run each module individually:
   ```
   python 1_Genetic_Background-2.py
   python 2_MSA_Alignment.py
   # etc.
   ```

## What to Expect

**Module 1: Genetic Sequence Retrieval**
- Downloads genetic sequences from NCBI for specified TaxIDs
- Will automatically select genes with sufficient data
- Creates FASTA files in the "1_" directory
- Typical runtime: 5-15 minutes depending on organism

**Subsequent Modules**
- Module 2: Aligns sequences using MAFFT
- Module 3: Generates consensus sequences
- Module 4: Compares consensus to NCBI database
- Module 5: Designs PCR primers
- Module 6: Tests primer specificity
- Module 7: Analyzes primer set specificity
- Module 8: Creates visual specificity report

## Troubleshooting

- **Long-running operations**: Some NCBI queries take time - the progress bar shows the application is working
- **Process seems stuck**: Use "Stop Current Process" and try again. The script includes status updates to show it's still working
- **Error messages**: Check logs directory for detailed error information
- **MAFFT errors**: Make sure MAFFT is installed and in your system PATH

## Output

Each module creates a numbered directory with its results:
- 1_: Contains FASTA files with gene sequences
- 2_: Contains multiple sequence alignments
- 3_: Contains consensus sequences
- ...
- 8_: Contains HTML dashboard of species specificity

## Need Help?

For more detailed information, refer to the README.md file or visit our GitHub repository.
