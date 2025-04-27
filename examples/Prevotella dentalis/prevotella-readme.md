# Prevotella dentalis Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Prevotella dentalis*, an anaerobic Gram-negative bacterium found in the human oral cavity and associated with periodontal disease.

## Overview

*Prevotella dentalis* (TaxID: 52227) is a member of the genus Prevotella that is implicated in periodontal infections and dental abscesses. This example demonstrates how to design specific primers that can identify *P. dentalis* in mixed samples containing multiple oral bacterial species or other clinical specimens.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Prevotella dentalis/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 52227_16S_RIBOSOMAL_RNA.fasta
│   ├── 52227_DNAJ.fasta
│   ├── 52227_GYRB.fasta
│   ├── 52227_HSP60.fasta
│   ├── 52227_RECA.fasta
│   ├── 52227_RPOB.fasta
│   ├── 52227_gene_stats.txt
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_0.txt
│   ├── BLAST_1.txt
│   ├── BLAST_2.txt
│   ├── BLAST_3.txt
│   ├── BLAST_4.txt
│   ├── BLAST_5.txt
│   ├── CONSENSUS_0.fasta
│   ├── CONSENSUS_1.fasta
│   ├── CONSENSUS_2.fasta
│   ├── CONSENSUS_3.fasta
│   ├── CONSENSUS_4.fasta
│   ├── CONSENSUS_5.fasta
│   ├── MAFFT_0.txt
│   ├── MAFFT_1.txt
│   ├── MAFFT_2.txt
│   ├── MAFFT_3.txt
│   ├── MAFFT_4.txt
│   └── MAFFT_5.txt
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
│   ├── CONSENSUS_1_primers.csv
│   ├── CONSENSUS_2_primers.csv
│   ├── CONSENSUS_3_primers.csv
│   ├── CONSENSUS_4_primers.csv
│   ├── CONSENSUS_5_primers.csv
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
│   ├── specificity_report_CONSENSUS_1_primers.html
│   ├── specificity_report_CONSENSUS_2_primers.html
│   ├── specificity_report_CONSENSUS_3_primers.html
│   ├── specificity_report_CONSENSUS_4_primers.html
│   ├── specificity_report_CONSENSUS_5_primers.html
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 52227 (*Prevotella dentalis*). The analysis identified 50 unique genes represented in the database as shown in the `52227_gene_stats.txt` file. The most abundant genes retrieved were:

1. RECA (20 sequences)
2. DNAJ (20 sequences)
3. RPOB (17 sequences)
4. GYRB (17 sequences)
5. HSP60 (12 sequences)
6. 16S_RIBOSOMAL_RNA (10 sequences)

For this example, six genes with sufficient sequence representation were selected for alignment:
- 16S_RIBOSOMAL_RNA
- DNAJ
- GYRB
- HSP60
- RECA
- RPOB

These genes represent a mix of conserved housekeeping genes and ribosomal RNA that are commonly used for bacterial identification and phylogenetic analysis.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. Six alignment groups were created, resulting in six consensus sequences:

1. `CONSENSUS_0.fasta`: A 1444 bp consensus from the 16S ribosomal RNA gene.
2. `CONSENSUS_1.fasta`: A 755 bp consensus sequence based on a conserved region.
3. `CONSENSUS_2.fasta`: A 1107 bp consensus sequence from another gene group.
4. `CONSENSUS_3.fasta`: A 1489 bp consensus sequence also from 16S ribosomal RNA.
5. `CONSENSUS_4.fasta`: A 600 bp consensus sequence.
6. `CONSENSUS_5.fasta`: A 2031 bp consensus sequence, the largest in this analysis.

The alignment files (MAFFT_*.txt) contain the full multiple sequence alignments, while the BLAST_*.txt files show the alignments in BLAST-like format for easier visual inspection.

### Module 3: PCR Primer Design

Using the consensus sequences, Module 3 designed primer-probe sets for qPCR. The log file shows that amplicon ranges of 70-120 bp were selected for all six consensus sequences, and 10 primer sets were generated for each consensus, totaling 60 candidate primer sets.

Each CSV file contains these fields for the 10 primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_3_primers.csv` (primer set #1):
- Forward primer: TGCGTCCGATTAGCTTGTTG
- Reverse primer: GTCATCCTGCACGCTACTTG
- Probe: GTTCTGAGAGGAAGGTCCCC
- Product size: 181 bp

These primers target different regions of the consensus sequences, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *P. dentalis* (TaxID: 52227) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- CONSENSUS_0 primers showed high match counts (ranging from 0 to 749 sequences)
- CONSENSUS_1 primers were more specific with most sets having fewer than 5 matches
- CONSENSUS_2 primers matched between 0-8 sequences
- CONSENSUS_3 primers showed similar match patterns to CONSENSUS_0
- CONSENSUS_4 primers consistently matched around 7 sequences
- CONSENSUS_5 primers were among the most specific, with several sets having 0 matches in GenBank

The specificity reports (HTML files) provide detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *P. dentalis* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed the best specificity for *Prevotella dentalis*:

1. From CONSENSUS_5, primer set #4:
   - Forward: CAAGCTCGACGACCTGAAAG
   - Reverse: TGCATCTTGTCCTCCACCAT
   - Probe: ATCTACGATGGTGAGACGGG
   - Product size: 162 bp
   - No cross-reactivity with non-target species

2. From CONSENSUS_1, primer set #1:
   - Forward: TCGAAATCAACATTCCCGCC
   - Reverse: CTTCTCTTCCTCCACCAGCA
   - Probe: GGCATCAACGGCGACATAC
   - Product size: 122 bp
   - Excellent specificity for *P. dentalis*

3. From CONSENSUS_4, primer set #10:
   - Forward: ACGCCCTGAAGTTCTATGCA
   - Reverse: GTGATTGCCCACCACATTGT
   - Probe: AGCGTTCGCGTGGACATC
   - Product size: 86 bp
   - Highly specific for the target species

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 52227

3. Run each module in sequence:
   - In Module 1, select the six genes mentioned above
   - In Module 2, create at least 6 alignment groups
   - In Module 3, use amplicon size range 70-120 bp
   - In Module 4, specify TaxID 52227 for specificity testing

4. Review the HTML reports in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- Pure cultures of *P. dentalis*
- Closely related species (especially other Prevotella species)
- Clinical samples from periodontal infections

For qPCR validation, recommended conditions:
- Annealing temperature: 60°C
- Cycling conditions: 95°C for 15s, 60°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *P. dentalis* demonstrates the challenges of designing species-specific primers within the Prevotella genus, which contains many closely related species in the oral microbiome. The 16S ribosomal RNA gene (CONSENSUS_0 and CONSENSUS_3) showed higher conservation across species, while some of the housekeeping genes (especially CONSENSUS_5) provided better specificity.

For applications requiring absolute specificity, consider using multiple primer sets in combination for confirmation, or conducting post-PCR melt curve analysis to distinguish closely related amplicons.
