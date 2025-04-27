# Puccinia recondita Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Puccinia recondita*, a fungal plant pathogen that causes leaf rust in wheat and other cereal crops.

## Overview

*Puccinia recondita* (TaxID: 27348) is an economically important biotrophic fungal pathogen responsible for leaf rust disease in wheat and other cereals worldwide. This example demonstrates how to design specific primers that can identify *P. recondita* in mixed samples containing other fungal species or plant material.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Puccinia recondita/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 27348_TEF1.fasta
│   ├── 27348_INTERNAL_TRANSCRIBED_SPACER_2.fasta
│   ├── 27348_28S_RIBOSOMAL_RNA.fasta
│   ├── 27348_5.8S_RIBOSOMAL_RNA.fasta
│   ├── 27348_INTERNAL_TRANSCRIBED_SPACER_1.fasta
│   ├── 27348_gene_stats.txt
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_0.txt
│   ├── CONSENSUS_0.fasta
│   ├── MAFFT_0.txt
│   └── ...
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 27348 (*Puccinia recondita*). The analysis identified 15 unique genes represented in the database as shown in the `27348_gene_stats.txt` file. The most abundant genes retrieved were:

1. INTERNAL_TRANSCRIBED_SPACER_2 (118 sequences)
2. 28S_RIBOSOMAL_RNA (68 sequences)
3. 5.8S_RIBOSOMAL_RNA (67 sequences)
4. LARGE_SUBUNIT_RIBOSOMAL_RNA (45 sequences)
5. INTERNAL_TRANSCRIBED_SPACER_1 (38 sequences)
6. TEF1 (32 sequences)

For this example, the Translation Elongation Factor 1-alpha (TEF1) gene was selected as the primary target, as it provides sufficient sequence variation for species differentiation within rust fungi while maintaining adequate conservation for reliable primer design.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were aligned using MAFFT to identify conserved regions. The alignment resulted in a consensus sequence:

`CONSENSUS_0.fasta`: A consensus sequence based on TEF1 gene alignment.

The alignment files (MAFFT_0.txt) contain the full multiple sequence alignment, while the BLAST_0.txt file shows the alignment in BLAST-like format for easier visual inspection.

The consensus sequence initially contained a high percentage of ambiguous bases (N), which presented challenges for primer design. After iterative refinement, the consensus sequence was improved to contain 15.7% ambiguous bases, which allowed successful primer design.

### Module 3: PCR Primer Design

Using the consensus sequence, Module 3 designed primer-probe sets for qPCR. The log file shows that amplicon ranges of 100-125 bp were selected, and 20 primer sets were generated.

The CSV file contains these fields for the 20 primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_0_primers.csv` (primer set #1):
- Forward primer: GCTATGGAAGTTCGAAACCCC
- Reverse primer: GATGAGGATAGCACAATCGGC
- Probe: ACATGATCACCGGTACCTCC
- Product size: 118 bp

These primers target conserved regions of the TEF1 gene, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *P. recondita* (TaxID: 27348) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- Primer pair 1 matched 242 sequences
- Primer pair 4 matched 540 sequences (least specific)
- Primer pair 13 matched 58 sequences (most specific)
- Primer pair 17 matched 175 sequences

The specificity report (HTML file) provides detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *P. recondita* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed the best specificity for *Puccinia recondita*:

1. Primer set #13:
   - Forward: ACCACTGGTCACTTGATCTACA
   - Reverse: GATCCTTTGCCCAACTCGG
   - Probe: GTGCGGTGGTATTGACAAGC
   - Product size: 92 bp
   - Lowest non-target cross-reactivity (58 matches)

2. Primer set #20:
   - Forward: AGCTATGGAAGTTCGAAACCC
   - Reverse: TGATGAGGATAGCACAATCGG
   - Probe: ACATGATCACCGGTACCTCC
   - Product size: 120 bp
   - Good specificity (128 matches)

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 27348

3. Run each module in sequence:
   - In Module 1, select the TEF1 gene
   - In Module 2, create alignment groups
   - In Module 3, use amplicon size range 100-125 bp and generate 20 primer sets
   - In Module 4, specify TaxID 27348 for specificity testing

4. Review the HTML report in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- Pure cultures of *P. recondita*
- Closely related species (especially other *Puccinia* species)
- Infected plant material containing mixed fungal species

For qPCR validation, recommended conditions:
- Annealing temperature: 60°C
- Cycling conditions: 95°C for 15s, 60°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *P. recondita* demonstrates the challenge of designing species-specific primers within the *Puccinia* genus, which contains many closely related rust species that share similar genomic regions. The TEF1 gene was selected for its balance of conservation and variation, providing regions that are conserved enough for primer binding but variable enough for species discrimination.

For applications requiring absolute specificity, consider using multiple primer sets in combination for confirmation, or conducting post-PCR melt curve analysis to distinguish closely related amplicons.
