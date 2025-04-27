# Blumeria graminis Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Blumeria graminis*, a fungal plant pathogen that causes powdery mildew disease in cereals and grasses.

## Overview

*Blumeria graminis* (TaxID: 34373) is an obligate biotrophic fungus that is highly host-specific, with distinct formae speciales that infect different grass species. This example demonstrates how to design specific primers that can identify *B. graminis* in mixed samples and potentially differentiate between its various specialized forms.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Blumeria graminis/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 34373_BTUB2.fasta
│   ├── 34373_CYP51.fasta
│   ├── 34373_TEF1.fasta
│   ├── 34373_TRANSLATION_ELONGATION_FACTOR_1_ALPHA.fasta
│   ├── 34373_gene_stats.txt
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_0.txt
│   ├── BLAST_1.txt
│   ├── BLAST_2.txt
│   ├── BLAST_3.txt
│   ├── CONSENSUS_0.fasta
│   ├── CONSENSUS_1.fasta
│   ├── CONSENSUS_2.fasta
│   ├── CONSENSUS_3.fasta
│   ├── MAFFT_0.txt
│   ├── MAFFT_1.txt
│   ├── MAFFT_2.txt
│   └── MAFFT_3.txt
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
│   ├── CONSENSUS_1_primers.csv
│   ├── CONSENSUS_2_primers.csv
│   └── CONSENSUS_3_primers.csv
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
│   ├── specificity_report_CONSENSUS_1_primers.html
│   ├── specificity_report_CONSENSUS_2_primers.html
│   └── specificity_report_CONSENSUS_3_primers.html
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 34373 (*Blumeria graminis*). The analysis identified 50 unique genes represented in the database as shown in the `34373_gene_stats.txt` file. The most abundant genes retrieved were:

1. INTERNAL_TRANSCRIBED_SPACER_1 (455 sequences)
2. INTERNAL_TRANSCRIBED_SPACER_2 (443 sequences)
3. TUB2 (410 sequences)
4. 5.8S_RIBOSOMAL_RNA (396 sequences)
5. ELONGATION_FACTOR_1_ALPHA (357 sequences)

For this example, genes with sufficient sequence representation were selected for alignment, including:
- BTUB2 (Beta-tubulin)
- CYP51 (Lanosterol 14-alpha demethylase)
- TEF1 (Translation elongation factor 1-alpha)
- TRANSLATION_ELONGATION_FACTOR_1_ALPHA

These genes represent a mix of housekeeping genes (BTUB2, TEF1) and fungicide target genes (CYP51) to provide multiple options for species-specific detection.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. Four alignment groups were created, resulting in four consensus sequences:

1. `CONSENSUS_0.fasta`: A 449 bp consensus based on beta-tubulin sequences.
2. `CONSENSUS_1.fasta`: A 317 bp consensus sequence derived from a different group of sequences.
3. `CONSENSUS_2.fasta`: A 443 bp consensus focused on elongation factor 1-alpha.
4. `CONSENSUS_3.fasta`: A 392 bp consensus also derived from elongation factor 1-alpha sequences.

The alignment files (MAFFT_*.txt) contain the full multiple sequence alignments, while the BLAST_*.txt files show the alignments in BLAST-like format for easier visual inspection.

### Module 3: PCR Primer Design

Using the consensus sequences, Module 3 designed primer-probe sets for qPCR. The log file shows that amplicon ranges of 70-250 bp were selected for all four consensus sequences, and 6 primer sets were generated for each consensus.

Each CSV file contains these fields for the primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

For example, from `CONSENSUS_0_primers.csv` (primer set #1):
- Forward primer: TGAGCCATACAACGCAACTC
- Reverse primer: ACCTGACATTACGGCAGACA
- Probe: AAACTCCGACGAGACGTTCT
- Product size: 166 bp

These primers target different regions of the consensus sequences, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *B. graminis* (TaxID: 34373) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- CONSENSUS_0 primers had moderate match counts (ranging from 56 to 228 sequences)
- CONSENSUS_1 primers consistently matched 69 sequences
- CONSENSUS_2 primers showed high specificity with 0 matching sequences
- CONSENSUS_3 primers also showed high specificity with 0 matching sequences

The specificity reports (HTML files) provide detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *B. graminis* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed the best specificity for *Blumeria graminis*:

1. From CONSENSUS_2, any primer set (sets #1-6):
   - Example (set #1):
     - Forward: TGGCCGTGTCCATCTTGTTA
     - Reverse: TGTGCCATTCTCATCATCGC
     - Probe: ACCTAATGTGAACGCAAGCA
     - Product size: 145 bp
     - No cross-reactivity with non-target species

2. From CONSENSUS_3, any primer set (sets #1-6):
   - Example (set #1): 
     - Forward: TGGCCGTGTCCATCTTGTTA
     - Reverse: TGTGCCATTCTCATCATCGC
     - Probe: ACCTAATGTGAACGCAAGCA
     - Product size: 145 bp
     - No cross-reactivity with non-target species

These primer sets targeting the elongation factor 1-alpha gene (TEF1) showed perfect specificity for *B. graminis*, with no matches to non-target organisms in the database.

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 34373

3. Run each module in sequence:
   - In Module 1, select the genes mentioned above
   - In Module 2, create at least 4 alignment groups
   - In Module 3, use amplicon size range 70-250 bp
   - In Module 4, specify TaxID 34373 for specificity testing

4. Review the HTML reports in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- Pure cultures of *B. graminis* (if available)
- Closely related fungi
- Plant samples infected with powdery mildew disease

For qPCR validation, recommended conditions:
- Annealing temperature: 59°C
- Cycling conditions: 95°C for 15s, 59°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *B. graminis* demonstrates the challenge of designing species-specific primers for an obligate biotroph that cannot be grown in axenic culture. The elongation factor 1-alpha gene (TEF1) provided the best specificity, showing no cross-reactivity with other organisms.

The primers generated from CONSENSUS_2 and CONSENSUS_3 could potentially be used to differentiate between different formae speciales of *B. graminis* that affect different host plants, though additional validation would be needed to confirm this capability.

For applications requiring identification of *B. graminis* directly from plant material, these highly specific primers offer a promising approach that avoids the need for culturing this obligate biotrophic pathogen.
