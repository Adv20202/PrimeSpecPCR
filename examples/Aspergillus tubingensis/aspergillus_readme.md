# Aspergillus tubingensis Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Aspergillus tubingensis*, a filamentous fungus that is a common food contaminant and potential producer of mycotoxins.

## Overview

*Aspergillus tubingensis* (TaxID: 5068) is a black Aspergillus species closely related to *A. niger* but distinct in its mycotoxin production profile and ecological niche. This example demonstrates how to design specific primers that can identify *A. tubingensis* in mixed samples containing multiple Aspergillus species or other fungal contaminants.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Aspergillus tubingensis/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 5068_ATWU_R58.fasta
│   ├── 5068_BETA_TUBULIN.fasta
│   ├── 5068_DNA_DEPENDENT_RNA_POLYMERASE_BETA.fasta
│   ├── 5068_RPB2.fasta
│   ├── 5068_STEROL_14_ALPHA_DEMETHYLASE.fasta
│   ├── 5068_gene_stats.txt
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_1.txt
│   ├── BLAST_2.txt
│   ├── BLAST_3.txt
│   ├── BLAST_4.txt
│   ├── CONSENSUS_1.fasta
│   ├── CONSENSUS_2.fasta
│   ├── CONSENSUS_3.fasta
│   ├── CONSENSUS_4.fasta
│   ├── MAFFT_1.txt
│   ├── MAFFT_2.txt
│   ├── MAFFT_3.txt
│   └── MAFFT_4.txt
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_1_primers.csv
│   ├── CONSENSUS_2_primers.csv
│   ├── CONSENSUS_3_primers.csv
│   └── CONSENSUS_4_primers.csv
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_1_primers.html
│   ├── specificity_report_CONSENSUS_2_primers.html
│   ├── specificity_report_CONSENSUS_3_primers.html
│   └── specificity_report_CONSENSUS_4_primers.html
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 5068 (*Aspergillus tubingensis*). The analysis identified 50 unique genes represented in the database as shown in the `5068_gene_stats.txt` file. The most abundant genes retrieved were:

1. BETA_TUBULIN (1306 sequences)
2. CALMODULIN (1244 sequences)
3. INTERNAL_TRANSCRIBED_SPACER_1 (1088 sequences)
4. INTERNAL_TRANSCRIBED_SPACER_2 (1073 sequences)
5. 5.8S_RIBOSOMAL_RNA (895 sequences)

For this example, five genes with sufficient sequence representation were selected for alignment:
- BETA_TUBULIN
- ATWU_R58
- DNA_DEPENDENT_RNA_POLYMERASE_BETA
- RPB2
- STEROL_14_ALPHA_DEMETHYLASE

These genes represent a mix of housekeeping genes (BETA_TUBULIN), RNA polymerase components (DNA_DEPENDENT_RNA_POLYMERASE_BETA, RPB2), and specialized metabolism genes (STEROL_14_ALPHA_DEMETHYLASE) to provide multiple options for species-specific detection.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. Four alignment groups were created, resulting in four consensus sequences:

1. `CONSENSUS_1.fasta`: Focused on the beta-tubulin gene, a 556 bp consensus was generated.
2. `CONSENSUS_2.fasta`: A 983 bp consensus sequence from a different group of sequences.
3. `CONSENSUS_3.fasta`: An 838 bp consensus sequence.
4. `CONSENSUS_4.fasta`: The largest consensus at 1551 bp, though contains 20.5% ambiguous bases.

The alignment files (MAFFT_*.txt) contain the full multiple sequence alignments, while the BLAST_*.txt files show the alignments in BLAST-like format for easier visual inspection.

### Module 3: PCR Primer Design

Using the consensus sequences, Module 3 designed primer-probe sets for qPCR. The log file shows that amplicon ranges of 70-120 bp were selected for all four consensus sequences, and 10 primer sets were generated for each consensus.

Each CSV file contains these fields for the 10 primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_1_primers.csv` (primer set #1):
- Forward primer: ATCACACCGTCCCTGAGTTT
- Reverse primer: TTACCGCTAGCCTGCTGAAG
- Probe: ACGACAATATCATCAATGTCCTGA
- Product size: 73 bp

These primers target different regions of the consensus sequences, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *A. tubingensis* (TaxID: 5068) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- CONSENSUS_1 primers had high match counts (ranging from 1 to 897 sequences)
- CONSENSUS_2 primers matched between 90-240 sequences
- CONSENSUS_3 primers showed similar match patterns to CONSENSUS_2
- CONSENSUS_4 primers were most specific, with most sets having no matches in GenBank

The specificity reports (HTML files) provide detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *A. tubingensis* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed the best specificity for *Aspergillus tubingensis*:

1. From CONSENSUS_4, primer set #5:
   - Forward: AGTGCCAGGAACCGATTACT
   - Reverse: GGAGTTGTCTTGGATGCGTC
   - Probe: TCTTCTCCGGTCCCTTGAAG
   - Product size: 88 bp
   - No cross-reactivity with non-target species

2. From CONSENSUS_2, primer set #3:
   - Forward: GTGGTTCACTGGTGCTCAAC
   - Reverse: CCATCCCAACCAAAGTAGCG
   - Probe: CGCAAACTGGAGCAAGACAA
   - Product size: 112 bp
   - Limited cross-reactivity with closely related species

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 5068

3. Run each module in sequence:
   - In Module 1, select the five genes mentioned above
   - In Module 2, create at least 4 alignment groups
   - In Module 3, use amplicon size range 70-120 bp
   - In Module 4, specify TaxID 5068 for specificity testing

4. Review the HTML reports in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- Pure cultures of *A. tubingensis*
- Closely related species (especially *A. niger*)
- Environmental samples containing mixed fungi

For qPCR validation, recommended conditions:
- Annealing temperature: 60°C
- Cycling conditions: 95°C for 15s, 60°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *A. tubingensis* demonstrates the challenge of designing species-specific primers within the *Aspergillus* genus, which contains many closely related black Aspergilli. The beta-tubulin gene (CONSENSUS_1) showed higher conservation across species, while some of the other targets (especially CONSENSUS_4) provided better specificity.

For applications requiring absolute specificity, consider using multiple primer sets in combination for confirmation, or conducting post-PCR melt curve analysis to distinguish closely related amplicons.
