# Treponema denticola Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Treponema denticola*, an oral spirochete bacterium associated with periodontal disease and a key component of the "red complex" of periodontal pathogens.

## Overview

*Treponema denticola* (TaxID: 158) is a motile, anaerobic spirochete found in the human oral cavity. This example demonstrates how to design specific primers that can identify *T. denticola* in dental plaque samples and clinical specimens where multiple oral bacterial species are present.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Treponema denticola/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 158_ERA.fasta
│   ├── 158_RADC.fasta 
│   ├── 158_RECA.fasta
│   ├── 158_ERA_blast_results.xml
│   ├── 158_RADC_blast_results.xml
│   ├── 158_RECA_blast_results.xml
│   ├── 158_gene_stats.txt
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_0.txt
│   ├── BLAST_1.txt
│   ├── BLAST_2.txt
│   ├── CONSENSUS_0.fasta
│   ├── CONSENSUS_1.fasta
│   ├── CONSENSUS_2.fasta
│   ├── MAFFT_0.txt
│   ├── MAFFT_1.txt
│   ├── MAFFT_2.txt
│   └── ...
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
│   ├── CONSENSUS_1_primers.csv
│   ├── CONSENSUS_2_primers.csv
│   └── ...
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
│   ├── specificity_report_CONSENSUS_1_primers.html
│   ├── specificity_report_CONSENSUS_2_primers.html
│   └── ...
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 158 (*Treponema denticola*). The analysis identified 50 unique genes represented in the database as shown in the `158_gene_stats.txt` file. The most abundant genes retrieved were:

1. PYRH (634 sequences)
2. 16S_RIBOSOMAL_RNA (284 sequences)
3. GAP (222 sequences)
4. RECA (98 sequences)
5. RADC (86 sequences)
6. ERA (86 sequences)

Several genes with sufficient sequence representation were selected for alignment, including ERA, RADC, and RECA. The files retrieved include both the sequence data (FASTA format) and BLAST result information.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. Three alignment groups were created, resulting in three consensus sequences:

1. `CONSENSUS_0.fasta`: A 885 bp consensus sequence with some ambiguous bases.
2. `CONSENSUS_1.fasta`: A 678 bp consensus sequence.
3. `CONSENSUS_2.fasta`: A 1245 bp consensus sequence, containing 10.9% ambiguous bases (N).

The alignment files (MAFFT_*.txt) contain the full multiple sequence alignments, while the BLAST_*.txt files show the alignments in BLAST-like format for easier visual inspection.

### Module 3: PCR Primer Design

Using the consensus sequences, Module 3 designed primer-probe sets for qPCR. According to the log file, an amplicon length range of 70-120 bp was specified for the design, and 10 primer sets were generated for each consensus sequence.

Each CSV file contains these fields for the 10 primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_0_primers.csv` (primer set #1):
- Forward primer: GGAGTTGTAACGATAATAGGCCG
- Reverse primer: GGATTTGTGGTAGCCGGGG
- Probe: TCGTAAACACTACCAAGGGACA
- Product size: 183 bp

These primers target different regions of the consensus sequences, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *T. denticola* (TaxID: 158) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- CONSENSUS_0 primers had 13 matching sequences across all primer sets
- CONSENSUS_1 primers had 12 matching sequences across all primer sets
- CONSENSUS_2 primers showed more variability with between 26-34 matching sequences

The specificity reports (HTML files) provide detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *T. denticola* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed good specificity for *Treponema denticola*:

1. From CONSENSUS_0, primer set #1:
   - Forward: GGAGTTGTAACGATAATAGGCCG
   - Reverse: GGATTTGTGGTAGCCGGGG
   - Probe: TCGTAAACACTACCAAGGGACA
   - Product size: 183 bp
   - Limited cross-reactivity with non-target species

2. From CONSENSUS_1, primer set #3:
   - Forward: AAGGCTTTTAGAGTACGGGC
   - Reverse: CATCGGCCAATTCTTTTACAGG
   - Probe: TTGGTGGCCATTCTTTTACG
   - Product size: 104 bp
   - Good specificity for *T. denticola*

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 158

3. Run each module in sequence:
   - In Module 1, select genes with high representation like PYRH, ERA, RADC, and RECA
   - In Module 2, create at least 3 alignment groups
   - In Module 3, use amplicon size range 70-120 bp
   - In Module 4, specify TaxID 158 for specificity testing

4. Review the HTML reports in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- Pure cultures of *T. denticola*
- Closely related oral treponemes
- Dental plaque samples containing mixed oral bacteria

For qPCR validation, recommended conditions:
- Annealing temperature: 59-60°C (based on primer Tm values)
- Cycling conditions: 95°C for 15s, 60°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *T. denticola* demonstrates the challenge of designing species-specific primers within the Treponema genus, which contains many closely related oral spirochetes. The specificity testing results suggest that the CONSENSUS_0 and CONSENSUS_1 primers provide better specificity than those from CONSENSUS_2, which had more non-specific matches.

For clinical applications, these primers could be valuable for quantifying *T. denticola* in periodontal disease specimens, monitoring treatment efficacy, or studying the microbial ecology of oral biofilms. The short amplicon sizes (70-120 bp) make these primers particularly suitable for qPCR applications where sample DNA may be partially degraded.
