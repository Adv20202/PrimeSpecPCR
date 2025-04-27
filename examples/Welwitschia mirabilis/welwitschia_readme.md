# Welwitschia mirabilis Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Welwitschia mirabilis*, a unique gymnosperm native to the Namib Desert that is the only living species in its family and genus.

## Overview

*Welwitschia mirabilis* (TaxID: 3377) is a remarkable long-lived plant with only two leaves that grow continuously throughout its lifespan, which can exceed 1,000 years. This example demonstrates how to design specific primers that can identify *W. mirabilis* in samples for taxonomic studies, conservation monitoring, or evolutionary analysis.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Welwitschia mirabilis/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 3377_ATPB.fasta
│   ├── 3377_MATK.fasta
│   ├── 3377_RBCL.fasta
│   ├── 3377_RPL16.fasta
│   ├── 3377_RPS12.fasta
│   ├── 3377_RPS3.fasta
│   ├── 3377_RPS4.fasta
│   ├── 3377_TRNA.fasta
│   ├── 3377_TRNH.fasta
│   ├── 3377_ATPB_blast_results.xml
│   ├── 3377_MATK_blast_results.xml
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_0.txt
│   ├── BLAST_1.txt
│   ├── BLAST_2.txt
│   ├── BLAST_3.txt
│   ├── BLAST_4.txt
│   ├── BLAST_5.txt
│   ├── BLAST_6.txt
│   ├── BLAST_7.txt
│   ├── BLAST_8.txt
│   ├── CONSENSUS_0.fasta
│   ├── CONSENSUS_1.fasta
│   ├── CONSENSUS_2.fasta
│   ├── CONSENSUS_3.fasta
│   ├── CONSENSUS_4.fasta
│   ├── CONSENSUS_5.fasta
│   ├── CONSENSUS_6.fasta
│   ├── CONSENSUS_7.fasta
│   ├── CONSENSUS_8.fasta
│   ├── MAFFT_0.txt
│   ├── MAFFT_1.txt
│   └── ...
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
│   ├── CONSENSUS_1_primers.csv
│   ├── CONSENSUS_2_primers.csv
│   ├── CONSENSUS_3_primers.csv
│   ├── CONSENSUS_4_primers.csv
│   ├── CONSENSUS_5_primers.csv
│   ├── CONSENSUS_6_primers.csv
│   ├── CONSENSUS_7_primers.csv
│   └── CONSENSUS_8_primers.csv
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
│   ├── specificity_report_CONSENSUS_1_primers.html
│   ├── specificity_report_CONSENSUS_2_primers.html
│   ├── specificity_report_CONSENSUS_3_primers.html
│   ├── specificity_report_CONSENSUS_4_primers.html
│   ├── specificity_report_CONSENSUS_5_primers.html
│   ├── specificity_report_CONSENSUS_6_primers.html
│   ├── specificity_report_CONSENSUS_7_primers.html
│   └── specificity_report_CONSENSUS_8_primers.html
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 3377 (*Welwitschia mirabilis*). The analysis identified 50 unique genes represented in the database as shown in the `3377_gene_stats.txt` file. The most abundant genes retrieved were:

1. TRNA (15 sequences)
2. RBCL (11 sequences)
3. RPS3 (9 sequences)
4. RPS12 (8 sequences)
5. RPS4 (8 sequences)
6. RRN5 (7 sequences)
7. MATK (7 sequences)
8. ATPB (7 sequences)
9. RPL16 (6 sequences)
10. TRNH (6 sequences)

For this example, several genes with sufficient sequence representation were selected for alignment, including ATPB, MATK, RBCL, RPL16, RPS12, RPS3, RPS4, TRNA, and TRNH. These genes represent a mix of chloroplast genes, ribosomal components, and transfer RNAs to provide multiple options for species-specific detection.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. Nine alignment groups were created, resulting in nine consensus sequences:

1. `CONSENSUS_0.fasta`: A 1137 bp consensus sequence with some ambiguous bases, identified as cox2/cox1 genes.
2. `CONSENSUS_1.fasta`: An 826 bp consensus sequence.
3. `CONSENSUS_2.fasta`: A 1137 bp consensus sequence similar to CONSENSUS_0, labeled as cox2/cox1 genes.
4. `CONSENSUS_3.fasta`: A 721 bp consensus sequence.
5. `CONSENSUS_4.fasta`: A 1209 bp consensus sequence identified as cox3 gene.
6. `CONSENSUS_5.fasta`: A shorter 293 bp consensus sequence labeled as cox2/cox1/cox3 genes.
7. `CONSENSUS_6.fasta`: A 575 bp consensus sequence.
8. `CONSENSUS_7.fasta`: A 1043 bp consensus sequence.
9. `CONSENSUS_8.fasta`: A 681 bp consensus sequence.

The alignment files (MAFFT_*.txt) contain the full multiple sequence alignments, while the BLAST_*.txt files show the alignments in BLAST-like format for easier visual inspection.

### Module 3: PCR Primer Design

Using the consensus sequences, Module 3 designed primer-probe sets for qPCR. According to the log file, 5 primer sets were generated for each consensus sequence. The amplicon length range varied for each consensus, but generally fell within 70-150 bp ranges.

Each CSV file contains these fields for the 5 primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_5_primers.csv` (primer set #1):
- Forward primer: CGGAGTTAGGGCAGGAGAAA
- Reverse primer: GGCCAGATCCTTACAAAGCG
- Probe: GTCAAGTCGACCGGGCTAC
- Product size: 78 bp

These primers target different regions of the consensus sequences, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *W. mirabilis* (TaxID: 3377) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- CONSENSUS_0 primers had 17 matching sequences for primer set #1, but 0 matches for sets #2-4, and 1 match for set #5
- CONSENSUS_1 primers had consistently 1 matching sequence across all primer sets
- CONSENSUS_4 primers were exceptionally specific, with 0 matching sequences for all primer sets
- CONSENSUS_5 primers had consistent 2 matches across all primer sets
- CONSENSUS_8 primers had 6 matching sequences for primer set #1, and 2 matches for sets #2-5

The specificity reports (HTML files) provide detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *W. mirabilis* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed the best specificity for *Welwitschia mirabilis*:

1. From CONSENSUS_4, any primer set:
   - Example set #1:
     - Forward: CGTCAACAAGGGCGTTCTAG
     - Reverse: ATAGCCTGCTCTAGTCCCCT
     - Probe: TGGCTAACCTCATAACAAACGT
     - Product size: 143 bp
     - No cross-reactivity with any non-target species

2. From CONSENSUS_1, primer set #1:
   - Forward: ACGTGTAGTTCGGAAAGAGTTC
   - Reverse: GCAAACCGCGAACAAAGAAA
   - Probe: TTCCAAAGGAAATCGAATTCGA
   - Product size: 96 bp
   - High specificity with only a single match

3. From CONSENSUS_3, primer set #3:
   - Forward: TCTGCAATCTTCGGGTAGGA
   - Reverse: TCACAGGCGAATATGGACTCT
   - Probe: TGGCGAATGAAACCAAAAGAGA
   - Product size: 117 bp
   - Good specificity with minimal cross-reactivity

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 3377

3. Run each module in sequence:
   - In Module 1, select the genes mentioned above with high representation
   - In Module 2, create at least 9 alignment groups
   - In Module 3, use appropriate amplicon size ranges (typically 70-150 bp)
   - In Module 4, specify TaxID 3377 for specificity testing

4. Review the HTML reports in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- DNA extracted from *W. mirabilis* plant tissues
- DNA from related gymnosperm species
- Environmental samples from habitats where the species occurs

For qPCR validation, recommended conditions:
- Annealing temperature: 58-60°C (based on primer Tm values)
- Cycling conditions: 95°C for 15s, 60°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *Welwitschia mirabilis* demonstrates the value of targeting unique mitochondrial and chloroplast genes for a phylogenetically distinct species. The CONSENSUS_4 primers (targeting the cox3 gene) showed perfect specificity, likely due to the evolutionary distinctiveness of *W. mirabilis* as the sole member of its family Welwitschiaceae.

For conservation genetics and biodiversity monitoring, these primers could be valuable for non-invasive detection of this protected species from environmental samples. The short amplicon sizes make these primers particularly suitable for applications involving degraded DNA samples or herbarium specimens, which is especially valuable for this rare and protected species.
