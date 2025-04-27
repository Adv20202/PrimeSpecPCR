# Erwinia psidii Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Erwinia psidii*, a gram-negative bacterial plant pathogen that causes bacterial blight in guava and other Myrtaceae species.

## Overview

*Erwinia psidii* (TaxID: 69224) is a phytopathogenic bacterium that affects economically important fruit crops, particularly guava. This example demonstrates how to design specific primers that can identify *E. psidii* in mixed samples containing multiple bacterial species for rapid diagnosis and monitoring of plant diseases.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Erwinia psidii/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 69224_16S_RIBOSOMAL_RNA.fasta
│   ├── 69224_GAPA.fasta
│   ├── 69224_GYRB.fasta
│   ├── 69224_RECA.fasta
│   ├── 69224_RPOB.fasta
│   ├── 69224_gene_stats.txt
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_0.txt
│   ├── BLAST_1.txt
│   ├── BLAST_2.txt
│   ├── BLAST_3.txt
│   ├── BLAST_4.txt
│   ├── CONSENSUS_0.fasta
│   ├── CONSENSUS_1.fasta
│   ├── CONSENSUS_2.fasta
│   ├── CONSENSUS_3.fasta
│   ├── CONSENSUS_4.fasta
│   ├── MAFFT_0.txt
│   ├── MAFFT_1.txt
│   ├── MAFFT_2.txt
│   ├── MAFFT_3.txt
│   └── MAFFT_4.txt
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
│   ├── CONSENSUS_1_primers.csv
│   ├── CONSENSUS_2_primers.csv
│   ├── CONSENSUS_3_primers.csv
│   └── CONSENSUS_4_primers.csv
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
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

Module 1 retrieved sequences from NCBI GenBank using TaxID 69224 (*Erwinia psidii*). The analysis identified 50 unique genes represented in the database as shown in the `69224_gene_stats.txt` file. The most abundant genes retrieved were:

1. RPOB (69 sequences)
2. RECA (65 sequences)
3. GAPA (54 sequences)
4. TNPA (35 sequences)
5. TSSK (27 sequences)

For this example, five genes with sufficient sequence representation were selected for alignment:
- 16S_RIBOSOMAL_RNA
- GAPA (Glyceraldehyde-3-phosphate dehydrogenase A)
- GYRB (DNA gyrase subunit B)
- RECA (Recombinase A)
- RPOB (RNA polymerase subunit beta)

These genes represent a mix of housekeeping genes (GAPA, GYRB, RECA, RPOB) and ribosomal RNA genes (16S_RIBOSOMAL_RNA) that are commonly used for bacterial identification and phylogenetic analysis, providing multiple options for species-specific detection.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. Five alignment groups were created, resulting in five consensus sequences:

1. `CONSENSUS_0.fasta`: Focused on the 16S ribosomal RNA gene, a 746 bp consensus sequence was generated.
2. `CONSENSUS_1.fasta`: A 410 bp consensus sequence derived from GAPA gene sequences.
3. `CONSENSUS_2.fasta`: A 742 bp consensus sequence based on recA gene sequences.
4. `CONSENSUS_3.fasta`: A 623 bp consensus sequence from another group of genes.
5. `CONSENSUS_4.fasta`: A 637 bp consensus sequence from the rpoB gene.

The alignment files (MAFFT_*.txt) contain the full multiple sequence alignments, while the BLAST_*.txt files show the alignments in BLAST-like format for easier visual inspection.

### Module 3: PCR Primer Design

Using the consensus sequences, Module 3 designed primer-probe sets for qPCR. The log file shows that 10 primer sets were generated for CONSENSUS_0, CONSENSUS_1, CONSENSUS_2, and CONSENSUS_3, while 20 primer sets were generated for CONSENSUS_4.

Each CSV file contains these fields for the primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_0_primers.csv` (primer set #1):
- Forward primer: CGGTGGAGCATGTGGTTTAA
- Reverse primer: ACAAAGGATAAGGGTTGCGC
- Probe: CGTCAGCTCGTGTTGTGAAA
- Product size: 188 bp

Example from `CONSENSUS_1_primers.csv` (primer set #1):
- Forward primer: TGTTGACGTGGTTGCTGAAG
- Reverse primer: GATTTCCTGGCCAGCGTATG
- Probe: ACGTAAGCACATTGAAGCCG
- Product size: 166 bp

These primers target different regions of the consensus sequences, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *E. psidii* (TaxID: 69224) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- CONSENSUS_0 primers: Primer pairs #1, #2, #3, #6, #7, #8, #9, and #10 had no matching sequences, while pairs #4 and #5 matched 220 and 219 sequences respectively.
- CONSENSUS_1 primers had moderate match counts (ranging from 16 to 22 sequences)
- CONSENSUS_2 primers were very specific, with most primer sets having only 1 or 4 matching sequences
- CONSENSUS_3 primers showed variable specificity (ranging from 1 to 141 matches)
- CONSENSUS_4 primers had varying match counts (ranging from 0 to 215 sequences)

The specificity reports (HTML files) provide detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *E. psidii* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed the best specificity for *Erwinia psidii*:

1. From CONSENSUS_2, primer sets #1, #2, #3, and #6:
   - Example (set #1):
     - Forward: TTGCGGTTGTGTCAGTGAAG
     - Reverse: ACGGCAGATTTCACCTCAGA
     - Probe: GCCGGACCCTAAATTCTCCT
     - Product size: 88 bp
     - Only 1 matching sequence, showing high specificity

2. From CONSENSUS_3, primer sets #1, #2, #4, and #9:
   - Example (set #1):
     - Forward: GTCCCTGGATATTGCGCTTG
     - Reverse: GCTCCGCATCAATAAAGGCA
     - Probe: ACCACGCTGACTTTACAGGT
     - Product size: 155 bp
     - Only 1 matching sequence, showing high specificity

3. From CONSENSUS_4, primer sets #8, #13, and #15:
   - Example (set #8):
     - Forward: CGCTGTCTGAGATTACGCAC
     - Reverse: TACTCGTTGGTCTGTGCGTA
     - Probe: GTGCAGGCTTTGAAGTTCGA
     - Product size: 193 bp
     - No matching sequences outside the target organism

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 69224

3. Run each module in sequence:
   - In Module 1, select the five genes mentioned above
   - In Module 2, create at least 5 alignment groups
   - In Module 3, use default amplicon size ranges
   - In Module 4, specify TaxID 69224 for specificity testing

4. Review the HTML reports in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- Pure cultures of *E. psidii*
- Closely related *Erwinia* species
- Plant samples infected with bacterial blight
- Environmental samples from agricultural settings

For qPCR validation, recommended conditions:
- Annealing temperature: 59°C
- Cycling conditions: 95°C for 15s, 59°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *Erwinia psidii* demonstrates the challenge of designing species-specific primers within the genus *Erwinia*, which contains many closely related plant pathogenic species. The housekeeping genes used in this example (particularly recA and rpoB) provided excellent specificity for distinguishing *E. psidii* from related species.

For field diagnostics, these primers could be incorporated into portable qPCR systems for on-site detection of the pathogen in plant material. The high specificity primers from CONSENSUS_2 and CONSENSUS_3 would be particularly valuable for regulatory purposes, where false positives must be minimized.

For applications requiring absolute specificity, consider using multiple primer sets in combination, targeting different genes, for confirmation of positive results.
