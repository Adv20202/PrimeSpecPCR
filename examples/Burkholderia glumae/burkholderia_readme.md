# Burkholderia glumae Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Burkholderia glumae*, a bacterial pathogen that causes bacterial panicle blight in rice and other cereal crops.

## Overview

*Burkholderia glumae* (TaxID: 337) is a gram-negative bacterium that produces toxoflavin, a major virulence factor contributing to its pathogenicity. This example demonstrates how to design specific primers that can identify *B. glumae* in mixed samples containing multiple bacterial species or environmental contaminants.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Burkholderia glumae/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 337_GYRB.fasta
│   ├── 337_TOXA.fasta
│   ├── 337_gene_stats.txt
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_0.txt
│   ├── BLAST_1.txt
│   ├── CONSENSUS_0.fasta
│   ├── CONSENSUS_1.fasta
│   ├── MAFFT_0.txt
│   ├── MAFFT_1.txt
│   └── ...
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
│   ├── CONSENSUS_1_primers.csv
│   └── ...
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
│   ├── specificity_report_CONSENSUS_1_primers.html
│   └── ...
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 337 (*Burkholderia glumae*). The analysis identified 50 unique genes represented in the database as shown in the `337_gene_stats.txt` file. The most abundant genes retrieved were:

1. GYRB (327 sequences)
2. ISTA (311 sequences)
3. ISTB (307 sequences)
4. TNPB (304 sequences)
5. RPOD (287 sequences)

For this example, genes with sufficient sequence representation were selected for alignment, including:
- GYRB (DNA gyrase subunit B)
- TOXA (Toxoflavin biosynthesis gene)

These genes represent a mix of housekeeping genes (GYRB) and virulence factors (TOXA) to provide multiple options for species-specific detection. The TOXA gene is particularly useful for identifying pathogenic strains of *B. glumae*.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. Two alignment groups were created, resulting in two consensus sequences:

1. `CONSENSUS_0.fasta`: A 616 bp consensus focused on the GYRB gene, which encodes the B subunit of DNA gyrase, an essential enzyme involved in DNA replication.
2. `CONSENSUS_1.fasta`: A 683 bp consensus sequence derived from a different group of sequences.

The alignment files (MAFFT_*.txt) contain the full multiple sequence alignments, while the BLAST_*.txt files show the alignments in BLAST-like format for easier visual inspection.

### Module 3: PCR Primer Design

Using the consensus sequences, Module 3 designed primer-probe sets for qPCR. The log file shows that amplicon ranges of 100-200 bp were selected for both consensus sequences, and 5 primer sets were generated for each consensus.

Each CSV file contains these fields for the primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_0_primers.csv` (primer set #1):
- Forward primer: AAGGGCTTCGTCGAGTACAT
- Reverse primer: CGTTCTCGTTGTAGCTGTCG
- Probe: ATCATCGGTCGAGAAGGACG
- Product size: 131 bp

Example from `CONSENSUS_1_primers.csv` (primer set #1):
- Forward primer: GTGAGGAGTCGCGCAAATAC
- Reverse primer: GTAATTGAAGAGCCAGGCGG
- Probe: GGCAATTCGACCTGGTCAAT
- Product size: 113 bp

These primers target different regions of the consensus sequences, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *B. glumae* (TaxID: 337) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- CONSENSUS_0 primers had variable match counts (ranging from 1 to 487 sequences)
- CONSENSUS_1 primers showed higher specificity with fewer matches (ranging from 10 to 14 sequences)

The specificity reports (HTML files) provide detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *B. glumae* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed the best specificity for *Burkholderia glumae*:

1. From CONSENSUS_0, primer set #1:
   - Forward: AAGGGCTTCGTCGAGTACAT
   - Reverse: CGTTCTCGTTGTAGCTGTCG
   - Probe: ATCATCGGTCGAGAAGGACG
   - Product size: 131 bp
   - Found only 1 matching sequence, showing high specificity

2. From CONSENSUS_0, primer set #3:
   - Forward: AGGGCTTCGTCGAGTACATC
   - Reverse: CGTTCTCGTTGTAGCTGTCG
   - Probe: ATCATCGGTCGAGAAGGACG
   - Product size: 130 bp
   - Found only 1 matching sequence, showing high specificity

3. From CONSENSUS_1, primer set #3:
   - Forward: GTCCGCACGATTTTCCACAT
   - Reverse: GTATTTGCGCGACTCCTCAC
   - Probe: GGCGTGGATATCTCCGAGAA
   - Product size: 174 bp
   - Limited cross-reactivity with 10 matching sequences

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 337

3. Run each module in sequence:
   - In Module 1, select the genes mentioned above (GYRB, TOXA)
   - In Module 2, create at least 2 alignment groups
   - In Module 3, use amplicon size range 100-200 bp
   - In Module 4, specify TaxID 337 for specificity testing

4. Review the HTML reports in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- Pure cultures of *B. glumae*
- Closely related species (especially other *Burkholderia* species)
- Environmental samples containing mixed bacteria from rice fields

For qPCR validation, recommended conditions:
- Annealing temperature: 59°C
- Cycling conditions: 95°C for 15s, 59°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *B. glumae* demonstrates the importance of selecting appropriate genetic targets for detection of plant pathogens. The GYRB gene provided highly specific primer sets (as seen in CONSENSUS_0 primers #1 and #3) that showed minimal cross-reactivity with non-target organisms.

For applications requiring detection of pathogenic strains, consider using primers targeting the TOXA gene, which is involved in toxoflavin production and is a key virulence factor. This can help distinguish pathogenic from non-pathogenic *Burkholderia* species.

For monitoring bacterial panicle blight in field conditions, these primers could be incorporated into a multiplex PCR assay that targets multiple *B. glumae* genes simultaneously, increasing both sensitivity and specificity for early detection of this important rice pathogen.
