# Daubentonia madagascariensis Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Daubentonia madagascariensis* (aye-aye), a nocturnal primate endemic to Madagascar with unique morphological features.

## Overview

*Daubentonia madagascariensis* (TaxID: 31869) is the only extant member of the family Daubentoniidae and genus Daubentonia. This example demonstrates how to design specific primers that can identify *D. madagascariensis* DNA in mixed samples for conservation research, evolutionary studies, or forensic applications.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Daubentonia madagascariensis/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 31869_V1R.fasta
│   ├── 31869_V1R_blast_results.xml
│   ├── 31869_gene_stats.txt
│   └── ...
├── 2_                         # Module 2 output (Multiple Sequence Alignment)
│   ├── BLAST_0.txt
│   ├── CONSENSUS_0.fasta
│   ├── MAFFT_0.txt
│   └── ...
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
│   └── ...
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
│   └── ...
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 31869 (*Daubentonia madagascariensis*). The analysis identified 50 unique genes represented in the database as shown in the `31869_gene_stats.txt` file. The most abundant genes retrieved were:

1. V1R (43 sequences) - Vomeronasal receptor type 1
2. Various genomic loci (WCI35_* genes) with 2 sequences each

For this example, the V1R gene was selected for alignment, which encodes vomeronasal receptors involved in pheromone detection. The V1R gene family has undergone extensive diversification in primates and provides a good molecular marker for species identification.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. One alignment group was created, resulting in a consensus sequence:

1. `CONSENSUS_0.fasta`: A consensus sequence derived from 41 different V1R gene sequences.

The alignment file (MAFFT_0.txt) contains the full multiple sequence alignment of V1R genes, showing all 41 different sequences in the alignment. The BLAST_0.txt file shows the alignment in BLAST-like format for easier visual inspection.

The MAFFT alignment reveals several conserved regions across the V1R genes that are suitable for primer design, despite the overall sequence variability typical of this gene family.

### Module 3: PCR Primer Design

Using the consensus sequence, Module 3 designed primer-probe sets for qPCR. The log file shows that 10 primer sets were generated for the consensus sequence. The sequence contained 15.9% ambiguous bases (N), which was noted in the log file but did not prevent successful primer design.

The CSV file contains these fields for the 10 primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_0_primers.csv` (primer set #1):
- Forward primer: GTTTCCCTCAGCACCACCT
- Reverse primer: CAGCAGAGGAAACAGCAGAA
- Probe: TCAGGCCATTAAGCTTTGCC
- Product size: 128 bp

These primers target different regions of the V1R gene consensus sequence, providing multiple options for species detection.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *D. madagascariensis* (TaxID: 31869) from related species.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows that all primer sets had high match counts:
- Primer set #1: 824 matching sequences
- Primer set #2: 830 matching sequences
- Primer set #3: 824 matching sequences
- Primer set #4: 830 matching sequences
- Primer set #5: 824 matching sequences
- Primer set #6: 829 matching sequences
- Primer set #7: 839 matching sequences
- Primer set #8: 839 matching sequences
- Primer set #9: 878 matching sequences
- Primer set #10: 832 matching sequences

The specificity report (HTML file) provides detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to *D. madagascariensis* vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets are recommended for *Daubentonia madagascariensis*:

1. Primer set #1:
   - Forward: GTTTCCCTCAGCACCACCT
   - Reverse: CAGCAGAGGAAACAGCAGAA
   - Probe: TCAGGCCATTAAGCTTTGCC
   - Product size: 128 bp

2. Primer set #3:
   - Forward: GTTTCCCTCAGCACCACCT
   - Reverse: CCAGCAGAGGAAACAGCAGA
   - Probe: TCAGGCCATTAAGCTTTGCC
   - Product size: 129 bp

While all primer sets had high match counts, these two offer good PCR properties with optimal melting temperatures, GC content, and minimal self-complementarity.

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 31869

3. Run each module in sequence:
   - In Module 1, select the V1R gene
   - In Module 2, create alignment groups for the V1R sequences
   - In Module 3, use default amplicon size range and generate 10 primer sets
   - In Module 4, specify TaxID 31869 for specificity testing

4. Review the HTML report in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated in the laboratory using:
- DNA samples from *D. madagascariensis*
- DNA from closely related lemur species
- Mixed environmental samples (when applicable)

For qPCR validation, recommended conditions:
- Annealing temperature: 59°C
- Cycling conditions: 95°C for 15s, 59°C for 60s
- Probe chemistry: TaqMan with FAM/BHQ1 labels
- Standard curve preparation: 10-fold serial dilutions of genomic DNA

## Additional Notes

The primer design for *Daubentonia madagascariensis* demonstrates the utility of vomeronasal receptor genes (V1R) as genetic markers for species identification. The V1R gene family has undergone significant expansion and diversification in primates, making it a valuable target for distinguishing closely related species.

Due to the high number of matches found during specificity testing, these primers might cross-react with V1R genes from other species, particularly other primates. For applications requiring absolute specificity, consider combining these primers with other genetic markers or using post-PCR melt curve analysis to differentiate closely related sequences.

For conservation genetic studies of the endangered aye-aye, these primers could be useful for non-invasive sampling methods, such as analyzing DNA from hair or fecal samples collected in the field.
