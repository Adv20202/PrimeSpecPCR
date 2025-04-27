# Aviadenovirus Example

This directory contains a complete example of using PrimeSpecPCR to design species-specific primers for *Aviadenovirus*, a genus of double-stranded DNA viruses that cause significant diseases in avian hosts.

## Overview

*Aviadenovirus* (TaxID: 10552) includes viral pathogens that affect both wild and domestic birds, causing diseases such as inclusion body hepatitis, hydropericardium syndrome, and egg drop syndrome. This example demonstrates the design of specific primers that can detect Aviadenovirus species in diagnostic samples, which is critical for disease surveillance and control programs in poultry.

## Directory Structure

The example directory follows the standard PrimeSpecPCR workflow with outputs from all four modules:

```
examples/Aviadenovirus/
├── 1_                         # Module 1 output (Genetic Sequence Retrieval)
│   ├── 10552_DNA_POLYMERASE.fasta
│   ├── 10552_FIBER_PROTEIN.fasta
│   ├── 10552_L3.fasta
│   ├── 10552_gene_stats.txt
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
│   └── MAFFT_2.txt
├── 3_                         # Module 3 output (PCR Primer Design)
│   ├── CONSENSUS_0_primers.csv
│   ├── CONSENSUS_1_primers.csv
│   └── CONSENSUS_2_primers.csv
├── 4_                         # Module 4 output (Primer Specificity Testing)
│   ├── specificity_report_CONSENSUS_0_primers.html
│   ├── specificity_report_CONSENSUS_1_primers.html
│   └── specificity_report_CONSENSUS_2_primers.html
├── 3_PCR_Primers_Design.log   # Log file from Module 3
├── 4_Primers_Specificity.log  # Log file from Module 4
├── PCR_primer_settings.txt    # Primer design parameters
└── primer_specificity_settings.txt  # Specificity test parameters
```

## Step-by-Step Workflow Results

### Module 1: Genetic Sequence Retrieval

Module 1 retrieved sequences from NCBI GenBank using TaxID 10552 (*Aviadenovirus*). The analysis identified 50 unique genes represented in the database as shown in the `10552_gene_stats.txt` file. The most abundant genes retrieved were:

1. HEXON (3664 sequences)
2. HEXON_PROTEIN (1507 sequences)
3. FIBER (645 sequences)
4. DNA_POLYMERASE (401 sequences)
5. PROTEASE (368 sequences)

For this example, three genes with sufficient sequence representation were selected for alignment:
- DNA_POLYMERASE
- FIBER_PROTEIN
- L3

These genes were selected due to their conservation across Aviadenovirus species and their importance in viral replication and structure. The DNA polymerase is essential for viral genome replication, the fiber protein is involved in host cell recognition and attachment, and L3 is involved in capsid formation.

### Module 2: Multiple Sequence Alignment

The sequences from Module 1 were grouped and aligned using MAFFT to identify conserved regions. Three alignment groups were created, resulting in three consensus sequences:

1. `CONSENSUS_0.fasta`: A 986 bp consensus sequence containing 15.1% ambiguous bases.
2. `CONSENSUS_1.fasta`: A 1197 bp consensus sequence with minimal ambiguity.
3. `CONSENSUS_2.fasta`: A 761 bp consensus sequence containing 10.2% ambiguous bases.

The alignment files (MAFFT_*.txt) contain the full multiple sequence alignments, while the BLAST_*.txt files show the alignments in BLAST-like format for easier visual inspection of conserved regions across viral strains.

### Module 3: PCR Primer Design

Using the consensus sequences, Module 3 designed primer-probe sets for qPCR detection. The log file shows that amplicon ranges of 70-200 bp were selected for all three consensus sequences, and 8 primer sets were generated for each consensus sequence.

Each CSV file contains these fields for the primer sets:
- Primer sequences (forward, reverse, probe)
- Melting temperatures (Tm)
- GC content percentages
- Self-complementarity scores
- Hairpin formation scores
- Expected product size

Example from `CONSENSUS_0_primers.csv` (primer set #1):
- Forward primer: CGGGCAAAACAGTACAACCA
- Reverse primer: GTCTCGTTGGGAGGAACGTG
- Probe: CGGTCGTCTGCCCCGAAAGC
- Product size: 78 bp

For the CONSENSUS_0 sequence, the initial primer design with standard settings failed to generate primers, likely due to the high number of ambiguous bases. The system automatically switched to a simplified approach, which successfully generated 8 primer sets.

### Module 4: Primer Specificity Testing

Module 4 evaluated the specificity of each primer set against GenBank using BLAST, with special focus on discriminating *Aviadenovirus* (TaxID: 10552) from related viruses and host genetic material.

Key specificity settings used:
- Maximum BLAST hits: 2000
- Minimum identity: 60.0%
- Word size: 5
- Maximum total mismatches: 3
- Maximum 3' end mismatches: 2

The log shows varying levels of specificity among the primer sets:
- CONSENSUS_0 primers showed excellent specificity with no matches found outside the target group
- CONSENSUS_1 primers showed variable specificity with primer sets #2 and #3 matching 27 sequences
- CONSENSUS_2 primers showed consistent matches (118 sequences) for most primer sets, with primer set #7 showing no matches

The specificity reports (HTML files) provide detailed analyses of each primer set including:
- Total number of matching sequences
- Matches to Aviadenovirus vs. non-target sequences
- Visualization of primer binding sites
- Mismatch patterns at critical positions

## Recommended Primer Sets

Based on the specificity testing, the following primer sets showed the best specificity for *Aviadenovirus*:

1. From CONSENSUS_0, primer set #1:
   - Forward: CGGGCAAAACAGTACAACCA
   - Reverse: GTCTCGTTGGGAGGAACGTG
   - Probe: CGGTCGTCTGCCCCGAAAGC
   - Product size: 78 bp
   - No cross-reactivity with non-target sequences

2. From CONSENSUS_2, primer set #7:
   - Forward: TCAGTCAACTAGCATCCCGG
   - Reverse: TGCTGTCCATGACCCAGTAA
   - Probe: GACGTTTCGCCAAGGTAGAC
   - Product size: 186 bp
   - No cross-reactivity with non-target sequences

These primer sets would be particularly useful for specific detection of Aviadenovirus species in mixed samples or clinical specimens.

## Usage Instructions

To replicate this example:

1. Start PrimeSpecPCR by running `python run.py`

2. Configure with your NCBI credentials and TaxID 10552

3. Run each module in sequence:
   - In Module 1, select the three genes mentioned above
   - In Module 2, create at least 3 alignment groups
   - In Module 3, use amplicon size range 70-200 bp
   - In Module 4, specify TaxID 10552 for specificity testing

4. Review the HTML reports in the 4_ directory to identify the most specific primer sets for your application

## Laboratory Validation

These primer sets should be experimentally validated using:

- **Viral samples**: Reference strains of Aviadenovirus propagated in appropriate cell cultures (e.g., chicken embryo liver cells, chicken embryo fibroblasts) or from embryonated chicken eggs
- **Clinical specimens**: Liver, spleen, kidney, and intestinal samples from birds with suspected aviadenovirus infection
- **Environmental samples**: Fecal samples, cloacal swabs, or drinking water from poultry houses
- **Cross-reactivity panel**: Other avian viruses (e.g., avian influenza virus, infectious bronchitis virus) to confirm specificity

For qPCR validation protocol:
- **Nucleic acid extraction**: Use commercial viral DNA extraction kits optimized for tissues and fecal samples
- **Reaction setup**: Prepare reactions in dedicated pre-PCR area to prevent contamination
- **Cycling conditions**: Initial denaturation at 95°C for 5 min, followed by 40 cycles of 95°C for 15s, 60°C for 60s
- **Probe chemistry**: TaqMan probes with FAM reporter and BHQ1 quencher
- **Controls**:
  - Positive controls: Plasmids containing target sequences or DNA extracted from reference virus
  - Negative controls: No-template controls to detect contamination
  - Internal amplification controls: To identify PCR inhibition in clinical samples
- **Validation parameters**: Determine analytical sensitivity (limit of detection), specificity, repeatability, and reproducibility
- **Field validation**: Test primers on samples from multiple geographic locations to account for viral strain variation

For comprehensive validation, compare results with established diagnostic methods like virus isolation or commercial PCR kits if available.

## Additional Notes

The primer design for *Aviadenovirus* demonstrates the challenge of designing primers for viruses with significant genetic variability across strains. The high number of ambiguous bases in the consensus sequences reflects this diversity but also highlights the importance of focusing on conserved regions for reliable detection.

The CONSENSUS_0 sequence required the simplified primer design approach, illustrating how PrimeSpecPCR can adapt to challenging templates. The resulting primers from this approach showed excellent specificity, demonstrating that the software can successfully generate useful primers even from highly variable viral sequences.

For diagnostic applications, combining multiple primer sets targeting different genes (e.g., one from CONSENSUS_0 and one from CONSENSUS_2) may provide the most reliable detection system for Aviadenovirus infections in field conditions.
