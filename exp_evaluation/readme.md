# PrimeSpecPCR Experimental Validation

This directory contains the experimental validation of PrimeSpecPCR-designed primers using biological samples.

> **Note**: The complete methodology and detailed results will be published in a scientific article. The reference will be added here once the publication is available.

## What Was Tested

Five different organisms were selected to test primers on real samples collected from agricultural fields near Lublin, Poland:

### Plant Pathogens (Fungi)
- **Blumeria graminis f. sp. tritici** - causes powdery mildew on wheat
- **Blumeria hordei** - causes powdery mildew on barley  
- **Zymoseptoria tritici** - causes septoria leaf blotch on wheat

### Plants
- **Capsella bursa-pastoris** - shepherd's purse (a common weed)
- **Equisetum arvense** - field horsetail (prehistoric-looking plant)

## Sample Collection

### Diseased Plant Material
Samples were collected from wheat and barley fields by looking for plants showing obvious disease symptoms. For powdery mildew, leaves covered in white, powdery fungal growth were collected. For septoria, wheat leaves with the characteristic brown spots and yellowing caused by this disease were selected, etc.

All samples were photographed to document exactly what was collected and ensure correct identification.

### Healthy Plant Samples
For the two plant species, fresh material was collected based on their distinctive features - shepherd's purse has small white flowers and a rosette growth pattern, while horsetail has segmented stems.

## Laboratory Workflow

### DNA Extraction
A commercial plant/fungi DNA extraction kit was used to obtain clean DNA from all samples. This kit works for both the fungal pathogens growing on plants and the plant tissue itself.

### Primer Selection
After running the full PrimeSpecPCR pipeline on each organism, the results were analyzed and the primer sets that scored best for specificity were selected. A total of 7 different primer sets were tested.

### Primer Comparison with Primer3
The same consensus sequences were also processed through Primer3 (the gold standard primer design software) using identical settings. This allowed comparison of whether PrimeSpecPCR found different or better primers than the established tool.

### PCR Testing
For each primer set, two different PCR reactions were performed:
1. **Standard PCR**: Using the forward and reverse primers as intended
2. **Alternative PCR**: Using the probe sequence as a forward primer with the reverse primer

All reactions were run in duplicate.

### Gel Electrophoresis
All PCR products were run on agarose gels to confirm bands of the expected size.

### DNA Sequencing
Getting a PCR band is good, but absolute confirmation required sequencing. All PCR products were sent for Sanger sequencing to read the actual DNA sequence.

### Sequence Analysis
The sequences were run through BLAST against GenBank to confirm they really came from the target organisms.

## Results Summary

### The Good News
**Every single primer set worked!** All 7 primer sets successfully amplified DNA from their target organisms and produced clean bands on gels at the expected sizes.

### Sequence Confirmation
When the PCR products were sequenced and searched against GenBank, they all came back as perfect matches to the target species. No contamination, no off-target amplification.

### PrimeSpecPCR vs Primer3
Here's where it gets interesting: 6 out of 7 primer sets that worked were identical between PrimeSpecPCR and Primer3. But one primer set (for Capsella bursa-pastoris) was unique to PrimeSpecPCR - Primer3 didn't find this particular combination. This shows that the software can sometimes find good primers that other tools miss.

### Reproducibility
All duplicate PCR reactions gave identical results.

## What This Means

This validation proves that PrimeSpecPCR doesn't just look good on paper - it actually designs primers that work in real laboratory conditions with biological samples. The software successfully:

- Identified conserved regions suitable for primer binding
- Designed primers with good thermodynamic properties
- Predicted specificity that held up under experimental testing
- Found at least one primer set that established tools missed

## Data Organization

Each organism's folder contains:

**Pipeline Outputs (1_ through 4_/)**
Complete results from running PrimeSpecPCR, including sequence retrieval, alignment, primer design, and specificity testing.

**PCR_verification/**
Photos of agarose gels showing successful PCR amplification. Each gel clearly shows bands at the expected molecular weights.

**Primer3_verification/**  
Results from running the same analysis with Primer3 software for comparison.

**Sanger_sequencing/**
Raw sequencing files (.ab1 format) and sequence data (.seq files) from all PCR products.

**Sequence_analysis/**
BLAST search results confirming the identity and specificity of sequenced PCR products.

## Important Notes

### What Wasn't Tested
This validation focused on basic PCR functionality and specificity. The following were not tested:
- Quantitative PCR performance (efficiency, sensitivity)
- Probe functionality for real-time PCR
- Cross-reactivity with closely related species in the same reaction
- Performance under different PCR conditions or cycling parameters

### Limitations
Success in this validation doesn't guarantee that every primer set designed by PrimeSpecPCR will work perfectly. PCR success depends on many factors including:
- DNA quality and concentration
- PCR reaction conditions
- Thermal cycling parameters
- Presence of inhibitors in samples
- and more!

### Future Validation
While these results are encouraging, more extensive validation with additional organisms and under different conditions would strengthen confidence in the software's reliability.


---

*For detailed methodology, statistical analysis, and complete experimental protocols, please refer to the associated scientific publication [reference will be added upon publication].*
