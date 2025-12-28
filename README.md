# A Comparison between CARLIN and DNA Typewriter in CRISPR-mediated Lineage Tracing

## Authors
Fengshuo Liu¹†, Xiang Zhang²†, Yipeng Yang³†*

† These authors contributed equally to this work  
* Corresponding author

---

## Affiliations
¹ Graduate Program in Cancer and Cell Biology,  
Baylor College of Medicine, Houston, TX, USA  

² Lester and Sue Smith Breast Center,  
Baylor College of Medicine, Houston, TX, USA  

³ Department of Mathematics and Statistics,  
University of Houston–Clear Lake, Houston, TX, USA  

---

## Corresponding Author
**Yipeng Yang**  
📧 yangy@uhcl.edu  

> This repository is maintained under **Xiang Zhang’s Lab**,  
> while correspondence regarding the manuscript should be directed to **Yipeng Yang**.

---

## Manuscript Overview
This repository accompanies the manuscript:

**“A Comparison between CARLIN and DNA Typewriter in CRISPR-mediated Lineage Tracing”**

In this study, we systematically compare two CRISPR-based lineage tracing strategies—**CARLIN** and **DNA Typewriter**—using a unified stochastic simulation framework with known ground-truth lineages. By explicitly modeling CRISPR editing dynamics, barcode evolution, and cell division processes, the framework enables quantitative benchmarking of lineage reconstruction accuracy across diverse experimental parameter regimes.

Both methods are evaluated using multiple accuracy metrics, including **Robinson–Foulds (RF) accuracy** and **triplet accuracy**, allowing a comprehensive assessment of lineage reconstruction performance under varying mutation rates, sampling depths, and lineage lengths.

---

## Major Conclusions
- **DNA Typewriter consistently outperforms CARLIN** in lineage reconstruction accuracy when sufficient numbers of recording targets are used, particularly in longer cell cycles.
- **Sequential, ordered recording** in DNA Typewriter substantially reduces ambiguity in lineage inference compared to unordered CRISPR barcode editing.
- **CARLIN’s lineage-recording potential exhausts rapidly** under continuous induction, limiting its effectiveness in long-term lineage tracing.
- **Target number is a dominant determinant of accuracy** in DNA Typewriter, whereas tandem length and insBC diversity play secondary roles.
- **Triplet accuracy provides a more permissive and informative metric** than strict RF accuracy, especially under partial sampling scenarios.
- Hybrid (dual) reconstruction approaches combining CARLIN and DNA Typewriter do **not** consistently outperform DNA Typewriter alone, due to competing lineage signals.

Together, this work provides **quantitative guidance for experimental design** and highlights the advantages of ordered molecular recording systems for high-resolution lineage tracing.

---

## Data and Code Availability
Simulation code and data supporting this study are provided in this repository and correspond to the analyses presented in the manuscript.

---

## Data Generation Code (MATLAB) — Script Guide

A separate folder (to be uploaded) contains the MATLAB scripts used to **simulate barcode evolution**, **generate datasets**, **reconstruct lineage trees**, and **compute accuracy metrics** (RF accuracy and triplet accuracy). Below is a brief guide explaining what each file does.

### Main entry point
- `MCnew.m`  
  Main entrance script that runs the simulation/reconstruction pipeline and produces the outputs used for downstream analyses.  
  **Note:** Lines calling `tripaccu` can be commented out to skip **triplet accuracy** calculation (useful for faster runs).

### Accuracy metrics
- `tripaccu.m`  
  Computes **triplet accuracy** for a reconstructed lineage tree, using the **simulated (ground-truth) lineage tree** as the benchmark.

- `findm.m`  
  Finds the number of matched entries between two matrices, used in the calculation of **RF accuracy**.

### Core simulation + reconstruction
- `funbarNBJNF_TW.m`  
  Core simulator/driver that:
  1) simulates **CARLIN barcode** evolution and **DNA Typewriter (DNA Tape)** evolution,  
  2) constructs lineage trees from simulated data (ground truth and/or reconstructed trees), and  
  3) computes **Robinson–Foulds (RF) accuracies** for benchmarking.

- `btree.m`  
  Builds the original/benchmark (ground-truth) lineage tree produced by the simulation.

- `valcheck.m`  
  Validates a rebuilt lineage tree by checking structural consistency, including the **number of cells per generation**.

### DNA Typewriter (DNA Tape) allele utilities
- `twalleles.m`  
  Generates a **child DNA Tape allele** from a **parent DNA Tape allele**, modeling the inheritance/editing process.

- `twscore.m`  
  Computes a **matching score** between two DNA Tape alleles (used for comparison/alignment/scoring during reconstruction).

- `twparent.m`  
  Rebuilds an inferred **parent DNA Tape** state from two children DNA Tape alleles (used in tree reconstruction steps).

### CARLIN barcode utilities
- `genchildcarlin.m`  
  Generates a **child CARLIN barcode** from a **parent barcode**, modeling editing/inheritance during cell division.

- `proend.m`  
  Processes a cut end on a barcode due to **DSB (double-strand break)**, implementing the barcode editing outcome logic.

- `CARLIN_raw`  
  Contains the **CARLIN reference barcode** used by the simulator (reference/template barcode data for CARLIN).

### Barcode alignment / vector comparison helpers
These helper functions support barcode comparison/alignment routines used by the reconstruction logic.

- `vcomp_cons.m`  
  Compares similarity between two vectors by **rewarding consecutive matches** and applying **less penalty for consecutive mismatches**.  
  Used for barcode alignment.

- `vcomp5.m`  
  Compares similarity between two vectors (including `0` entries).  
  Used for barcode alignment.

- `trim.m`  
  Trims leading and trailing zeros from a vector.  
  Used as a preprocessing step for barcode alignment.

### Visualization
- `redrawHeatmap.m`  
  Plots a heatmap of **DNA Tape alleles** (used for visualizing allele distributions/patterns).

---

## Citation
Manuscript under preparation

---

## Competing Interests
The authors declare no competing interests.
