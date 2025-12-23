# A Comparison between CARLIN and DNA Typewriter in CRISPR-mediated Lineage Tracing

## Authors
Fengshuo Liu¹†, Xiang Zhang²†, Yipeng Yang³†*

† These authors contributed equally to this work  
** Corresponding author

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

## Citation
Manuscript under preparation

---

## Competing Interests
The authors declare no competing interests.
