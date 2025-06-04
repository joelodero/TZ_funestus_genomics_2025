## Genomic Evidence of Spatially Structured Gene Flow and Divergent Insecticide Resistance Backgrounds of the Malaria Vector Anopheles funestus in Tanzania

Population genomic analysis of 334 Tanzanian An. funestus samples.
The code notebooks and scripts in this repo are in support of this analysis: <Publication link>

Anopheles funestus is the dominant malaria vector in Tanzania and most parts of East and Southern Africa. However, studies on its population structure are limited. Such studies would help understand the distribution of insecticide resistance alleles, devise sustainable insecticide-based vector control approaches, and determine how malaria vector populations are structured in space. 

We analysed the whole genome sequences of 334 An. funestus mosquitoes sampled from 11 regions with varying ecologies and malaria burdens across mainland Tanzania. Signs of reduced migration between western and eastern cohorts across the semi-arid central region containing the Rift Valley suggest a partial barrier to gene flow between these populations. This was evidenced by population structure between the eastern and western cohorts, as well as asynchronous selective sweeps and copy number variant profiles at the Cyp9k1 gene, and Cyp6p gene cluster. Eastern cohorts, despite having less diversity and greater inbreeding, also share genetic histories characterised by low genome-wide Fst values with those in the west. This suggests that the barrier to gene flow is porous and likely represents a continuous spatial structure rather than a complete barrier to migration. 

This repository contains the code used to analyse the data generated in this analysis. The R-markdown contains:
1. Script for generating the ADMIXTURE plots.
2. Script for generating the runs of homozygosity (ROH) plot.

The Jupyter notebook contains:
1. The PCA plots on all chromosome X and 2RL, with 2RL revealing inland and coastal clusters
2. Effective migration surfaces using fEEMS
3. Genome-wide population diversity - nucleotide diversity (π) and Tajima’s D
4. Changes in effective population size (Ne) over time by using Stairwayplot2.
5. Between-population differentiation (Fst)
6. Between-sample relatedness using PC-Relate implemented in sgkit.
7. Genome-wide selection scans (GWSS) with the H12 statistic 
8. Diplotype clustering dendrograms for the genes Cyp6p9 and Cyp9k1
9. Amino acid variation around the Cyp6p9 and Cyp9k1 genes
10. Copy number variation around the Cyp6p gene cluster (Rp1) and Cyp9k1 genes

If you'd like to reproduce the analysis, please clone this repo and run the Jupyter notebook first. Make sure you can access malariagen_data](https://malariagen.github.io/vector-data/vobs/vobs-data-access.html). The notebook should run and download the data required to plot the rest of the analyses in the RMarkdown.

Joel Odero, June 2025
