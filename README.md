## Distinct Genetic Populations and Resistance Backgrounds of the Malaria Vector Anopheles funestus in Tanzania

Population genomic analysis of 334 Tanzanian An. funestus samples.
The code notebooks and scripts in this repo are in support of this analysis: <BIORXiv link>

Anopheles funestus is the dominant malaria vector in Tanzania and most parts of East and Southern Africa. However, studies on its population structure are limited. Such studies would help understand the distribution of insecticide resistance alleles, devise sustainable insecticide-based vector control approaches, and determine how malaria vector populations are structured in space. 

We analysed the whole genome sequences of 334 An. funestus mosquitoes sampled from 11 administrative regions with varying ecologies and malaria burdens across mainland Tanzania. We found not only marked population structure within and between populations of An. funestus, potentially determined by the Rift Valley, but also that the distribution of metabolic resistance determinants varied according to population, suggesting divergent genomic architectures underlying a consistent resistance phenotype across Tanzania.

This repository contains the code used to analyse the data generated in this analysis. The R-markdown contains:
1. Script for generating the ADMIXTURE plots.
2. Script for generating the runs of homozygosity (ROH) plot.

The jupyter notebook contains:
1. The PCA plots on all chroms arm with 2RL revealing inland and coastal clusters
2. Genome-wide population diversity - nucleotide diversity (π) and Tajima’s D
3. Between-population differentiation (Fst)
4. Genome-wide selection scans (GWSS) with the H12 statistic 
5. Haplotype clustering dendrograms for the genes Cyp6p9 and Cyp9k1
6. Amino acid variation around the Cyp6p9 and Cyp9k1 genes
7. Copy number variation around the Cyp6p gene cluster (Rp1) and Cyp9k1 genes

If you wish to replicate the analysis, clone this repo, and run the jupyter notebook first. Make sure [you can access malariagen_data](https://malariagen.github.io/vector-data/vobs/vobs-data-access.html). The notebook should run and download the data required to plot the rest of the analyses in the RMarkdown.

Joel Odero, Jan 2025
