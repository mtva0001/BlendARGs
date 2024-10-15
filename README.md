This repository contains scripts that have been used during the BlendARGs project with the aim to investigate the dissemination of antibiotic resistance genes (ARGs) by horizontal gene transfers (HGTs) in both marine and freshwater environments. 

Publicly available metagenomic datasets were processed by nf-core/mag Nextflow pipeline (for details, see multiQC summary reports). The resulted bins (>= 40 % completeness) were used as input data for the following analyses:

1. [MetaCHIP](https://github.com/songweizhi/MetaCHIP) to detect putative HGT events and the involved genes across sampled water depths,
2. Functional annotation of genes using [DeepNOG](https://github.com/univieCUBE/deepnog) (db: eggNOG5),
3. The detection of ARGs by:
   (a) using [RGI](https://github.com/arpcard/rgi) with the CARD database,
   (b) using BLASTp against the [ResFinderFG2.0 database](https://github.com/RemiGSC/ResFinder_FG_Construction), and
   (c) using latent ARG database compiled by [Inda-Díaz et al. (2023)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-023-01479-0).

The results of the project are presented in our preprint: ....
