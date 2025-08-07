This repository contains scripts that have been used during the BlendARGs project with the aim to investigate the vertical dissemination of antibiotic resistance genes (ARGs) by horizontal gene transfers (HGTs) in both various marine and freshwater environments. 

Publicly available metagenomic data (Freshwater: [Buck et al. 2021](https://www.nature.com/articles/s41597-021-00910-1), and Marine: [Biller et al. 2018](https://www.nature.com/articles/sdata2018176)) were processed by using [nf-core/mag](https://nf-co.re/mag/2.5.1/) Nextflow pipeline (for details, see multiQC summary reports). The resulted bins (>= 40 % completeness) were used as input data for the following analyses:

1. [MetaCHIP](https://github.com/songweizhi/MetaCHIP) to detect putative HGT events and the involved genes across sampled water depths,
2. Functional annotation of genes using [DeepNOG](https://github.com/univieCUBE/deepnog) (db: eggNOG5),
3. The identification of ARGs by using:
   (a) [RGI](https://github.com/arpcard/rgi) tool with the CARD database,
   (b) BLASTp against the [ResFinderFG2.0 database](https://github.com/RemiGSC/ResFinder_FG_Construction), and
   (c) latent ARG database compiled by [Salvador et al. (2023)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-023-01479-0).
4. Identification of mobile genetic elements (MGEs) with [geNomad](https://portal.nersc.gov/genomad/index.html) and calculation of their nucleotide-level genomic similarities with [FastANI](https://github.com/ParBLiSS/FastANI).
5. Virus lifestyle predictions was done with the PhaTYP module of [PhaBOX](https://phage.ee.cityu.edu.hk/).
6. Site-wise linking of MGEs with MAG-ARGs using BLASTp.



<img width="4656" height="3895" alt="Flowchart_new" src="https://github.com/user-attachments/assets/ec74219b-1e53-45b6-a3fe-4ec6047cf9e3" />




The results of the project are presented in our preprint: ....
