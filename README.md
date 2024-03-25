# Pharmacogenomic-eQTL
Codes and sample data to perform PGx-eQTL analysis and visualization

This is the shared code to perform Pharmacogenomic (PGx)-eQTL analysis and visualization using ligand/drug-treated LCL panel. Technically, these codes can be used to perform any PGx-eQTL analysis. However, due to limited time, these codes were not extensively optimized for running options, computation speed, or memory efficiency. 

For purpose of manuscript review, here is a brief summary of resources in this project:

**Source Code in R/R-studio:**
PGx_eQTL_Main.Rmd: The major R studio command file. Each section is designated to perform a certain analysis or generate a plot used in the manuscript. The functionality of each titled section should be relatively intuitive.

Functions.R: All functions defined for and behind the running of PGx_eQTL_Main.Rmd. This file was "sourced" at the beginning of most of the sections in PGx_eQTL_Main.Rmd. 

**Data:**
Exp_data_ARhigh.tsv: Sample RNAseq data
Gene_loc.tsv: Sample gene annotation data
Exon_sub.tsv: Sample exon annotation data used in regional plot
SNP_data_ARhigh.tsv: Sample genotyping data
SNP_loc.tsv: Sample SNP annotation data with SNP position and allelle genotypes
snp_CHM_annotatation.tsv: Sample SNP annotation data with ChromHMM predicted chromotin status
ChIP_data.tsv: Sample ChIPseq data
GWAS_data.tsv: Sample GWAS data


**Concept of PGx-eQTL:**
![image](https://github.com/huanyaogao/Pharmacogenomic-eQTL/assets/88346194/2847d2bc-9d85-4c30-b294-0a5fb388b981)


If you used the code, please cite:
Thanh Thanh L Nguyen, Huanyao Gao, Duan Liu, Trudy Janice Philips, Zhenqing Ye, Jeong-Heon Lee, Geng-xian Shi, Kaleigh Copenhaver, Lingxin Zhang, Lixuan Wei, Jia Yu, Huan Zhang, Abhijeet Barath, Maggie Luong, Cheng Zhang, Alexandre Gaspar-Maia, Hu Li, Liewei Wang, Tamas Ordog, Richard M Weinshilboum, Glucocorticoids unmask silent non-coding genetic risk variants for common diseases, Nucleic Acids Research, Volume 50, Issue 20, 11 November 2022, Pages 11635–11653, https://doi.org/10.1093/nar/gkac1045
