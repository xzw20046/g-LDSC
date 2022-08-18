# g-LDSC ```V1.0.0```
g-LDSC is a tool for estimating heritability and functional enrichment from GWAS summary statistics. g-LDSC is written under R-4.1.0. In this tutorial, we would give a tutorial
# Tutorial
g-LDSC is written under R-4.1.0. In this tutorial, we would give a demonstration of how g-LDSC is run under a Linux-based system.
## Installation
We need to install all the packages we need in R before we run the tool.
```
Rscript Install.gLDSC.R
```
## Input file
3 files are required to run g-LDSC:
- GWAS summary statistics
- Pre-calculated (partitioned) LD score matrix
- g-LDSC function file
### GWAS summary statistics
For GWAS summmary statistics, the required formate is shown as follows:
```
SNP A1 A2 N Z
rs1000000 G A 361194 1.5397
rs10000010 T C 361194 -0.850433
rs1000002 C T 361194 0.672368
```
To convert your GWAS result in such format, you could use ```munge_sumstats.py``` from [ldsc](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability).
### Pre-calculated LD score matrix
This file contains the information of LD matrix and annotation information, to get this file, see tutorial of **Calculate LD score matrix** in the futher section.
### g-LDSC function file
The R file ```functions.R``` that contain all g-LDSC functions.
## Usage
