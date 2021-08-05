# RNAseqR
# RNA-Seq based transcriptomic analysis for differential expressed genes (DEGs)
RNA-Seq analysis pipeline is developed in R language and is a R wrapper for multiple functions; 

In a terminal, conda create a "DEGs" environment and install three software fastp(v0.20.1) and rsem (v1.3.3). 
```
codna create –n DEGs –c bioconda fastp rsem

```
Go to work directory containing raw reads;
Go to R environment and install the following required packages, which will be needed to be installed once.  


in R, install the required packages and it only needs to install once

```
R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("tximportData")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("csaw")
BiocManager::install("mixOmics")
install.packages("data.table")
install.packages("devtools")
install.packages("readr")
```
```
R
```

After installing the required software and R packages, install the final RNASEQ package
``` 
library(devtools)
devtools::install_github("susanruimingao/RNAseqR")
```
Load the required package each time before running your job
```
library(tximportData)
library(tximport)
library(limma)
library(edgeR)
library(RColorBrewer)
library(mixOmics)
library(magrittr)
library("data.table")
library("devtools")
library("readr")
```

To run the RNA-Seq pipeline, there are five main steps are involved, using general function: including trim reads, gene reads mapping and expression, statistics (mdsplot), counts table and differential expressed genes (DEGs)

The default RawReads format suffixNameR1 = "_R1.fastq.gz",suffixNameR2 = "_R2.fastq.gz". Otherwise, specifiy as following:

The prepared metadata file ID names need to match the ID of raw reads.


```
RNASEQ::trimExpressionStatisticsCtsTableDEGs(inputDir = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice", 
                         suffixNameR1 = "_R1.fastq.gz", 
                         suffixNameR2 = "_R2.fastq.gz",
                         transcriptome = "Bm_atcc23344cDNA.fa", 
                         output = "Bm_atcc23344cDNA",
                         reference = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/Bm_atcc23344cDNA",
                         metaData = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/Bmallei_meta_practice.txt")
```
