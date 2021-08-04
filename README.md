# RNAseqR
# RNA-Seq based transcriptomic analysis for differential expressed genes (DEGs)
RNA-Seq analysis pipeline is developed in R language and is a R wrapper for multiple functions; 

In a terminal, conda create a "DEGs" environment and install three software fastp(v0.20.1) and rsem (v1.3.3). 
```
codna create –n DEGs –c bioconda fastp rsem

```
Go to work directory containing raw reads, and the reads name is preferred in "_R1.fastq.gz"/"_R2.fastq.gz"; otherwise specify in the command line.
Go to R environment and install three required packages, which will be needed to be installed once.  

Download RNA-SeqR codes from https://github.com/susanruimingao/RNAseqR.git

in R, install the required packages and load necessary packages

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
```
```
library(tximportData)
library(tximport)
library(limma)
library(edgeR)
library(RColorBrewer)
library(mixOmics)
library(magrittr)
library("data.table")
```

Every time for enterring the R envrionment, the above downloaded RNA-SeqR_0.1.0.tar.gz package need to be installed
``` 
install.packages("/home/CFIA-ACIA/gaoru/R_package/RNA-SeqR_0.1.0.tar.gz", repos = NULL, type="source")
```

To run the RNA-Seq pipeline, there are five main steps are involved, using general function: including trim reads, gene reads mapping and expression, statistics (mdsplot), counts table and differential expressed genes (DEGs)

The default RawReads format suffixNameR1 = "_R1.fastq.gz",suffixNameR2 = "_R2.fastq.gz". Otherwise, specifiy as following:

The example format of the prepared meta file
isolate	repeat	description
23344-1	23344-Bio1	WTstd_23344
23344-2	23344-Bio2	WTstd_23344
23344-3	23344-Bio3	WTstd_23344
DD3008-1	DD3008-Bio1	DD3008
DD3008-2	DD3008-Bio2	DD3008
DD3008-3	DD3008-Bio3	DD3008

```
RNA-SeqR::trimExpressionStatisticsCtsTableDEGs(inputDir = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice", 
                         suffixNameR1 = "_R1.fastq.gz", 
                         suffixNameR2 = "_R2.fastq.gz",
                         transcriptome = "Bm_atcc23344cDNA.fa", 
                         output = "Bm_atcc23344cDNA",
                         reference = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/Bm_atcc23344cDNA",
                         metaData = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/Bmallei_meta_practice.txt")
```
