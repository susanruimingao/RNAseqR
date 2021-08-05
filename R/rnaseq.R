#setwd("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice")

# in R studio, install the required packages and load necessary packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("csaw")
# BiocManager::install("mixOmics")

library(tximportData)
library(tximport)
library(limma)
library(edgeR)
library(RColorBrewer)
library(mixOmics)
library(magrittr)
library(data.table)
library(devtools)

trimReads = function(inputDir,
                     outputDir = "trimmedReads",
                     suffixNameR1 = "_R1.fastq.gz", #must be these forms; _R1.fastq.gz, _1.fastq.gz, _R1.fq.gz, _1.fq.gz;,
                     suffixNameR2 = "_R2.fastq.gz",
                     nThreads = 16, ...){
  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse);
  library(parallel); library(crayon); library(stringr);
  
  #check number of threads;
  nCores = detectCores();
  if(nCores < nThreads){
    nThreads = nCores;
    cat(format(Sys.time(), usetz = TRUE), yellow(" using ", nThreads), green("threads"), "\n");
  }
  
  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }
  
  #get file infor
  fileNames = list.files(inputDir, pattern = ".fastq.gz$");
  fileNamesLeft = fileNames[str_detect(fileNames, suffixNameR1)] %>% str_remove(suffixNameR1);
  fileNamesRight = fileNames[str_detect(fileNames, suffixNameR2)] %>% str_remove(suffixNameR2);
  nSamples = 0;
  
  if(all(fileNamesLeft == fileNamesRight)){
    nSamples = length(fileNamesLeft);
    cat(format(Sys.time(), usetz = TRUE), yellow(" all read files (", nSamples, ") are paired\n"));
    fileLeftInput = file.path(inputDir, paste0(fileNamesLeft, suffixNameR1));
    fileRightInput = file.path(inputDir, paste0(fileNamesLeft, suffixNameR2));
    fileLeftOutput = file.path(outputDir, paste0(fileNamesLeft, "_clean_R1.fastq.gz"));
    fileRightOutput = file.path(outputDir, paste0(fileNamesRight, "_clean_R2.fastq.gz"));
    fileHtml = file.path(outputDir, paste0(fileNamesLeft, ".html"));
  } else {
    cat(format(Sys.time(), usetz = TRUE));
    stop(red(" pairend reads are not paired, please check and ensure pairedend reads"), "\n");
  }
  
  #using fastp to do QC and trim;
  for(i in 1:nSamples){
    cmd = paste("fastp -V",
                "-i", fileLeftInput[i],
                "-I", fileRightInput[i],
                "--correction",
                "-o", fileLeftOutput[i],
                "-O", fileRightOutput[i],
                "-w", nThreads,
                "--html", fileHtml[i],
                sep = " ");
    cat(format(Sys.time(), usetz = TRUE), 
        yellow(" begin to process quality check using fastp for sample ", 
               fileNamesLeft[i]),
        " ", 
        green(i, " out of ", nSamples, "\n"));
    system(cmd);
    cat(format(Sys.time(), usetz = TRUE), 
        yellow(" finish to process quality check using fastp for sample ", fileNamesLeft[i]), 
        " ", 
        green(i, " out of ", nSamples, "\n"));
  }
}

#inputDir = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice";
#trimReads("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice" )
#trimReads(inputDir = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice" )

prepareIndex <- function(transcriptome,
                         output = "rsem_index",
                         ...){
  cmd = paste("rsem-prepare-reference",
              "--bowtie2",
              transcriptome,
              output,
              sep = " ");
  print(cmd);
  #cat(format(Sys.time(), usetz = TRUE), yellow(" begin to prepare reference \n"));
  system(cmd);
}

#prepareIndex("Bm_atcc23344cDNA.fa", "Bm_atcc23344cDNA")


expressionValue = function(trimmedReadsDir = "trimmedReads",
                           outputDir = "rsemValueOut",
                           nMemory = 30, #memory in GB;
                           reference = "rsem_index",
                           nThreads = 12, ...){
  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(parallel); library(crayon); library(stringr); library(seqinr);
  
  #check number of threads;
  nCores = detectCores();
  if(nCores < nThreads){
    nThreads = nCores;
    cat(format(Sys.time(), usetz = TRUE), yellow(" using ", nThreads), green("threads"), "\n");
  }
  
  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }
  
  #get file infor
  fileNames = list.files(trimmedReadsDir, pattern = "_clean_R1.fastq.gz") %>%
    str_remove("_clean_R1.fastq.gz");
  
  nSamples = length(fileNames);
  
  if(nSamples < 1){
    cat(format(Sys.time(), usetz = TRUE));
    stop(red(" detect ", nSamples, " , please check your ", trimmedReadsDir), "\n");
  } else {
    cat(format(Sys.time(), usetz = TRUE), yellow(" detect", nSamples, ") for rsem expression values\n"));
    #fileMergedInput = file.path(trimmedReadsDir, paste0(fileNames, "_clean_merged.fastq.gz"));
    fileLeftUnMerged = file.path(trimmedReadsDir, paste0(fileNames, "_clean_R1.fastq.gz"));
    fileRightUnMerged = file.path(trimmedReadsDir, paste0(fileNames, "_clean_R2.fastq.gz"));
  }
  
  for(i in 1 : nSamples){
    cmd = paste("rsem-calculate-expression",
                "-p", nThreads,
                "--bowtie2",
                "-paired-end",
                fileLeftUnMerged[i],
                fileRightUnMerged[i],
                reference,
                file.path(outputDir, fileNames[i]));
    print(cmd);
    system(cmd);
    cat(format(Sys.time(), usetz = TRUE), yellow(" begin to caculate gene expression value with rsem ", fileNames[i]), " ", green(i, " out of ", nSamples, "\n"));
  }
}

#expressionValue(trimmedReadsDir = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/trimmedReads", reference = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/Bm_atcc23344cDNA")


statistics_rsem = function(inputDir = "rsemValueOut",
                           suffix = ".genes.results",
                           outputDir = "rsemStatistics",
                           ...){
  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(crayon);
  
  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }
  
  
  #get file infor
  fileNames = list.files(inputDir, pattern = ".genes.results") %>%
    str_remove(".genes.results");
  
  # if(nSamples < 1){
  # cat(format(Sys.time(), usetz = TRUE));
  # stop(red(" detect ", nSamples, " , please check your ", inputDir), "\n");
  # } else {
  # cat(format(Sys.time(), usetz = TRUE), yellow(" detect", nSamples, ") for rsem statistics\n"));
  #fileMergedInput = file.path(trimmedReadsDir, paste0(fileNames, "_clean_merged.fastq.gz"));
  #fileLeftUnMerged = file.path(trimmedReadsDir, paste0(fileNames, "_clean_R1.fastq.gz"));
  #fileRightUnMerged = file.path(trimmedReadsDir, paste0(fileNames, "_clean_R2.fastq.gz"));
  #}
  
  for(i in 1: length(fileNames)){
    fileName = fileNames[i];
    #geneFileName = paste0(fileName, ".genes.results");
    inputFile = file.path(inputDir, fileName);
    print(inputFile);
    outputFile = file.path(outputDir, paste0(fileName, ".pdf"));
    print(outputFile);
    cmd = paste("rsem-plot-model",
                inputFile,
                outputFile);
    #file.path(inputDir, fileNames[i]));
    #fileRightUnMerged[i]
    
    print(cmd);
    cat(format(Sys.time(), usetz = TRUE), yellow(" begin rsem statistics analysis ", fileName), " ", green(i, " out of ", length(fileNames), "\n"));
    system(cmd);
  }
}


#statistics_rsem("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/rsemValueOut")

countsTableDEGs = function(inputDir = "rsemValueOut",
                           suffix = ".genes.results",
                           #outputDir = "CtsDEGOut",
                           metaData,
                           ...){
  # #create output file;
  # if(!dir.exists(outputDir)){
  #   cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
  #   dir.create(outputDir);
  # } else {
  #   cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
  #   unlink(outputDir, recursive = TRUE);
  #   dir.create(outputDir);
  # }
  
  samples <- fread(metaData, header = TRUE)
  files <- file.path(inputDir, 
                     paste0(samples$isolate, ".genes.results"))
  
  # oldNames = list.files("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/rsemValueOut",
  #                       pattern = ".genes.results");
  # 
  # for(i in oldNames){
  #   fp = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/rsemValueOut";
  #   oldFile = file.path(fp, i);
  #   newName = i %>%  
  #     str_remove("HI.*_i5_") %>% 
  #     str_remove(".{2}\\.");
  #   
  #   newFile = file.path(fp, newName);
  #   file.rename(oldFile, newFile)
  # }
  
  names(files) <- samples$isolate
  txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
  #head(txi.rsem$counts)
  cts <- txi.rsem$counts
  #head(cts)
  write.table(cts, file = "counts.txt", sep = "\t")
  
  counts <- txi.rsem$counts
  
  d0 <- DGEList(counts)
  d0 <- calcNormFactors(d0)
  
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  #dim(d) # number of genes left
  
  group <- as.factor(samples$description)
  
  pdf(file="mdsplot.pdf", width = 10, height = 6);
  plotMDS(d, col = as.numeric(group))
  dev.off()
  
  mm <- model.matrix(~0 + group)
  pdf(file="voom.pdf", width = 10, height = 6);
  y <- voom(d, mm, plot = T)
  dev.off()
  
  pdf(file="voom_d0.pdf", width = 10, height = 6);
  tmp <- voom(d0, mm, plot = T)
  dev.off()
  
  fit <- lmFit(y, mm)
  #head(coef(fit))
  
  contr <- makeContrasts(groupDD3008 - groupWTstd_23344, levels = colnames(coef(fit)))
  
  tmp <- contrasts.fit(fit, contr)
  
  tmp <- eBayes(tmp)
  
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  #head(top.table, 20)
  
  length(which(top.table$adj.P.Val < 0.05))
  filDEG <- length(which(top.table$adj.P.Val < 0.05))
  #view(filDEG)
  
  top.table %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    write.table(file = "filDEG_DD3008_WT23344.txt", row.names = F, sep = "\t", quote = F)
}

#countsTableDEGs(inputDir = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/rsemValueOut",
#metaData = "/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/RNA-Seq_practice/Bmallei_meta_practice.txt")


trimExpressionStatisticsCtsTableDEGs <- function(inputDir,
                                                 suffixNameR1 = "_R1.fastq.gz", #must be these forms; _R1.fastq.gz, _1.fastq.gz, _R1.fq.gz, _1.fq.gz;,
                                                 suffixNameR2 = "_R2.fastq.gz",
                                                 nThreads = 16,
                                                 nMemory = 30,
                                                 transcriptome,
                                                 metaData,
                                                 output = "rsem_index",
                                                 reference = "rsem_index",
                                                 ...){
  RNASEQ::trimReads(inputDir = inputDir,
            outputDir = "trimmedReads",
            suffixNameR1 = suffixNameR1, #must be these forms; _R1.fastq.gz, _1.fastq.gz, _R1.fq.gz, _1.fq.gz;,
            suffixNameR2 = suffixNameR2,
            nThreads = nThreads);
  
  RNASEQ::prepareIndex(transcriptome = transcriptome,
               output = output)  
  
  RNASEQ::expressionValue(trimmedReadsDir = "trimmedReads",
                  outputDir = "rsemValueOut",                    
                  nMemory = nMemory, #memory in GB;
                  reference = output,
                  nThreads = nThreads);
  
  RNASEQ::statistics_rsem(inputDir = "rsemValueOut",
                  suffix = ".genes.results");
  
  RNASEQ::countsTableDEGs(inputDir = "rsemValueOut",
                  suffix = ".genes.results",
                  outputDir = "CtsDEGOut",
                  metaData = metaData);
}
