# RNA-seq Analysis: Effect of Sodium Butyrate in Diabetic Kidney Disease Model

## Project Overview
This project implements an RNA-seq analysis workflow using publicly available data, expecifically a dataset to study the effect of sodium butyrate (NaB) on gene expression in a cellular model of diabetic kidney disease (DKD). The aim is to identify transcriptional changes associated with NaB treatment under glucolipotoxic (GLT) conditions using a reproducible RNA-seq workflow. 

## Dataset
**Source:** NCBI Gene Expression Omnibus (GEO)\
**Accession:** GSE266108\
**Organism:** Homo sapiens\
**Cell line:** HK-2 (human renal proximal tubular epithelial cells)\
**Sequencing type:** RNA-seq\
**Raw data:** FASTQ files available via NCBI SRA

**Summary:** This study investigates the role of sodium butyrate (NaB) in diabetic kidney disease (DKD), a condition characterized by impaired lipid metabolism and defective lipophagy. Using glucolipotoxicity (GLT)-stimulated human renal tubular epithelial cells (HK-2), RNA-seq was performed to compare gene expression profiles with and without NaB treatment. The results show that NaB improves lipid accumulation and renal dysfunction by restoring lipophagy, likely through activation of the PP2A–TFEB signaling pathway. This dataset enables the analysis of differential gene expression associated with NaB treatment in a controlled cellular model of DKD.

**Experimental design:** RNA-seq was performed on GLT-stimulated HK-2 cells, comparing three biological replicates treated with sodium nutyrate (NaB_GLT) against three untreated controls (BSA_GLT).

**Sample information**: A clean sample annotation table ['samples.csv'](samples.csv) was generated form the original SRA metadata and contains the run accesions an experimental conditions used in this analysis.

## Analysis plan
1. Quality control of raw reads.
2. Read alignment / quantification.
3. Gene-level expression analysis.
4. Differential expression analysis.
5. Visualization of results and interpretation.

## Goal
This project is designed as a learning and portfolio project to demonstrate practicall skills in RNA-seq data handling, experimental design interpretation, and differential expression analysis.


## Status
Analysing the dataset.

## Research Question
Which genes are differentially expressed between two experimental conditions in an RNA-seq experiemnt using public data?

This project focuses on applying a standard RNA-seq workflow to answer a biologically relevant question in a reproducible way.

## VS Studio
This last line is edited and updated exclusively from VS Studio. It´s a try out.

