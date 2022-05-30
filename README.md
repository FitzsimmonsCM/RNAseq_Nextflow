# RNAseq_Nextflow
RNAseq Nextflow pilot for Batista Lab

The goal of this project is to develop a working pipeline to process RNA-seq data from the lab. 
This pipeline will take input data from the commonly used Illumina TruSeq library kit and perform adapter trimming, genome mapping, deduplication, and counting. 
This pipeline is written using the nextflow pipeline framework, and designed to be a cluster-agnostic pipeline by making use of docker and singularity containers. 
