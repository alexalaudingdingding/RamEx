# RamEx

![Version](https://img.shields.io/badge/Version-1.0-brightgreen)
![Release date](https://img.shields.io/badge/Release%20date-Jan.%2010%2C%202022-brightgreen)



## Contents

- [Introduction](#introduction)
- [System Requirement and dependency](#system-requirement-and-dependency)
- [Installation guide](#installation-guide)
- [Usage](#usage)
- [Example dataset](#example-dataset)
- [Supplementary](#supplementary)

# Introduction

A ramanome represents a single-cell-resolution metabolic phenome that is information-rich, revealing functional heterogeneity among cells, and universally applicable to all cell types.Ramanome Explorer (RamEx, Fig.1) is a toolkit for comprehensive and efficient analysis and comparison of ramanomes. Results from the multidimensional analysis are visualized via intuitive graphics. Implemented via R and C++, RamEx is fully extendable and supports cross-platform use.By providing simple-to-use modules for computational tasks that otherwise would require substantial programming experience and algorithmic skills, RamEx should greatly facilitate the computational mining of ramanomes.
![Fig. 1. The software architecture of RamEx.](http://bioinfo.single-cell.cn/images/RamEx.png)  
<div align=center>Fig.1. The software architecture of RamEx. RamEx integrates five core modules (data governance and preprocessing, real-time ramanome depth analysis, microbiome heterogeneity index analysis, IRCA, and RBCS) to support a complete ramanome-analysis pipeline.</div>



# System Requirement and dependency

## Hardware Requirements

RamEx only requires a standard computer with >1GB RAM to support the operations defined by a user.

## Software Requirements

RamEx only requires a C++ compiler (e.g. g++) to build the source code. Most Linux releases have g++ already been installed in the system. In Mac OS, to install the compiler, we recommend to install the Xcode application from the App store.

# Installation guide

## Automatic Installation (recommended)

At present, RamEx provides a fully automatic installer for easy installation.

**a. Download the package**
```
git clone https://github.com/qibebt-bioinfo/RamEx.git	
```

**b. Install by installer**
```
cd RamEx
source install.sh
```

The package should take less than 1 minute to install on a computer with the specifications recommended above.

The example dataset could be found at “examples” folder. Check the “examples/Readme” for details about the demo run.

## Manual Installation

If the automatic installer fails, RamEx can still be installed manually.

**a. Download the package**
```
git clone https://github.com/qibebt-bioinfo/RamEx.git	
```

**b. Configure the environment variables (the default environment variable configuration file is “~/.bashrc”)**
```
export RamEX=Path to RamEx
export PATH="$PATH:$RamEX/bin/"
source ~/.bashrc
```
**c. Compile the source code**
```
cd RamEx
make
```
# Usage
The RamEx consists of three modules: QC（pretreatment）IRCA and RBCS. Currently the RamEx requires all CAST-R, Horiba and Renishaw files as txt object.

**I. QC（pretreatment）** 

Ramanome analysis starts with data management and preprocessing of SCRS
In the QC step, RamEx provides several functions: initialization, read data,  smooth, quanlity control, baseline, normalization and visual variance analysis, including PCA,t-SNE and PLS-DA. 

Currently the RamEx requires files as txt object. All spectrum data stored in few subfolders of the working path by group. 

**Initialization**

The initialization function makesure the working directory and creat the ouput directry.

**Read data**

For CAST-R and Horiba spectrum data, RamEx reads data from the subfolders of the working path by group. And for Renishaw mapping data, RamEx can read mapping data from the working directory and creat some subfolders by group. Then, RamEx can read them like CAST-R and Horiba spectrum data.

**Preprocessing(smooth; quanlity control; baseline; normalization)**

All data can be preprocessed by combining the provided functions as desired.

**Visual variance analysis(PCA,t-SNE and PLS-DA)**

Low-dimensional spatial visualization provides us with an intuitive understanding of ramanome. 


**II. Calibration for KO abundance**

In the calibration step, RamEx calibrates the predicted functional profiles of more amplicon samples using the model built in the training step. The predicted gene profiles of amplicons for calibration should be processed in the same way as the training amplicon samples (e.g. by PICRUSt 2).  

For calibration, the RamEx also accepts gene profiles of amplicon samples in two formats as well as the training step.  

**a. Abundance table**  

```
RamEx-calibrate -t 16s.ko.abd -m RamEx.model -o 16s.ko.calibrated.abd
```
The output file “16s.ko.calibrated.abd” is the calibrated gene abundance table. 

**b. Sample list**  
```
RamEx-calibrate -l 16s.ko.list -m RamEx.model -o 16s.ko.calibrated.out
```
The output folder “16s.ko.calibrated.out” contains the calibrated gene profiles of each single input sample, and the file list of calibrated samples is also output to “16s.ko.calibrated.out.list”.  

# Example dataset

Here we provide a demo dataset in "examples" folder. This dataset contains 15 WGS-amplicon training pairs, and 622 amplicon samples for calibration. All samples were produced by Human Microbiome Project stage I. Gene profiles of WGS samples were directly analyzed by HuMANn 2. Amplicon samples were sequenced by V3-V5 region 16S rRNA sequencing, and gene profiles were inferred by PICRUSt 2.   

In this dataset, " training.wgs.ko.abd" is the gene relative abundance table of 15 training WGS samples, and "training.16s.ko.abd" is that of 15 matched training amplicon samples. The gene profiles of samples for calibration is in “16s.ko.abd”, and their meta-data is in “meta.txt”.  


To run the demo, you can either run the script “Readme”:  

```
cd example
sh Readme
```
or type the following command:
```
RamEx-train -T training.wgs.ko.abd -t training.16s.ko.abd -o RamEx.model
RamEx-calibrate -t 16s.ko.abd -m RamEx.model -o 16s.ko.calibrated.abd
```
Then the output file "16s.ko.calibrated.abd " is the calibrated relative abundance of the amplicon microbiomes.  

This demo run should take less than 2 minutes on a recommended computer.  



# Supplementary
We also provide more datasets for test run, which is used in the manuscript “RamEx improves accuracy of 16S-amplicon-based prediction of microbiome function”. Gene profiles of WGS samples were directly analyzed by HuMANn 2. and gene profiles of amplicons were inferred by PICRUSt 2. 

[Dataset 1](http://bioinfo.single-cell.cn/Released_Software/RamEx/dataset1.tar.gz) (~109 MB) contains 622 paired samples of WGS and V3-V5 region 16S rRNA amplicon samples (also integrated in the source code package as example dataset for demo run).  
[Dataset 2](http://bioinfo.single-cell.cn/Released_Software/RamEx/dataset2.tar.gz) (~47 MB) contains 295 paired samples of WGS and V1-V3 region 16S rRNA amplicon samples.   
[Dataset 3](http://bioinfo.single-cell.cn/Released_Software/RamEx/dataset3.tar.gz) (~827 MB) contains unpaired 2,354 WGS samples and 5,350 V3-V5 16S rRNA amplicon samples.  
[Dataset 4](http://bioinfo.single-cell.cn/Released_Software/RamEx/dataset4.tar.gz) (~325 MB) contains unpaired 2,045 WGS samples and 2,186 V1-V3 16S rRNA amplicon samples.  
[Dataset 5](http://bioinfo.single-cell.cn/Released_Software/RamEx/dataset5.tar.gz) (~25 MB) contains partially paired 18 WGS samples and 150 V1-V3 16S rRNA amplicon samples.

# Reference
1. Abubucker, S., et al., Metabolic Reconstruction for Metagenomic Data and Its Application to the Human Microbiome. _PLoS Computational Biology_, 2012, 8(6): e1002358.
2. Douglas, G.M., et al., PICRUSt2 for prediction of metagenome functions. _Nat Biotechnol_, 2020, 38(6): 685-688.

