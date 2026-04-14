# GWAS_Tutorial
## Overview
This repository provides an overview and tutorials for analyzing genetic data. These tutorials lead users through the processing and quality control (QC) of genotyping array data as it is done at the WUSTL Neurogenomics and Informatics Center. These tutorials include:
•	QC of array-based data
•	PseudoGWAS to remove batch effects
•	Imputation
•	Merging batches of data
•	Calculating Principle Components (PCs) for the data
•	Performing association studies
While example data is included for practice/testing purposes, these tutorials can also be used as a template for processing real data. Familiarity with a Unix/Linux environment is assumed and the tutorials rely on running in that environment.
## Setup
Setting up requires cloning this GitHub repository, which will download all scripts and data and keep them organized. This can be done by opening a shell on your computer or server, navigating to the location you want to keep the data, and entering “git clone https://github.com/NeuroGenomicsAndInformatics/GWAS_Tutorial.git” and hitting enter.
## Organization
Each step of the tutorial builds on the previous steps, so they must be completed in order. Each step has a Step#_name.txt and a Step#_name_automated.txt file to lead you through the tutorial. The automated files utilize production scripts to process data efficiently, while the non-automated files provide step by step walkthroughs. The code in the files should be copied and pasted into your Unix environment to be run.
