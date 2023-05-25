# PRRDetect R package

## Table of content

-   [Introduction](#intro)
-   [Install the package](#install)
-   [Package documentation](#PackageDoc)

<a name="intro"/>

## Introduction

PRRDetect is a R package that computes the probability for a sample to be Post-Replicative Repair (PRR) deficient.

<a name="install"/>

## Installing the package

To install PRRDetect in first place the github repo has to be cloned and you have to enter tin hte downloaded folder.

``` bash
git clone https://github.com/Nik-Zainal-Group/PRRDetect.git
cd PRRDetect
```

Successively open **R** and make sure you have the *devtools* package installed.

``` r
devtools::install()
```

<a name="PackageDoc"/>

## Package Documentation

The package consists mainly in two functions: **prepare_PRRDetect** and **PRRDetect**. As the name states, the first one prepares your input for PRRDetect. 

### **prepare_PRRDetect** has several parameters:

- **InDel_Catalogs**: Indel Catalogs dataframe.
- **SNV_Catalogs**: SNV catalogs dataframe.
- **Indel_VCF_path**: Vector of paths or list of data frames. If it is a vector of paths the file will be rearranged by the function, otherwise it has to be a list of dataframe having the following structure: "Sample", "chr", "position", "REF", "ALT".
- **total_InDel**: Total number of indels, it can be defined previously if wanted.
- **SNV_VCF_path**: Vector of paths or list of data frames.
- **genome.v**: Either "hg19" or "hg38", which is the default option.
- **nparallel**: Number of threads for signature fit.
- **is.filtered_SNV**: if yes, SNV VCF get filtered. Default equal to False. SNV information gets filtered only if provided as a vector of paths.
- **is.filtered_InDel**: if yes InDel VCF get filtered. Default equal to False. InDel information gets filtered only if provided as a vector of paths.
- **sample_name**: Vector of sample names, it has to be in the same order of Indel VCF vector or list.
- **organ**: Organ for SBS signature fit.
- **setseed**: Seed for reproducibility purposes.

It outputs a list with the following objects:
- **Indel fits**
- **SNV fits**
- **total number of SNV**
- **total number of InDels**

#### Important Notes:
  - **Make sure that the colnames of the catalogs or the vector of path and/or the sample names are ordered accordingly.**
  - **If no sample name is given and one of the inputs have sample names, especially is catalogs are provided, those names are going to be used**
  - **If a list of VCF inputs is provided, it is assumed that the input is already been filtered.**
  - **If a VCF file is given as input, then the total number of InDels is the total number of mutations, otherwise the sum of the mutations in the catalog will be given**
  

### **PRRDetect** has several parameters:

- **InDel_fits**: Indel Signatures Fits
- **SNV_fits**: SNV Signatures Fits
- **total_SNV**: Total number of SNVs of the sample
- **total_InD**: Total number of InDel mutation in the sample
- **prediction_type**: The type of model used. Default to "prob", which is the model published. Otherwise "abs" can be chosen.

It outputs a table defining the probabilities of each sample to belong to the four classes: *Negative*,*MMRd*,*MMRd+Poly-dys*,*Poly-dys*. 

#### Important Notes:
  * *The fits have to defined as*:
    - The InDel fits have as rows the samples and as columns the signatures as "RefSig.InD"
    - The SNV fits have as rows the samples and as columns the signatures as "SBS"
  
  

  
