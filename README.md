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
git clone https://github.com/GRinaldi97/synthetic.catalogs.git
cd synthetic.catalogs
```

Successively open **R** and make sure you have the *devtools* package installed.

``` r
devtools::install()
```

<a name="PackageDoc"/>

The package consists mainly in two functions: **prepare_PRRDetect** and **PRRDetect**. As the name states, the first one prepares your input for PRRDetect. The function has several parameter:
- *InDel_Catalogs*: Indel Catalogs dataframe
- *SNV_Catalogs*: SNV catalogs daframe
- *Indel_VCF_path*: Vector of paths or list of data frames. If it is a vector of paths the file will be rearranged by the function, otherwise it has to be a list of dataframe having the following structure: "Sample", "chr", "position", "REF", "ALT"
- *total_InDel*: Total number of indels, it can be defined previously
- *SNV_VCF_path*: Vector of paths or list of data frames.
- *genome.v*: Either "hg19" or "hg38", which is the default option.
- *nparallel*: Number of threads for signature fit
- *is.filtered_SNV*: if yes, SNV VCF get filtered. Default equal to False. SNV information gets filtered only if provided as a vector of paths.
- *is.filtered_InDel*: if yes InDel VCF get filtered. Default equal to False. InDel information gets filtered only if provided as a vector of paths.
- *sample_name*: Vector of sample names, it has to be in the same order of Indel VCF vector or list.
- *organ*: Organ for SBS signature fit
- *setseed*: Seed for reproducibility purposes
