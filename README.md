# PRRDetect R package

## Table of content

-   [Introduction](#intro)
-   [Install the package](#install)
-   [Package documentation](#PackageDoc)

<a name="intro"/>

## Introduction

PRRDetect is a R package that computes the probability of a sample being MMRd, Polymerses Dysfunctional (Poly-dys), MMRd+Poly-dys or Negative based on mutational signature analysis.

<a name="install"/>

## Installing the package

To install PRRDetect in first place the github repo has to be cloned and you have to enter in the downloaded folder.

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
  
  
## Using PRRDetect


  
### Loading PRRDetect
```{r setup}
library(PRRDetect)
```

### Launching PRRDetect

The following examples are done using experimental cell line data from a _Mixed_ cell line (Δ MLH1+POLD1) and a wild-type cell line, mimicking a _Negative_ sample.
In this case, for signature fitting I chose as organ "Colorectal" as seemed appropriate even though FitMS - the signature fit algorithm embedded in *signature.tools.lib* and used by PRRDetect - is designed for in-vivo cancer data.   

#### Single sample analysis from VCF files

```{r, eval = FALSE}
examples <- prepare_PRRDetect(Indel_VCF_path = "~/Cam_home/Giuseppe_ecDNA/PRRDetect/examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.pindel-filtered.vcf.gz", SNV_VCF_path = "~/Cam_home/Giuseppe_ecDNA/PRRDetect/examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.caveman-filtered.vcf.gz",nparallel = 4,organ = "Colorectal", sample_name = "Mixed_Samples")
PRRDetect_result <- PRRDetect(InDel_fits = examples$InDel, SNV_fits = examples$SNV, total_SNV = examples$total_SNV, total_InD = examples$total_InDel)

# > PRRDetect_result
#                   Negative        MMRd MMRd+Poly-dys     Poly-dys    Prediction
# Mixed_Samples 0.0002255015 5.08031e-06     0.9997611 8.313859e-06 MMRd+Poly-dys
```

#### Multiple samples analysis from VCF files

```{r, eval=FALSE}
examples <- prepare_PRRDetect(Indel_VCF_path = c("~/Cam_home/Giuseppe_ecDNA/PRRDetect/examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.pindel-filtered.vcf.gz", "~/Cam_home/Giuseppe_ecDNA/PRRDetect/examples/RPE1_WT_1_vs_RPE1_WT_0-re-pl2.03-pu1.pindel-filtered.vcf.gz"), SNV_VCF_path = c("~/Cam_home/Giuseppe_ecDNA/PRRDetect//examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.caveman-filtered.vcf.gz", "~/Cam_home/Giuseppe_ecDNA/PRRDetect/examples/RPE1_WT_1_vs_RPE1_WT_0-re-pl2.03-pu1.caveman-filtered.vcf.gz"),nparallel = 4,organ = "Colorectal", sample_name = c("Mixed_Samples","Negative_examples"))
PRRDetect_result <- PRRDetect(InDel_fits = examples$InDel, SNV_fits = examples$SNV, total_SNV = examples$total_SNV, total_InD = examples$total_InDel)

# > PRRDetect_result
#                       Negative         MMRd MMRd+Poly-dys     Poly-dys    Prediction
# Mixed_Samples     0.0002204726 5.011765e-06  0.9997662429 8.272818e-06 MMRd+Poly-dys
# Negative_examples 0.9984848933 8.905789e-04  0.0003948388 2.296891e-04      Negative
```


#### Multiple samples analysis from catalogues


If the user starts from catalogues then, by default, the SNV and the InDel burden are approximated as the number of mutations in the catalogue. Otherwise the user can specify these parameters from the _PRRDetect_ function.


```{r, eval = F}

SBS_VCF_PATHs <- c("examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.caveman-filtered.vcf.gz", "examples/RPE1_WT_1_vs_RPE1_WT_0-re-pl2.03-pu1.caveman-filtered.vcf.gz")

SBS_catalogues <- sapply(SBS_VCF_PATHs, simplify = F,  function(x) {signature.tools.lib::vcfToSNVcatalogue(x, genome.v = "hg38")$catalogue})
SBS_catalogues <- do.call("cbind", SBS_catalogues)
colnames(SBS_catalogues) <- c("Mixed_sample", "Negative_sample")

InDel_Catalogue <- sapply(c("examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.pindel-filtered.vcf.gz", "examples/RPE1_WT_1_vs_RPE1_WT_0-re-pl2.03-pu1.pindel-filtered.vcf.gz"), simplify = F, FUN = function(x){
  w <- read.table(x)[c(3,1,2,4,5)]
  w$V3 <- basename(x)
  colnames(w) <- c("Sample", "chr", "position", "REF", "ALT")
  muts <- indelsig.tools.lib::indel_classifier89(w, genome.v = "hg38")
  muts <- indelsig.tools.lib::indel_highspecific(muts)
  muts <- indelsig.tools.lib::gen_catalogue89(muts, sample_col = 1)
  muts <- muts[PRRDetect:::InDel_sigs_order,,drop = F]
  return(muts)
})

InDel_Catalogue <- do.call("cbind", InDel_Catalogue)
colnames(InDel_Catalogue) <- c("Mixed_sample", "Negative_sample")

examples <- prepare_PRRDetect(InDel_Catalogs = InDel_Catalogue, SNV_Catalogs = SBS_catalogues,nparallel = 4,organ = "Colorectal", sample_name = c("Mixed_Samples","Negative_examples"))

PRRDetect_result <- PRRDetect(InDel_fits = examples$InDel, SNV_fits = examples$SNV, total_SNV = examples$total_SNV, total_InD = examples$total_InDel)


# > PRRDetect_result
#                       Negative         MMRd MMRd+Poly-dys     Poly-dys    Prediction
# Mixed_Samples     0.0001924083 4.551485e-06  0.9997956859 7.354270e-06 MMRd+Poly-dys
# Negative_examples 0.9985891513 8.027154e-04  0.0003618234 2.463098e-04      Negative

```
