library(PRRDetect)


### Single

examples <- prepare_PRRDetect(Indel_VCF_path = "examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.pindel-filtered.vcf.gz", SNV_VCF_path = "examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.caveman-filtered.vcf.gz",nparallel = 4,organ = "Colorectal", sample_name = "Mixed_Samples")
PRRDetect_result <- PRRDetect(InDel_fits = examples$InDel, SNV_fits = examples$SNV, total_SNV = examples$total_SNV, total_InD = examples$total_InDel)

### Multiple samples

examples <- prepare_PRRDetect(Indel_VCF_path = c("examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.pindel-filtered.vcf.gz", "examples/RPE1_WT_1_vs_RPE1_WT_0-re-pl2.03-pu1.pindel-filtered.vcf.gz"), SNV_VCF_path = c("examples/DheML6d1_vs_DheML6-re-pl2.03-pu1.caveman-filtered.vcf.gz", "examples/RPE1_WT_1_vs_RPE1_WT_0-re-pl2.03-pu1.caveman-filtered.vcf.gz"),nparallel = 4,organ = "Colorectal", sample_name = c("Mixed_Samples","Negative_examples"))
PRRDetect_result <- PRRDetect(InDel_fits = examples$InDel, SNV_fits = examples$SNV, total_SNV = examples$total_SNV, total_InD = examples$total_InDel)


#############################
#############################
## Starting from Catalogues

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


