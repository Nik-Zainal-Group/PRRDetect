prepare_PRRDetect <- function(Indel_VCF_path, SNV_VCF_path, genome.v, nparallel, is.filtered_SNV=F, is.filtered_InDel=F, sample_name, organ, setseed = 42){
  if(is.filtered_SNV == F){
      SNV_Filter(SNV_VCF_path, vcfout = paste0(SNV_VCF_path, ".filtered"), genomev = genome.v,sample_name = sample_name )
      SNV_VCF_path <- paste0(SNV_VCF_path, ".filtered.bgz")
   }
  if(is.filtered_InDel == F){
    InDelFilter(vcffile = Indel_VCF_path, vcfout = paste0(Indel_VCF_path, ".filtered"), genomev = genome.v, sample_name = sample_name )
    Indel_VCF_path <- paste0(Indel_VCF_path, ".filtered.bgz")
   }

  catalogs  <- generate_catalogs_from_vcf(Indel_VCF_path, SNV_VCF_path, genome.v, sample_name)
  signature_fits <- Signature_fit(SNV_Catalog = catalogs$SNV_Catalog, InDel_Catalog = catalogs$Indel_Catalog, organ, nparallel, setseed)
  SNV_fit <- signature_fits$SNV$exposures
  InDel_fit <- signature_fits$InDel$exposures
  colnames(SNV_fit) <- unlist(lapply(strsplit(colnames(SNV_fit), split = "_", fixed = T), function(x){tail(x,1)}))
  return(list("SNV"=SNV_fit, "InDel"=InDel_fit, "total_SNV"=catalogs$total_SNV, "total_InDel"=catalogs$total_InD))
}



# TODO : check names of the fits

# TODO : test it


#' Compute PRRDetect probability and label to file
#'
#' @param ind_fits : "Indel Signatures Fits"
#' @param snv_fits : "SNV Signatures Fits"
#' @return PRRDetect table
#' @export
PRRDetect <- function(InDel_fits, SNV_fits, total_SNV, total_InD){

  pure_MMR_sigs <- paste0("SBS", c(6,15,26,44,97))
  pure_POL_sigs <- paste0("SBS", c("10a", "10d"))
  mixed_sigs <- paste0("SBS", c(14,20))

  pure_MMR_indS <- paste0("RefSig.InD", c(7,19))
  pure_POL_indS <- paste0("RefSig.InD", c(14,15))
  mixed_indS <- paste0("RefSig.InD", c("16a", "16b", 20, 21))

  ## if the name of the signature is not present, which returns "integer(0)", then the sum is equal to 0, otherwise is equal to the rowsum

  ## SNV
  ifelse(identical(which(colnames(SNV_fits) %in% pure_MMR_sigs), integer(0)), MMR_SBS <- 0 , MMR_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% pure_MMR_sigs), drop=F])/ total_SNV)

  ifelse(identical(which(colnames(SNV_fits) %in% pure_POL_sigs), integer(0)), POL_SBS <- 0 , POL_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% pure_POL_sigs), drop=F])/ total_SNV)

  ifelse(identical(which(colnames(SNV_fits) %in% mixed_sigs), integer(0)), MIX_SBS <- 0 , MIX_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% mixed_sigs), drop=F])/ total_SNV)


  ## INDEL
  ifelse(identical(which(colnames(InDel_fits) %in% pure_MMR_indS), integer(0)), MMR_IND <- 0 , MMR_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% pure_MMR_indS), drop=F])/total_InD)

  ifelse(identical(which(colnames(InDel_fits) %in% pure_POL_indS), integer(0)), POL_IND <- 0 , POL_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% pure_POL_indS), drop=F])/total_InD)

  ifelse(identical(which(colnames(InDel_fits) %in% mixed_indS), integer(0)), MIX_IND <- 0 , MIX_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% mixed_indS), drop=F])/total_InD)


  Ratio <-  c(total_InD + 1)/(total_SNV + 1)

  final_df <- as.data.frame(cbind(MMR_SBS, MIX_SBS, POL_SBS, MMR_IND, MIX_IND, POL_IND, Ratio))

  final_df[,1:6] <- log2(final_df[,1:6]+1)
  final_df[,7] <- log2(final_df[,7])

  final_df <- prediction_function(final_df)
  colnames(final_df)[1:4] <- c("MMRd", "MMRd+Poly-dys", "Negative", "Poly-dys")
  final_df$Prediction = prediction_label(final_df)
  final_df <- final_df[,c(3,1,2,4,5)]
  return(final_df)

}


