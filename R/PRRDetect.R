prepare_PRRDetect <- function(Indel_VCF_path, SNV_VCF_path, genome.v, nparallel, is.filtered_SNV=F, is.filtered_InDel=F, sample_name){
  if(is.filtered_SNV == F){
      SNV_Filter(SNV_VCF_path, vcfout = paste0(SNV_VCF_path, ".filtered"), genomev = genome.v,sample_name = sample_name )
      SNV_VCF_path <- paste0(SNV_VCF_path, ".filtered")
   }
  if(is.filtered_InDel == F){
    InDelFilter(vcffile = Indel_VCF_path, vcfout = paste0(Indel_VCF_path, ".filtered"), genomev = genome.v, sample_name = sample_name )
    Indel_VCF_path <- paste0(Indel_VCF_path, ".filtered")
   }

  catalogs  <- generate_catalogs_from_vcf(Indel_VCF_path, SNV_VCF_path, genome.v, sample_name)
  signature_fits <- Signature_fit(catalogs$SNV_catalog, catalogs$Indel_Catalog, organ, nparallel)
  return(signature_fits)
}



# TODO : check names of the fits

# TODO : test it


#' Compute PRRDetect probability and label to file
#'
#' @param ind_fits : "Indel Signatures Fits"
#' @param snv_fits : "SNV Signatures Fits"
#' @return PRRDetect table
#' @export
PRRDetect <- function(InDel_fits, SNV_fits){

  pure_MMR_sigs <- paste0("SBS", c(6,15,26,44,97))
  pure_POL_sigs <- paste0("SBS", c("10a", "10d"))
  mixed_sigs <- paste0("SBS", c(14,20))

  pure_MMR_indS <- paste0("RefSig.InD", c(7,19))
  pure_POL_indS <- paste0("RefSig.InD", c(14,15))
  mixed_indS <- paste0("RefSig.InD", c("16a", "16b", 20, 21))

  #############################################

  ## add prediction and and label

  #############################################
  return(sample_out)
}
