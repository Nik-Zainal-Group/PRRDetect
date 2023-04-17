prepare_PRRDetect <- function(Indel_VCF_path, SNV_VCF_path, genome.v, nparallel){
  catalogs  <- generate_catalogs_from_vcf(Indel_VCF_path, SNV_VCF_path, genome.v)
  signature_fits <- Signature_fit(catalogs$SNV_catalog, catalogs$Indel_Catalog, organ, nparallel)
  return(signature_fits)
}





#' Compute PRRDetect probability and label to file
#'
#' @param ind_fits : "Indel Signatures Fits"
#' @param snv_fits : "SNV Signatures Fits"
#' @param nparallel : Number of cores for fitting
#' @param organ : Organ fit
#' @return PRRDetect table
#' @export
PRRDetect <- function(ind_fits, snv_fits, organ=NULL){




}
