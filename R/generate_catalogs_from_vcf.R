# TODO: check if InDel VCF has to be prepared for generate_catalogs_from_vcf
# Do we assume that the Variants are already processed? Yes ??



#' Generate catalogs from VCF files
#'
#' @param Indel_VCF_path : "Path of the InDel VCF file"
#' @param SNV_VCF_path : "Path of the SNV VCF file"
#' @param genome.v : "Genome version"
#' @return List containing the catalogs
#' @export
generate_catalogs_from_vcf <- function(Indel_VCF_path, SNV_VCF_path, genome.v){
  InDel_catalog <- indelsig.tools.lib::indel_classifier89(indels = read.table(Indel_VCF_path, sep = "\t", header = T), genome.v = genome.v)
  SNV_catalog <- signature.tools.lib::vcfToSNVcatalogue(vcfFilename = SNV_VCF_path, genome.v = genome.v)
  return(list("Indel_Catalog"=InDel_catalog, "SNV_Catalog"=SNV_catalog))
}
