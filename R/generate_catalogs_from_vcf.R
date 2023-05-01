# TODO: check if InDel VCF has to be prepared for generate_catalogs_from_vcf
# Do we assume that the Variants are already processed? Yes ??


InDel_sigs_order <- c("A[Ins(C):R0]A", "A[Ins(C):R0]T", "Ins(C):R(0,3)", "Ins(C):R(4,6)",
                      "Ins(C):R(7,9)", "A[Ins(T):R(0,4)]A", "A[Ins(T):R(0,4)]C", "A[Ins(T):R(0,4)]G",
                      "C[Ins(T):R(0,4)]A", "C[Ins(T):R(0,4)]C", "C[Ins(T):R(0,4)]G",
                      "G[Ins(T):R(0,4)]A", "G[Ins(T):R(0,4)]C", "G[Ins(T):R(0,4)]G",
                      "A[Ins(T):R(5,7)]A", "A[Ins(T):R(5,7)]C", "A[Ins(T):R(5,7)]G",
                      "C[Ins(T):R(5,7)]A", "C[Ins(T):R(5,7)]C", "C[Ins(T):R(5,7)]G",
                      "G[Ins(T):R(5,7)]A", "G[Ins(T):R(5,7)]C", "G[Ins(T):R(5,7)]G",
                      "A[Ins(T):R(8,9)]A", "A[Ins(T):R(8,9)]C", "A[Ins(T):R(8,9)]G",
                      "C[Ins(T):R(8,9)]A", "C[Ins(T):R(8,9)]C", "C[Ins(T):R(8,9)]G",
                      "G[Ins(T):R(8,9)]A", "G[Ins(T):R(8,9)]C", "G[Ins(T):R(8,9)]G",
                      "Ins(2,4):R0", "Ins(5,):R0", "Ins(2,4):R1", "Ins(5,):R1", "Ins(2,):R(2,4)",
                      "Ins(2,):R(5,9)", "[Del(C):R1]A", "[Del(C):R1]T", "[Del(C):R2]A",
                      "[Del(C):R2]T", "[Del(C):R3]A", "[Del(C):R3]T", "[Del(C):R(4,5)]A",
                      "[Del(C):R(4,5)]T", "[Del(C):R(1,5)]G", "Del(C):R(6,9)", "A[Del(T):R(1,4)]A",
                      "A[Del(T):R(1,4)]C", "A[Del(T):R(1,4)]G", "C[Del(T):R(1,4)]A",
                      "C[Del(T):R(1,4)]C", "C[Del(T):R(1,4)]G", "G[Del(T):R(1,4)]A",
                      "G[Del(T):R(1,4)]C", "G[Del(T):R(1,4)]G", "A[Del(T):R(5,7)]A",
                      "A[Del(T):R(5,7)]C", "A[Del(T):R(5,7)]G", "C[Del(T):R(5,7)]A",
                      "C[Del(T):R(5,7)]C", "C[Del(T):R(5,7)]G", "G[Del(T):R(5,7)]A",
                      "G[Del(T):R(5,7)]C", "G[Del(T):R(5,7)]G", "A[Del(T):R(8,9)]A",
                      "A[Del(T):R(8,9)]C", "A[Del(T):R(8,9)]G", "C[Del(T):R(8,9)]A",
                      "C[Del(T):R(8,9)]C", "C[Del(T):R(8,9)]G", "G[Del(T):R(8,9)]A",
                      "G[Del(T):R(8,9)]C", "G[Del(T):R(8,9)]G", "Del(2,4):R1", "Del(5,):R1",
                      "Del(2,8):U(1,2):R(2,4)", "Del(2,):U(1,2):R(5,9)", "Del(3,):U(3,):R2",
                      "Del(3,):U(3,):R(3,9)", "Del(2,5):M1", "Del(3,5):M2", "Del(4,5):M(3,4)",
                      "Del(6,):M1", "Del(6,):M2", "Del(6,):M3", "Del(6,):M(4,)", "Complex"
)


#' Generate catalogs from VCF files
#'
#' @param Indel_VCF_path : "Path of the InDel VCF file"
#' @param SNV_VCF_path : "Path of the SNV VCF file"
#' @param genome.v : "Genome version"
#' @return List containing the catalogs
#' @export
generate_catalogs_from_vcf <- function(Indel_VCF_path, SNV_VCF_path, genome.v, sample_name){
  InDel_VCF <- read.table(Indel_VCF_path)[c(3,1,2,4,5)]
  InDel_VCF$V3 <- sample_name
  colnames(InDel_VCF) <- c("Sample", "chr", "position", "REF", "ALT")

  InDel_catalog <- indelsig.tools.lib::indel_classifier89(indels = InDel_VCF, genome.v = genome.v)
  InDel_catalog <- indelsig.tools.lib::gen_catalogue89(InDel_catalog, sample_col = "Sample")
  InDel_catalog <- InDel_catalog[InDel_sigs_order,, drop=F]

  colnames(InDel_catalog) <- sample_name
  SNV_catalog <- signature.tools.lib::vcfToSNVcatalogue(vcfFilename = SNV_VCF_path, genome.v = genome.v)$catalogue
  colnames(SNV_catalog) <- sample_name
  return(list("Indel_Catalog"=InDel_catalog, "SNV_Catalog"=SNV_catalog, "total_SNV"=colSums(SNV_catalog), "total_InD"=nrow(InDel_VCF)))
}
