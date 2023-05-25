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

#' Generate catalogs from VCF data
#'
#' This function returns the catalogs and the total number of mutations from VCF files. You can give just InDel VCF param or SNV VCF param.
#'
#' @param Indel_VCF: Vector of paths or list of data frames. If it is a vector of paths the file will be rearranged by the function, otherwise it has to be a list of dataframe having the following structure: "Sample", "chr", "position", "REF", "ALT"
#' @param SNV_VCF:  Vector of paths or list of data frames.
#' @param genome.v: Either "hg19" or "hg38", which is the default option.
#' @param sample_name: Vector of sample names, it has to be in the same order of Indel VCF vector or list.
#' @return Returns a list containing the Indel Catalos, total InD and/or SNV catalogs, total SNV
#' @export
generate_catalogs_from_mutations <- function(Indel_VCF=NULL, SNV_VCF=NULL, genome.v="hg38", sample_name=NULL){
  ## if it is a list then it is a list of dataframes
  ## if it is a vector is anything else
  outlist = list()

  library(indelsig.tools.lib)
  class_and_catalog <- function(x, genome.v){
    InDel_catalog <- indelsig.tools.lib::indel_classifier89(indels = x, genome.v = genome.v)
    InDel_catalog <- indelsig.tools.lib::gen_catalogue89(InDel_catalog, sample_col = "Sample")
    InDel_catalog <- InDel_catalog[InDel_sigs_order,, drop=F]
    return(InDel_catalog)
  }

  if(is.null(Indel_VCF)==F){

    if(is.atomic(Indel_VCF)){
      Indel_VCF <- apply(array(Indel_VCF),1, function(y){
        w <- read.table(y)[c(3,1,2,4,5)]
        w$V3 <- basename(y)
        colnames(w) <- c("Sample", "chr", "position", "REF", "ALT")
        return(w)
      })
    }

    nrow_InDel_catalogs <- unlist(lapply(Indel_VCF,function(y){
      return(nrow(y))
    }))

    Indel_catalogs <- lapply(Indel_VCF, function(y){class_and_catalog(y, genome.v = genome.v)})

    Indel_catalogs <- as.data.frame(do.call("cbind",Indel_catalogs))

    if(is.null(sample_name)==F){
      colnames(Indel_catalogs) <- sample_name
    }
    names(nrow_InDel_catalogs) <- colnames(Indel_catalogs)
    outlist$Indel_Catalog = Indel_catalogs
    outlist$total_InD = nrow_InDel_catalogs
  }


  if(is.null(SNV_VCF)==F){
    if(is.atomic(SNV_VCF)){
      SNV_catalogs <- lapply(SNV_VCF, function(x){signature.tools.lib::vcfToSNVcatalogue(vcfFilename = x, genome.v = genome.v)$catalogue})
      SNV_catalogs <- do.call("cbind", SNV_catalogs)
      if(is.null(sample_name)==F){
        colnames(SNV_catalogs) <- sample_name
      }
    }else if(is(SNV_VCF, "list")){
      SNV_catalogs <- lapply(SNV_VCF, function(x) {
        signature.tools.lib::tabToSNVcatalogue(vcfFilename = x, genome.v = genome.v)$catalogue})
      SNV_catalogs <- do.call("cbind", SNV_catalogs)
      if(is.null(sample_name)==F){
        colnames(SNV_catalogs) <- sample_name
      }
    }

    outlist$SNV_Catalog=SNV_catalogs

    outlist$total_SNV=colSums(SNV_catalogs)
  }

  return(outlist)
}
