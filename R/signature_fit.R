#' Signature fit
#' @param SNV_Catalog: SNV Catalog
#' @param InDel_Catalog: Indel Catalog
#' @param organ: Organ necessary for SNV signature fit
#' @param nparallel: Number of threads
#' @param setseed: seed for reproducibility
#' @return List containing SNV fit and Indel fit
Signature_fit <- function(SNV_Catalog, InDel_Catalog, organ, nparallel, setseed){

  set.seed(setseed)

  outlist <- list()
  if( is.null(SNV_Catalog)==F ){

      SNV_fit <- signature.tools.lib::FitMS(catalogues = SNV_Catalog, organ = organ,useBootstrap = T, nparallel = nparallel, commonSignatureTier = "T2"  )
      outlist$SNV=SNV_fit

  }

  if( is.null(InDel_Catalog)==F ){

    InDel_fit <- signature.tools.lib::Fit(catalogues = InDel_Catalog, signatures = read.table("data/InDel_Consensus_Sigs.tsv", sep = "\t"), useBootstrap = T, nboot = 200,  nparallel = nparallel)
    outlist$InDel=InDel_fit
  }

  return(outlist)

}
