Signature_fit <- function(SNV_Catalog, InDel_Catalog, organ, nparallel){

  SNV_fit <- signature.tools.lib::FitMS(catalogues = SNV_Catalog, organ = organ,useBootstrap = T, nparallel = nparallel )

  InDel_fit <- signature.tools.lib::Fit(catalogues = InDel_Catalog, signatures = read.table("data/InDel_Consensus_Sigs.tsv", sep = "\t"), useBootstrap = T, nparallel= nparallel)

  return(list("SNV"=SNV_fit, "InDel"=InDel_fit))
}
