SNV_Filter <- function(vcffile,vcfout,genomev,sample_name,clpm_max=0,asmd_min=140,filterPASS=T){
  tempvcf <- tempfile()
  system(paste0("zgrep -v vcfProcessLog ",vcffile," > ",tempvcf))

  subs_VCF <- VariantAnnotation::readVcf(tempvcf,genome = genomev)
  e.snv <- VariantAnnotation::expand(subs_VCF)

  if(filterPASS){
    selected_snv <- VariantAnnotation::fixed(e.snv)[,"FILTER"]=="PASS"
    e.snv <- e.snv[selected_snv,]
  }

  selected_snv <- VariantAnnotation::info(e.snv)$CLPM<=clpm_max & VariantAnnotation::info(e.snv)$ASMD>=asmd_min
  e.snv <- e.snv[selected_snv,]
  VariantAnnotation::writeVcf(e.snv,vcfout,index=TRUE)
}

InDelFilter <- function(vcffile,vcfout,genomev,sample_name,qualmin=250,repmax=9,filterPASS=T){
  indels_VCF <- VariantAnnotation::readVcf(vcffile,genome = genomev)
  e.indels <- VariantAnnotation::expand(indels_VCF)
  if(filterPASS){
    selected_indels <- VariantAnnotation::fixed(e.indels)[,"FILTER"]=="PASS"
    e.indels <- e.indels[selected_indels,]
  }
  selected_indels <- VariantAnnotation::fixed(e.indels)[,"QUAL"]>=qualmin & VariantAnnotation::info(e.indels)$REP<=repmax
  e.indels <- e.indels[selected_indels,]
  VariantAnnotation::writeVcf(e.indels,vcfout,index=TRUE)
}
