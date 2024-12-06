#' prepare_PRRDetect
#'
#' Prepares data for PRRDetect input
#'
#' @param InDel_Catalogs: Indel Catalogs dataframe
#' @param SNV_Catalogs: SNV catalogs daframe
#' @param Indel_VCF_path: Vector of paths or list of data frames. If it is a vector of paths the file will be rearranged by the function, otherwise it has to be a list of dataframe having the following structure: "Sample", "chr", "position", "REF", "ALT"
#' @param total_InDel: Total number of indels, it can be defined previously
#' @param SNV_VCF_path: Vector of paths or list of data frames.
#' @param genome.v: Either "hg19" or "hg38", which is the default option.
#' @param nparallel: Number of threads for signature fit
#' @param is.filtered_SNV: if yes, SNV VCF get filtered. Default equal to False. SNV information gets filtered only if provided as a vector of paths.
#' @param is.filtered_InDel: if yes InDel VCF get filtered. Default equal to False. InDel information gets filtered only if provided as a vector of paths.
#' @param sample_name: Vector of sample names, it has to be in the same order of Indel VCF vector or list.
#' @param organ: Organ for SBS signature fit
#' @param setseed: Seed for reproducibility purposes
#' @return A list containing the indel and snv fits, the total number of SNVs, total number of InDels. If InDel VCF paths are provided, the total number of InDels is the number of InDels in the VCF file, otherwise it is the sum of the mutations found in the catalog.
#' @export
prepare_PRRDetect <- function(InDel_Catalogs=NULL,SNV_Catalogs=NULL, Indel_VCF_path=NULL, total_InDel=NULL, SNV_VCF_path=NULL, genome.v="hg38", nparallel=NULL, is.filtered_SNV=T, is.filtered_InDel=T, sample_name=NULL, organ=NULL, setseed = 42){


  if(is.null(InDel_Catalogs)==F & is.null(Indel_VCF_path)==F){

    message("Either InDel Catalogs or VCF information has to be provided")
    exit()

  }

  if(is.null(SNV_Catalogs)==F & is.null(SNV_VCF_path)==F){

    message("Either SNV Catalogs or VCF information has to be provided")
    exit()

  }

  if(is.null(Indel_VCF_path)==F){

    if(is.atomic(Indel_VCF_path) & is.filtered_InDel == F){

      InDelFilter(vcffile = Indel_VCF_path, vcfout = paste0(Indel_VCF_path, ".filtered"), genomev = genome.v, sample_name = sample_name )
      Indel_VCF_path <- paste0(Indel_VCF_path, ".filtered.bgz")

    }

  }


  if(is.null(SNV_VCF_path)==F){

    if(is.atomic(SNV_VCF_path) & is.filtered_SNV == F){

      SNV_Filter(SNV_VCF_path, vcfout = paste0(SNV_VCF_path, ".filtered"), genomev = genome.v )
      SNV_VCF_path <- paste0(SNV_VCF_path, ".filtered.bgz")

    }

  }


  if( is.null(Indel_VCF_path) == F | is.null(SNV_VCF_path)==F ){

     catalogs  <- generate_catalogs_from_mutations(Indel_VCF = Indel_VCF_path, SNV_VCF =  SNV_VCF_path,genome.v =  genome.v, sample_name = sample_name)
  }


  if(is.null(InDel_Catalogs)==F & is.null(sample_name)){

    sample_name = colnames(InDel_Catalogs)
  }else if(!is.null(sample_name)){
    colnames(InDel_Catalogs) <- sample_name
  }


  if(is.null(SNV_Catalogs)==F & is.null(sample_name)){
    sample_name= colnames(SNV_Catalogs)
  }else if(!is.null(sample_name)){
    colnames(SNV_Catalogs) <- sample_name
  }





  if(is.null(SNV_VCF_path) == F){

    SNVC <- catalogs$SNV_Catalog

  }else{

    SNVC <- SNV_Catalogs

  }

  if( is.null(Indel_VCF_path)==F){

    INDC <- catalogs$Indel_Catalog

  }else{

    INDC <- InDel_Catalogs

  }

  signature_fits <- Signature_fit(SNV_Catalog = SNVC , InDel_Catalog = INDC, organ, nparallel, setseed)

  SNV_fit <- signature_fits$SNV$exposures

  InDel_fit <- signature_fits$InDel$exposures

  colnames(SNV_fit) <- unlist(lapply(strsplit(colnames(SNV_fit), split = "_" , fixed = T), function(x){tail(x,1)}))

  outlist <- list(SNV=SNV_fit, InDel=InDel_fit)

  if(exists("catalogs$total_SNV")){
    outlist$total_SNV = catalogs$total_SNV
  }else{
    outlist$total_SNV = colSums(SNV_Catalogs)
  }

  if(exists("catalogs$total_InD")){
    outlist$total_InDel = catalogs$total_InD
  }else{
    outlist$total_InDel = colSums(InDel_Catalogs)
  }

  return(outlist)

}



#' Compute PRRDetect probability and label
#'
#' @param InDel_fits: Indel Signatures Fits
#' @param SNV_fits: SNV Signatures Fits
#' @param total_SNV: total number of SNVs of the sample
#' @param total_InD: total number of InDel mutation in the sample
#' @param prediction_type: Models to use, "prop" is the one published. Otherwise "abs" for other model
#' @return PRRDetect table
#' @export
PRRDetect <- function(InDel_fits, SNV_fits, total_SNV, total_InD, prediction_type="prop"){

  pure_MMR_sigs <- paste0("SBS", c(6,15,26,44,97))
  pure_POL_sigs <- paste0("SBS", c("10a", "10d"))
  mixed_sigs <- paste0("SBS", c(14,20))

  pure_MMR_indS <- paste0("RefSig.InD", c(7,19))
  pure_POL_indS <- paste0("RefSig.InD", c(14,15))
  mixed_indS <- paste0("RefSig.InD", c("16a", "16b", 20, 21))

  ## if the name of the signature is not present, which returns "integer(0)", then the sum is equal to 0, otherwise is equal to the rowsum
  if(prediction_type == "prop"){

  ## SNV
      if(identical(which(colnames(SNV_fits) %in% pure_MMR_sigs), integer(0))){

          MMR_SBS <- rep(0,nrow(SNV_fits))

      }else{

          MMR_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% pure_MMR_sigs), drop=F])/total_SNV

      }

      if(identical(which(colnames(SNV_fits) %in% pure_POL_sigs), integer(0))){

          POL_SBS <- rep(0,nrow(SNV_fits))

      }else{

          POL_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% pure_POL_sigs), drop=F])/total_SNV

      }

      if(identical(which(colnames(SNV_fits) %in% mixed_sigs), integer(0))){

          MIX_SBS <- rep(0,nrow(SNV_fits))

      }else{

          MIX_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% mixed_sigs), drop=F])/total_SNV

      }

      if(identical(which(colnames(InDel_fits) %in% pure_MMR_indS), integer(0))){

        MMR_IND <- rep(0,nrow(SNV_fits))

      }else{

        MMR_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% pure_MMR_indS), drop=F])/total_InD

      }

      if(identical(which(colnames(InDel_fits) %in% pure_POL_indS), integer(0))){

        POL_IND <- rep(0,nrow(SNV_fits))

      }else{

        POL_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% pure_POL_indS), drop=F])/total_InD

      }

      if(identical(which(colnames(InDel_fits) %in% mixed_indS), integer(0))){

        MIX_IND <- rep(0,nrow(SNV_fits))

      }else{

        MIX_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% mixed_indS), drop=F])/total_InD

      }

  }else if(prediction_type == "abs"){

    ## SNV
    if(identical(which(colnames(SNV_fits) %in% pure_MMR_sigs), integer(0))){

      MMR_SBS <- rep(0,nrow(SNV_fits))

    }else{

      MMR_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% pure_MMR_sigs), drop=F])

    }

    if(identical(which(colnames(SNV_fits) %in% pure_POL_sigs), integer(0))){

      POL_SBS <- rep(0,nrow(SNV_fits))

    }else{

      POL_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% pure_POL_sigs), drop=F])

    }

    if(identical(which(colnames(SNV_fits) %in% mixed_sigs), integer(0))){

      MIX_SBS <- rep(0,nrow(SNV_fits))

    }else{

      MIX_SBS <- rowSums(SNV_fits[,which(colnames(SNV_fits) %in% mixed_sigs), drop=F])

    }

    if(identical(which(colnames(InDel_fits) %in% pure_MMR_indS), integer(0))){

      MMR_IND <- rep(0,nrow(SNV_fits))

    }else{

      MMR_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% pure_MMR_indS), drop=F])

    }

    if(identical(which(colnames(InDel_fits) %in% pure_POL_indS), integer(0))){

      POL_IND <- rep(0,nrow(SNV_fits))

    }else{

      POL_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% pure_POL_indS), drop=F])

    }

    if(identical(which(colnames(InDel_fits) %in% mixed_indS), integer(0))){

      MIX_IND <- rep(0,nrow(SNV_fits))

    }else{

      MIX_IND <- rowSums(InDel_fits[,which(colnames(InDel_fits) %in% mixed_indS), drop=F])

    }

  }

  Ratio <-  c(total_InD + 1)/(total_SNV + 1)

  final_df <- as.data.frame(cbind(MMR_SBS, MIX_SBS, POL_SBS, MMR_IND, MIX_IND, POL_IND, Ratio))

  final_df[,1:6] <- log2(final_df[,1:6]+1)
  final_df[,7] <- log2(final_df[,7])

  final_df <- prediction_function(final_df, prediction_type = prediction_type)

  colnames(final_df)[1:4] <- c("MMRd", "MMRd+Poly-dys", "Negative", "Poly-dys")
  final_df$Prediction = prediction_label(final_df)
  final_df <- final_df[,c(3,1,2,4,5)]
  return(final_df)

}


