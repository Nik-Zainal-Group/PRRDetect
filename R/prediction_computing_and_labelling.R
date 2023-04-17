prediction_label <- function(w){
  apply(w[,1:4], 1, function(x){
    if(x["Negative"] > 0.5){
      return("Negative")
    }else{
      y = x[-which(names(x)=="Negative")]
      return(names(y)[which.max(y)])
    }
  }
  )
}


prediction_function <- function(x){
  coeffs <- list(MMRd = structure(c(2.55352163025099, 1.56823408939836, -0.346000966115436,
                                    -0.322362286532779, 1.23234257002924, 0.0453901226090131, -0.450877294206874,
                                    0.939723567976693), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)",
                                                                                             "MMR_SBS", "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                 MMRd_Poly_dys = structure(c(-2.53798408516766, -0.301979436643312, 1.05348660645481, 0.544253405127099,
                                             -0.0282833655889912,   2.27409090468062, 0.149011924311783, 0.493827923651865), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS", "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                 Negative = structure(c(1.2633033174444, -1.12627529102184,
                                        -0.202422944197919, -1.20872809953331, -1.10450089320956,
                                        -0.926012003184389, -0.872800758244794, -0.404926411020721), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS",
                                                                                                                                          "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                 Poly_dys = structure(c(-1.27884086252774, -0.0345918527521169,
                                        -0.399675187160369, 0.986836980938991, -0.001795005825622,
                                        -1.39346902410524, 1.17466612813989, -1.02862508060784), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS", "MIX_SBS","POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")))


  means = c(MMR_SBS = 7.94627449171039, MIX_SBS = 0.638655735499494, POL_SBS = 2.89858120355565,
            MMR_IND = 8.37113346307278, MIX_IND = 2.65950427832639, POL_IND = 2.02803595358004,
            Ratio = -2.13509082353101)
  st_d = c(MMR_SBS = 7.90359191873798, MIX_SBS = 3.46257726866702, POL_SBS = 6.95923315361267,
           MMR_IND = 8.01359221128261, MIX_IND = 5.40720196788962, POL_IND = 4.72299044617658,
           Ratio = 2.61021885269447)


  ####################################################################################

  value_df <- as.data.frame(apply(array(colnames(x)),1, function(y){
    return((x[y] - means[y])/st_d[y])
  }))
  colnames(value_df) <- lapply(strsplit(colnames(value_df),split = ".", fixed=T), function(x){x[1]})

  if(any(rownames(coeffs[["MMRd"]])[2:8] != colnames(value_df))){
    message("different rownames")
    exit
  }

  prediction_formula <- function(value_df, coeffs){
    values = c()
    for(class in names(coeffs)){
      values <- c(values , exp(coeffs[[class]]["(Intercept)",] +  as.vector(t(coeffs[[class]][2:8,]) %*% as.numeric(value_df))))
    }
    names(values) <- names(coeffs)
    return(values/sum(values))
  }

  out_obj =  as.data.frame(t(apply(value_df, 1, function(z){prediction_formula(z, coeffs)})))
  colnames(out_obj) <- c("MMRd", "MMRd_Poly_dys", "Negative", "Poly_dys")
  return(out_obj)
}
