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


prediction_function <- function(x, prediction_type="prop"){

  if(prediction_type == "prop"){
      coeffs <- list(MMRd = structure(c(1.65933982816886, 1.72241915581865, -0.405701704733388,
                                          -0.248442572841702, 1.38782121603953, -0.02114157396552, -0.237053140646573,
                                          0.913754960613305), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)",
                                                                                                 "MMR_SBS", "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                       MMRd_Poly_dys = structure(c(-1.08028241, -0.29327753, 0.73960286, -0.10178405,
                                                   -0.30871994,   1.72304055, 0.52953112, 0.70818278), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS", "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                       Negative = structure(c(2.30263517, -1.35968913,
                                              -0.25330441, -1.12322343, -1.42603060,
                                              -0.85715729, -0.85727534, -0.37772872), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS",
                                                                                                                         "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                       Poly_dys = structure(c(-2.88169258, 0,
                                              0, 1.57880550, 0.34692932,
                                              -0.73938623, 0.56479736, -1.24420902), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS", "MIX_SBS","POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")))


      means = c(MMR_SBS = 0.25911078, MIX_SBS = 0.02319080, POL_SBS = 0.10029001,
                MMR_IND = 0.27304617, MIX_IND = 0.06701659, POL_IND = 0.05891483,
                Ratio = -2.36772891)
      st_d = c(MMR_SBS = 0.3394553, MIX_SBS = 0.1464743, POL_SBS = 0.2864869,
               MMR_IND =  0.3523247, MIX_IND = 0.1911812, POL_IND = 0.1786840,
               Ratio = 2.3154349)
  }else if(prediction_type == "abs"){
    coeffs <- list(MMRd = structure(c(1.99458292805284, 1.67401150028469, -0.349681897486237,
                                      -0.301216297643909, 1.34323277067465, 0.101439059773288, -0.450608683176479,
                                      0.976324255733506), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)",
                                                                                               "MMR_SBS", "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                   MMRd_Poly_dys = structure(c(-3.58203452997986, -0.294668565197675, 1.19150438545967, 0.543033331452076,
                                               -0.15971132167818,   2.5302184073533, 0.0456880625396468, 0.543220163717026), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS", "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                   Negative = structure(c(2.87848513447241, -1.26969990514736,
                                          -0.220391823317334, -1.24204398705193, -1.25516183900502,
                                          -0.946533228324055, -0.931586385430911, -0.400972904291052), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS",
                                                                                                                                            "MIX_SBS", "POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")),
                   Poly_dys = structure(c(-1.29103353254539, -0.00428757640430442,
                                          -0.516075211120755, 1.00022695324376, 0.0716403900085551,
                                          -1.68512423880253, 1.33650700606774, -1.11857151515948), .Dim = c(8L, 1L), .Dimnames = list(c("(Intercept)", "MMR_SBS", "MIX_SBS","POL_SBS", "MMR_IND", "MIX_IND", "POL_IND", "Ratio"), "1")))


    means = c(MMR_SBS = 5.9005611, MIX_SBS = 0.4742382, POL_SBS = 2.1523615,
              MMR_IND = 6.2160431, MIX_IND = 1.9748333, POL_IND = 1.5059321,
              Ratio = -2.3677289)
    st_d = c(MMR_SBS = 7.645191, MIX_SBS = 2.995920, POL_SBS = 6.127789,
             MMR_IND =  7.815107, MIX_IND = 4.801253, POL_IND = 4.164313,
             Ratio = 2.315435)

  }

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

