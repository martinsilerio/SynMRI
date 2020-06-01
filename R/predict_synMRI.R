
bloch_tranform <- function(params, TE, TR) {

  rho <- params[1]
  T1 <- params[2]
  T2 <- params[3]

  nu <- rho*exp(-TE/T2)*(1 - exp(-TR/T1))

  return(nu)

}

predict_SynMRI <- function(model_fit, TE, TR) {

  M <- length(TE)

  dims <- dim(model_fit$rho)
  n1 <- dims[1]
  n2 <- dims[2]
  n3 <- dims[3]
  n <- n1*n2*n3

  pred <- mapply(function(TE_i, TR_i){

    nu <- array(0, dim = c(n1, n2, n3))

    for(i1 in 1:n1)
      for(i2 in 1:n2)
        for(i3 in 1:n3)
           nu[i1, i2, i3] <- bloch_tranform(c(model_fit$rho[i1, i2, i3],
                                              model_fit$T1[i1, i2, i3],
                                              model_fit$T2[i1, i2, i3]), TE_i, TR_i)

    return(nu)

  }, TE, TR, SIMPLIFY = FALSE)

  return(pred)
}
