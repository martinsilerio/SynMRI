
#' @export
#' @useDynLib SynMRI fit_LSE_
fit_model_SynMRI <- function(R, TE, TR){

  # Check format is correct

  M <- length(R)
  dims <- dim(R[[1]])
  n1 <- dims[1]
  n2 <- dims[2]
  n3 <- dims[3]
  n <- n1*n2*n3
  # Change format
  Rvect <- lapply(R, c)
  r_ij <- do.call(cbind, Rvect)
  r <- as.vector(r_ij)

  # Here we set a seed to match the one we used in C
  RNGkind("Marsaglia-Multicarry")
  .Random.seed<-c(401L,111L,222L)

  tot_zeros <- sum(r==0)

  # Here we set the intensity values equal to zero to some very small value
  r[r==0] <- runif(tot_zeros, 0, 0.1)

  # Here we set a seed to match the one we used in C
  #RNGkind("Marsaglia-Multicarry")
  #.Random.seed<-c(401L,111L,222L)

  res <- .Call(fit_LSE_, M, n, n1, n2, n3, r, TE, TR)

  W1 <- array(res[[1]][1:n], c(n1, n2, n3))
  W2 <- array(res[[1]][(n+1):(2*n)], c(n1, n2, n3))
  W3 <- array(res[[1]][(2*n+1):(3*n)], c(n1, n2, n3))

  rho <- array(res[[2]][1:n], c(n1, n2, n3))
  T1 <- array(res[[2]][(n+1):(2*n)], c(n1, n2, n3))
  T2 <- array(res[[2]][(2*n+1):(3*n)], c(n1, n2, n3))

  list(W1 = W1, W2 = W2, W3 = W3, rho = rho, T1 = T1, T2 = T2)

}
