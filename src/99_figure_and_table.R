# this is a recovered version for previous lost file

# functions 

irnt_f <- function(x){
  res <- x
  res[!is.na(res)] <- qnorm(rank(res[!is.na(res)])/(max(rank(res[!is.na(res)]))+0.5))
}

irnt_df_f <- function(x){
  apply(x,2,irnt_f)
}


r2_nagelkerke_f <- function(model, l_base) {
  L.full <- insight::get_loglikelihood(model)
  D.full <- -2 * L.full
  
  D.base <- -2 * l_base
  # Is it still necessary?
  if (inherits(model, c("vglm", "vgam", "clm2"))) {
    n <- insight::n_obs(model)
  } else {
    n <- attr(L.full, "nobs")
    if (is.null(n)) n <- insight::n_obs(model)
  }
  
  r2_nagelkerke <- as.vector((1 - exp((D.full - D.base) / n)) / (1 - exp(-D.base / n)))
  
  names(r2_nagelkerke) <- "Nagelkerke's R2"
  r2_nagelkerke
}

