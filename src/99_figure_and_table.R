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

test_r2_conti_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- lm(y~x,d_r2_tmp)
  cd <- cooks.distance(fit_tmp)
  Ncd <- sum(cd>4)
  if(Ncd>0){
    d_r2_tmp <- d_r2_tmp[cd<=4,]
  }
  fit_tmp <- lm(y ~ x,d_r2_tmp)
  fit_tmp_0 <- lm(y ~ 1,d_r2_tmp)
  
  shapiro_test_res <- shapiro.test(resid(fit_tmp))$p.value
  cook_weisberg_res <- ncvTest(fit_tmp)$p
  #adj_r2 <- summary(fit_tmp)$adj.r.squared
  adj_r2 <- summary(fit_tmp)$r.squared
  df <- attr(logLik(fit_tmp),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  adj_r2_bootstrap <- sapply(1:B,function(i){
    fit_tmp <- lm(y ~ x+1,d_r2_tmp[sample(1:NROW(d_r2_tmp),NROW(d_r2_tmp),replace = TRUE),])
    #summary(fit_tmp)$adj.r.squared
    summary(fit_tmp)$r.squared
  })
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  
  c(shapiro_test_res,cook_weisberg_res,adj_r2,adj_r2_bootstrap_se,lrtest_p,Nt,Ncd)
}


## create dataset

load("../RLS_META_MG/to_stepwise_conditional.RData")





