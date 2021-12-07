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



test_r2_disc_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- lm(y~x,d_r2_tmp)
  cd <- cooks.distance(fit_tmp)
  Ncd <- sum(cd>4 | is.na(cd))
  if(Ncd>0){
    d_r2_tmp <- d_r2_tmp[which(cd<=4),]
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
    balanced_resampling <- tapply(1:NROW(d_r2_tmp),d_r2_tmp$x,function(xin){sample(xin,NROW(xin),replace = T)})
    balanced_resampling <- unlist(balanced_resampling)
    fit_tmp <- lm(y ~ x+1,d_r2_tmp[balanced_resampling,])
    #summary(fit_tmp)$adj.r.squared
    summary(fit_tmp)$r.squared
  })
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  
  c(shapiro_test_res,cook_weisberg_res,adj_r2,adj_r2_bootstrap_se,lrtest_p,Nt,Ncd)
}

test_r2_discBIN_outcome_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- glm(x~y,d_r2_tmp,family = binomial(link="logit"))
  cd <- cooks.distance(fit_tmp)
  Ncd <- sum(cd>4 | is.na(cd))
  if(Ncd>0){
    d_r2_tmp <- d_r2_tmp[which(cd<=4),]
  }
  fit_tmp <- glm(x ~ y,d_r2_tmp,family = binomial(link="logit"))
  fit_tmp_0 <- glm(x ~ 1,d_r2_tmp,family = binomial(link="logit"))
  
  df <- attr(logLik(fit_tmp),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  adj_r2 <- r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  
  converged <- fit_tmp$converged
  B=50
  cl <- makePSOCKcluster(rep("localhost",2))
  adj_r2_bootstrap <- parSapply(cl,1:B,function(i,
                                                d_r2_tmp=d_r2_tmp,
                                                balanced_resampling=balanced_resampling,
                                                r2_nagelkerke=r2_nagelkerke){
    balanced_resampling <- tapply(1:NROW(d_r2_tmp),d_r2_tmp$x,function(xin){sample(xin,NROW(xin),replace = T)})
    balanced_resampling <- unlist(balanced_resampling)
    fit_tmp <- glm(x ~ y+1,d_r2_tmp[balanced_resampling,],family = binomial(link="logit"))
    fit_tmp_0 <- glm(x ~ 1,d_r2_tmp[balanced_resampling,],family = binomial(link="logit"))
    #summary(fit_tmp)$adj.r.squared
    r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  },
  d_r2_tmp=d_r2_tmp,
  balanced_resampling=balanced_resampling,
  r2_nagelkerke=r2_nagelkerke)
  stopCluster(cl)
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  c(converged,NA,adj_r2,adj_r2_bootstrap_se,lrtest_p,Nt,Ncd)
}

test_r2_disccateg_outcome_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- multinom(x~y,d_r2_tmp,trace=FALSE)
  
  # fit_tmp_res <- rowSums(resid(fit_tmp)*(model.matrix(y0~as.factor(x0)+0)))
  # 
  # cd <- outliers::scores(fit_tmp_res,type = "mad",prob = TRUE)
  # Ncd <- sum(cd<0.005 | is.na(cd))
  # if(Ncd>0){
  #   d_r2_tmp <- d_r2_tmp[which(cd>0.005),]
  # }
  Ncd <- 0
  fit_tmp <- multinom(x ~ y,d_r2_tmp,trace=FALSE)
  fit_tmp_0 <- multinom(x ~ 1,d_r2_tmp,trace=FALSE)
  
  df <- attr(logLik(fit_tmp),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  adj_r2 <- r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  
  converged <- fit_tmp$convergence==0
  
  adj_r2_bootstrap <- sapply(1:B,function(i){
    balanced_resampling <- tapply(1:NROW(d_r2_tmp),d_r2_tmp$x,function(xin){sample(xin,NROW(xin),replace = T)})
    balanced_resampling <- unlist(balanced_resampling)
    fit_tmp <- multinom(x ~ y+1,d_r2_tmp[balanced_resampling,],trace=FALSE)
    fit_tmp_0 <- multinom(x ~ 1,d_r2_tmp[balanced_resampling,],trace=FALSE)
    #summary(fit_tmp)$adj.r.squared
    r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  })
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  c(converged,NA,adj_r2,adj_r2_bootstrap_se,lrtest_p,Nt,Ncd)
}


test_r2_discordinal_outcome_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- polr(x~y,d_r2_tmp,method = c("logistic"))
  x_surrogate <- surrogate(fit_tmp)
  fit_tmp_l <- lm(x_surrogate~d_r2_tmp$y)
  cd <- cooks.distance(fit_tmp_l)
  Ncd <- sum(cd>4 | is.na(cd))
  if(Ncd>0){
    d_r2_tmp <- d_r2_tmp[which(cd<=4),]
  }
  fit_tmp <- polr(x ~ y,d_r2_tmp,method = c("logistic"))
  fit_tmp_0 <- polr(x ~ 1,d_r2_tmp,method = c("logistic"))
  
  df <- attr(logLik(fit_tmp),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  adj_r2 <- r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  
  converged <- fit_tmp$convergence==0
  
  adj_r2_bootstrap <- sapply(1:B,function(i){
    balanced_resampling <- tapply(1:NROW(d_r2_tmp),d_r2_tmp$x,function(xin){sample(xin,NROW(xin),replace = T)})
    balanced_resampling <- unlist(balanced_resampling)
    fit_tmp <- polr(x ~ y+1,d_r2_tmp[balanced_resampling,],method = c("logistic"))
    fit_tmp_0 <- polr(x ~ 1,d_r2_tmp[balanced_resampling,],method = c("logistic"))
    #summary(fit_tmp)$adj.r.squared
    r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  })
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  c(converged,NA,adj_r2,adj_r2_bootstrap_se,lrtest_p,Nt,Ncd)
}


## load old / create dataset

load("../RLS_META_MG/to_stepwise_conditional.RData")

## test MG+GO+DG3 test as exposure





