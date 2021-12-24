# this is a recovered version for previous lost file

# functions 

irnt_f <- function(x){
  res <- x
  res[!is.na(res)] <- qnorm(rank(res[!is.na(res)])/(max(rank(res[!is.na(res)]))+0.5))
  res
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


test_r2_conti_3_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- lm(y~.,d_r2_tmp)
  cd <- cooks.distance(fit_tmp)
  Ncd <- sum(cd>4)
  if(Ncd>0){
    d_r2_tmp <- d_r2_tmp[cd<=4,]
  }
  fit_tmp <- lm(y ~ .,d_r2_tmp)
  fit_tmp_0 <- lm(y ~ 1,d_r2_tmp)
  
  shapiro_test_res <- shapiro.test(resid(fit_tmp))$p.value
  cook_weisberg_res <- ncvTest(fit_tmp)$p
  #adj_r2 <- summary(fit_tmp)$adj.r.squared
  adj_r2 <- summary(fit_tmp)$r.squared
  df <- attr(logLik(fit_tmp),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  adj_r2_bootstrap <- sapply(1:B,function(i){
    fit_tmp <- lm(y ~ .,d_r2_tmp[sample(1:NROW(d_r2_tmp),NROW(d_r2_tmp),replace = TRUE),])
    #summary(fit_tmp)$adj.r.squared
    summary(fit_tmp)$r.squared
  })
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  
  c(shapiro_test_res,cook_weisberg_res,adj_r2,adj_r2_bootstrap_se,lrtest_p,Nt,Ncd)
}

test_r2_discordinal_outcome_3_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- polr(x~.,d_r2_tmp,method = c("logistic"))
  x_surrogate <- surrogate(fit_tmp)
  fit_tmp_l <- lm(x_surrogate~.,data.frame(x_surrogate,d_r2_tmp[,1:3]))
  cd <- cooks.distance(fit_tmp_l)
  Ncd <- sum(cd>4 | is.na(cd))
  if(Ncd>0){
    d_r2_tmp <- d_r2_tmp[which(cd<=4),]
  }
  fit_tmp <- polr(x ~ .,d_r2_tmp,method = c("logistic"))
  fit_tmp_0 <- polr(x ~ 1,d_r2_tmp,method = c("logistic"))
  
  df <- attr(logLik(fit_tmp),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  adj_r2 <- r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  
  converged <- fit_tmp$convergence==0
  
  adj_r2_bootstrap <- sapply(1:B,function(i){
    balanced_resampling <- tapply(1:NROW(d_r2_tmp),d_r2_tmp$x,function(xin){sample(xin,NROW(xin),replace = T)})
    balanced_resampling <- unlist(balanced_resampling)
    fit_tmp <- polr(x ~ .,d_r2_tmp[balanced_resampling,],method = c("logistic"))
    fit_tmp_0 <- polr(x ~ 1,d_r2_tmp[balanced_resampling,],method = c("logistic"))
    #summary(fit_tmp)$adj.r.squared
    r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  })
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  c(converged,NA,adj_r2,adj_r2_bootstrap_se,lrtest_p,Nt,Ncd)
}


test_r2_disccateg_outcome_3_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- multinom(x~.,d_r2_tmp,trace=FALSE)
  
  # fit_tmp_res <- rowSums(resid(fit_tmp)*(model.matrix(y0~as.factor(x0)+0)))
  # 
  # cd <- outliers::scores(fit_tmp_res,type = "mad",prob = TRUE)
  # Ncd <- sum(cd<0.005 | is.na(cd))
  # if(Ncd>0){
  #   d_r2_tmp <- d_r2_tmp[which(cd>0.005),]
  # }
  Ncd <- 0
  fit_tmp <- multinom(x ~ .,d_r2_tmp,trace=FALSE)
  fit_tmp_0 <- multinom(x ~ 1,d_r2_tmp,trace=FALSE)
  
  df <- attr(logLik(fit_tmp),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  adj_r2 <- r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  
  converged <- fit_tmp$convergence==0
  
  adj_r2_bootstrap <- sapply(1:B,function(i){
    balanced_resampling <- tapply(1:NROW(d_r2_tmp),d_r2_tmp$x,function(xin){sample(xin,NROW(xin),replace = T)})
    balanced_resampling <- unlist(balanced_resampling)
    fit_tmp <- multinom(x ~ .,d_r2_tmp[balanced_resampling,],trace=FALSE)
    fit_tmp_0 <- multinom(x ~ 1,d_r2_tmp[balanced_resampling,],trace=FALSE)
    #summary(fit_tmp)$adj.r.squared
    r2_nagelkerke(fit_tmp,logLik(fit_tmp_0))
  })
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  c(converged,NA,adj_r2,adj_r2_bootstrap_se,lrtest_p,Nt,Ncd)
}

test_r2_discBIN_outcome_3_f <- function(d_r2_tmp,B=200){
  d_r2_tmp <- d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
  Nt <- nrow(d_r2_tmp)
  fit_tmp <- glm(x~.,d_r2_tmp,family = binomial(link="logit"))
  cd <- cooks.distance(fit_tmp)
  Ncd <- sum(cd>4 | is.na(cd))
  if(Ncd>0){
    d_r2_tmp <- d_r2_tmp[which(cd<=4),]
  }
  fit_tmp <- glm(x ~ . ,d_r2_tmp,family = binomial(link="logit"))
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
    fit_tmp <- glm(x ~ .,d_r2_tmp[balanced_resampling,],family = binomial(link="logit"))
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

f_get_x <- function(a0,feature_order=feature_order,res_r2_mg=res_r2_mg){
  fe <- feature_order[a0,2]
  
  x0 <- d_kora_analysis[,fe]
  
  type <- feature_order[a0,4]
  trans <- res_r2_mg[1,a0]
  
  
  if(grepl("_log",trans)){
    if(min(x0,na.rm=T)<=0){
      x0 <- x0 + 1
    }
    x0 <- log2(x0)
  }else if(grepl("_irn",trans)){
    x0 <- irnt_f(x0)
  }else if(type!="continous"){
    x0 <- as.factor(x0)
  }
  x0
}



test_r2_cpls_3_f <- function(X,Y,B=200){
  
  NX = rank.condition(X)$condition
  NY = rank.condition(Y)$condition
  
  cc_res <- cppls(Y~X,1,data=list(X=X,Y=Y))
  X_proj <- cc_res$scores[,1]
  cc_res_comp1 <- cancor(X_proj,Y)
  Y_proj <- Y%*%cc_res_comp1$ycoef[,1]
  fit_tmp_1 <- lm(Y_proj~X_proj)
  cd_outlier <- cooks.distance(fit_tmp_1)>2
  
  X <- X[!(cd_outlier ),,drop=FALSE]
  Y <- Y[!(cd_outlier ),]
  
  NX = rank.condition(X)$condition
  NY = rank.condition(Y)$condition
  
  cc_res <- cppls(Y~X,1,data=list(X=X,Y=Y))
  X_proj <- cc_res$scores[,1]
  cc_res_comp1 <- cancor(X_proj,Y)
  Y_proj <- Y%*%cc_res_comp1$ycoef[,1]
  adj_r2 <- cc_res_comp1$cor^2
  
  lrtest_p <- cor.test(X_proj,Y_proj,method="spearman",exact = FALSE)$p.value
  N <- NROW(X)
  
  yi_cut <- apply(Y,2,function(xi){
    cut(xi,breaks = c(-Inf,unique(quantile(xi))))
  })
  yi_cut_q <- apply(yi_cut,1,paste0,collapse="_")
  
  adj_r2_bootstrap <- sapply(1:B,function(i){
    balanced_resampling <- tapply(1:NROW(X),yi_cut_q,function(xin){sample(xin,NROW(xin),replace = T)})
    balanced_resampling <- unlist(balanced_resampling)
    X_boot=X[balanced_resampling,,drop=FALSE]
    Y_boot=Y[balanced_resampling,,drop=FALSE]
    cc_res <- cppls(Y~X,1,data=list(X=X_boot,Y=Y_boot))
    X_proj <- cc_res$scores[,1]
    cc_res_comp1 <- cancor(X_proj,Y_boot)
    Y_proj <- Y_boot%*%cc_res_comp1$ycoef[,1]
    adj_r2 <- cc_res_comp1$cor^2
    adj_r2
    
  })
  
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  c(NX,NY,adj_r2,adj_r2_bootstrap_se,lrtest_p,N,sum(cd_outlier))
}

test_r2_rda_3_f <- function(X,Y,B=200,nperm=10000,parallel = NULL){
  
  as.mlm.rda <- function (x){
      X <- as.data.frame(qr.X(x$CCA$QR))
      WA <- x$CCA$wa
      list(model_1 = lm(WA ~ . , data = X),model_0 = lm(WA ~ 1 , data = X))
    }
  
  na_idx <- which(colSums(apply(cbind(X,Y),1,is.na))==0)
  
  X <- X[na_idx,,drop=FALSE]
  Y <- Y[na_idx,,drop=FALSE]
  
  if(NCOL(X)>1){
    X_rm_sing <- apply(X,2,function(x){min(table(x))})>1
    X <- X[,X_rm_sing,drop=FALSE]
  }
  
  rda_res <- rda(Y~X)
  options(warn = -1)
  na_idx <- which(colSums(apply(cooks.distance(rda_res),1,is.na))>0)
  
  while(length(na_idx)>0){
    X <- X[-na_idx,,drop=FALSE]
    Y <- Y[-na_idx,,drop=FALSE]
    rda_res <- rda(Y~X)
    na_idx <- which(colSums(apply(cooks.distance(rda_res),1,is.na))>0)
    if(sum(class(apply(X,2,unique))=="matrix")==1){
      Xcolidx <- apply(apply(X,2,unique),2,length)>1
    }else{
      Xcolidx <- sapply(apply(X,2,unique),length)>1
    }
    X <- X[,Xcolidx,drop=FALSE]
  }
  options(warn = 1)
  rda_res <- rda(Y~X)
  
  cd_outlier <- rowSums(apply(cooks.distance(rda_res),2,function(cd){cd>4 | outliers::scores(cd,type="chisq",prob=TRUE)>0.99999}))>0
  
  X_tmp <- X[!(cd_outlier ),,drop=FALSE]
  if(length(unique(X_tmp))>1){
    X <- X[!(cd_outlier ),,drop=FALSE]
    Y <- Y[!(cd_outlier ),,drop=FALSE]
  }
  
  
  if(sum(class(apply(X,2,unique))=="matrix")==1){
    Xcolidx <- apply(apply(X,2,unique),2,length)>1
  }else{
    Xcolidx <- sapply(apply(X,2,unique),length)>1
  }
  X <- X[,Xcolidx,drop=FALSE]
  
  NX = rank.condition(X)$condition
  NY = rank.condition(Y)$condition
  
  rda_res <- rda(Y~X)
  adj_r2 <- RsquareAdj(rda_res)$adj.r.squared
  
  rda_res_mlm <- as.mlm.rda(rda_res)
  # p <- linearHypothesis(rda_res_mlm[[1]],names(Anova(rda_res_mlm[[1]])[[3]]))
  if(NCOL(X)==1){
    p <- Anova(rda_res_mlm[[1]])['X','Pr(>F)']
  }else{
   lhtest <- linearHypothesis(rda_res_mlm[[1]],names(Anova(rda_res_mlm[[1]])[[3]]))
   p <- getresult.linearHypothesis.mlm(lhtest )['Pillai',"Pr(>F)"]
  }

  anova_test <- anova.cca(rda_res,permutations = how(nperm = nperm),model="direct",parallel=parallel)
  #p <- pf(anova_test$F[1],anova_test$Df[1],anova_test$Df[2],lower.tail = F)
  perm_p <- anova_test $`Pr(>F)`[1]
  N <- NROW(X)
  
  xi_cut <- apply(X,2,function(xi){
    cut(xi,breaks = c(-Inf,unique(quantile(xi))))
  })
  xi_cut_q <- apply(xi_cut,1,paste0,collapse="_")
  
  adj_r2_bootstrap <- sapply(1:B,function(i){
    balanced_resampling <- tapply(1:NROW(X),xi_cut_q,function(xin){sample(xin,NROW(xin),replace = T)})
    balanced_resampling <- unlist(balanced_resampling)
    X_boot=X[balanced_resampling,,drop=FALSE]
    Y_boot=Y[balanced_resampling,,drop=FALSE]
    rda_res <- rda(Y_boot~X_boot)
    adj_r2 <- RsquareAdj(rda_res)$adj.r.squared
    adj_r2
  })
  
  adj_r2_bootstrap_se = sqrt(sum((adj_r2_bootstrap-mean(adj_r2_bootstrap))^2)/(B-1))
  
  c(NX,NY,adj_r2,adj_r2_bootstrap_se,p,perm_p,N)
}

Pillai <- function (eig, q, df.res) {
  test <- sum(eig/(1 + eig))
  p <- length(eig)
  s <- min(p, q)
  n <- 0.5 * (df.res - p - 1)
  m <- 0.5 * (abs(p - q) - 1)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * n + s + 1
  c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s * tmp2)
}

Wilks <- function (eig, q, df.res) {
  test <- prod(1/(1 + eig))
  p <- length(eig)
  tmp1 <- df.res - 0.5 * (p - q + 1)
  tmp2 <- (p * q - 2)/4
  tmp3 <- p^2 + q^2 - 5
  tmp3 <- if (tmp3 > 0) 
    sqrt(((p * q)^2 - 4)/tmp3)
  else 1
  c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q, 
    p * q, tmp1 * tmp3 - 2 * tmp2)
}

HL <- function (eig, q, df.res) {
  test <- sum(eig)
  p <- length(eig)
  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (df.res - p - 1)
  s <- min(p, q)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * (s * n + 1)
  c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
}

Roy <- function (eig, q, df.res) {
  p <- length(eig)
  test <- max(eig)
  tmp1 <- max(p, q)
  tmp2 <- df.res - tmp1 + q
  c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
}

getresult.linearHypothesis.mlm <- function(x, SSP=TRUE, SSPE=SSP,
                                       digits=10, ...){
  test <- x$test
  if (!is.null(x$P) && SSP){
    P <- x$P
  }
  if (SSP){
  }
  if (SSPE){
  }
  if ((!is.null(x$singular)) && x$singular){
    warning("the error SSP matrix is singular; multivariate tests are unavailable")
    return(invisible(x))
  }
  SSPE.qr <- qr(x$SSPE)
  # the following code is adapted from summary.manova
  eigs <- Re(eigen(qr.coef(SSPE.qr, x$SSPH), symmetric = FALSE)$values)
  tests <- matrix(NA, 4, 4)
  rownames(tests) <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
  if ("Pillai" %in% test)
    tests[1, 1:4] <- Pillai(eigs, x$df, x$df.residual)
  if ("Wilks" %in% test)
    tests[2, 1:4] <- Wilks(eigs, x$df, x$df.residual)
  if ("Hotelling-Lawley" %in% test)
    tests[3, 1:4] <- HL(eigs, x$df, x$df.residual)
  if ("Roy" %in% test)
    tests[4, 1:4] <- Roy(eigs, x$df, x$df.residual)
  tests <- na.omit(tests)
  ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
  ok <- !is.na(ok) & ok
  tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4],
                                 lower.tail = FALSE))
  colnames(tests) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")
  tests <- structure(as.data.frame(tests),
                     heading = paste("\nMultivariate Test",
                                     if (nrow(tests) > 1) "s", ": ", x$title, sep=""),
                     class = c("anova", "data.frame"))
  tests
}


plot_biplot <- function (x, choices = c(1, 2), display = c("sp", "wa", "cn"), 
                         scaling = "species", type, xlim, ylim, const, correlation = FALSE, 
                         hill = FALSE, col_species,col_biplot,col_points_site_fill = rgb(0.5,0.5,0.5,0.5),
                         xlim_force,
                         ylim_force,
                         mul_rescale=1,
                         mul_rescale_biplot_label=1,...) {
  TYPES <- c("text", "points", "none")
  g <- scores(x, choices, display, scaling, const, correlation = correlation, 
              hill = hill)
  if (length(g) == 0 || all(is.na(g))) 
    stop("nothing to plot: requested scores do not exist")
  if (!is.list(g)) 
    g <- list(default = g)
  for (i in seq_along(g)) {
    if (length(dim(g[[i]])) > 1) 
      rownames(g[[i]]) <- rownames(g[[i]], do.NULL = FALSE, 
                                   prefix = substr(names(g)[i], 1, 3))
  }
  if (!is.null(g$centroids)) {
    if (is.null(g$biplot)) 
      g$biplot <- scores(x, choices, "bp", scaling)
    if (!is.na(g$centroids)[1]) {
      bipnam <- rownames(g$biplot)
      cntnam <- rownames(g$centroids)
      g$biplot <- g$biplot[!(bipnam %in% cntnam), , drop = FALSE]
      if (nrow(g$biplot) == 0) 
        g$biplot <- NULL
    }
  }
  if (missing(type)) {
    nitlimit <- 80
    nit <- max(nrow(g$spe), nrow(g$sit), nrow(g$con), nrow(g$def))
    if (nit > nitlimit) 
      type <- "points"
    else type <- "text"
  }
  else type <- match.arg(type, TYPES)
  if (length(choices) == 1) {
    if (length(g) == 1) 
      pl <- linestack(g[[1]], ...)
    else {
      hasSpec <- names(g)[1] == "species"
      ylim <- range(c(g[[1]], g[[2]]), na.rm = TRUE)
      pl <- linestack(g[[1]], ylim = ylim, side = ifelse(hasSpec, 
                                                         "left", "right"), ...)
      linestack(g[[2]], ylim = ylim, side = ifelse(hasSpec, 
                                                   "right", "left"), add = TRUE, ...)
    }
    return(invisible(pl))
  }
  if (missing(xlim)) {
    xlim <- range(g$species[, 1], g$sites[, 1], g$constraints[, 
                                                              1], g$biplot[, 1], if (length(g$centroids) > 0 && 
                                                                                     all(is.na(g$centroids))) NA else g$centroids[, 1], 
                  g$default[, 1], na.rm = TRUE)
  }
  if (!any(is.finite(xlim))) 
    stop("no finite scores to plot")
  if (missing(ylim)) {
    ylim <- range(g$species[, 2], g$sites[, 2], g$constraints[, 
                                                              2], g$biplot[, 2], if (length(g$centroids) > 0 && 
                                                                                     all(is.na(g$centroids))) NA else g$centroids[, 2], 
                  g$default[, 2], na.rm = TRUE)
  }
  if(!missing(xlim_force)){
    xlim = xlim_force
  }
  if(!missing(ylim_force)){
    ylim = ylim_force
  }
  print(xlim)
  print(ylim)
  
  plot(g[[1]], xlim = xlim, ylim = ylim, type = "n", asp = 1, 
       ...)
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
  if (!is.null(g$default) && type != "none") {
    if (type == "text") 
      text(g$default, rownames(g$default), cex = 0.5)
    else if (type == "points") 
      points(g$default, pch = 1, cex = 0.5)
  }
  if (!is.null(g$sites)) {
    if (type == "text") 
      text(g$sites, rownames(g$sites), cex = 0.5)
    else if (type == "points"){
      points(g$sites, pch=21,type="p",cex=0.7,bg=col_points_site_fill,col=rgb(0,0,0,0.5))
    }
  }
  
  if (!is.null(g$species)) {
    if (type == "text") 
      text(g$species, rownames(g$species), col = col_species, 
           cex = 0.7)
    else if (type == "points") 
      points(g$species, pch = "+",  cex = 1, col=col_species)
  }
  
  if (!is.null(g$constraints)) {
    if (type == "text") 
      text(g$constraints, rownames(g$constraints), cex = 0.7, 
           col = "darkgreen")
    else if (type == "points") 
      points(g$constraints, pch = 2, cex = 0.7, col = "darkgreen")
  }
  if (!is.null(g$biplot) && nrow(g$biplot) > 0 && type != 
      "none") {
    if (length(display) > 1) {
      mul <- ordiArrowMul(g$biplot)
    }
    else mul <- 1
    mul <- mul*mul_rescale
    attr(g$biplot, "arrow.mul") <- mul
    arrows(0, 0, mul * g$biplot[, 1], mul * g$biplot[, 2], 
           length = 0.05, col = col_biplot)
    biplabs <- ordiArrowTextXY(mul * g$biplot , rownames(g$biplot))
    text(biplabs*mul_rescale_biplot_label, rownames(g$biplot), col = col_biplot)
  }
  if (!is.null(g$regression) && nrow(g$regression > 0) && 
      type != "none") {
    rcol <- "purple4"
    if (length(display) > 1) {
      mul <- ordiArrowMul(g$regression)
    }
    else mul <- 1
    attr(g$regression, "arrow.mul") <- mul
    arrows(0, 0, mul * g$regression[, 1], mul * g$regression[, 
                                                             2], length = 0.05, col = rcol)
    biplabs <- ordiArrowTextXY(mul * g$regression, rownames(g$regression))
    text(biplabs*mul_rescale_biplot_label, rownames(g$regression), col = rcol)
  }
  
  if (!is.null(g$centroids) && all(!is.na(g$centroids)) && 
      type != "none") {
    if (type == "text") 
      text(g$centroids, rownames(g$centroids), col = col_biplot)
    else if (type == "points"){
      #points(g$centroids, pch = 4, col = col_biplot)
    }
      
    
  }
  
  class(g) <- "ordiplot"
  invisible(g)
}


cb_d <- function(range_res){
  write.table(range_res, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
}

## required packages

library(car)
library(quantreg)
library(rcompanion)
library(MASS)
library(sure)
library(nnet)
library(CCP)
library(pls)
library(corpcor)
library(lmtest)
library(outliers)
library(parallel)
library(magrittr)
library(SIBER)
library(vegan)
library(candisc)

pkgs <- c("car","quantreg","rcompanion","MASS","sure","nnet","CCP","pls","corpcor","lmtest","outliers","parallel","vegan",
          "DescTools","multcompView","EMT","coin")
require(pkgs)
# install.packages("https://cran.r-project.org/src/contrib/Archive/rcompanion/rcompanion_2.3.26.tar.gz",repo=NULL,type="src")

## load old / create dataset

load("../RLS_META_MG/to_stepwise_conditional.RData")

## test MG+GO+DG3 test as exposure

tr <- cbind("MG","GO","DG3")
x <- as.character(feature_order[1,])
batch <- as.factor(d_kora_analysis$batch)

res_r2_merged_outcome <- apply(feature_order,1,function(x){
  print(x)
  x0 <- d_kora_analysis[,x[2]]
  type <- x[4]
  y0 <- d_kora_analysis[,tr]
  y0 <- apply(y0,2,function(yi){ yi-tapply(yi,batch,median)[batch] + median(yi)}) # median center
  
  if(type=="continous"){
    x0_min <- min(x0,na.rm=T)
    if(min(x0_min)<=0){
      x0 <- x0+1
    }
    res <- rbind(test_r2_conti_3_f(d_r2_tmp = data.frame(x=log2(y0),y=x0)),
                 test_r2_conti_3_f(d_r2_tmp = data.frame(x=log2(y0),y=log2(x0))),
                 test_r2_conti_3_f(d_r2_tmp = data.frame(x=log2(y0),y=irnt_f(x0))),
                 test_r2_conti_3_f(d_r2_tmp = data.frame(x=irnt_df_f(y0),y=x0)),
                 test_r2_conti_3_f(d_r2_tmp = data.frame(x=irnt_df_f(y0),y=log2(x0))),
                 test_r2_conti_3_f(d_r2_tmp = data.frame(x=irnt_df_f(y0),y=irnt_f(x0)))
                 )
    rownames(res) <- c("log_raw","log_log","log_irn","irn_raw","irn_log","irn_irn")
    colnames(res) <- c("shapiro_test_res","cook_weisberg_res","r2","r2_bootstrap_se","lrtest_p","Nt","Ncd")
    r2_model <- names(which.max(res[c(1,2,4,5),3]))
    r2_model_r2 <- res[r2_model,3]
    r2_model_r2_se <- res[r2_model,4]
    if(sum(res[r2_model,1:2]>0.01)==2){
      p_model <- r2_model; p <- res[r2_model,5]; N <- res[r2_model,6]-res[r2_model,7]
    }else if(sum(res["irn_raw",1:2]>0.01)==2){
      p_model <- "irn_raw"; p <- res["irn_raw",5]; N <- res["irn_raw",6]-res["irn_raw",7]
    }else if(sum(res["irn_log",1:2]>0.01)==2){
      p_model <- "irn_log"; p <- res["irn_log",5]; N <- res["irn_log",6]-res["irn_log",7]
    }else if(sum(res["irn_irn",1:2]>0.01)==2){
      p_model <- "irn_irn"; p <- res["irn_irn",5]; N <- res["irn_irn",6]-res["irn_irn",7]
    }else{
      options(warn = -1)
      d_r2_tmp = data.frame(x=irnt_df_f(y0),y=irnt_f(x0))
      d_r2_tmp = d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
      res_rq <- lm(y~.,data = d_r2_tmp )
      p_model <- "irnt_spearman";
      p <- cor.test(d_r2_tmp$y,predict(res_rq))$p.value
      N <- NROW(d_r2_tmp)
      options(warn = 1)
    }
    c(r2_model,r2_model_r2,r2_model_r2_se,p_model,p,N)
  }else if(type=="ordinal"){
    x0 <- as.factor(x0)
    res <- rbind(test_r2_discordinal_outcome_3_f(d_r2_tmp = data.frame(y=log2(y0),x=x0)),
                 test_r2_discordinal_outcome_3_f(d_r2_tmp = data.frame(y=irnt_df_f(y0),x=x0)),
                 test_r2_disccateg_outcome_3_f(d_r2_tmp = data.frame(y=log2(y0),x=x0)),
                 test_r2_disccateg_outcome_3_f(d_r2_tmp = data.frame(y=irnt_df_f(y0),x=x0))
    )
    rownames(res) <- c("log_ord","irnt_ord","log_cat","irn_cat")
    colnames(res) <- c("shapiro_test_res","cook_weisberg_res","r2","r2_bootstrap_se","lrtest_p","Nt","Ncd")
    r2_model <- names(which.max(res[res[,1]==1,3]))
    r2_model_r2 <- res[r2_model,3]
    r2_model_r2_se <- res[r2_model,4]
    if(res[r2_model,1]==1){
      p_model <- r2_model; p <- res[r2_model,5]; N <- res[r2_model,6]-res[r2_model,7]
      c(r2_model,r2_model_r2,r2_model_r2_se,p_model,p,N)
    }else{
      d_r2_tmp = data.frame(y=log2(y0),x=x0)
      d_r2_tmp = d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
      X = as.matrix(d_r2_tmp[,1:3])
      Y = as.matrix(model.matrix(~d_r2_tmp$x)[,-1])
      NX = rank.condition(X)$condition
      NY = rank.condition(Y)$condition
      
      cc_res <- cppls(Y~X,1,data=list(X=X,Y=Y))
      X_proj <- cc_res$scores[,1]
      cc_res_comp1 <- cancor(X_proj,Y)
      Y_proj <- Y%*%cc_res_comp1$ycoef[,1]
      fit_tmp_1 <- lm(Y_proj~X_proj)
      cd_outlier <- cooks.distance(fit_tmp_1)>2
      
      X <- X[!(cd_outlier ),]
      Y <- Y[!(cd_outlier ),]
      
      NX = rank.condition(X)$condition
      NY = rank.condition(Y)$condition
      
      cc_res <- cppls(Y~X,1,data=list(X=X,Y=Y))
      X_proj <- cc_res$scores[,1]
      # kruskal.test(X_proj,as.factor(rowSums(Y)))
      cc_res_comp1 <- cancor(X_proj,Y)
      Y_proj <- Y%*%cc_res_comp1$ycoef[,1]
      
      p <- cor.test(X_proj,Y_proj,method="spearman",exact = FALSE)$p.value
      p_model <- "cpls_spearman"
      N <- NROW(X)
      c(r2_model,r2_model_r2,r2_model_r2_se,paste0(p_model,"_nx:",round(NX,2),"_ny:",round(NY,2)),p,N)
    }
  }else if(type=="catological"){
    unique_rm = is.na(match(x0,names(which(table(x0)>1))))
    x0 <- as.factor((x0[!unique_rm]))
    y0 <- y0[!unique_rm,]
    res <- rbind(test_r2_disccateg_outcome_3_f(d_r2_tmp = data.frame(y=log2(y0),x=x0)),
                 test_r2_disccateg_outcome_3_f(d_r2_tmp = data.frame(y=irnt_df_f(y0),x=x0))
    )
    rownames(res) <- c("log_cat","irn_cat")
    colnames(res) <- c("shapiro_test_res","cook_weisberg_res","r2","r2_bootstrap_se","lrtest_p","Nt","Ncd")
    r2_model <- names(which.max(res[,3]))
    
    r2_model_r2 <- res[r2_model,3]
    r2_model_r2_se <- res[r2_model,4]
    if(res[r2_model,1]==1){
      p_model <- r2_model; p <- res[r2_model,5]; N <- res[r2_model,6]-res[r2_model,7]
      c(r2_model,r2_model_r2,r2_model_r2_se,p_model,p,N)
    }else{
      d_r2_tmp = data.frame(y=log2(y0),x=x0)
      d_r2_tmp = d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
      X = as.matrix(d_r2_tmp[,1:3])
      Y = as.matrix(model.matrix(~d_r2_tmp$x)[,-1])
      NX = rank.condition(X)$condition
      NY = rank.condition(Y)$condition
      
      cc_res <- cppls(Y~X,1,data=list(X=X,Y=Y))
      X_proj <- cc_res$scores[,1]
      cc_res_comp1 <- cancor(X_proj,Y)
      Y_proj <- Y%*%cc_res_comp1$ycoef[,1]
      fit_tmp_1 <- lm(Y_proj~X_proj)
      cd_outlier <- cooks.distance(fit_tmp_1)>2
      
      X <- X[!(cd_outlier ),]
      Y <- Y[!(cd_outlier ),]
      
      NX = rank.condition(X)$condition
      NY = rank.condition(Y)$condition
      
      cc_res <- cppls(Y~X,1,data=list(X=X,Y=Y))
      X_proj <- cc_res$scores[,1]
      # kruskal.test(X_proj,as.factor(rowSums(Y)))
      cc_res_comp1 <- cancor(X_proj,Y)
      Y_proj <- Y%*%cc_res_comp1$ycoef[,1]
      
      p <- cor.test(X_proj,Y_proj,method="spearman",exact = FALSE)$p.value
      p_model <- "cpls_spearman"
      N <- NROW(X)
      c(r2_model,r2_model_r2,r2_model_r2_se,paste0(p_model,"_nx:",round(NX,2),"_ny:",round(NY,2)),p,N)
    }
  }else if(type=="binary"){
    x0 <- as.factor(x0)
    res <- rbind(test_r2_discBIN_outcome_3_f(d_r2_tmp = data.frame(y=log2(y0),x=x0)),
                 test_r2_discBIN_outcome_3_f(d_r2_tmp = data.frame(y=irnt_df_f(y0),x=x0))
    )
    rownames(res) <- c("log_cat","irn_cat")
    colnames(res) <- c("shapiro_test_res","cook_weisberg_res","r2","r2_bootstrap_se","lrtest_p","Nt","Ncd")
    r2_model <- names(which.max(res[res[,1]==1,3]))
    r2_model_r2 <- res[r2_model,3]
    r2_model_r2_se <- res[r2_model,4]
    if(res[r2_model,1]==1){
      p_model <- r2_model; p <- res[r2_model,5]; N <- res[r2_model,6]-res[r2_model,7]
      c(r2_model,r2_model_r2,r2_model_r2_se,p_model,p,N)
    }else{
      d_r2_tmp = data.frame(y=log2(y0),x=x0)
      d_r2_tmp = d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
      X = as.matrix(d_r2_tmp[,1:3])
      Y = as.matrix(model.matrix(~d_r2_tmp$x)[,-1])
      NX = rank.condition(X)$condition
      NY = 1
      
      cc_res <- cppls(Y~X,1,data=list(X=X,Y=Y))
      X_proj <- cc_res$scores[,1]
      cc_res_comp1 <- cancor(X_proj,Y)
      Y_proj <- Y%*%cc_res_comp1$ycoef[,1]
      fit_tmp_1 <- lm(Y_proj~X_proj)
      cd_outlier <- cooks.distance(fit_tmp_1)>2
      
      X <- X[!(cd_outlier ),]
      Y <- Y[!(cd_outlier ),]
      
      NX = rank.condition(X)$condition
      NY = 1
      
      cc_res <- cppls(Y~X,1,data=list(X=X,Y=Y))
      X_proj <- cc_res$scores[,1]
      
      
      
      p <- kruskal.test(X_proj,as.factor((Y)))$p.value
      p_model <- "cpls_kruskal"
      N <- NROW(X)
      c(r2_model,r2_model_r2,r2_model_r2_se,paste0(p_model,"_nx:",round(NX,2),"_ny:",round(NY,2)),p,N)
    }
  }
})


## test MG+GO+DG3 test as outcome

tr <- cbind("MG","GO","DG3")
x <- as.character(feature_order[1,])
batch <- as.factor(d_kora_analysis$batch)
res_r2_merged <- matrix(NA,85,7)
ii <- 0
j <- ii

for(ii in (j+1):NROW(feature_order)){
  
  print(ii)
  a0 <- ii
  x <- as.character(feature_order[ii,])
  x_list_tmp <- d_kora_analysis[,x[2]]
  
  type <- x[4]
  y0 <- d_kora_analysis[,tr]
  y0 <- apply(y0,2,function(yi){ yi-tapply(yi,batch,median)[batch] + median(yi)}) # median center
  
  d_r2_tmp_0 <- data.frame(y=log2(y0),as.data.frame(x_list_tmp))
  idxNA <- colSums(apply(d_r2_tmp_0,1,is.na))==0
  
  Y <- as.matrix(log2(y0)[idxNA,])
  if(type=="ordinal" || type=="catological"){
    x_list_tmp <- as.factor(x_list_tmp[which(idxNA)])
    X_tmp <- model.matrix(~x_list_tmp)[,-1]
  }else{
    X_tmp <- as.data.frame(x_list_tmp)[idxNA,]
  }
  X <- as.matrix(X_tmp)
  
  if(type=="continous"){
    X_min <- min(X,na.rm=T)
    if(min(X_min)<=0){
      X <- X + 1
    }
    res <- rbind(test_r2_rda_3_f(X=X,Y=Y,nperm = 1000000,parallel = 18),
                 test_r2_rda_3_f(X=log2(X),Y=Y,nperm = 1000000,parallel = 18),
                 test_r2_rda_3_f(X=irnt_df_f(X),Y=Y,nperm = 1000000,parallel = 18)
                 )
    
    rownames(res) <- c("log_raw","log_log","log_irn")
    colnames(res) <- c("shapiro_test_res","cook_weisberg_res","r2","r2_bootstrap_se","p","p.perm","N")
    r2_model <- names(which.max(res[,3]))
    r2_model_r2 <- res[r2_model,3]
    r2_model_r2_se <- res[r2_model,4]
    if(sum(res[r2_model,1:2]<100)==2){
      p_model <- r2_model; NX <- res[r2_model,1];  NY <- res[r2_model,2]; 
      p <- res[r2_model,5] ; perm_p <- res[r2_model,6]; N <- res[r2_model,7];
    }else{
      p_model <- "ill_matrix"; NX <- res[r2_model,1];  NY <- res[r2_model,2]; 
      p <- res[r2_model,5] ; perm_p <- res[r2_model,6]; N <- res[r2_model,7];
    }
    
    res_r2_merged_tmp <- c(r2_model,r2_model_r2,r2_model_r2_se,paste0(p_model,"_nx:",round(NX,2),"_ny:",round(NY,2)),p,perm_p,N)
  }else{
    res <- rbind(test_r2_rda_3_f(X=X,Y=Y,nperm = 1000000,parallel = 18))
    rownames(res) <- c("log_raw")
    colnames(res) <- c("shapiro_test_res","cook_weisberg_res","r2","r2_bootstrap_se","p","p.perm","N")
    r2_model <- "log_raw"
    r2_model_r2 <- res[r2_model,3]
    r2_model_r2_se <- res[r2_model,4]
    if(sum(res[r2_model,1:2]<100)==2){
      p_model <- r2_model; p <- res[r2_model,5]; perm_p <- res[r2_model,6]; N <- res[r2_model,7];
    }else{
      p_model <- "ill_matrix"; NX <- res[r2_model,1];  NY <- res[r2_model,2]; 
      p <- res[r2_model,5] ; perm_p <- res[r2_model,6]; N <- res[r2_model,7];
    }
    
    res_r2_merged_tmp <- c(r2_model,r2_model_r2,r2_model_r2_se,paste0(p_model,"_nx:",round(NX,2),"_ny:",round(NY,2)),p,perm_p,N)
  }
  res_r2_merged[ii,] <- res_r2_merged_tmp
  print(res_r2_merged_tmp)
}

cb_d(res_r2_merged)

### conditional 

res_r2_merged

get_x_rdc_f <- function(a0,feature_order=feature_order,res_r2_merged=res_r2_merged){
  fe <- feature_order[a0,2]
  
  x0 <- d_kora_analysis[,fe]
  
  type <- feature_order[a0,4]
  trans <- res_r2_merged[1,a0]
  
  
  if(grepl("_log",trans)){
    if(min(x0,na.rm=T)<=0){
      x0 <- x0 + 1
    }
    x0 <- log2(x0)
  }else if(grepl("_irn",trans)){
    x0 <- irnt_f(x0)
  }else if(type!="continous"){
    x0 <- as.factor(x0)
    x0 <- model.matrix(~x0)[,-1,drop=FALSE]
  }
  x0 <- as.matrix(x0)
  colnames(x0) <- paste0(fe,".",1:NCOL(x0))
  x0
}

#res_r2_merged <- t(res_r2_merged)

a0s <- which.min(res_r2_merged[5,])

x_list <- lapply(a0s,function(a0){
  get_x_rdc_f(a0=a0,feature_order=feature_order,res_r2_merged=res_r2_merged)
})

if(NROW(res_r2_merged)==85){
  res_r2_merged <- t(res_r2_merged)
}

tr <- c("MG","GO","DG3")
batch <- d_kora_analysis[,"batch"]
y0 <- d_kora_analysis[,tr]
y0 <- apply(y0,2,function(yi){ yi-tapply(yi,batch,median)[batch] + median(yi)}) # median center

a0s_totest <- setdiff(1:NCOL(res_r2_merged),a0s)

res_r2_merged_semipart <- matrix(NA,5,85)
#for(ii in 2:85){
for(ii in c(64:66,69)){
    print(ii)
    a0 <- ii
    x_list_tmp <- get_x_rdc_f(a0=a0,feature_order=feature_order,res_r2_merged=res_r2_merged)
    x_list_tmp_level <- f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_merged)
    
    if(length(rownames(x_list_tmp))==0){
      rownames(x_list_tmp) <- 1:NROW(x_list_tmp)
    }
    d_r2_tmp_0 <- data.frame(as.data.frame(x_list))
    d_r2_tmp_0 <- d_r2_tmp_0[as.numeric(rownames(x_list_tmp)),,drop=FALSE]
    
    if(identical(batch,x_list_tmp)){
      d_r2_tmp_1 <- d_r2_tmp_0
    }else{
      d_r2_tmp_1 <- data.frame(as.data.frame(d_r2_tmp_0),as.data.frame(x_list_tmp))
    }
    
    d_r2_tmp_0$y <- log2(y0)[as.numeric(rownames(x_list_tmp)),,drop=FALSE]
    d_r2_tmp_1$y <- log2(y0)[as.numeric(rownames(x_list_tmp)),,drop=FALSE]
    
    rownames(d_r2_tmp_0$y) <- rownames(x_list_tmp)
    rownames(d_r2_tmp_1$y) <- rownames(x_list_tmp)
    
    idxna <- complete.cases(d_r2_tmp_1)
    d_r2_tmp_0 <- d_r2_tmp_0[idxna,,drop=FALSE]
    d_r2_tmp_1 <- d_r2_tmp_1[idxna,,drop=FALSE]
    x_list_tmp_level <- x_list_tmp_level[idxna]
    
    y <- d_r2_tmp_0$y
    fit_tmp_0 <- rda(y~.,data=d_r2_tmp_0)
    fit_tmp_1 <- rda(y~.,data=d_r2_tmp_1)
    
    cd_0 <- cooks.distance(fit_tmp_0)
    cd_1 <- cooks.distance(fit_tmp_1)
    
    cd_0_in <- apply(as.matrix(cd_0),1,max)<4
    cd_1_in <- apply(as.matrix(cd_1),1,max)<4
    cd_1_sc <- apply(as.matrix(cd_1),1,max)
    
    ot_0_in <- apply(apply(as.matrix(cd_0),2,outliers::scores,type="chisq",prob=TRUE)<=0.9999,1,min)==1
    ot_1_in <- apply(apply(as.matrix(cd_1),2,outliers::scores,type="chisq",prob=TRUE)<=0.9999,1,min)==1
    
    if(sum(is.na(ot_0_in))>length(ot_0_in)/2){
      ot_0_in[is.na(ot_0_in)] <- TRUE
    }
    if(sum(is.na(ot_1_in))>length(ot_1_in)/2){
      ot_1_in[is.na(ot_1_in)] <- TRUE
    }
    ot_1_sc <- apply(apply(as.matrix(cd_1),2,outliers::scores,type="chisq",prob=TRUE),1,max)
    
    cd_in <- which((cd_0_in+cd_1_in+ot_0_in+ot_1_in)==4)
    
    #idx_disc_n1 <- apply(d_r2_tmp_1[cd_in,],2,function(di){if(length(table(di))<=2){min(table(di))<2 || length(table(di))==1}else{FALSE}})
    if(class(x_list_tmp_level)=="factor"){
      idx_disc_n1 <- sum(table(x_list_tmp_level[cd_in])>1)<2
      levels_missed <- as.numeric(names(which(table(x_list_tmp_level[cd_in])<2)))
      idx_disc_missed <- 1-!is.na(match(x_list_tmp_level,levels_missed)>0)
      
    }else{
      idx_disc_n1 <- FALSE
    }
    
    if(idx_disc_n1){
      
      idx_disc_n1_add <- do.call("order",data.frame(idx_disc_missed,ot_1_sc,cd_1_sc))[1:2]
      cd_in <- sort(unique(c(cd_in,unlist(idx_disc_n1_add)))) # patch remove no enought data; passed return NA NA p=1
      
    }
    
    if(length(cd_in)>0){
      d_r2_tmp_0 <- d_r2_tmp_0[cd_in,]
      d_r2_tmp_1 <- d_r2_tmp_1[cd_in,]
      x_list_tmp <- x_list_tmp[cd_in]
      y <- y[cd_in,]
    }
    fit_tmp_0 <- rda(y~.,data=d_r2_tmp_0)
    fit_tmp_1 <- rda(y~.,data=d_r2_tmp_1)
    
    r2_0 <- RsquareAdj(fit_tmp_0)$r.squared
    r2_1 <- RsquareAdj(fit_tmp_1)$r.squared
    
    parallel = 18
    model = "direct"
    nperm=1000000
    
    perm_p <- vegan::anova.cca(fit_tmp_0,fit_tmp_1,permutations = how(nperm=nperm),model=model,parallel=parallel)
    
    m <- perm_p$Df[2]
    n <- NROW(d_r2_tmp_1)
    
    R2 <- r2_1-r2_0
    r2_semipart = 1-(1-R2)*(n-1)/(n-m-1)
    
    perm_p <- perm_p $`Pr(>F)`[2]
    
    rda_res_mlm_1 <- as.mlm.rda(fit_tmp_1)[[1]]
    
    var_nms <- names(rda_res_mlm_1$model)
    var_fll <- gsub("\\.[0-9]+$","",var_nms)
    var_all <- unique(var_fll)
    var_1 <- var_all[2:(length(var_all))]
    var_0 <- var_all[2:(length(var_all)-1)]
    
    formula_1 <- paste0("WA ~ ",paste(var_nms[which(match(var_fll,var_1)>0)],collapse = "+"))
    formula_0 <- paste0("WA ~ ",paste(var_nms[which(match(var_fll,var_0)>0)],collapse = "+"))
    
    rda_res_mlm_lm_1 <- lm(formula_1,rda_res_mlm_1$model)
    rda_res_mlm_lm_0 <- lm(formula_0,rda_res_mlm_1$model)
    
    p <- stats::anova(rda_res_mlm_lm_1,rda_res_mlm_lm_0,test="Pillai")
    lrtest_p <- p$`Pr(>F)`[2]
    
    B=200
    
    r2_bootstrap <- sapply(1:B,function(i){
      dattype <- apply(d_r2_tmp_1,2,function(di){length(unique(di))==2})
      if(length(which(dattype)>0)){
        disc_strings <- apply(d_r2_tmp_1[,dattype,drop=FALSE],1,paste0,collapse="__")
        balanced_resampling <- 1:NROW(disc_strings)
        balanced_resampling <- tapply(balanced_resampling,disc_strings,function(iii){sample(as.character(iii),NROW(iii),replace = TRUE)})
        balanced_resampling <- as.numeric(unlist(balanced_resampling))
      }else{
        balanced_resampling <- sample(1:NROW(d_r2_tmp_1),NROW(d_r2_tmp_1),replace = TRUE)
      }
      
      fit_tmp_0_boot <- rda(y~.,data=d_r2_tmp_0[balanced_resampling,])
      fit_tmp_1_boot <- rda(y~.,data=d_r2_tmp_1[balanced_resampling,])
      r2_0_boot <- RsquareAdj(fit_tmp_0_boot)$r.squared
      r2_1_boot <- RsquareAdj(fit_tmp_1_boot)$r.squared
      r2_1_boot - r2_0_boot
    })
    r2_bootstrap <- 1-(1-r2_bootstrap)*(n-1)/(n-m-1)
    r2_semipart_bootstrap_se = sqrt(sum((r2_bootstrap-mean(r2_bootstrap))^2)/(B-1))
    
    res_cd <- c(r2_semipart,r2_semipart_bootstrap_se,lrtest_p,perm_p,n)
    
    print(res_cd)
    res_r2_merged_semipart[,ii] <- res_cd
}

### stepwise forward



cutoff <- 0.05/85
min_p_less_bon5 <- TRUE
a0s <- which.min(as.numeric(res_r2_merged[5,]))

tr <- c("MG","GO","DG3")
batch <- d_kora_analysis[,"batch"]
y0 <- d_kora_analysis[,tr]
y0 <- apply(y0,2,function(yi){ yi-tapply(yi,batch,median)[batch] + median(yi)}) # median center

res_r2_merged_semipart_FS_res <- matrix(NA,85,6)


while(min_p_less_bon5){
  x_list <- lapply(a0s,function(a0){
    get_x_rdc_f(a0=a0,feature_order=feature_order,res_r2_merged=res_r2_merged)
  })
  
  a0s_totest <- setdiff(1:NCOL(res_r2_go),a0s)
  
  res_r2_merged_semipart_FS_tmp <- mclapply(a0s_totest,function(ii){
    print(ii)
    a0 <- ii
    x_list_tmp <- get_x_rdc_f(a0=a0,feature_order=feature_order,res_r2_merged=res_r2_merged)
    x_list_tmp_level <- f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_merged)
    
    if(length(rownames(x_list_tmp))==0){
      rownames(x_list_tmp) <- 1:NROW(x_list_tmp)
    }
    d_r2_tmp_0 <- data.frame(as.data.frame(x_list))
    d_r2_tmp_0 <- d_r2_tmp_0[as.numeric(rownames(x_list_tmp)),,drop=FALSE]
    
    if(identical(batch,x_list_tmp)){
      d_r2_tmp_1 <- d_r2_tmp_0
    }else{
      d_r2_tmp_1 <- data.frame(as.data.frame(d_r2_tmp_0),as.data.frame(x_list_tmp))
    }
    
    d_r2_tmp_0$y <- log2(y0)[as.numeric(rownames(x_list_tmp)),,drop=FALSE]
    d_r2_tmp_1$y <- log2(y0)[as.numeric(rownames(x_list_tmp)),,drop=FALSE]
    
    rownames(d_r2_tmp_0$y) <- rownames(x_list_tmp)
    rownames(d_r2_tmp_1$y) <- rownames(x_list_tmp)
    
    idxna <- complete.cases(d_r2_tmp_1)
    d_r2_tmp_0 <- d_r2_tmp_0[idxna,,drop=FALSE]
    d_r2_tmp_1 <- d_r2_tmp_1[idxna,,drop=FALSE]
    x_list_tmp_level <- x_list_tmp_level[idxna]
    
    y <- d_r2_tmp_0$y
    fit_tmp_0 <- rda(y~.,data=d_r2_tmp_0)
    fit_tmp_1 <- rda(y~.,data=d_r2_tmp_1)
    
    cd_0 <- cooks.distance(fit_tmp_0)
    cd_1 <- cooks.distance(fit_tmp_1)
    
    cd_0_in <- apply(as.matrix(cd_0),1,max)<4
    cd_1_in <- apply(as.matrix(cd_1),1,max)<4
    cd_1_sc <- apply(as.matrix(cd_1),1,max)
    
    ot_0_in <- apply(apply(as.matrix(cd_0),2,outliers::scores,type="chisq",prob=TRUE)<=0.9999,1,min)==1
    ot_1_in <- apply(apply(as.matrix(cd_1),2,outliers::scores,type="chisq",prob=TRUE)<=0.9999,1,min)==1
    
    if(sum(is.na(ot_0_in))>length(ot_0_in)/2){
      ot_0_in[is.na(ot_0_in)] <- TRUE
    }
    if(sum(is.na(ot_1_in))>length(ot_1_in)/2){
      ot_1_in[is.na(ot_1_in)] <- TRUE
    }
    ot_1_sc <- apply(apply(as.matrix(cd_1),2,outliers::scores,type="chisq",prob=TRUE),1,max)
    
    cd_in <- which((cd_0_in+cd_1_in+ot_0_in+ot_1_in)==4)
    
    #idx_disc_n1 <- apply(d_r2_tmp_1[cd_in,],2,function(di){if(length(table(di))<=2){min(table(di))<2 || length(table(di))==1}else{FALSE}})
    if(class(x_list_tmp_level)=="factor"){
      idx_disc_n1 <- sum(table(x_list_tmp_level[cd_in])>1)<2
      levels_missed <- as.numeric(names(which(table(x_list_tmp_level[cd_in])<2)))
      idx_disc_missed <- 1-!is.na(match(x_list_tmp_level,levels_missed)>0)
      
    }else{
      idx_disc_n1 <- FALSE
    }
    
    if(idx_disc_n1){
      
      idx_disc_n1_add <- do.call("order",data.frame(idx_disc_missed,ot_1_sc,cd_1_sc))[1:2]
      cd_in <- sort(unique(c(cd_in,unlist(idx_disc_n1_add)))) # patch remove no enought data; passed return NA NA p=1
      
    }
    
    if(length(cd_in)>0){
      d_r2_tmp_0 <- d_r2_tmp_0[cd_in,]
      d_r2_tmp_1 <- d_r2_tmp_1[cd_in,]
      x_list_tmp <- x_list_tmp[cd_in]
      y <- y[cd_in,]
    }
    fit_tmp_0 <- rda(y~.,data=d_r2_tmp_0)
    fit_tmp_1 <- rda(y~.,data=d_r2_tmp_1)
    
    r2_0 <- RsquareAdj(fit_tmp_0)$r.squared
    r2_1 <- RsquareAdj(fit_tmp_1)$r.squared
    
    parallel = 18
    model = "direct"
    nperm=1000000
    
    perm_p <- vegan::anova.cca(fit_tmp_0,fit_tmp_1,permutations = how(nperm=nperm),model=model,parallel=parallel)
    
    m <- perm_p$Df[2]
    n <- NROW(d_r2_tmp_1)
    
    R2 <- r2_1-r2_0
    r2_semipart = 1-(1-R2)*(n-1)/(n-m-1)
    
    perm_p <- perm_p $`Pr(>F)`[2]
    
    rda_res_mlm_1 <- as.mlm.rda(fit_tmp_1)[[1]]
    
    var_nms <- names(rda_res_mlm_1$model)
    var_fll <- gsub("\\.[0-9]+$","",var_nms)
    var_all <- unique(var_fll)
    var_1 <- var_all[2:(length(var_all))]
    var_0 <- var_all[2:(length(var_all)-1)]
    
    formula_1 <- paste0("WA ~ ",paste(var_nms[which(match(var_fll,var_1)>0)],collapse = "+"))
    formula_0 <- paste0("WA ~ ",paste(var_nms[which(match(var_fll,var_0)>0)],collapse = "+"))
    
    rda_res_mlm_lm_1 <- lm(formula_1,rda_res_mlm_1$model)
    rda_res_mlm_lm_0 <- lm(formula_0,rda_res_mlm_1$model)
    
    p <- stats::anova(rda_res_mlm_lm_1,rda_res_mlm_lm_0,test="Pillai")
    lrtest_p <- p$`Pr(>F)`[2]
    
    B=200
    
    r2_bootstrap <- sapply(1:B,function(i){
      dattype <- apply(d_r2_tmp_1,2,function(di){length(unique(di))==2})
      if(length(which(dattype)>0)){
        disc_strings <- apply(d_r2_tmp_1[,dattype,drop=FALSE],1,paste0,collapse="__")
        balanced_resampling <- 1:NROW(disc_strings)
        balanced_resampling <- tapply(balanced_resampling,disc_strings,function(iii){sample(as.character(iii),NROW(iii),replace = TRUE)})
        balanced_resampling <- as.numeric(unlist(balanced_resampling))
      }else{
        balanced_resampling <- sample(1:NROW(d_r2_tmp_1),NROW(d_r2_tmp_1),replace = TRUE)
      }
      
      fit_tmp_0_boot <- rda(y~.,data=d_r2_tmp_0[balanced_resampling,])
      fit_tmp_1_boot <- rda(y~.,data=d_r2_tmp_1[balanced_resampling,])
      r2_0_boot <- RsquareAdj(fit_tmp_0_boot)$r.squared
      r2_1_boot <- RsquareAdj(fit_tmp_1_boot)$r.squared
      r2_1_boot - r2_0_boot
    })
    r2_bootstrap <- 1-(1-r2_bootstrap)*(n-1)/(n-m-1)
    r2_semipart_bootstrap_se = sqrt(sum((r2_bootstrap-mean(r2_bootstrap))^2)/(B-1))
    
    res_cd <- c(r2_semipart,r2_semipart_bootstrap_se,lrtest_p,perm_p,n)
    
    print(res_cd)
    res_cd
  },mc.cores=2)
  
  res_r2_merged_semipart_FS_tmp <- sapply(res_r2_merged_semipart_FS_tmp,function(x)x)
  res_r2_merged_semipart_FS_tmp[3:4,][which(is.na(res_r2_merged_semipart_FS_tmp[3:4,]))] <- 1
  res_r2_merged_semipart_FS_res[a0s_totest,] <- cbind(t(res_r2_merged_semipart_FS_tmp),length(a0s))
  a0s <- c(a0s,a0s_totest[which.min(res_r2_merged_semipart_FS_tmp[3,])[1]])
  min_p_less_bon5 <- min(res_r2_merged_semipart_FS_tmp[3,])<0.05/85
  
}

# stepforward RDA

feature_order2 <- read.table("../RLS_META_MG/feature_order.txt",header=TRUE,stringsAsFactors = F,sep="\t")

tr <- c("MG","GO","DG3")
batch <- d_kora_analysis[,"batch"]
y0 <- d_kora_analysis[,tr]
y0 <- apply(y0,2,function(yi){ yi-tapply(yi,batch,median)[batch] + median(yi)}) # median center
y0 <- log2(y0)

get_x_fillna_f <- function(a0,feature_order=feature_order,res_r2_mg=res_r2_mg){
  fe <- feature_order[a0,2]
  
  x0 <- d_kora_analysis[,fe]
  
  type <- feature_order[a0,4]
  trans <- res_r2_mg[1,a0]
  
  
  if(grepl("_log",trans)){
    if(min(x0,na.rm=T)<=0){
      x0 <- x0 + min(x0,na.rm=T)+1
    }
    x0 <- log2(x0)
    if(sum(is.na(x0))>0){x0[is.na(x0)] <- median(x0,na.rm=TRUE)}
  }else if(grepl("_irn",trans)){
    x0 <- irnt_f(x0)
    if(sum(is.na(x0))>0){x0[is.na(x0)] <- median(x0,na.rm=TRUE)}
  }else if(grepl("_raw",trans)){
    if(sum(is.na(x0))>0){x0[is.na(x0)] <- median(x0,na.rm=TRUE)}
  }
  if(type!="continous"){
    x0 <- as.numeric(x0)+2
    min_1 <- min(x0,na.rm=TRUE)-1
    if(sum(is.na(x0))>0){x0[is.na(x0)] <- min_1}
    x0 <- as.factor(paste0("f_",x0))
  }
  x0
}

cutoff_at_FDR005 <- p.adjust(res_r2_merged_semipart_FS_res[-1,4],method="fdr")
cutoff_at_FDR005 <- max(res_r2_merged_semipart_FS_res[-1,4][which(cutoff_at_FDR005<0.05)])

# cutoff_at_FDR005 <- p.adjust(as.numeric(res_r2_merged[6,]),method="fdr")
# cutoff_at_FDR005 <- max(as.numeric(res_r2_merged[6,])[which(cutoff_at_FDR005<0.05)])
# a0s_fw_all <- c(which(as.numeric(res_r2_merged[6,])<cutoff_at_FDR005))

# in pc
# a0s_fw_all <-  c(which(as.numeric(res_r2_merged[6,])<0.05))

a0s_fw_cond <- c(1,which(res_r2_merged_semipart_FS_res[,3]<0.1))
a0s_fw_all <- c(1,which(res_r2_merged_semipart_FS_res[,4]<cutoff_at_FDR005))

x_list_fillna_cond <- lapply(a0s_fw_cond,function(a0){
  get_x_fillna_f(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_merged)
})

x_list_fillna_all <- lapply(a0s_fw_all,function(a0){
  get_x_fillna_f(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_merged)
})

names(x_list_fillna_all) <- feature_order$raw_name[a0s_fw_all]
d_kora_analysis_features <- data.frame(x_list_fillna_all)
colnames(d_kora_analysis_features) <- feature_order$raw_name[a0s_fw_all]

complete.cases(d_kora_analysis_features)

d_kora_analysis_features$y0 <- y0


rda_0 <- rda(y0 ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia, data=d_kora_analysis_features)
rda_full <- rda(y0 ~ ., data=d_kora_analysis_features)

anova(rda_full,by="margin")

cd_full_0_sc <- cooks.distance(rda_full)
cd_full_0_sc[is.na(cd_full_0_sc)] <- 0
cd_full_0 <- apply(cd_full_0_sc,2,function(di){
  otsc <- outliers::scores(di,type="mad",prob=TRUE)
  di>4 | di==0 | outliers::scores(di,type="mad",prob=TRUE)>quantile(otsc,0.95) | is.na(di)
})
cd_full_0 <- rowSums(cd_full_0)==0

y0_cd_0 <- y0[cd_full_0,]
rda_0 <- rda(y0_cd_0 ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia, data=d_kora_analysis_features[cd_full_0,])
rda_full <- rda(y0_cd_0 ~ ., data=d_kora_analysis_features[cd_full_0,])

complete.cases(d_kora_analysis_features[cd_full_0,])


step_both_P_full_2 <- ordistep(rda_0, scope = formula(rda_full), direction = "both", 
                             Pin = 0.05,Pout = 0.1,permutations = how(nperm = 999),trace=TRUE,
                             parallel = 18, model = "direct")

step_both_P_full_2_P <- anova(step_both_P_full_2,by="margin",permutations = how(nperm = 49999),parallel = 18)


### rda diagnosis

step_both_P_full_2

inertcomp(step_both_P_full_2,display = "sp",proportional = TRUE)


# VIF<3
sort(vif.cca(step_both_P_full_2))

## plot biplot

full_2_model <- "y0_cd_0 ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia + utrauchp + uthyact + ul_ggt"
full_2_model_0 <- full_2_model

model_both_P_full_0 <- rda(formula(full_2_model_0),d_kora_analysis_features[cd_full_0,])
cd_full_1_sc <- cooks.distance(model_both_P_full_0)
cd_full_1_sc[is.na(cd_full_1_sc)] <- 0
cd_full_1 <- apply(cd_full_1_sc,2,function(di){
  otsc <- outliers::scores(di,type="mad",prob=TRUE)
  di>4 | di==0 | otsc>quantile(otsc,probs=0.975) | is.na(di)
})
cd_full_1 <- rowSums(cd_full_1)==0
cd_full_1 <- !is.na(match(1:length(cd_full_0),which(cd_full_0)[cd_full_1]))
sum(cd_full_1)

y0_cd <- y0[cd_full_1,]
full_2_model_0 <- "y0_cd ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia + utrauchp + uthyact + ul_ggt"
y0_cd <- y0[cd_full_1,]
model_both_P_full_2 <- rda(formula(full_2_model_0),d_kora_analysis_features[cd_full_1,])
range(cooks.distance(model_both_P_full_2))

pl <- plot(model_both_P_full_2,scaling = 2,xlim=c(-3,3),ylim=c(-3,3))
pl_ell <- ordiellipse(model_both_P_full_2, rep(1,NROW(model_both_P_full_2$CCA$u)), label = FALSE,conf=0.9998,lty=3,scaling = 2 )
pl_scores <- vegan::scores(model_both_P_full_2,display = "site",scaling = 2)
Z_samp <- pointsToEllipsoid(pl_scores, mu = pl_ell$`1`$center, Sigma=pl_ell$`1`$cov)
pl_ell_ot_idx <- which(!ellipseInOut(Z_samp, p = 1-pl_ell$`1`$wt[1]))

sum(cd_full_1)
cd_full_1[pl_ell_ot_idx] <- FALSE
points(pl_scores[as.numeric(pl_ell_ot_idx),],cex=2,pch=19)
sum(cd_full_1)
cd_full_1[which(d_kora_analysis_features$utwhoish=="f_8")] <- FALSE

# detect a extreme value for utwhoishf_8

# loop_cd2 <- 0
# while(length(pl_ell_ot_idx)>0){
#   loop_cd2 <- loop_cd2 + 1
#   print(loop_cd2)
#   print(sum(cd_full_1))
#   y0_cd <- y0[cd_full_1,]
#   full_2_model_0 <- "y0_cd_0 ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia + utrauchp + ul_got + utwhoish + utmgicht + utdiamm + utmrelax + uv77a_2_1 + ua3_11 + utht"
#   model_both_P_full_2 <- rda(formula(full_2_model),d_kora_analysis_features[cd_full_1,])
#   pl_ell <- ordiellipse(model_both_P_full_2, rep(1,NROW(model_both_P_full_2$CCA$u)), label = FALSE,conf=0.9995,lty=3,scaling = 2 )
#   pl_scores <- vegan::scores(model_both_P_full_2,display = "site",scaling = 2)
#   Z_samp <- pointsToEllipsoid(pl_scores, mu = pl_ell$`1`$center, Sigma=pl_ell$`1`$cov)
#   pl_ell_ot_idx <- which(!ellipseInOut(Z_samp, p = 1-pl_ell$`1`$wt[1]))
#   cd_full_1[pl_ell_ot_idx] <- FALSE
# }

full_2_model_0 <- "y0_cd ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia + utrauchp + uthyact + ul_ggt"

y0_cd <- y0[cd_full_1,]
model_both_P_full_2 <- rda(formula(full_2_model_0),d_kora_analysis_features[cd_full_1,])
range(cooks.distance(model_both_P_full_2))

conf = 1-0.05/400

pl <- plot(model_both_P_full_2,scaling = 2,xlim=c(-3,3),ylim=c(-3,3))
pl_ell <- ordiellipse(model_both_P_full_2, rep(1,NROW(model_both_P_full_2$CCA$u)), label = FALSE,conf=conf,lty=3,scaling = 2 )
pl_scores <- vegan::scores(model_both_P_full_2,display = "site",scaling = 2)
Z_samp <- pointsToEllipsoid(pl_scores, mu = pl_ell$`1`$center, Sigma=pl_ell$`1`$cov)
pl_ell_ot_idx <- which(!ellipseInOut(Z_samp, p = 1-pl_ell$`1`$wt[1]))



loop_cd2 <- 0
while(max(cooks.distance(model_both_P_full_2))>4 | loop_cd2==0){
  print(loop_cd2)
  loop_cd2 <- loop_cd2+1
  cd_full_2_sc <- cooks.distance(model_both_P_full_2)
  cd_full_2_sc[is.na(cd_full_2_sc)] <- 0
  cd_full_2 <- apply(cd_full_2_sc,2,function(di){
    otsc <- outliers::scores(di,type="mad",prob=TRUE)
    di>4 | di==0 | otsc>max(quantile(otsc,probs=0.975),0.9999) | is.na(di)
  })
  cd_full_2 <- rowSums(cd_full_2)==0
  cd_full_2 <- !is.na(match(1:length(cd_full_0),which(cd_full_1)[cd_full_2]))
  sum(cd_full_2)
  
  y0_cd <- y0[cd_full_2,]
  full_2_model_0 <- "y0_cd ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia + utrauchp + uthyact + ul_ggt"
  y0_cd <- y0[cd_full_2,]
  model_both_P_full_2 <- rda(formula(full_2_model_0),d_kora_analysis_features[cd_full_2,])
  cd_full_1 <- cd_full_2
}

sum(cd_full_2)

pl <- plot(model_both_P_full_2,scaling = 2,xlim=c(-3,3),ylim=c(-3,3))
pl_ell <- ordiellipse(model_both_P_full_2, rep(1,NROW(model_both_P_full_2$CCA$u)), label = FALSE,conf=conf,lty=3,scaling = 2 )
pl_scores <- vegan::scores(model_both_P_full_2,display = "site",scaling = 2)
Z_samp <- pointsToEllipsoid(pl_scores, mu = pl_ell$`1`$center, Sigma=pl_ell$`1`$cov)
pl_ell_ot_idx <- which(!ellipseInOut(Z_samp, p = 1-pl_ell$`1`$wt[1]))
points(pl_scores[as.numeric(pl_ell_ot_idx),],cex=2,pch=19)


model_both_P_full_P <- anova(model_both_P_full_2,by="margin",permutations = how(maxperm = 99999,nperm=3999))
model_both_P_full_P

cd_full_2_max <- apply(cooks.distance(model_both_P_full_2),1,max)

fixInNamespace(plot.cca,ns = "vegan")

anova_cont_p <- model_both_P_full_P$`Pr(>F)`[1:(NROW(model_both_P_full_P)-1)] 
anova_cont_names <- rownames(model_both_P_full_P)[1:(NROW(model_both_P_full_P)-1)]
col_biplot <- as.numeric(cut(anova_cont_p,breaks = c(0,1e-4,1e-3,1e-2,0.05)))
col_biplot <- colorRampPalette(c("darkblue","lightblue"),alpha=1)(4)[col_biplot]
names(col_biplot) <- anova_cont_names

library(shape)

pl <- plot(model_both_P_full_2)

pdf("rda_biplot_raw_scale2.pdf",height=5,width=5)

scaling = 2
pl <- plot_biplot(model_both_P_full_2, scaling=scaling, xlim=c(-2,2),ylim=c(-2,2),
                  col_specie=rgb(0.5,0.5,0.5,0),
                  col_biplot = col_biplot[gsub("f_[0-9]+$","",rownames(pl$biplot))],
                  col_points_site_fill = rgb(0.5,0.5,0.5,0.5))

#pl <- plot(model_both_P_full_2, scaling=2, xlim=c(-2,2),ylim=c(-1.5,1.5),type="n")
#pl <- plot(model_both_P_full_2, scaling=2, type="n")
#
spe2.sc <- scores(model_both_P_full_2, scaling = scaling, choices=1:2, display="sp")
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')

text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")

spe2.ce <- pl$centroids
spe2.ce_names <- gsub("f_[0-9]+$","",rownames(pl$centroids))

spe2.ce_col <- col_biplot[match(spe2.ce_names,anova_cont_names)]
points(spe2.ce[,1], spe2.ce[,2],pch=4,cex=1.5,col=spe2.ce_col,lwd=1.5)
spe2.ce_p <- anova_cont_p[match(spe2.ce_names,anova_cont_names)]
text(spe2.ce[spe2.ce_p<0.001,1], spe2.ce[spe2.ce_p<0.001,2],labels = rownames(pl$centroids)[spe2.ce_p<0.001],
     col=spe2.ce_col[spe2.ce_p<0.001])
spe2.ce_w <- apply(spe2.ce,1,function(di)sum(di^2))
spe2.ce_w_show <- spe2.ce_w>0.1
text(spe2.ce[spe2.ce_w_show,1], spe2.ce[spe2.ce_w_show,2],
     labels = rownames(pl$centroids)[spe2.ce_w_show],
     col=spe2.ce_col[spe2.ce_w_show])

ordiellipse(model_both_P_full_2, rep(1,sum(cd_full_2)), label = FALSE,conf=conf,lty=3,scaling = scaling )

dev.off()

pdf("rda_biplot_raw_scale1.pdf",height=5,width=5)

scaling = 1
pl <- plot_biplot(model_both_P_full_2, scaling=scaling, 
                  xlim_force=c(-0.8,0.8),
                  ylim_force=c(-0.8,0.8),
                  col_specie=rgb(0.5,0.5,0.5,0),
                  col_biplot = col_biplot[gsub("f_[0-9]+$","",rownames(pl$biplot))],
                  col_points_site_fill = rgb(0.5,0.5,0.5,0.5),
                  mul_rescale = 1,
                  mul_rescale_biplot_label = 1)
#pl <- plot(model_both_P_full_2, scaling=2, xlim=c(-2,2),ylim=c(-1.5,1.5),type="n")
#pl <- plot(model_both_P_full_2, scaling=2, type="n")
#
spe2.sc <- scores(model_both_P_full_2, scaling = scaling, choices=1:2, display="sp")*0.25
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')

text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col=col_biplot)


spe2.ce <- pl$centroids
spe2.ce_names <- gsub("f_[0-9]+$","",rownames(pl$centroids))

spe2.ce_col <- col_biplot[match(spe2.ce_names,anova_cont_names)]
points(spe2.ce[,1], spe2.ce[,2],pch=4,cex=1.5,col=spe2.ce_col,lwd=1.5)
# spe2.ce_p <- anova_cont_p[match(spe2.ce_names,anova_cont_names)]
# text(spe2.ce[spe2.ce_p<0.001,1], spe2.ce[spe2.ce_p<0.001,2],labels = rownames(pl$centroids)[spe2.ce_p<0.001],
#      col=spe2.ce_col[spe2.ce_p<0.001])
# spe2.ce_w <- apply(spe2.ce,1,function(di)sum(di^2))
# spe2.ce_w_show <- spe2.ce_w>0.1
# text(spe2.ce[spe2.ce_w_show,1], spe2.ce[spe2.ce_w_show,2],
#      labels = rownames(pl$centroids)[spe2.ce_w_show],
#      col=spe2.ce_col[spe2.ce_w_show])

ordiellipse(model_both_P_full_2, rep(1,sum(cd_full_2)), label = FALSE,conf=conf,lty=3,scaling = scaling )

dev.off()


# bip <- scores(model_both_P_full_2, choices = 1:2, scale= 2, display = "bp")
# (mul <- ordiArrowMul(bip, fill = 0.75))
# bip.scl <- bip * mul                    # Scale the biplot scores
# labs <- rownames(bip)
# 
# (bip.lab <- ordiArrowTextXY(bip.scl, rescale = FALSE, labels = labs))
# 
# ## draw arrows and text labels
# arrows(0, 0, bip.scl[,1], bip.scl[,2], length = 0.1,color="red")
# text(bip.lab, labels = labs)
# 
# 
# arrows(0,0,pl$biplot[,1],pl$biplot[,2],length=0, lty=1, col='darkblue')

RsquareAdj(model_both_P_full_2)

library(ggord)
ggord(model_both_P_full_2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

model_both_P_full_2_raw <- model_both_P_full_2
anova_raw <- anova(model_both_P_full_2_raw,by="margin",permutations = how(maxperm = 99999,nperm=9999))
pl_2_raw <- pl

# figure in continous trasformed 


full_2_model<- "y0_cd ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia + utrauchp + uthyact + ul_ggt"
full_2_model_0 <- "y0_cd ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia + utrauchp + uthyact + ul_ggt"

d_kora_analysis_features_numeric <- data.frame(lapply(d_kora_analysis_features[,1:82],function(x)scale(as.numeric(gsub("f_","",x)))))
d_kora_analysis_features_numeric$utoc <- d_kora_analysis_features$utoc
d_kora_analysis_features_numeric$utoc <- as.character(d_kora_analysis_features$utoc)
d_kora_analysis_features_numeric$utoc[d_kora_analysis_features_numeric$utoc=="f_4"] <- "f_1"
d_kora_analysis_features_numeric$utoc <- scale(as.numeric(gsub("f_","",d_kora_analysis_features_numeric$utoc)))
model_both_P_full_2 <- rda(formula(full_2_model),d_kora_analysis_features_numeric[cd_full_2,],scale=FALSE)
model_both_P_full_0 <- rda(y0_cd ~ 1, d_kora_analysis_features_numeric[cd_full_2,],scale=FALSE)

step_both_P_full_cont <- ordistep(model_both_P_full_0, scope = formula(model_both_P_full_2), direction = "both", 
                                  Pin = 0.001,Pout = 0.1,permutations = how(nperm = 5000),trace=TRUE, model = "direct")

anova_cont <- anova(step_both_P_full_cont,by="margin",permutations = how(maxperm = 99999,nperm=50000))
cd_full_2_max <- apply(cooks.distance(model_both_P_full_2),1,max)



anova_cont <- anova(model_both_P_full_2,by="margin",permutations = how(maxperm = 99999,nperm=50000))
anova_cont_p <- anova_cont$`Pr(>F)`[-NROW(anova_cont)]

full_2_model<- "y0_cd ~ utgfr_ckd_cr + utoc + utglukfast_a + utbmi + utmeddia + utrauchp + uthyact + ul_ggt"

d_kora_analysis_features_numeric <- data.frame(lapply(d_kora_analysis_features[,1:82],function(x)scale(as.numeric(gsub("f_","",x)))))
d_kora_analysis_features_numeric$utoc <- d_kora_analysis_features$utoc
d_kora_analysis_features_numeric$utoc <- as.character(d_kora_analysis_features$utoc)
d_kora_analysis_features_numeric$utoc[d_kora_analysis_features_numeric$utoc=="f_4"] <- "f_1"
d_kora_analysis_features_numeric$utoc <- scale(as.numeric(gsub("f_","",d_kora_analysis_features_numeric$utoc)))

model_both_P_full_2 <- rda(formula(full_2_model),d_kora_analysis_features_numeric[cd_full_2,],scale=FALSE)

anova_cont <- anova(model_both_P_full_2,by="margin",permutations = how(maxperm = 99999,nperm=50000))
anova_cont_p <- anova_cont$`Pr(>F)`[-NROW(anova_cont)]
model_both_P_full_P <- anova_cont
model_both_P_full_P

anova_cont_p <- model_both_P_full_P$`Pr(>F)`[1:(NROW(model_both_P_full_P)-1)] 
anova_cont_names <- rownames(model_both_P_full_P)[1:(NROW(model_both_P_full_P)-1)]
col_biplot <- as.numeric(cut(anova_cont_p,breaks = c(0,1e-4,1e-3,1e-2,0.05)))
col_biplot <- colorRampPalette(c("darkblue","lightblue"),alpha=1)(4)[col_biplot]
names(col_biplot) <- anova_cont_names

library(shape)

pdf("rda_biplot_con_scale2.pdf",height=5,width=5)

#scores <- function(...){-scores(...)}

scaling = 2
pl <- plot_biplot(model_both_P_full_2, scaling=scaling, xlim = c(-3,3),ylim=c(-2.2,2.2),
                  col_specie=rgb(0.5,0.5,0.5,0),
                  col_biplot = col_biplot[gsub("f_[0-9]+$","",rownames(pl$biplot))],
                  col_points_site_fill = rgb(0.5,0.5,0.5,0.5),
                  mul_rescale = 0.8,
                  mul_rescale_biplot_label = 0.8
                  )
#pl <- plot(model_both_P_full_2, scaling=2, xlim=c(-2,2),ylim=c(-1.5,1.5),type="n")
#pl <- plot(model_both_P_full_2, scaling=2, type="n")
#
spe2.sc <- scores(model_both_P_full_2, scaling = scaling, choices=1:2, display="sp")
spe2.sc <- spe2.sc*2
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')

text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")

ordiellipse(model_both_P_full_2, rep(1,sum(cd_full_2)), label = FALSE,conf=conf,lty=3,scaling = scaling )

dev.off()

pdf("rda_biplot_con_scale1.pdf",height=5,width=5)

scaling = 1
pl <- plot_biplot(model_both_P_full_2, scaling=scaling, 
                  xlim_force=c(-0.8,0.8),
                  ylim_force=c(-0.8,0.8),
                  col_specie=rgb(0.5,0.5,0.5,0),
                  col_biplot = col_biplot[gsub("f_[0-9]+$","",rownames(pl$biplot))],
                  col_points_site_fill = rgb(0.5,0.5,0.5,0.5),
                  mul_rescale = 1,
                  mul_rescale_biplot_label = 1)
#pl <- plot(model_both_P_full_2, scaling=2, xlim=c(-2,2),ylim=c(-1.5,1.5),type="n")
#pl <- plot(model_both_P_full_2, scaling=2, type="n")
#
spe2.sc <- scores(model_both_P_full_2, scaling = scaling, choices=1:2, display="sp")*0.25
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')

text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")

ordiellipse(model_both_P_full_2, rep(1,sum(cd_full_2)), label = FALSE,conf=conf,lty=3,scaling = scaling )

dev.off()


RsquareAdj(model_both_P_full_2)

library(ggord)
ggord(model_both_P_full_2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

vif.cca(model_both_P_full_2)

###

# individual figure

boxplot(d_kora_analysis_features$y0[,"MG"]~d_kora_analysis_features$utoc)
anova(lm(d_kora_analysis_features$y0[,"MG"]~d_kora_analysis_features$utoc))

# association plot

feature_order2 <- read.table("../RLS_META_MG/feature_order.txt",header=TRUE,stringsAsFactors = F,sep="\t")

exclude <- "Age"
idx_exclude_in <- which(!is.na(match(feature_order2$name,exclude)))

feature_order2_in <- feature_order2[-idx_exclude_in,]

pdf("assoc_plot_mg.pdf",height=9,width=6)
mg_name_loc <- f_assoc_plot(feature_order2_in,t(res_r2_mg_outcome[,-idx_exclude_in]))
dev.off()

pdf("assoc_plot_go.pdf",height=9,width=6)
go_name_loc <- f_assoc_plot(feature_order2_in,t(res_r2_go_outcome[,-idx_exclude_in]))
dev.off()


pdf("assoc_plot_dg3.pdf",height=9,width=6)
dg3_name_loc <- f_assoc_plot(feature_order2_in,t(res_r2_dg3_outcome[,-idx_exclude_in]))
dev.off()

pdf("assoc_plot_merged.pdf",height=9,width=6)
merge_name_loc <- f_assoc_plot(feature_order2_in=feature_order2_in,
                               data_in=t(res_r2_merged_outcome[,-idx_exclude_in]),legend_p=TRUE)
dev.off()

assoc_plot_data_FS_mg <- read.table("assoc_plot_mg_data_semipartial_FS.txt",header=T,sep="\t")
pdf("assoc_plot_mg_fs.pdf",height=9,width=6)
mg_name_loc_fs <- f_assoc_plot(feature_order2_in=feature_order2_in,
                               data_in = assoc_plot_data_FS_mg[-idx_exclude_in,c(1,1:2,3,3:5)],
                               legend_p=TRUE,
                               r2adj=FALSE)
dev.off()

assoc_plot_data_FS_go <- read.table("assoc_plot_go_data_semipartial_FS.txt",header=T,sep="\t")
pdf("assoc_plot_go_fs.pdf",height=9,width=6)
go_name_loc_fs <- f_assoc_plot(feature_order2_in=feature_order2_in,
                               data_in = assoc_plot_data_FS_go[-idx_exclude_in,c(1,1:2,3,3:5)],
                               legend_p=TRUE,
                               r2adj=FALSE)
dev.off()

assoc_plot_data_FS_dg3 <- read.table("assoc_plot_dg3_data_semipartial_FS.txt",header=T,sep="\t")
pdf("assoc_plot_dg3_fs.pdf",height=9,width=6)
dg3_name_loc_fs <- f_assoc_plot(feature_order2_in=feature_order2_in,
                               data_in = assoc_plot_data_FS_dg3[-idx_exclude_in,c(1,1:2,3,3:5)],
                               legend_p=TRUE,
                               r2adj=FALSE)
dev.off()

assoc_plot_data_FS_merge <- read.table("assoc_plot_merge_data_semipartial_FS.txt",header=T,sep="\t")
pdf("assoc_plot_merged_fs.pdf",height=9,width=6)
merge_name_loc_fs <- f_assoc_plot(feature_order2_in=feature_order2_in,
                               data_in=assoc_plot_data_FS_merge[-idx_exclude_in,c(1,1:2,3,3:5)],
                               legend_p=TRUE,
                               r2adj=FALSE)
dev.off()

f_assoc_plot <- function(feature_order2_in,data_in,legend_p=FALSE,r2adj=TRUE){
  
  idx_class <- do.call("order",feature_order2_in[,c("cat1","raw_name")])
  x_loc <- 1:NROW(feature_order2_in)
  x_loc_group <- as.numeric(as.factor(feature_order2_in[idx_class,c("cat1")]))
  x_loc <- x_loc + x_loc_group*5
  x_loc_col <- c("darkblue","lightblue")[x_loc_group%%2+1]
  
  x_loc_cat_at <- tapply(x_loc,x_loc_group,median)
  x_loc_cat_label <- unique(feature_order2_in[idx_class,c("cat1")])
  x_loc_cat_name <- unique(feature_order2_in[idx_class,c("name")])
  
  x_loc_cat_label[grep("Neurode",x_loc_cat_label)] <- "Neurological disorder"
  x_loc_cat_label[grep("Cholesterol and triglycerides",x_loc_cat_label)] <- "Cholesterol and \n triglycerides"
  
  p <- as.numeric(data_in[idx_class,5])
  #p_adj005 <- max(p[p.adjust(p,method = "bonferroni")<0.05])
  p_adj005 <- 0.05/85
  
  x_loc_names_idx <- tapply(1:NROW(p),x_loc_group,function(ti){
    ti[which.min(p[ti])]
  })
  
  x_loc_name_p <- p[x_loc_names_idx]
  x_loc_name_at <- x_loc[x_loc_names_idx][x_loc_name_p<p_adj005]
  x_loc_name_label <- x_loc_cat_name[x_loc_names_idx][x_loc_name_p<p_adj005]
  
  r2 <- as.numeric(data_in[idx_class,2])
  if(r2adj){
    r2_n <- as.numeric(data_in[idx_class,6])
    adj_r2 <- 1-(1-r2)*(r2_n-1)/(r2_n-1-1)
  }else{
    adj_r2 <- r2
  }
  
  p_r2_size <- as.numeric(cut(r2, breaks = c(-1,0.05,0.1,0.2,0.3,0.6,1)))/4+0.5
  
  par(mar=c(4,10,4,10))
  plot(-log10(p),
       x_loc,
       ylim = range(x_loc)[2:1],
       pch = 21,
       cex = p_r2_size,
       bg = x_loc_col,
       axes=FALSE,xlab="-log10 p-value",ylab="")
  
  axis(2,x_loc_cat_at,x_loc_cat_label,las=1)
  axis(4,x_loc_name_at,x_loc_name_label,las=1)
  
  axis(1)
  box()
  abline(v=-log10(p_adj005),lty=3)
  if(legend_p){
    legend("topright",legend = c("(0,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.6]","(0.6,1]"),
           pch=21,
           pt.cex=1:6/4+0.5,
           bg="white",
           bty = "n")
  }
  x_loc_name_at
}

## individual trait plot

res_r2_list <- list(res_r2_mg,res_r2_go,res_r2_dg3)

fi_to_plot <- "utgfr_ckd_cr"
jitter_mod = FALSE
fi_to_plot <- "ul_hbava"
jitter_mod = TRUE
jitter_factor = 15
fi_to_plot <- "utwhrat"
jitter_mod = FALSE
fi_to_plot <- "utglukfast_a"
jitter_mod = TRUE
jitter_factor = 15
fi_to_plot <- "utsysmm"
jitter_mod = FALSE
fi_to_plot <- "utdiamm"
jitter_mod = FALSE
fi_to_plot <- "uh_eisen"
jitter_mod = FALSE
fi_to_plot <- "uh_ferri"
jitter_mod = FALSE
fi_to_plot <- "ul_gpt"
jitter_mod = FALSE
fi_to_plot <- "ul_got"
jitter_mod = TRUE
jitter_factor = 15
fi_to_plot <- "ul_ggt"
jitter_mod = TRUE
jitter_factor = 15
fi_to_plot <- "ul_hdla"
jitter_mod = TRUE
jitter_factor = 15


{
a0 <- which(feature_order$raw_name==fi_to_plot)
a0_name <- feature_order2[a0,"name"]
fi_dtype <- rbind(res_r2_mg[,a0],res_r2_go[,a0],res_r2_dg3[,a0])
fi_dtype_min <- which.min(as.numeric(fi_dtype[,5]))
di_to_plot <- f_get_x(a0 = a0,feature_order = feature_order,res_r2_mg = res_r2_list[[fi_dtype_min]])
if(fi_to_plot=="utglukfast_a"){
  di_to_plot <- irnt_f(di_to_plot)
}
fi_dtype

fi_col <- rep("grey",3)
fi_col[as.numeric(fi_dtype[,5])<0.05] <- "lightblue"
fi_col[as.numeric(fi_dtype[,5])<0.05/485] <- "darkblue"

di_mg <- log2(d_kora_analysis$MG)
di_go <- log2(d_kora_analysis$GO)
di_dg3 <- log2(d_kora_analysis$DG3)

di_major <- list(di_mg,di_go,di_dg3)[[fi_dtype_min]]

ylab = fi_dtype[fi_dtype_min,1]
ylab = gsub("log_","",ylab)
if(fi_to_plot=="utglukfast_a"){
  ylab = "irn"
}

idx_rm_na <- complete.cases(cbind(di_major,di_to_plot,di_mg,di_go,di_dg3))

di_mg <- di_mg[idx_rm_na]
di_go <- di_go[idx_rm_na]
di_dg3 <- di_dg3[idx_rm_na]
di_major <- di_major[idx_rm_na]

di_to_plot <- di_to_plot[idx_rm_na]
par(mfrow=c(1,1))
plot(di_major,di_to_plot)

idx_rm_1 <- f_model_cd_residul_outlier_conti(di_major,di_to_plot)
idx_rm_2 <- f_model_cd_residul_outlier_conti(di_to_plot,di_major)
idx_rm_3 <- outliers::scores(di_major,type = "chisq",prob = 1-0.05/483)
idx_rm_4 <- outliers::scores(di_to_plot,type = "chisq",prob = 1-0.05/483)
idx_rm_5 <- outliers::scores(di_mg,type = "chisq",prob = 1-0.05/483)
idx_rm_6 <- outliers::scores(di_go,type = "chisq",prob = 1-0.05/483)
idx_rm_7 <- outliers::scores(di_dg3,type = "chisq",prob = 1-0.05/483)

idx_rm <- idx_rm_1 | idx_rm_2 | idx_rm_3 | idx_rm_4 | idx_rm_5 | idx_rm_6 | idx_rm_7
sum(idx_rm)

li_major <- lm(di_major[!idx_rm]~di_to_plot[!idx_rm])
summary(li_major)

li_mg <- lm(di_mg[!idx_rm]~di_to_plot[!idx_rm])
summary(li_mg)
li_go <- lm(di_go[!idx_rm]~di_to_plot[!idx_rm])
summary(li_go)
li_dg3 <- lm(di_dg3[!idx_rm]~di_to_plot[!idx_rm])
summary(li_dg3)

if(jitter_mod){
  di_to_plot <- jitter(di_to_plot,factor = jitter_factor)
}

pdf(paste0("scatter_plot_",fi_to_plot,".pdf"),width=8,height = 3)

par(mfrow=c(1,3))

par(mar=c(4,4,4,0))
plot(di_mg[!idx_rm],di_to_plot[!idx_rm],pch=19,cex=.7,col=fi_col[1],ylab="",xlab="MG(log2)")
mtext(paste0(fi_to_plot,"(",ylab,")"),2,line=2.5)
par(mar=c(4,2,4,2))
plot(di_go[!idx_rm],di_to_plot[!idx_rm],pch=19,cex=.7,col=fi_col[2],ylab="",axes=F,xlab="GO(log2)")
axis(1)
box()
par(mar=c(4,0,4,4))
plot(di_dg3[!idx_rm],di_to_plot[!idx_rm],pch=19,cex=.7,col=fi_col[3],ylab="",axes=F,xlab="3DG(log2)")
axis(1)
box()
mtext(a0_name,4,line = 1)
dev.off()
}


#lines(loess(di_to_plot[!idx_rm]~di_mg[!idx_rm],span = 1.5))


f_model_cd_residul_outlier_conti <- function(di_mg,di_to_plot){
  lm_tmp <- lm(di_mg~di_to_plot)
  lm_cd <- cooks.distance(lm_tmp)
  lm_rs <- resid(lm_tmp)
  lm_cd_idx <- outliers::scores(lm_cd,type = "chisq",prob = 1-0.05/483)
  lm_rs_idx <- outliers::scores(lm_rs,type = "chisq",prob = 1-0.05/483)
  lm_cd_idx | lm_rs_idx
}



plot(di_mg,di_to_plot)

##


step_both_P <- ordistep(rda_0, scope = formula(rda_1), direction = "both", 
                        Pin = 0.05,Pout = 0.1,permutations = how(nperm = 999),trace=TRUE,
                        parallel = 18, model = "direct")


step_both_R2 <- ordiR2step(rda_0, scope = formula(rda_1), direction = "both", R2scope = TRUE,
                        Pin = 0.05,Pout = 0.1,
                        R2permutations = 1000,
                        permutations = how(nperm = 999),
                        trace=TRUE,
                        parallel = 18, model = "direct")

# group traits

feature_order2 <- read.table("../RLS_META_MG/feature_order.txt",header=TRUE,stringsAsFactors = F,sep="\t")

table(feature_order2$cat1)
g1 <- "Kidney function"

a0s <- feature_order2[feature_order2$cat1==g1,1]




#


### test cppls.fit


Rcal(X, Y, Yprim, weights) 

d_r2_tmp = data.frame(y=log2(y0),x=x0)
d_r2_tmp = d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
X = as.matrix(d_r2_tmp[,1:3])
Y = as.matrix(model.matrix(~d_r2_tmp$x)[,-1])

n <- NROW(X)
cx <- colMeans(X)
cy <- colMeans(Y)
X <- X - rep(cx, each = n)
Y <- Y - rep(cy, each = n)

cpls_res = cppls(Y~X,3,data=list(X,Y), validation = "CV")
cpls_res_rev = cppls(X~Y,3,data=list(X,Y), validation = "CV")

cpls_res_pred <- predict(cpls_res,X)

cpls_res_pred_cancor <- cancor(cpls_res_pred[,,1],Y)
Y_scores <- Y%*%cpls_res$Yloadings
X_scores <- X%*%cpls_res$projection
X_scores_2 <- X%*%cpls_res$loadings

p_perm <- p.perm(X,Y,nboot = 5000)


cancor(cpls_res$Yscores,cpls_res$scores)$cor
cancor(X,Y)$cor

Y_score <- t(t(Y%*%cpls_res$Yloadings)/((crossprod(cpls_res$Yloadings)[1])))
Y_score <- cpls_res$Yscores
Y_score <- Y%*%cpls_res_rev$projection


cancor(X%*%X_proj[,1],Y_score)$cor

cancor(Y_score,cpls_res$scores)$cor
cancor(Y,cpls_res$scores)$cor
#cancor(Y,X)$cor

#cor(cpls_res$scores[,1],cpls_res$Yscores[,1])


x_list <- lapply(a0s,function(a0){
  f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_mg)
})
a0s_totest <- setdiff(1:NCOL(res_r2_mg),a0s)

res_r2_mg_semipart_FS_tmp <- mclapply(a0s_totest,function(ii){
  print(ii)
  a0 <- ii
  x_list_tmp <- f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_mg)
  
  d_r2_tmp_0 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list))
  fit_tmp_0 <- lm(y~.,data=d_r2_tmp_0)
  if(identical(batch,x_list_tmp)){
    d_r2_tmp_1 <- d_r2_tmp_0
  }else{
    d_r2_tmp_1 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list),as.data.frame(x_list_tmp))
  }
  fit_tmp_1 <- lm(y~.,data=d_r2_tmp_1)
  cd_0 <- cooks.distance(fit_tmp_0)
  cd_1 <- cooks.distance(fit_tmp_1)
  
  cd_in <- sort(as.numeric(intersect(names(which(cd_0<=4)), names(which(cd_1<=4)))))
  if(length(cd_in)>0){
    d_r2_tmp_0 <- d_r2_tmp_0[cd_in,]
    d_r2_tmp_1 <- d_r2_tmp_1[cd_in,]
    x_list_tmp <- x_list_tmp[cd_in]
  }
  fit_tmp_0 <- lm(y~.,data=d_r2_tmp_0)
  fit_tmp_1 <- lm(y~.,data=d_r2_tmp_1)
  df <- attr(logLik(fit_tmp_1),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp_1)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  
  r2_semipart <- summary(fit_tmp_1)$r.squared-summary(fit_tmp_0)$r.squared
  B=200
  r2_bootstrap <- sapply(1:B,function(i){
    if(length(which(sapply(d_r2_tmp_1,class)=="factor"))>0){
      disc_strings <- apply(d_r2_tmp_1[,which(sapply(d_r2_tmp_1,class)=="factor"),drop=FALSE],1,paste0,collapse="__")
      balanced_resampling <- 1:NROW(disc_strings)
      balanced_resampling <- tapply(balanced_resampling,disc_strings,function(iii){sample(as.character(iii),NROW(iii),replace = TRUE)})
      balanced_resampling <- as.numeric(unlist(balanced_resampling))
    }else{
      balanced_resampling <- sample(1:NROW(d_r2_tmp_1),NROW(d_r2_tmp_1),replace = TRUE)
    }
    
    fit_tmp_0_boot <- lm(y~.,data=d_r2_tmp_0[balanced_resampling,])
    fit_tmp_1_boot <- lm(y~.,data=d_r2_tmp_1[balanced_resampling,])
    summary(fit_tmp_1_boot)$r.squared-summary(fit_tmp_0_boot)$r.squared
  })
  r2_semipart_bootstrap_se = sqrt(sum((r2_bootstrap-mean(r2_bootstrap))^2)/(B-1))
  
  
  permutation_test <- sapply(1:5e3,function(i){
    d_r2_tmp_0_perm <- d_r2_tmp_0
    d_r2_tmp_1_perm <- d_r2_tmp_1
    set.seed(i)
    perm_idx <- sample(1:NROW(d_r2_tmp_0),NROW(d_r2_tmp_0))
    d_r2_tmp_1_perm$y <- d_r2_tmp_0_perm$y <- d_r2_tmp_0[perm_idx,1]
    fit_tmp_0_perm <- lm(y~.,data=d_r2_tmp_0_perm)
    fit_tmp_1_perm <- lm(y~.,data=d_r2_tmp_1_perm)
    summary(fit_tmp_1_perm)$r.squared-summary(fit_tmp_0_perm)$r.squared
  })
  
  perm_p <- 1-sum(r2_semipart>permutation_test)/5000
  c(r2_semipart,r2_semipart_bootstrap_se,lrtest_p,perm_p)
},mc.cores=12)

cancor(cpls_res$Yscores,cpls_res$scores)

x_list <- lapply(a0s,function(a0){
  f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_mg)
})
a0s_totest <- setdiff(1:NCOL(res_r2_mg),a0s)

res_r2_mg_semipart_FS_tmp <- mclapply(a0s_totest,function(ii){
  print(ii)
  a0 <- ii
  x_list_tmp <- f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_mg)
  
  d_r2_tmp_0 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list))
  fit_tmp_0 <- lm(y~.,data=d_r2_tmp_0)
  if(identical(batch,x_list_tmp)){
    d_r2_tmp_1 <- d_r2_tmp_0
  }else{
    d_r2_tmp_1 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list),as.data.frame(x_list_tmp))
  }
  fit_tmp_1 <- lm(y~.,data=d_r2_tmp_1)
  cd_0 <- cooks.distance(fit_tmp_0)
  cd_1 <- cooks.distance(fit_tmp_1)
  
  cd_in <- sort(as.numeric(intersect(names(which(cd_0<=4)), names(which(cd_1<=4)))))
  if(length(cd_in)>0){
    d_r2_tmp_0 <- d_r2_tmp_0[cd_in,]
    d_r2_tmp_1 <- d_r2_tmp_1[cd_in,]
    x_list_tmp <- x_list_tmp[cd_in]
  }
  fit_tmp_0 <- lm(y~.,data=d_r2_tmp_0)
  fit_tmp_1 <- lm(y~.,data=d_r2_tmp_1)
  df <- attr(logLik(fit_tmp_1),"df")
  df_0 <- attr(logLik(fit_tmp_0),"df")
  lrtest_p <- pchisq(2*abs(logLik(fit_tmp_1)-logLik(fit_tmp_0)),df=df-df_0,lower.tail = F)
  
  
  r2_semipart <- summary(fit_tmp_1)$r.squared-summary(fit_tmp_0)$r.squared
  B=200
  r2_bootstrap <- sapply(1:B,function(i){
    if(length(which(sapply(d_r2_tmp_1,class)=="factor"))>0){
      disc_strings <- apply(d_r2_tmp_1[,which(sapply(d_r2_tmp_1,class)=="factor"),drop=FALSE],1,paste0,collapse="__")
      balanced_resampling <- 1:NROW(disc_strings)
      balanced_resampling <- tapply(balanced_resampling,disc_strings,function(iii){sample(as.character(iii),NROW(iii),replace = TRUE)})
      balanced_resampling <- as.numeric(unlist(balanced_resampling))
    }else{
      balanced_resampling <- sample(1:NROW(d_r2_tmp_1),NROW(d_r2_tmp_1),replace = TRUE)
    }
    
    fit_tmp_0_boot <- lm(y~.,data=d_r2_tmp_0[balanced_resampling,])
    fit_tmp_1_boot <- lm(y~.,data=d_r2_tmp_1[balanced_resampling,])
    summary(fit_tmp_1_boot)$r.squared-summary(fit_tmp_0_boot)$r.squared
  })
  r2_semipart_bootstrap_se = sqrt(sum((r2_bootstrap-mean(r2_bootstrap))^2)/(B-1))
  
  
  permutation_test <- sapply(1:5e3,function(i){
    d_r2_tmp_0_perm <- d_r2_tmp_0
    d_r2_tmp_1_perm <- d_r2_tmp_1
    set.seed(i)
    perm_idx <- sample(1:NROW(d_r2_tmp_0),NROW(d_r2_tmp_0))
    d_r2_tmp_1_perm$y <- d_r2_tmp_0_perm$y <- d_r2_tmp_0[perm_idx,1]
    fit_tmp_0_perm <- lm(y~.,data=d_r2_tmp_0_perm)
    fit_tmp_1_perm <- lm(y~.,data=d_r2_tmp_1_perm)
    summary(fit_tmp_1_perm)$r.squared-summary(fit_tmp_0_perm)$r.squared
  })
  
  perm_p <- 1-sum(r2_semipart>permutation_test)/5000
  c(r2_semipart,r2_semipart_bootstrap_se,lrtest_p,perm_p)
},mc.cores=12)

### vegan cca rda

test_data <- d_kora_analysis[,c("MG","GO","DG3","utgfr_ckd_cr","utglukfast_a","utwhoish","uh_eisen")]
test_data <- test_data[colSums(apply(test_data,1,is.na))==0,]
test_data_metabolic <- log2(test_data[,1:3])
test_data_env <- test_data[,-(1:3)]

rda_res <- rda(test_data_metabolic~utgfr_ckd_cr+utglukfast_a+utwhoish+uh_eisen,test_data_env)

plot(hatvalues(rda_res))
pairs(cooks.distance(rda_res))
cd_rm <- rowSums(apply(cooks.distance(rda_res),2,function(cd){cd>4 | outliers::scores(cd,type="chisq",prob=TRUE)>0.9999}))>0

test_data <- d_kora_analysis[,c("MG","GO","DG3","utgfr_ckd_cr","utglukfast_a","utwhoish","uh_eisen")]
test_data <- test_data[colSums(apply(test_data,1,is.na))==0,]
test_data <- test_data[!cd_rm,]
test_data_metabolic <- log2(test_data[,1:3])
test_data_env <- test_data[,-(1:3)]
rda_res <- rda(test_data_metabolic~utgfr_ckd_cr+utglukfast_a+utwhoish+uh_eisen,test_data_env)
cca_res <- cca(test_data_metabolic~utgfr_ckd_cr+utglukfast_a+utwhoish+uh_eisen,test_data_env)

goodness(rda_res)
goodness(cca_res)

inertcomp(rda_res)
inertcomp(cca_res)

vif.cca(rda_res)
vif.cca(cca_res)

plot(rda_res, dis=c("wa","lc"), type="p")
ordispider(rda_res)

spenvcor(rda_res)
spenvcor(cca_res)


varpart_res <- varpart(test_data_metabolic,log2(test_data[,4,drop=FALSE]),test_data[,5,drop=FALSE],test_data[,6:7,drop=FALSE],
                       chisquare = FALSE)
varpart_res

showvarparts(3, bg=2:4)
plot(varpart_res, bg=2:4)

anova(rda_res)
anova(cca_res)
anova(rda_res,rda_res)
anova(cca_res,cca_res)


plot(rda_res, type = "n")
ordispider(rda_res) # segment from LC to WA scores
points(rda_res, dis="si", cex=5*hatvalues(rda_res), pch=21, bg=2) # WA scores
text(rda_res, dis="bp", col=4)

plot(rda_res)
rda_res_sum <- summary(rda_res)
R2adj <- RsquareAdj(rda_res)$adj.r.squared
R2adj

rda_res_sum$cont

plot(rda_res)
spe2.sc <- scores(rda_res, choices=1:2, display="sp")
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')
text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8)

rda_res_step <- ordistep(rda_res, direction="forward", pstep=1000, R2scop=TRUE)
anova.cca(rda_res, step=1000)

anova_cca_res <- anova.cca(rda_res,by='margin', permutations = how(nperm=9999))
anova_cca_res <- anova.cca(rda_res,by='axis', permutations = how(nperm=9999))
anova.cca(rda_res,step=10000)

pf(29.369,df1 = 4,df2 = 463,lower.tail = F)

rda_res <- metaMDS(test_data[,-6])

plot(rda_res, type="p")
spe2.sc <- scores(rda_res, choices=1:2, display="sp")
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')
text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8)

rda_res <- cca(test_data_metabolic~utgfr_ckd_cr+utglukfast_a+utwhoish+uh_eisen,test_data_env)
rda_res <- capscale(test_data_metabolic~utgfr_ckd_cr+utglukfast_a+utwhoish+uh_eisen,test_data_env)

plot(rda_res)
summary(rda_res)
R2adj <- RsquareAdj(rda_res)$adj.r.squared
R2adj

plot(rda_res)
spe2.sc <- scores(rda_res, choices=1:2, display="sp")
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')
text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8)

rda_res_step <- ordistep(rda_res, direction="forward", pstep=1000, R2scop=TRUE)
anova.cca(rda_res, step=1000)

anova.cca(rda_res,by='margin', step=1000)

rda_res <- rda(test_data_metabolic~utglukfast_a+utwhoish+uh_eisen+Condition(utgfr_ckd_cr),test_data_env)

plot(rda_res)
summary(rda_res)
R2adj <- RsquareAdj(rda_res)$adj.r.squared
R2adj

plot(rda_res)
spe2.sc <- scores(rda_res, choices=1:2, display="sp")
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')
text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8)

rda_res_step <- ordistep(rda_res, direction="forward", pstep=1000, R2scop=TRUE)
anova.cca(rda_res, step=1000)

anova.cca(rda_res,by='margin', step=1000)


# options(repos = c(
#   fawda123 = 'https://fawda123.r-universe.dev',
#   CRAN = 'https://cloud.r-project.org'))

# Install ggord
# install.packages('ggord')
# library(ggord)
# library(ggplot2)
# ggord(rda_res) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# save.image(file="results_20211208.combined_assoc.RData")

### rda condition

Z <- X+rnorm(NROW(X))*5


rda_res_cond <- rda(Y~Z+Condition(X))

rda_res_full <- rda(Y~X)
rda_res_full_1 <- rda(Y~X+Z)
rda_res_full_2 <- rda(Y~Z)


anova(rda_res_full)
anova(rda_res_cond)

adj_r2 <- RsquareAdj(rda_res_full_1)$r.squared-RsquareAdj(rda_res_full_2)$r.squared
adj_r2

rda_res_cond_mlm <- as.mlm.rda(rda_res_cond)
p <- Anova(rda_res_cond_mlm[[1]])['Z','Pr(>F)']

rda_res_mlm_1 <- as.mlm.rda(rda_res_full_1)
rda_res_mlm_2 <- as.mlm.rda(rda_res_full_2)

rda_res_mlm_lm_1 <- lm(rda_res_mlm_1[[1]]$model$WA~X+Z,rda_res_mlm_1[[1]]$model)
rda_res_mlm_lm_0 <- lm(rda_res_mlm_1[[1]]$model$WA~X,rda_res_mlm_1[[1]]$model)
Pillai_test <- anova(rda_res_mlm_lm_1,rda_res_mlm_lm_0,test="Pillai")
Pillai_test
anova(rda_res_mlm_lm_1,rda_res_mlm_lm_0,test="Spherical")
anova(rda_res_full_1,rda_res_full_2)


# test.statistic=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
Anova(rda_res_mlm_lm_1,test.statistic = c("Pillai"))
Anova(rda_res_mlm_lm_1,test.statistic = c("Wilks"))

compareCoefs(rda_res_mlm_lm_1,rda_res_mlm_lm_0)

# p <- linearHypothesis(rda_res_mlm[[1]],names(Anova(rda_res_mlm[[1]])[[3]]))
if(NCOL(X)==1){
  p <- Anova(rda_res_mlm_1[[1]])['X','Pr(>F)']
}else{
  lhtest <- linearHypothesis(rda_res_mlm[[1]],names(Anova(rda_res_mlm[[1]])[[3]]))
  p <- getresult.linearHypothesis.mlm(lhtest )['Pillai',"Pr(>F)"]
}

anova_test <- anova(rda_res_full_2,rda_res_full_1,permutations = how(nperm = nperm),model="direct",parallel=parallel)
# p <- pf(anova_test$F[2],anova_test$ResDf[1],anova_test$ResDf[2],lower.tail = F)
perm_p <- anova_test $`Pr(>F)`[1]
N <- NROW(X)

### add rda pca mds

library(vegan)

#pdf("pca_rda_biplot_con_scale2.pdf",height=5,width=5)

#scores <- function(...){-scores(...)}

y0_outliers_idx <- apply(y0,2,function(di){outliers::scores(di,type="chisq",prob=TRUE)>1-0.05/483})
y0_outliers_idx <- rowSums(y0_outliers_idx)==0

y0_outliers <- y0[y0_outliers_idx,]

scaling = 1

rda_pca_1 <- rda(y0_outliers~1,scale=TRUE)
#0.6521 0.2293 0.1186
rda_pca_1 <- rda(y0_outliers~1,scale=FALSE)
#0.6467 0.24898 0.10430

pl <- plot_biplot(rda_pca_1, scaling=scaling, xlim = c(-0.6,0.6),ylim=c(-0.6,0.6),
                  col_specie=rgb(0.5,0.5,0.5,0),
                  col_biplot = col_biplot[gsub("f_[0-9]+$","",rownames(pl$biplot))],
                  col_points_site_fill = rgb(0.5,0.5,0.5,0.5),
                  mul_rescale = 0.8,
                  mul_rescale_biplot_label = 0.8
)
#pl <- plot(model_both_P_full_2, scaling=2, xlim=c(-2,2),ylim=c(-1.5,1.5),type="n")
#pl <- plot(model_both_P_full_2, scaling=2, type="n")
#
spe2.sc <- scores(rda_pca_1, scaling = scaling, choices=1:2, display="sp")
spe2.sc <- spe2.sc*0.2
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')

text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")

ordiellipse(model_both_P_full_2, rep(1,sum(cd_full_2)), label = FALSE,conf=1-0.05/485,lty=3,scaling = scaling )

dev.off()

pl <- plot_biplot(rda_pca_1, choices = c(1,3),scaling=scaling, xlim = c(-0.6,0.6),ylim=c(-0.6,0.6),
                  col_specie=rgb(0.5,0.5,0.5,0),
                  col_biplot = col_biplot[gsub("f_[0-9]+$","",rownames(pl$biplot))],
                  col_points_site_fill = rgb(0.5,0.5,0.5,0.5),
                  mul_rescale = 0.8,
                  mul_rescale_biplot_label = 0.8
)

spe2.sc <- scores(rda_pca_1, choices = c(1,3),scaling = scaling,  display="sp")
spe2.sc <- spe2.sc*0.2
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')

text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")

#pca_comp <- prcomp(y0_outliers,retx = TRUE)

y0_dist <- dist(y0_outliers)

y0_mds <- wcmdscale(y0_dist,3,eig=TRUE,x.ret=TRUE,w=rep(1,NROW(y0_dist)))
plot(y0_mds)


y0_mds <- wcmdscale(y0_dist,2,eig=TRUE,x.ret=TRUE,w=rep(1,NROW(y0_dist)))
plot(y0_mds)

y0_mds$eig/sum(y0_mds$eig)
plot(y0_mds$points[,1],-y0_mds$points[,2])
y0_mds_points <- y0_mds$points

rda_mds_0 <- rda(y0_mds_points~y0_outliers)
plot(rda_mds_0,scaling=1)

#
y0_outliers <- scale(y0_outliers,center = TRUE,scale=FALSE)
scaling = 1
rda_pca_1 <- rda(y0_outliers~1)
choices = c(1,2)
pl <- plot_biplot(rda_pca_1, choices = choices,scaling=scaling, xlim = c(-0.6,0.6),ylim=c(-0.6,0.6),
                  col_specie=rgb(0.5,0.5,0.5,0),
                  col_biplot = col_biplot[gsub("f_[0-9]+$","",rownames(pl$biplot))],
                  col_points_site_fill = rgb(0.5,0.5,0.5,0.5),
                  mul_rescale = 0.8,
                  mul_rescale_biplot_label = 0.8
)

spe2.sc <- scores(rda_pca_1, choices = choices,scaling = scaling,  display="sp")
spe2.sc <- spe2.sc*0.2
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')
text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")

rda_pca_1$CA$v
# PC1        PC2         PC3
# MG  0.4801000 -0.3905538  0.78547546
# GO  0.6289263 -0.4709612 -0.61858485
# DG3 0.6115191  0.7909888  0.01952109

# p1 shared direction 
# p2 DG3-(MG+GO)
# p3 MG-GO

# compare mds and cca

head(scores(rda_pca_1)$sites)
head(y0_mds$points)

head(y0_mds$points)/head(scores(rda_pca_1,choices = 1:3,scaling = 1)$sites)

y0_mds$GOF

cumsum(rda_pca_1$CA$eig)/sum(rda_pca_1$CA$eig)
#0.6467124 0.8956970 1.0000000

## score and loadings

choices=1:2
scores(rda_pca_1, choices = choices,scaling = scaling,  display="sp")

scores(rda_pca_1, choices = choices,scaling = scaling,  display="sp")/3.670813
rda_pca_1$CA$v

## loadings dim 1
w_v <- (head(y0_mds$points)/head(scores(rda_pca_1)$sites))[1,]
spe2.sc <- t(rda_pca_1$CA$v[,1:2])*w_v

par(mfrow=c(1,2))
plot(y0_mds$points[,1],y0_mds$points[,2],xlim=c(-2,2),ylim=c(-2,2))
spe2.sc <- scores(rda_pca_1, choices = choices,scaling = scaling,  display="sp")
w_v <- (head(y0_mds$points)/head(scores(rda_pca_1)$sites))[1,]
spe2.sc <- t(rda_pca_1$CA$v[,1:3])*w_v
spe2.sc <- spe2.sc
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')
text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")


pl <- plot_biplot(rda_pca_1, choices = choices,scaling=scaling, xlim = c(-0.6,0.6),ylim=c(-0.6,0.6),
                  col_specie=rgb(0.5,0.5,0.5,0),
                  col_biplot = col_biplot[gsub("f_[0-9]+$","",rownames(pl$biplot))],
                  col_points_site_fill = rgb(0.5,0.5,0.5,0.5),
                  mul_rescale = 0.8,
                  mul_rescale_biplot_label = 0.8
)

### map to mds 2
# varpart 

y0_mds <- wcmdscale(y0_dist,2,eig=TRUE,x.ret=TRUE,w=rep(1,NROW(y0_dist)))
plot(y0_mds)

qr_w <- qr.solve(y0_outliers,y0_mds$points)
# qr_w = rda_pca_1$CA$v

par(mfrow=c(1,1),mar=c(4,4,4,4))
plot(y0_mds$points[,1],-y0_mds$points[,2],xlim=c(-2,2),ylim=c(-2,2),type="n")
spe2.sc <- rda_pca_1$CA$v[,1:2]
spe2.sc <- spe2.sc*2
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')
text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")
gfr_col <- (d_kora_analysis$utgfr_ckd_cr[y0_outliers_idx]<60)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","red")[gfr_col])
#glu_col <- (d_kora_analysis$utglukfast_a[y0_outliers_idx]>100)+1
#points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","blue")[glu_col])

gfr_col <- (d_kora_analysis$utgfr_ckd_cr[y0_outliers_idx]<60)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","red")[gfr_col])

gfr_col <- (d_kora_analysis$utgfr_ckd_cr[y0_outliers_idx]>120)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","blue")[gfr_col])

cors <- apply(d_kora_analysis[y0_outliers_idx,feature_order$raw_name],2,function(x){
  cor(as.numeric(x),y0_mds$points[,1],use="pairwise.complete.obs")
})

uthyact_col <- (d_kora_analysis$uthyact[y0_outliers_idx]>1)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","red")[uthyact_col])

cor(d_kora_analysis$uthyact,d_kora_analysis$utdiamm)

utsysmm_col <- (d_kora_analysis$utsysmm[y0_outliers_idx]<90)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","red")[utsysmm_col])

utdiamm_col <- (d_kora_analysis$utdiamm[y0_outliers_idx]<60)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","red")[utdiamm_col])


ul_ggt_col <- (d_kora_analysis$ul_ggt[y0_outliers_idx]>50)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","red")[ul_ggt_col])

ul_ggt_col <- (d_kora_analysis$ul_ggt[y0_outliers_idx]<50)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","red")[ul_ggt_col])

ul_ggt_col <- (d_kora_analysis$ul_ggt[y0_outliers_idx]==45)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","red")[ul_ggt_col],cex=2)


glu_col <- (d_kora_analysis$utglukfast_a[y0_outliers_idx]<80)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","blue")[glu_col])

which.min(as.numeric(y0_mds$points[,1]))
points(y0_mds$points[387,1],y0_mds$points[387,2],cex=2,col="blue")
d_kora_analysis[y0_outliers_idx,][387,]
2^(y0_outliers[387,])

d_kora_analysis[y0_outliers_idx,][387,c("ul_ggt","utgfr_ckd_cr")]

gfr_ggt_col <- (d_kora_analysis$ul_ggt[y0_outliers_idx]<50 & d_kora_analysis$utgfr_ckd_cr[y0_outliers_idx]>110)+1
points(y0_mds$points[,1],-y0_mds$points[,2],col=c("black","blue")[gfr_ggt_col])

## plot mds 

library(mdatools)

rda_pca_1$CA$eig/sum(rda_pca_1$CA$eig)
y0_mds$eig/sum(y0_mds$eig)
mdsvars = y0_mds$eig/sum(y0_mds$eig)

pdf("mds_plot.pdf",height=5,width=5)
plot(y0_mds$points[,1],-y0_mds$points[,2],xlim=c(-2,2)*1.5,ylim=c(-2,2)*1.5,type="n",
     xlab=paste0("Dim.1(",round(mdsvars[1],4),")"),
     ylab=paste0("Dim.2(",round(mdsvars[2],4),")")
     )
spe2.sc <- rda_pca_1$CA$v[,1:2]
spe2.sc <- spe2.sc*3
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='darkred')
points(spe2.sc[,1], spe2.sc[,2],pch=3,cex=0.7,col='darkred')
points(y0_mds$points[,1],-y0_mds$points[,2],pch=21,bg=rgb(0.5,0.5,0.5,0.5),col=rgb(0,0,0,0.5))
text(spe2.sc[,1], spe2.sc[,2], labels = rownames(spe2.sc),adj = c(-0.2, NA),cex=0.8,col="darkred")
t1 <- as.numeric(y0_mds$points[,1])
t2 <- as.numeric(y0_mds$points[,2])
s1 <- sd(t1)^2
s2 <- sd(t2)^2
nobj <- length(t1)
conf.lim=1-0.05/482
T2lim <- (2 * (nobj - 1)/(nobj - 2)) * qf(conf.lim, 2, (nobj - 2))
a <- sqrt(T2lim * s1)
b <- sqrt(T2lim * s2)
ellipse(a = a, b = b, col = "darkgrey", lty = 3)
abline(h=0,v=0,col="darkgrey",lty=3)
dev.off()


## boot strap 




## output 

table(dat_MG_age_sex_nonoutlier$age_group,dat_MG_age_sex_nonoutlier$RLS )

rls_full_data_idx_in_ttest <- which(!is.na(match(d_rls$ID,
      dat_MG_age_sex_nonoutlier[which(dat_MG_age_sex_nonoutlier$RLS==1 & dat_MG_age_sex_nonoutlier$age_group=="(55,65]"),]$"ID")))
d_rls[rls_full_data_idx_in_ttest,c("ID","Serum_ID","IHG_Nr","Status")]


write.table(dat_MG_age_sex_nonoutlier,"dat_MG_age_sex_nonoutlier.txt",col.names = TRUE,row.names = F,sep="\t",quote=F)







