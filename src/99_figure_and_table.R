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
    if(sum(res[r2_model,1:2]>0.05)==2){
      p_model <- r2_model; p <- res[r2_model,5]; N <- res[r2_model,6]-res[r2_model,7]
    }else if(sum(res["irn_irn",1:2]>0.05)==2){
      p_model <- "irn_irn"; p <- res["irn_irn",5]; N <- res["irn_irn",6]-res["irn_irn",7]
    }else{
      options(warn = -1)
      d_r2_tmp = data.frame(x=log2(y0),y=x0)
      d_r2_tmp = d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
      res_rq <- rq(y~.,data = d_r2_tmp )
      p_model <- "median_reg";
      p <- nagelkerke(res_rq)$Likelihood.ratio.test[1,"p.value"]
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

res_r2_merged <- sapply(1:NroW(feature_order),function(ii){
  
  print(ii)
  a0 <- ii
  x <- as.character(feature_order[ii,])
  x_list_tmp <- d_kora_analysis[,x[2]]
  
  type <- x[4]
  y0 <- d_kora_analysis[,tr]
  y0 <- apply(y0,2,function(yi){ yi-tapply(yi,batch,median)[batch] + median(yi)}) # median center
  
  d_r2_tmp_0 <- data.frame(y=log2(y0),as.data.frame(x_list))
  idxNA <- colSums(apply(d_r2_tmp_0,1,is.na))==0
  
  Y <- as.matrix(log2(y0)[idxNA,])
  if(type=="ordinal" || type=="catological"){
    x_list_tmp <- as.factor(x_list_tmp)
    X_tmp <- model.matrix(~x_list_tmp)[idxNA,-1]
  }else{
    X_tmp <- as.data.frame(x_list_tmp)[idxNA,]
  }
  X <- as.matrix(X_tmp)
  
  if(type=="continous"){
    X_min <- min(X,na.rm=T)
    if(min(X_min)<=0){
      X <- X + 1
    }
    res <- rbind(test_r2_cpls_3_f(X=X,Y=Y),
                 test_r2_cpls_3_f(X=log2(X),Y=Y),
                 test_r2_cpls_3_f(X=irnt_df_f(X),Y=Y)
                 )
    
    rownames(res) <- c("log_raw","log_log","log_irn")
    colnames(res) <- c("shapiro_test_res","cook_weisberg_res","r2","r2_bootstrap_se","lrtest_p","Nt","Ncd")
    r2_model <- names(which.max(res[,3]))
    r2_model_r2 <- res[r2_model,3]
    r2_model_r2_se <- res[r2_model,4]
    if(sum(res[r2_model,1:2]<100)==2){
      p_model <- r2_model; p <- res[r2_model,5]; N <- res[r2_model,6]
    }else{
      p_model <- "ill_matrix"; p <- "NA"; N <- NA
    }
    
    c(r2_model,r2_model_r2,r2_model_r2_se,paste0(p_model,"_nx:",round(NX,2),"_ny:",round(NY,2)),p,N)
  }else{
    res <- rbind(test_r2_cpls_3_f(X=Y,Y=X))
    
  }
  
  fit_tmp_0 <- lm(y~.,data=d_r2_tmp_0)
  if(identical(batch,x_list_tmp)){
    d_r2_tmp_1 <- d_r2_tmp_0
  }else{
    d_r2_tmp_1 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list),as.data.frame(x_list_tmp))
  }
  
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
    if(sum(res[r2_model,1:2]>0.05)==2){
      p_model <- r2_model; p <- res[r2_model,5]; N <- res[r2_model,6]-res[r2_model,7]
    }else if(sum(res["irn_irn",1:2]>0.05)==2){
      p_model <- "irn_irn"; p <- res["irn_irn",5]; N <- res["irn_irn",6]-res["irn_irn",7]
    }else{
      options(warn = -1)
      d_r2_tmp = data.frame(x=log2(y0),y=x0)
      d_r2_tmp = d_r2_tmp[colSums(apply(d_r2_tmp,1,is.na))==0,]
      res_rq <- rq(y~.,data = d_r2_tmp )
      p_model <- "median_reg";
      p <- nagelkerke(res_rq)$Likelihood.ratio.test[1,"p.value"]
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



