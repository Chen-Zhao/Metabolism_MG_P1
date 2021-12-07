
tr <- "MG"
cutoff <- 0.05/85
min_p_less_bon5 <- TRUE
a0s <- which.min(res_r2_mg[5,])
y0 <- d_kora_analysis[,tr]
batch <- as.factor(d_kora_analysis[,"batch"])
res_r2_mg_semipart_FS_res <- matrix(NA,85,5)

while(min_p_less_bon5){
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
  
  res_r2_mg_semipart_FS_tmp <- sapply(res_r2_mg_semipart_FS_tmp,function(x)x)
  res_r2_mg_semipart_FS_res[a0s_totest,] <- cbind(t(res_r2_mg_semipart_FS_tmp),length(a0s))
  a0s <- c(a0s,a0s_totest[which.min(res_r2_mg_semipart_FS_tmp[3,])[1]])
  min_p_less_bon5 <- min(res_r2_mg_semipart_FS_tmp[3,])<0.05/85
}

cb_d(res_r2_mg_semipart_FS_res)



tr <- "GO"
cutoff <- 0.05/85
min_p_less_bon5 <- TRUE
a0s <- which.min(res_r2_go[5,])
y0 <- d_kora_analysis[,tr]
batch <- as.factor(d_kora_analysis[,"batch"])
res_r2_go_semipart_FS_res <- matrix(NA,85,5)

while(min_p_less_bon5){
  x_list <- lapply(a0s,function(a0){
    f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_go)
  })
  a0s_totest <- setdiff(1:NCOL(res_r2_go),a0s)
  
  res_r2_go_semipart_FS_tmp <- mclapply(a0s_totest,function(ii){
    print(ii)
    a0 <- ii
    x_list_tmp <- f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_go)
    
    d_r2_tmp_0 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list))
    fit_tmp_0 <- lm(y~.,data=d_r2_tmp_0)
    if(identical(batch,x_list_tmp)){
      d_r2_tmp_1 <- d_r2_tmp_0
    }else{
      d_r2_tmp_1 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list),as.data.frame(x_list_tmp))
    }
    fit_tmp_1 <- tryCatch({lm(y~.,data=d_r2_tmp_1)},error = function(e) {e})
    if(grepl("with 2 or more levels",fit_tmp_1)[1]){
      c(NA,NA,1,1)
    }else{
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
      fit_tmp_1 <- tryCatch({lm(y~.,data=d_r2_tmp_1)},error = function(e) {e})
      if(grepl("with 2 or more levels",fit_tmp_1)[1]){
        c(NA,NA,1,1)
      }else{
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
      }
    }
  },mc.cores=12)
  
  res_r2_go_semipart_FS_tmp <- sapply(res_r2_go_semipart_FS_tmp,function(x)x)
  res_r2_go_semipart_FS_res[a0s_totest,] <- cbind(t(res_r2_go_semipart_FS_tmp),length(a0s))
  a0s <- c(a0s,a0s_totest[which.min(res_r2_go_semipart_FS_tmp[3,])[1]])
  min_p_less_bon5 <- min(res_r2_go_semipart_FS_tmp[3,])<0.05/85
}

cb_d(res_r2_go_semipart_FS_res)


tr <- "DG3"
cutoff <- 0.05/85
min_p_less_bon5 <- TRUE
a0s <- which.min(res_r2_dg3[5,])
y0 <- d_kora_analysis[,tr]
batch <- as.factor(d_kora_analysis[,"batch"])
res_r2_dg3_semipart_FS_res <- matrix(NA,85,5)

while(min_p_less_bon5){
  x_list <- lapply(a0s,function(a0){
    f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_dg3)
  })
  a0s_totest <- setdiff(1:NCOL(res_r2_dg3),a0s)
  
  res_r2_dg3_semipart_FS_tmp <- mclapply(a0s_totest,function(ii){
    print(ii)
    a0 <- ii
    x_list_tmp <- f_get_x(a0=a0,feature_order=feature_order,res_r2_mg=res_r2_dg3)
    
    d_r2_tmp_0 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list))
    fit_tmp_0 <- lm(y~.,data=d_r2_tmp_0)
    if(identical(batch,x_list_tmp)){
      d_r2_tmp_1 <- d_r2_tmp_0
    }else{
      d_r2_tmp_1 <- data.frame(y=log2(y0),batch=batch,as.data.frame(x_list),as.data.frame(x_list_tmp))
    }
    fit_tmp_1 <- tryCatch({lm(y~.,data=d_r2_tmp_1)},error = function(e) {e})
    if(grepl("with 2 or more levels",fit_tmp_1)){
      c(NA,NA,1,1)
    }else{
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
    }
  },mc.cores=12)
  
  res_r2_dg3_semipart_FS_tmp <- sapply(res_r2_dg3_semipart_FS_tmp,function(x)x)
  res_r2_dg3_semipart_FS_res[a0s_totest,] <- cbind(t(res_r2_dg3_semipart_FS_tmp),length(a0s))
  a0s <- c(a0s,a0s_totest[which.min(res_r2_dg3_semipart_FS_tmp[3,])[1]])
  min_p_less_bon5 <- min(res_r2_dg3_semipart_FS_tmp[3,])<0.05/85
}

cb_d(res_r2_dg3_semipart_FS_res)

