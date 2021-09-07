tab_clin_cont <- function(clinical,groups,contrast=F){
  clinvar <- c()
  mean <- c()
  pvalanova <- c()
  pvalkruskal <- c()
  pvalshap <- c()
  pvalPairwise<-c()
  
  if(contrast){
    contrasts_col <- c()
  }
  
  for (j in names(clinical)){
    if (is.numeric(clinical[,j])){
      mean1 <- sapply(levels(groups),function(l){signif(mean(clinical[which(groups==l),j],na.rm=TRUE),3)})        
      sd <- sapply(levels(groups),function(l){signif(sd(clinical[which(groups==l),j],na.rm=TRUE),3)})
      median1 <- sapply(levels(groups),function(l){signif(median(clinical[which(groups==l),j],na.rm = T),3)})
      quantile25_1 <- sapply(levels(groups), function(l){signif(quantile(clinical[which(groups==l),j],c(0.25),na.rm = T),3)})
      quantile75_1 <- sapply(levels(groups), function(l){signif(quantile(clinical[which(groups==l),j],c(0.75),na.rm = T),3)})
      NAS <- sapply(levels(groups),function(l){sum(is.na(clinical[which(groups==l),j]))})
            
      dat <- data.frame(groups)
      mod <- lm(formula=clinical[,j] ~ .,dat)
      tmp <- tryCatch(anova(mod)['groups','Pr(>F)'],error=function(e){return(NA)})
      # TukeyHSD test for pairwise comparisons in normality conditions
	  tuk <- TukeyHSD(x=aov(clinical[,j] ~ groups), 'groups', conf.levels=0.95)
	  contr<-rbind(tuk$groups[,4])
      names_pair<-row.names(tuk$groups)	
 	
      clinvar <- c(clinvar,j)
      mean <- rbind(mean,paste(mean1," +/- ",sd, " or ", median1, " (",quantile25_1," - ",quantile75_1,"); NA's:",NAS ,sep=""))
      pvalshap <- c(pvalshap,tryCatch(shapiro.test(clinical[,j])$p.value, error=function(e){return(NA)}))
      pvalanova <- c(pvalanova,signif(tmp,3))
      pvalkruskal <- c(pvalkruskal,tryCatch(kruskal.test(formula= clinical[,j] ~  groups)$p.value,error=function(e){return(NA)}))
	  pvalPairwise <- rbind(pvalPairwise, as.numeric(contr[1,]))
	  
      if(contrast){
        contrasts_col <- c(contrasts_col,NA)
      }
    }
	else{
      clinical_data <- as.factor(as.character(clinical[,j]))
      count <- c()
      pval <- data.frame()
      for (l in levels(groups)){
        if(any(is.na(clinical_data[which(groups==l)]))){
          count <- c(count,paste(" ", c(levels(clinical_data),"NA"),":",
            summary(clinical_data[which(groups==l)]),sep="",collapse=";"))
        }else{
          count <- c(count,paste(" ", levels(clinical_data),":",
             summary(clinical_data[which(groups==l)]),sep="",collapse=";"))
        }
      }
      mean <- rbind(mean,count)
      if(contrast){
        for (l in 1:(length(levels(groups))-1)){
          grouph <- relevel(groups,ref = levels(groups)[l])
          dat <- data.frame(grouph)
              
          for (k in 1:(length(levels(clinical_data))-1)){
            clin <- relevel(clinical_data,ref=levels(clinical_data)[k])
            logi <-  tryCatch(multinom(clin ~  .,data = dat,trace=FALSE),error=function(e){return(NA)})
            summ <- tryCatch(summary(logi),error=function(e){return(NA)})
            z <- tryCatch(summ$coefficients/summ$standard.errors,error=function(e){return(NA)})
            p <- tryCatch((1 - pnorm(abs(z), 0, 1)) * 2,error=function(e){return(NULL)})
              
            if(is.null(dim(p))){
              if(!is.null(p)){
                p2 <- data.frame(matrix(data = NA,ncol=length(levels(groups))+1,nrow=1))
                p2[,(l+2):(length(levels(groups))+1)] <- p[(l+1):length(p)]
                p2[,1] <- paste(levels(clinical_data)[2],'vs',levels(clinical_data)[k])
                pval <- rbind(pval,p2)
              }
            }else{
              p2 <- data.frame(matrix(data = NA,ncol=length(levels(groups))+1,nrow=(dim(p)[1])+1-k))
              p2[,(l+2):(length(levels(groups))+1)] <- p[k:dim(p)[1],(l+1):dim(p)[2]]
              p2[,1] <- paste(row.names(p)[k:dim(p)[1]],'vs',levels(clinical_data)[k])
              pval <- rbind(pval,p2)
            }
          }
        }
        mean <- rbind(mean,apply(pval[,-1],2,function(x){signif(x,3)}))
      }
            
      tmp <- NULL
      try(tmp <- kruskal.test(formula=clinical_data~groups)$p.value,TRUE)
      # Dunn test for pairwise comparisons in non normality comparisons
	  dunn<-dunnTest(as.numeric(clinical_data) ~as.factor(groups), method="bh")
	  pval_dunn<-dunn$res[,4]
	  names_pair<-as.character(dunn$res[,1])
	  if(any(dim(pval)==0)){
        pvalanova <- c(pvalanova,NA)
        pvalshap <-c(pvalshap,NA)
        clinvar <- c(clinvar,j)
        if(!is.null(tmp)){
          pvalkruskal <- c(pvalkruskal,tmp)
		  pvalPairwise<-rbind(pvalPairwise, as.numeric(pval_dunn))
        }else{
          pvalkruskal <- c(pvalkruskal,NA)
		  pvalPairwise<-rbind(pvalPairwise, as.numeric(pval_dunn))
        }
      }else{
        pvalanova <- c(pvalanova,NA,rep(NA,dim(pval)[1]))
        pvalshap <-c(pvalshap,NA,rep(NA,dim(pval)[1]))
        clinvar <- c(clinvar,j,rep(NA,dim(pval)[1]))
        if(contrast){
          contrasts_col <- c(contrasts_col,NA,pval[,1])
        }
        if(!is.null(tmp)){
          pvalkruskal <- c(pvalkruskal,tmp,rep(NA,dim(pval)[1]))
		  pvalPairwise <- rbind(pvalPairwise, as.numeric(pval_dunn))
        }else{
          pvalkruskal <- c(pvalkruskal,NA,rep(NA,dim(pval)[1]))
		  pvalPairwise <- rbind(pvalPairwise, as.numeric(pval_dunn))
        }
      }
      tmp <- NULL
    }
  }
        
  if(contrast){
#     table1 <- data.frame(clinvar,contrasts,mean,signif(pvalshap,3),signif(pvalkruskal,3),pvalanova)
    table1 <- data.frame(clinvar,contrasts_col,mean,signif(pvalshap,3),signif(pvalkruskal,3),pvalanova)
    names(table1) <- c("Variables","contrasts",paste(levels(groups),'(n=',summary(groups),')'),"pvalshapiro","pvalkruskal","pvalanova")
  }
  else{
    table1 <- data.frame(clinvar,mean,signif(pvalshap,3),signif(pvalkruskal,3),pvalanova, pvalPairwise)
    names(table1) <- c("Variables",paste(levels(groups),'(n=',summary(groups),')'),"pvalshapiro","pvalkruskal","pvalanova",names_pair)
  }
  row.names(table1) <- NULL
  return(table1)
}