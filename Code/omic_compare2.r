omic_compare2 <- function(data, # numeric matrix or data frame of omics. Patients must correspond to rows.
                         group, # column number or name of 'clin_data' corresponding to group definition.
                         covariate=NULL, # numeric or character vector of the column number(s) or name(s) of 'clin_data' corresponding to covariate(s).
                         clin_data, # matrix or data frame of clinical metadata, it must include group and covariate(s). Patients must correspond to rows.
                         digit=5 # number of significant digits to display.
){
  ## test if input parameters are as expected
  if(!all(sapply(c('data','group','clin_data'),function(x){hasArg(name = x)}))){
    stop("You must at least specify 'data','group', and 'clin_data'")
  }
  if(!all(apply(data,2,is.numeric))){
    stop("'data' must only contain numeric values")
  }
  if(is.null(clin_data)){
    stop("'clin_data' must at least contain a column for corresponding to the 'group' argument")
  }else if (is.matrix(clin_data)){
    clin_data <- data.frame(clin_data)
  }
  if( !is.numeric(group) & !is.character(group)){
    stop("Invalid specification of 'group'. 'group' must correspond to either a column name or number of 'clin_data'.")
  }else if(is.character(group) & ! group %in% colnames(clin_data)){
    stop(paste("Invalid specification of 'group'. ", group, " is not a column of 'clin_data"))
  }else if(is.numeric(group) & !group<=ncol(clin_data)){
    stop(paste("Invalid specification of 'group'. 'group'=", group, " but 'clin_data' has only ", ncol(clin_data)," column(s)"))
  }
  
  if( !is.null(covariate)){
    if( !is.numeric(covariate) & !is.character(covariate)){
      stop("Invalid specification of 'covariate'. 'covariate' must correspond to either column name(s) or number(s) of 'clin_data'.")
    }else if(is.character(covariate) & !all(covariate %in% colnames(clin_data))){
      stop(paste("Invalid specification of 'covariate'. ", covariate[!(covariate %in% colnames(clin_data))], " is not a column of 'clin_data"))
    }else if( is.numeric(covariate) & !all(covariate<=ncol(clin_data))){
      stop(paste("Invalid specification of 'covariate'. 'covariate'=", covariate, " but 'clin_data' has only ", ncol(clin_data)," column(s)"))
    }
  }
  
  if(!is.null(row.names(data))&!is.null(row.names(clin_data))){
    patients <- intersect(row.names(data),row.names(clin_data))
    if(length(patients)==0){
      if((nrow(data)==nrow(clin_data))){
        cat("There is no common row names in 'clin_data' and 'data'\n")
        cat("However both have same number of rows.\n") 
        #print(interactive())
        tmp <- readline(prompt = "Are you sure your that 'data' and 'clin_data' are in the same order? (y/n) ")
        if(!(tolower(tmp) %in% c('yes','y'))){
          #print(str(tmp))
          #print(paste('1',tmp))
          stop("'data' and 'clin_data' must be ordered the same way or contain patients ID as row names")
        }else{
          patients <- 1:nrow(data)
        }
      }else{
        stop("row names of 'data' and 'clin_data' must contains the patients ID or have same number of rows")
      }
    }
  } else if(is.null(row.names(data))|is.null(row.names(clin_data))){
    if((nrow(data)==nrow(clin_data))){
      cat(paste("Are you sure your that 'data' and 'clin_data' are in the same order? (y/n)\n"))
      tmp1 <- readline()
      if(!(tolower(tmp1) %in% c('yes','y'))){
        #print(paste('2',tmp1))
        stop("'data' and 'clin_data' must be ordered the same way or contain patients ID as row names")
      }else{
        patients <- 1:nrow(data)
      }
    }else{
      stop("row names of 'data' and 'clin_data' must contains the patients ID or have same number of rows")
    }
  }
  
  ## calculate table since everything is ok
  if(exists('patients')){
    ## install needed packages
    if(! all(c('car','multcomp') %in% installed.packages())){
      cat("You need 'car' and 'multcomp' packages. \n Do you want to install it now? (y/n)")
      tmp2 <- readline()
      if(tolower(tmp2) %in% c('yes','y')){
        #print(paste('3',tmp2))
        install.packages(pkgs = c('car','multcomp')[!(c('car','multcomp') %in% installed.packages())])
      }else{
        stop()
      }
    }
    ## load needed package
    sapply(c('car','multcomp'),function(x){require(x,character.only = T)})
    options(show.error.messages = FALSE,warn=-1)
    omics <- data[patients,]
    clin <- data.frame(clin_data[patients,c(group,covariate)])
    names(clin)[1] <- 'group'
    clin[,'group'] <- as.factor(as.character(clin[,'group']))
    res <- c()
    meansd <- sapply(levels(clin[,'group']),function(y){
      sapply(1:ncol(omics),function(x){
        paste(signif(mean(omics[clin[,'group']==y,x],na.rm=T),digit),' +/- ',signif(sd(omics[clin[,'group']==y,x],na.rm=T),digit),sep="")
      })
    })
    medianiqr <- sapply(levels(clin[,'group']),function(y){
      sapply(1:ncol(omics),function(x){
        paste(signif(median(omics[clin[,'group']==y,x],na.rm=T),digit),' (',paste(signif(quantile(omics[clin[,'group']==y,x],probs = c(0.25,0.75),na.rm=T),digit),collapse=";"),")",sep="")
      })
    })
    nas <- sapply(levels(clin[,'group']),function(y){
      sapply(1:ncol(omics),function(x){
      sum(is.na(omics[clin[,'group']==y,x]))})
    })
    meansd2 <- t(sapply(1:nrow(meansd),function(x){ifelse(nas[x,]==0,meansd[x,],paste( meansd[x,]," (Na's:",nas[x,],")",sep=""))}))
    medianiqr2 <- t(sapply(1:nrow(medianiqr),function(x){ifelse(nas[x,]==0,medianiqr[x,],paste( medianiqr[x,]," (Na's:",nas[x,],")",sep=""))}))
    pval_shapiro <- signif(sapply(1:ncol(omics),function(x){tryCatch(shapiro.test(omics[,x])$p.value,error=function(e){NaN})}),digit)
    pval_bartlett <- signif(sapply(1:ncol(omics),function(x){tryCatch(bartlett.test(omics[,x]~clin[,'group'])$p.value,error=function(e){NaN}) }),digit)
    pval_levene <- signif(sapply(1:ncol(omics),function(x){tryCatch(leveneTest(omics[,x]~clin[,'group'])$'Pr(>F)'[1],error=function(e){NaN})}),digit)
    #interaction <- (sapply(1:ncol(data),function(x){anov <- tryCatch(Anova(lm(formula=eval(parse(text=paste('data[,',x,']~',paste(rep('.',n),collapse = ":"),sep=""))),clin)),error=function(e){NaN}); min(anov[grep(':',row.names(anov)),'Pr(>F)'])  })<0.05)    
    pval_chisq <- signif(sapply(1:ncol(omics),function(x){tryCatch(chisq.test(as.character(is.na(omics[,x])),y=clin[,'group'])$p.value,error=function(e){NaN})}),digit)
    pval_anova <-signif(sapply(1:ncol(omics),function(x){tryCatch(Anova(lm(formula=omics[,x]~. ,clin))['group','Pr(>F)'],error=function(e){NaN}) }),digit)
    pval_fdr<-p.adjust(pval_anova, method='fdr')
    pval_Tukey_tmp <- lapply(1:ncol(omics),function(x){tryCatch({summary(glht(tryCatch(lm(formula=omics[,x]~.,clin)), linfct=mcp(group="Tukey", interaction_average = TRUE, covariate_average = TRUE)))},error=function(e){NaN})})
    pval_Tukey_tmp2 <- sapply(pval_Tukey_tmp,function(x){is.list(x)})
    if(any(pval_Tukey_tmp2)){
      contrast_tmp <- lapply(pval_Tukey_tmp[pval_Tukey_tmp2],function(x){row.names(x$linfct)})
      contrast <- sort(unique(unlist(contrast_tmp)))
      pval_tmp <- lapply(pval_Tukey_tmp[pval_Tukey_tmp2],function(x){x$test$pvalues})
      pval_Tukey <- matrix(NaN,nrow=length(pval_Tukey_tmp),ncol=length(contrast))
      for(i in 1:length(contrast)){
        feat_in_contrast <- sapply(contrast_tmp,function(x){contrast[i] %in% x})
        pval_Tukey[pval_Tukey_tmp2,i][feat_in_contrast] <- sapply(1:sum(feat_in_contrast),function(x){signif(pval_tmp[feat_in_contrast][[x]][which(contrast_tmp[feat_in_contrast][[x]]==contrast[i])],digit)})
      }
    }else{
      pval_Tukey <- pval_Tukey_tmp
      contrast <- 'Pval'
    }
    if(all(is.nan(pval_chisq))){
      res <- data.frame(meansd2, medianiqr2, pval_shapiro, pval_levene, pval_bartlett, pval_anova, pval_fdr, pval_Tukey)
      names(res) <- c(paste(levels(clin[,'group']),' (n=',summary(clin[,'group']), ') Mean +/- SD',sep=""),
                      paste(levels(clin[,'group']),' (n=',summary(clin[,'group']), ') Median (IQR)',sep=""),
                      'Pval_Shapiro','Pval_Levene','Pval_Bartlett','Pval_ANOVA*', 'Pval_FDR',
                      paste(contrast," (Tukey*)",sep=""))
      
      row.names(res) <- colnames(omics)
      cat('Pval_Shapiro tests for normality \nPval_Levene and Pval_Bartlett test for homogeneity of variance \nPval_ANOVA* tests if at least one group mean is different than the rest, adjusted for covariates (linear model without interactions) \nTukey* is a post-hoc pairwise test between groups, adjusted for covariates (linear model without interactions)\n')
      
    }else{
      res <- data.frame(meansd2, medianiqr2, pval_shapiro, pval_levene, pval_bartlett, pval_chisq, pval_anova, pval_fdr, pval_Tukey)
      names(res) <- c(paste(levels(clin[,'group']),' (n=',summary(clin[,'group']), ') Mean +/- SD',sep=""),
                      paste(levels(clin[,'group']),' (n=',summary(clin[,'group']), ') Median (IQR)',sep=""),
                      'Pval_Shapiro','Pval_Levene','Pval_Bartlett','Pval_Chisq','Pval_ANOVA*', 'Pval_FDR',
                      paste(contrast," (Tukey*)",sep=""))
      
      row.names(res) <- colnames(omics)
      cat("Pval_Shapiro tests for normality \nPval_Levene and Pval_Bartlett test for homogeneity of variance \nPval_Chisq tests if the frequency of Na's differs in at least one group than what you'd expect by chance \nPval_ANOVA* tests if at least one group mean is different than the rest, adjusted for covariates (linear model without interactions), Pvalue FDR for Benjamini-Hochberg multiple testing correction \nTukey* is a post-hoc pairwise test between groups, adjusted for covariates (linear model without interactions)\n")
      
    }
      
    options(show.error.messages = TRUE,warn=0)
    return(res)
  }
}
