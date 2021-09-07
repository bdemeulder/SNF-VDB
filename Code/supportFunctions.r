############################
##### Support functions ####
############################

keepSamples <- function(dataSet, samplesVector){
  ## Returns the subset containing only the specified samples.
  # @params <dataSet> the matrix from which the subset is extracted; rownames must have the same format as the values from <samplesVector
  # @params <samplesVector> names of the kept samples; must be of the same format as the rownames of <dataSet>
  # @return subset of <dataSet> containing only samples from <sampleVector>
  return(dataSet[match(samplesVector,rownames(dataSet)),])
}

overlapDatasets <- function(datasets){  ## Screen provided matrices to return the corresponding subset reduced to rows shared by all of them.
  #  All datasets must have similar rownames format.
  # @params <datasets> the list of matrices to reduce
  # @call keepSamples
  # @return overlap the list of reduced matrices
  commonRows = rownames(datasets[[1]])
  for(d in 2:length(datasets)){  commonRows = intersect( commonRows, rownames(datasets[[d]]) )  }
  overlap = lapply(datasets, function(d) keepSamples(d, commonRows))
  names(overlap) = names(datasets)
  cat("Datasets have been reduced to",length(commonRows),"rows. Their new dimensions hence are:\n")
  print(sapply(overlap,dim))
  return(overlap)
}

checkList <- function(dataList){
  ## Make sure a list has more than one dataset, and that all datasets have the same number of rows.
  # @params <dataList> list of datasets
  # @return TRUE if everything is alright; calls stop() otherwise

  if(!is.list(dataList) || !length(dataList)>1){  stop("dataList must be a list of more than one matrices.")  }
  nrowsEquals = sapply(dataList, function(d) nrow(d)==nrow(dataList[[1]]))

  if(sum(nrowsEquals)!=length(dataList)){  
    cat("Row number in each dataset: ", sapply(dataList,nrow),"\n")
    stop("All datasets must have the same number of rows.")
  }

}

checkTitleDirectory <- function(title){
  ## Checks if a folder that has the provided <title> exists in working directory. If not, creates it.
  # @params <title> the name of the directory
  # @return void
  if(dir.exists(title)){ # if a directory with this name exists, it means that the title has been used before
    cat("WARNING: a directory named '",title,"' already exists!\nThe existing files using this title may be overwritten if used again.\n", sep="")
    if(!askConfirmation()){ stop("Please make the necessary modifications.")  } # if the user doesn't agree to continue, the execution will stop
  }else{  # otherwise create it along with the subfolders
    dir.create(title)
    dir.create(paste0(title,"/graphs/"))
    dir.create(paste0(title,"/tmpData/"))
    dir.create(paste0(title,"/outData/"))
  } 
}

readInteger <- function(desc="Enter an integer: "){
  ## Asks user for an int value. <Desc> should be a message indicating the meaning of the value.
  # @params <desc> a string describing the int for which the value is asked
  # @return <n> the user-inputed int
  n <- readline(prompt=desc)
  n <- as.integer(n)
  if (is.na(n)){
    n <- readInteger(desc)
  }
  return(n)
}

readDouble <- function(desc="Enter a double: "){
  ## Asks user for a double value. <Desc> should be a message indicating the meaning of the value.
  # @params <desc> a string describing the double for which the value is asked
  # @return <d> the user-inputed double
  d <- readline(prompt=desc)
  d <- as.double(d)
  if (is.na(d)){
    d <- readDouble(desc)
  }
  return(d)
}

readString <- function(desc){
  ## Asks user for a string. <Desc> should be a message indicating the meaning of the string.
  # @params <desc> a string describing the asked string's meaning
  # @return <str> the user-inputed string
  str <- readline(prompt=desc)
  if (is.na(str)){
    str <- readString(desc)
  }
  return(str)
}

runSNF <- function(dataList, dataType=NULL, neighboursN=20, alpha=0.5, iterationsN=50){
  {
  ## Performs the SNF algorithm on provided data.
  # @params <dataList> a list where each element is one of the datasets that are to be merged together;
  #                    it can be of different types (see <dataType>), but all datasets should already have the same number of samples (rows), of the same name and in the same order
  # @params <dataType> vector indicating what kind of data each element of the <dataList> list is: 
  #                    (1) is for "raw" data that will be normalised before computing distances
  #                    (2) and (3) respectively stand for distance and affinity matrices
  #                    default value is NULL when all elements of the list are raw datasets. Otherwise, a vector of the same length as <dataList> must be provided
  # @params <neighboursN> number of neighbours, must be positive, usually between 10 and 30
  # @params <alpha> hyperparameter for the method, usually between 0.3 and 0.8
  # @params <iterationsN> number of iterations, usually between 20 and 50
  #
  # @return <mergedNetwork> the network resulting from merging elements form <dataList> using the SNF method
  } # runSNF header
  
  cat("=====  Merging method ===== \n")
  
  ## Make some verifications.
  checkList(dataList)
  if(!is.null(dataType) && length(dataType)!=length(dataList)){  stop("When provided, dataType must be of the same length as dataList.")  }

  ## Separate <dataList> elements according to their type.
  cat("Preparing datasets for merging...\n")
  if(is.null(dataType)){  dataType = rep(1, length(dataList))  } # create a vector with only ones if all are "raw" datasets
  raws = which(dataType==1)
  distances = which(dataType==2)
  affinities = which(dataType==3)
  
  # First, normalise and compute distances for raw datasets.
  if(length(raws)>0){ # if there are raw datasets
    rawDatasets = lapply(raws, function(raw) dataList[[raw]]) # assemble them in a separate list
    rawsNormed = lapply(rawDatasets, function(raw) standardNormalization(as.matrix(raw))) # to get a (0,1) distribution for continuous values 
    distancesMatrices = lapply(  rawsNormed, function(mat) 1-cor(t(mat))  ) # (sqrt( dist2(as.matrix(mat),as.matrix(mat)) ))/sqrt(ncol(mat))  ) # distance matrices
  }
  # Then add provided distance matrices and compute affinity matrices.
  if(length(distances)>0){ # if distances matrices were provided to the function
    if(exists("distancesMatrices")){ # if distance matrices were already calculated from raw datasets
      for(dist in 1:length(distances)){  distancesMatrices[[length(raws)+dist]] = dataList[[ distances[dist] ]]  } # add the provided ones to the list
    }else{  distancesMatrices = lapply(distances, function(dist) dataList[[dist]])  } # otherwise, create the list with only the provided ones
  }
  if(exists("distancesMatrices")){ # if distance matrices exist, either directly provided or computed from datasets
    affinityMatrices = lapply(distancesMatrices, function(d) affinityMatrix(d, neighboursN, alpha) ) # compute the affinity matrices
  }
  # Add affinity matrices provided straight to the runSNF function to the list.
  if(length(affinities)>0){
    if(exists("affinityMatrices")){
      for(aff in 1:length(affinities)){  affinityMatrices[[aff+length(raws)+length(distances)]] = dataList[[ affinities[aff] ]]   }
    }else{  affinityMatrices = lapply(dataList, function(d) d)  }
  }
  
  ## Perform the actual SNF method and return.
  cat("Start merging...\n")
  mergedNetwork = SNF(affinityMatrices, neighboursN, iterationsN)
  return(mergedNetwork)
  
} # end of runSNF()

spectralClustering_hook <<- function(this_dist,k){
  # Use the algorithm form the SNFtool package in external functions.
  res <- SNFtool::spectralClustering(as.matrix(this_dist),k)
  return(res)
}
consensClustAUC <- function(consensusMatrix){
  mlow <- consensusMatrix[lower.tri(consensusMatrix)]
  beta <- sort(unique(c(0,mlow,1)))
  CDF <- sapply(  beta, function(i) sum(consensusMatrix<=i)/(length(consensusMatrix))  )
  CDF_0.5 <- sum(mlow<=0.5)/(length(mlow))
  length_interval <- beta[-1]-beta[-length(beta)]
  integral <- sum(sapply(  1:(length(beta)-1), function(j) abs(CDF[j]-CDF_0.5)*length_interval[j]  ))
  return(integral)
}

