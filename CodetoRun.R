library(SNFtool)
library(ConsensusClusterPlus)
library(FSA)

source('code/supportFunctions.r')
source('code/tab_clin.r')
source('code/tab_clin2.r')
source('code/omic_compare2.r')

############################
#### Handprint analysis ####
############################

#### Data preparation ####
# The script expects the omics data to be patients in row and features in columns
# The lists of patients do not have to be the same across platforms, the function will create the overlap of patients automatically
# Transcriptomics data should be log2 transformed

# dSets = list(Data1, Data2, Data3)
# Replace DataX with your formatted data matrices
# names(dSets) = c("Transcriptomics", "Proteomics", "Metabolomics")
# Replace names with relevant names

# dataType = rep(1, length(dSets))
# reduced = overlapDatasets(dSets)
# mergedNetwork = runSNF (reduced, dataType=dataType)

prot<-as.matrix(read.csv2("Proteomics.csv", sep=";", row.names=1))
prot2<-matrix(ncol=ncol(prot), nrow=nrow(prot))
for (i in 1:nrow(prot2))
{
  for (j in 1:ncol(prot2))
  {
    prot2[i, j]<-as.numeric(prot[i,j])
  }
}
row.names(prot2)<-row.names(prot)
colnames(prot2)<-colnames(prot)

metabo<-as.matrix(read.csv2("Metabolomics.csv", sep=";", row.names=1))
metabo2<-matrix(ncol=ncol(metabo), nrow=nrow(metabo))
for (i in 1:nrow(metabo2))
{
  for (j in 1:nrow(metabo2))
  {
    metabo2[i,j]<-as.numeric(metabo[i,j])
  }
}
row.names(metabo2)<-row.names(metabo)
colnames(metabo2)<-colnames(metabo)

clinical<-as.matrix(read.csv2("clinical.csv", sep=";", row.names=1))
clin.names<-clinical[,1]
clin<-clinical[, 2:ncol(clinical)]

clinical2<-matrix(ncol=ncol(clin), nrow=nrow(clin))
for (i in 1:nrow(clin))
{
  for (j in 1:ncol(clin))
  {
    clinical2[i, j]<-as.numeric(clin[i, j])
  }
}

row.names(clinical2)<-clin.names
colnames(clinical2)<-colnames(clin)
clinical<-as.data.frame(clinical2)

###########################
#### Linear correction ####
###########################

# Add any variable from the clinical data which effects we want to remove.
# ! Names in covariates have to be the same as in the clinical data.
covariates=c("age.highest", "sex")

# Example for the proteomics dataset
# Don't forget to replace the corrected datasets in the SNF below
prot2_Corr<-regressOutCovsEffect(dataSet=prot2, clinicalSet=clinical, covariates=covariates, signThreshold = 0.05)


################
## Manual SNF ##
################

dist1<-as.matrix(dist(prot2))
dist2<-as.matrix(dist(metabo2))

W1<-affinityMatrix(dist1, K=20, sigma=0.5)
W2<-affinityMatrix(dist2, K=20, sigma=0.5)
mergedNetwork<-SNF(list(W1, W2), K=20, t=20)

# We run SNF with default parameters, as recommended by the authors. 
# Parameters are: K (N of neighbours) (default = 20), hyperparameter sigma (default = 0.5) and t (number of iterations, default = 50). 

nk<-estimateNumberOfClustersGivenGraph(mergedNetwork, NUMC=2:20)
# In nk is the estimated number of clusters given by SNF without the consensus clustering step. 
# It gives four results based on different metrics; usually there will be agreement between nk and the consensus clustering step below.

############################
#### Consensus Clustering ##
############################

# This function performes the consensus clustering step. 
# MaxK = maximum of clusters to test. 
# reps = number of bootstrap iterations (usually 50-100).
cluster<-ConsensusClusterPlus(as.dist(mergedNetwork), maxK=20, reps=50, clusterAlg="spectralClustering_hook")

# Computing the deviation from ideal stability

DIS<-c()
for (i in 2:length(cluster))
	{
	temp<-consensClustAUC(cluster[[i]]$consensusMatrix)
	DIS[i-1]<-temp
	}
	
plot(DIS, type="b", xlab="Number of clusters", ylab="Deviation from Ideal Stability", main="Deviation from Ideal Stability for cluster numbers")

# You can then choose the rigt number of clusters by looking at the DIS plot (look for a local minimum in the curve), 
# combined with the consensus CDF plot (look for the flattest curve).

# Generating PDFs for Consensus clustering and DIS plot.

pdf("DISplot.pdf")
plot(DIS, type="b", xlab="Number of clusters", ylab="Deviation from Ideal Stability", main="Deviation from Ideal Stability for cluster numbers")
dev.off()

pdf("ConsensusClustering.pdf")
cluster<-ConsensusClusterPlus(as.dist(mergedNetwork), maxK=20, reps=50, clusterAlg="spectralClustering_hook")
dev.off()

save(mergedNetwork, file="mergedNetwork.Rdata", compress=T)
save(cluster, file="clustering.Rdata", compress=T)

#################################
#### Clusters characterisation ##
#################################

# This section will produce the statistical comparison tables between clusters for the clinical and omics data
# This assumes you have your clinical data as a data frame with patients in lines and features in rows. 
# Put the number of clusters you want to characterise in k
k = 3 
group<-cluster[[k]]$consensusClass
test<-rownames(mergedNetwork)
names(group)<-test
clinical$group<-group

# Producing the clinical tables

clin_table<-tab_clin(group=as.factor(clinical$group), clin=clinical)
# Wilks-Shapiro test: normality test. If p<0.05, data is not normal and you should use non-parametric test (Kruskal)
# Kruskal-Wallis test: statistical test on medians (if data is not normal)
# Anova test: statistical test on means (if data is normal)

# Chi-square test (for not ordered categorical variables, e.g. sex), with 3 significant figures
signif(chisq.test(x=clinical$group, y=clinical$Sex)$p.value, 3)
write.table(clin_table, file="clin_table.csv", sep=";")

# Producing the clinical table with all pairwise contrasts
# TukeyHSD test if data is normal, Dunn test if data is not normal, as per Wilks-Shapiro test

clin_table2<-tab_clin_cont(clinical=clinical, group=as.factor(clinical$group))
write.table(clin_table2, file="clin_table2.csv", sep=";")

# Producing the omics tables

prot_table<-omic_compare2(prot2, clin_data=clinical, group="group")
prot_table$FDR<-p.adjust(prot_table$'Pval_ANOVA*', method='fdr')
write.table(prot_table, file="prot_table.csv", sep=";")

########################################################################################################
### Characterise one cluster against the rest: example for K=4, comparing cluster #1 against the rest ##
########################################################################################################

group41<-cluster[[4]]$consensusClass
group41[which(group41!="1")]<-"4"
clinical$group41<-group41
clin_table41<-tab_clin(group=as.factor(clinical$group41), clin=clinical)
# Chi-square test (for not ordered categorical variables, e.g. sex), with 3 significant figures
signif(chisq.test(x=clinical$group41, y=clinical$Sex)$p.value, 3)
write.table(clin_table41, file="clin_table41.csv", sep=";")

# Omic tables
prot_table41<-omic_compare2(prot2, clin_data=clinical, group="group41")
prot_table41$FDR<-p.adjust(prot_table41$'Pval_ANOVA*', method='fdr')
write.table(prot_table41, file="prot_table41.csv", sep=";")

####################
## Alluvial plots ##
####################

library(ggplot2)
library(ggalluvial)
library(plyr)

.prePdata <- function(pcluster_list){
        temp <- pcluster_list
        names(temp) <- NULL
        sample <- unique(names(unlist(temp)))
        rm(temp)
        pdata_raw <- matrix(NA,nrow=length(sample),ncol=length(pcluster_list))
        rownames(pdata_raw) = sample
        for (i in 1:ncol(pdata_raw)){
                sample2 <- intersect(sample,names(pcluster_list[[i]]))
                pdata_raw[sample2,i] = pcluster_list[[i]]
                }
        pdata_raw <- pdata_raw[rowSums(is.na(pdata_raw))==0,]
        pdata_raw <- as.data.frame(pdata_raw)
        colnames(pdata_raw) <- names(pcluster_list)
        
        pdata.res <- count(pdata_raw)
        res <- pdata.res
        for (j in 1:(ncol(res)-1)){res[,j]<-as.factor(res[,j])}
        return(res)
}

# For alluvial plots for more than two clusters: (below is an example for 4 clusters with K=3, 6, 8 and 10)
# - Place the interesting clusters in clusterx. Place them in ascending number of clusters. 
# - place them all in pcluster
# - change the numbers in pdata to reflect the correct number of clusterings to plot
# - change the names(pdata) to correspond to the correct length. There should always be 'freq' added last
# - In the ggplot, add the correct list of clusterings in aes(y=freq, axis1=Cluster1, axis2=Cluster2...). The names here have to be the same as the names of pdata.
# - The colors of the alluvial are given by geom_flow(aes(fill=Cluster2)... 
# - Change the limit of breaks=1:4 for the number of clusterings 
# - Change the labels to reflect your choice of clusterings
# - Change the title of the figure in ggtitle

cluster1<-cluster[[3]]$consensusClass
cluster2<-cluster[[6]]$consensusClass
cluster3<-cluster[[8]]$consensusClass
cluster4<-cluster[[10]]$consensusClass

pcluster<-list(cluster1, cluster2, cluster3, cluster4)
pdata_all<-.prePdata(pcluster)
pdata<-.prePdata(pcluster[c(1, 2, 3, 4)])
# You can replace the Cluster 1 and Cluster2 names below to more relevant names
names(pdata)<-c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "freq")
cbbPalette<-c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# In the plot function, change ClusterX with the same names as in names(pdata)
# geom_flow(aes(fill=Cluster1)) decides which clustering should be used to color the plot. 
# labels=c("K=3", "K=6") also needs to be updated with your choice of number of clusters.
ggplot(pdata,
  aes(y = freq, axis1 = Cluster1, axis2 = Cluster2, axis3 = Cluster3, axis4 = Cluster4)) +
  geom_flow(aes(fill=Cluster2), aes.bind =TRUE, reverse = TRUE, width= 1/10) +
  geom_stratum(width= 1/10, reverse = TRUE) +
  geom_text(stat = "stratum", size= 3.5, label.strata = TRUE, reverse = TRUE) +
  scale_x_continuous(breaks=1:4, labels=c("K = 3", "K = 6", "K=8", "K=10")) +
  scale_fill_manual(values=cbbPalette) +
  theme(panel.background = element_rect(fill=NA), legend.position = "bottom", 
        axis.text.y = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.ticks.y = element_line(colour="white")) + 
  ggtitle("Stable patient allocations for K=3, K=6, K=8 and K=10 clusterings") 
  
