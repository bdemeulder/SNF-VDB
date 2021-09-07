library(SNFtool)
library(ConsensusClusterPlus)
library(FSA)

source('supportFunctions.r')
source('tab_clin.r')
source('tab_clin2.r')
source('omic_compare2.r')

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

# Manual SNF

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
# Put the number of clusters you want to characterise in K
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
write.table(prot_table, file="prot_table.csv", sep=";")

# Characterise one cluster against the rest: example for K=4, comparing cluster #1 against the rest

group41<-cluster[[4]]$consensusClass
group41[which(group41!="1")]<-"4"
clinical$group41<-group41
clin_table41<-tab_clin(group=as.factor(clinical$group41), clin=clinical)
# Chi-square test (for not ordered categorical variables, e.g. sex), with 3 significant figures
signif(chisq.test(x=clinical$group41, y=clinical$Sex)$p.value, 3)
write.table(clin_table41, file="clin_table41.csv", sep=";")

# Omic tables
prot_table41<-omic_compare2(prot2, clin_data=clinical, group="group41")
write.table(prot_table41, file="prot_table41.csv", sep=";")
