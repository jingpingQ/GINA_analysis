
library(corrplot)
library(wrapr)
library(jsonlite)
library(lsa)
library(rlang)
library(betapart)
library(dplyr)
library(ggpubr)
library(xlsx)
library(purrr)
library(reshape)
library(car)

############################################################
# This script analyze the major hypothesis of the GINA manuscript
# Part I analyze the linearity and truncality of phylogenetic trees between TN and other clinical subtypes
# Part II analyze the correlation between pairwise heterogeneity (pair_beta) and physical distance
############################################################


##Part I

setwd("~/GINA_analysis/PartI")


subtype=readxl::read_xlsx("gina_meta.xlsx",sheet = 1)
#get the TN labels
subtype$tn="no"
subtype$tn[subtype$clinical_subtypes=="TN"]="yes"
#calculate linearity and truncality
subtype$linear=subtype$linear_nodes/subtype$all_nodes
subtype$truncality=subtype$trunk1/subtype$all_nodes

#store focality and subtype for future use
focality=subtype$focality
names(focality)=subtype$patient
clinical_subtypes=subtype$tn
names(clinical_subtypes)=subtype$patient


#boxplot to see the diff
ggpubr::ggboxplot(subtype,x="tn", y="linear",
                  color="tn")
ggpubr::ggboxplot(subtype,x="tn", y="truncality",
                  color="tn")

#calculate mean
mean(subtype$linear[subtype$tn=="no"])
mean(subtype$linear[subtype$tn=="yes"])
mean(subtype$truncality[subtype$tn=="no"])
mean(subtype$truncality[subtype$tn=="yes"])




#shapiro test for normality
shapiro.test(subtype$linear) #normal
leveneTest(linear ~ tn, data = subtype) #No significant difference in variances

shapiro.test(subtype$truncality) #not normal

p_linearity=t.test(subtype$linear~subtype$tn)$p.value
p_truncality=wilcox.test(subtype$truncality~subtype$tn)$p.value




##Part II
setwd("~/GINA_analysis/PartII")

#read spatial coordinates
dist <- read.csv("dist.csv",header = T,row.names = 1)
dist <- dist[ , apply(dist, 2, function(x) !any(is.na(x)))]


common <-intersect(sample_names,colnames(dist))
#some of the foci are speculated on pathological level but not sequenced due to lack of tissue
#remove them.
dist <- dist[,colnames(dist)%in% common]



#read Pairtree output in json format, this will include eta (cellular prevelance) and sample name
files <- list.files(pattern = "*.json")

#calculate pair-wise beta using eta
beta_pair_eta=list()
sample_names=c()
for (i in 1:length(files)) {
  
  
  
  tmp <- jsonlite::fromJSON(files[i])
  file_name <- gsub("\\..*","",files[i])
  tmp_names <- tmp[["samples"]]
  ##in common
  if(length(tmp_names[tmp_names%in% common])>0){
    
    #remove the one not in list
    p=which(!tmp_names%in% common)
    
  }
  tmp_names <- tmp_names[tmp_names%in% common]

  sample_names <- c(sample_names,tmp_names)
  
  tmp_eta <- tmp[["eta"]]
  
  if(ncol(tmp_eta)!=length(tmp_names)){
    tmp_eta=tmp_eta[,-p]
    
  }

  colnames(tmp_eta) <- tmp_names
  
  tmp_eta[tmp_eta>0.03]=1
  tmp_eta[tmp_eta!=1]=0

  

  input=tmp_eta
  betatmp=beta.pair(t(input),index.family = "sorensen")$beta.sor
  

  beta_pair_eta[[i]] <- betatmp
  names(beta_pair_eta)[[i]] <- file_name
}

beta_pair_eta = beta_pair_eta[lapply(beta_pair_eta,length)>0]



#unique(gsub("B.*","",colnames(dist)))
colnames(dist)=paste(gsub("B.*","",colnames(dist)),"_B",gsub(".*B","",colnames(dist)),sep="")
dist_list <-  map(set_names(unique(gsub("B.*","",colnames(dist)))),~dplyr::select(dist,starts_with(.x)))

dist_list <- lapply(dist_list, t)
dist_pair = lapply(dist_list,dist)


dist_pair_melt=list()
for (i in 1:length(dist_pair)) {
  
  
  dist_tmp = dist_pair[[i]]
  dist_melt <- reshape::melt(as.matrix(dist_tmp), varnames = c("site1", "site2"))
  
  #remove the 0 distance when its self-paired
  #remove duplicated pairs
  dist_melt <- dist_melt[!duplicated(dist_melt$value),]
  dist_melt <- dist_melt[!dist_melt$value==0,]
  colnames(dist_melt)[3] = "dist"
  dist_pair_melt[[i]]=dist_melt
  names(dist_pair_melt)[[i]] <- gsub("\\_","",names(dist_pair))[[i]]
}

eta_pair_melt=list()
for (i in 1:length(beta_pair_eta)) {
  
  
  eta_tmp = beta_pair_eta[[i]]
  eta_melt <- reshape::melt(as.matrix(eta_tmp), varnames = c("site1", "site2"))
  eta_melt <- eta_melt[eta_melt$site1 != eta_melt$site2,]
  cols = c(1,2)
  
  
  newdf = eta_melt[,cols]
  
  
  tmp_sorted=t(apply(newdf,1,sort))
  colnames(tmp_sorted)=c('site1','site2')
  
  eta_melt=eta_melt[!duplicated(tmp_sorted),]
  
  colnames(eta_melt)[3] = "sorensen_eta"
  
  eta_pair_melt[[i]]=eta_melt
  names(eta_pair_melt)[[i]] <- names(beta_pair_eta)[[i]]
}


#combine physical distance and paired beta
beta_eta_df <- do.call(rbind.data.frame, eta_pair_melt)

#data manipulation to prepare combining the physical distance and the paired beta
for(i in 1:length(dist_pair_melt)){
  names(dist_pair_melt)[[i]]=gsub("\\_","",names(dist_pair_melt)[[i]])
  tmp=as.data.frame(dist_pair_melt[[i]],stringsAsFactors = F)
  tmp$site1=as.character(tmp$site1)
  tmp$site2=as.character(tmp$site2)
  
  for(j in 1:nrow(tmp)){
    
    tmp[j,1]<-gsub("\\_","",tmp[j,1])
    tmp[j,2]=gsub("\\_","",tmp[j,2])
    
  }
  dist_pair_melt[[i]]=tmp
  
}

dist_df <- do.call(rbind.data.frame, dist_pair_melt)
#check rownames
#rownames(beta_eta_df) == rownames(dist_df)
dist_df$sor_eta = beta_eta_df$sorensen_eta
dist_df$patient= gsub("B.*","",dist_df$site1)



#correlating the physical distance and the paired beta, spearman=0.18
cor(dist_df$dist,dist_df$sor_eta,method="spearman")
cor(dist_df$dist[dist_df$focality=="multi"],dist_df$sor_eta[dist_df$focality=="multi"],method="spearman")
cor(dist_df$dist[dist_df$focality=="uni"],dist_df$sor_eta[dist_df$focality=="uni"],method="spearman")


#merge previous info from meta data
dist_df$subtype=clinical_subtypes[dist_df$patient]
dist_df$focality=focality[dist_df$patient]

#write.csv(dist_df,file = "dist.diff.csv")

#test correlation for focality and subtype, noticed they are not all normally distributed
#use spearman instead of pearson.
cor(dist_df$sor_eta[dist_df$subtype=="yes"],dist_df$dist[dist_df$subtype=="yes"],method = "spearman")
cor(dist_df$sor_eta[dist_df$subtype=="no"],dist_df$dist[dist_df$subtype=="no"],method = "spearman")
#spearman shows strong positive correlation (r=0.80) in TN and weak correlation (r=0.11) in non-TN
#adjust by sample size bias in the two groups

# Function to compute Fisher Z-transformation
fisher_z_transform <- function(correlation) {
  return(0.5 * log((1 + correlation) / (1 - correlation)))
}

# Function to compute the Z-score difference
fisher_z_adjustment <- function(corr1, n1, corr2, n2) {
  z1 <- fisher_z_transform(corr1)
  z2 <- fisher_z_transform(corr2)
  
  se1 <- 1 / sqrt(n1 - 3)
  se2 <- 1 / sqrt(n2 - 3)
  
  z_diff <- (z1 - z2) / sqrt(se1^2 + se2^2)
  
  return(z_diff)
}

# Example usage with given correlation coefficients and sample sizes
corr_no <- 0.1131  # Spearman correlation for group "non_TN"
n_no <-  length(subset(dist_df, subtype == "no")$dist)  # Sample size for "non_TN"

corr_yes <- 0.7967  # Spearman correlation for group "TN"
n_yes <- length(subset(dist_df, subtype == "yes")$dist)  # Sample size for "TN"

# Compute Fisher Z-score difference
z_score_diff <- fisher_z_adjustment(corr_no, n_no, corr_yes, n_yes)

# Print the result
print(z_score_diff)
#z-score=-2.98, p=0.00291
p_cor=0.00291

#adjust all p
pval=c(p_linearity,p_truncality,p_cor)
p.adjust(pval)



