# Quantifying the Growth of Glioblastoma Tumors (GBM) using Multimodal MRI Brain Images

## Abstract:
Predicting the eventual volume of tumor cells, that might proliferate from a given tumor, can help in cancer early detection and medical procedures planning to prevent their migration to other organs. In this work, a new statistical framework is proposed using Bayesian techniques for detecting the eventual volume of cells expected to proliferate from a Glioblastoma (GBM) tumor. Specifically, the tumor region is first extracted using a parallel image segmentation algorithm. Once the tumor region is determined, we are interested in the number of cells that can proliferate from this tumor until its survival time. For this, we construct the posterior distribution of the tumor cell numbers based on proposed likelihood function and a certain prior. We also determine the corresponding probability that no tumor cell goes undetected when we find the ultimate eventual volume. The model so developed gives excellent results on our dataset. It is expected that the model will work on any dataset where the changes are not measured with regards time. This is particularly crucial for predicting tumor volumes for a deadly tumor like glioblastoma. Furthermore, we extend the detection model and conduct a Bayesian regression analysis by incorporating radiomic features to discover that their inclusion enhances the chances of those non-tumor cells remaining undetected. The main focus of the study was to develop a time-independent prediction model that can reliably predict the ultimate volume a malignant tumor of the fourth-grade severity can attain, and can also determine if incorporation of radiomic properties of the tumor enhances the chances of no malignant cells remaining undetected.

## Code Usage:
### (a) Simulation Run:

Total number of voxels tested = 10 
Frequency of Malignant cells detected at each voxel:
>
> ```ruby
> aj=c(61,9,34,83,14,21,5,11,8,5,4) ## malignant cells detected at each voxel
> ```

Frequency of Overall tumor cells detected at each voxel:
>
> ```ruby
> nj=c(226,107,119,307,74,21,33,225,268,125) ## overall tumor cells detected at each voxel
> ```
>
Choosing a sample of tests conducted until a malignant cell gets detected:
>
> ```ruby
> sub_aj = c(sample(aj[1],1),sample(aj[2],1),sample(aj[3],1),sample(aj[4],1),sample(aj[5],1),sample(aj[6],1),            >   
>             sample(aj[7],1),sample(aj[8],1),sample(aj[9],1),sample(aj[10],1))
> ```

> ```ruby
> if1 = sample(aj[1],sub_aj[1],replace=FALSE,prob=NULL) ## number of tests conducted in 1st voxel until a malignant cell gets detected
> if2 = sample(aj[2],sub_aj[2],replace=FALSE,prob=NULL) ## number of tests conducted in 2nd voxel until a malignant cell gets detected
> if3 = sample(aj[3],sub_aj[3],replace=FALSE,prob=NULL) ## number of tests conducted in 3rd voxel until a malignant cell gets detected
> if4 = sample(aj[4],sub_aj[4],replace=FALSE,prob=NULL) ## number of tests conducted in 4th voxel until a malignant cell gets detected
> if5 = sample(aj[5],sub_aj[5],replace=FALSE,prob=NULL) ## number of tests conducted in 5th voxel until a malignant cell gets detected
> if6 = sample(aj[6],sub_aj[6],replace=FALSE,prob=NULL) ## number of tests conducted in 6th voxel until a malignant cell gets detected
> if7 = sample(aj[7],sub_aj[7],replace=FALSE,prob=NULL) ## number of tests conducted in 7th voxel until a malignant cell gets detected
> if8 = sample(aj[8],sub_aj[8],replace=FALSE,prob=NULL) ## number of tests conducted in 8th voxel until a malignant cell gets detected
> if9 = sample(aj[9],sub_aj[9],replace=FALSE,prob=NULL) ## number of tests conducted in 9th voxel until a malignant cell gets detected
> if10 = sample(aj[10],sub_aj[10],replace=FALSE,prob=NULL) ## number of tests conducted in 10th voxel until a malignant cell gets detected
> Sij = list(if1,if2,if3,if4,if5,if6,if7,if8,if9,if10)
> Sij_sum = lapply(Sij,sum)
> p = mapply(FUN='/',Sij,Sij_sum)
> pij_sum = lapply(p,sum)
> ```

Evaluation of our model:
> ```ruby
> f_Si1 = dbinom(if1,sum(if1),p[[1]])
> f_Si2 = dbinom(if2,sum(if2),p[[2]])
> f_Si3 = dbinom(if3,sum(if3),p[[3]])
> f_Si4 = dbinom(if4,sum(if4),p[[4]])
> f_Si5 = dbinom(if5,sum(if5),p[[5]])
> f_Si6 = dbinom(if6,sum(if6),p[[6]])
> f_Si7 = dbinom(if7,sum(if7),p[[7]])
> f_Si8 = dbinom(if8,sum(if8),p[[8]])
> f_Si9 = dbinom(if9,sum(if9),p[[9]])
> f_Si10 = dbinom(if10,sum(if10),p[[10]])
> ```

Deriving the final distribution:
> ```ruby
> f_Sij = list(f_Si1,f_Si2,f_Si3,f_Si4,f_Si5,f_Si6,f_Si7, f_Si8, f_Si9, f_Si10)
> f_Sij_sum = lapply(f_Sij,sum)
> mean_Sij = lapply(f_Sij,mean)
> sd_Sij = lapply(f_Sij,sd)
> cv_Sij = mapply('/',sd_Sij,mean_Sij)
> prior = mapply(FUN='*',f_Sij,cv_Sij)
> posterior = lapply(prior,sum)
> g_Sij = mapply('/',prior,posterior)
> g_Sij_sum = lapply(g_Sij,sum)
> ```


### (b) Running the proposed method on our data:
We have considered 5 regions of interest or 5 voxel sets - ET, NET, ED, TC, WT
> ```ruby
> Si1 = data[[5]] ## VOLUME_ET
> Si2 = data[[6]] ## VOLUME_NET
> Si3 = data[[7]] ## VOLUME_ED
> Si4 = data[[8]] ## VOLUME_TC
> Si5 = data[[9]] ## VOLUME_WT
> Sij=list(Si1,Si2,Si3,Si4,Si5)
> Sij_sum=lapply(Sij,sum)
> p=mapply(FUN='/',Sij,Sij_sum)
> sum(p[,1])
> pij_sum=lapply(p,sum)
>
> f_Si1=dbinom(Si1,sum(Si1),p[,1])
> f_Si2=dbinom(Si2,sum(Si2),p[,2])
> f_Si3=dbinom(Si3,sum(Si3),p[,3])
> f_Si4=dbinom(Si4,sum(Si4),p[,4])
> f_Si5=dbinom(Si5,sum(Si5),p[,5])
> f_Sij=list(f_Si1,f_Si2,f_Si3,f_Si4,f_Si5)
> f_Sij_sum=lapply(f_Sij,sum)
> mean_Sij=lapply(f_Sij,mean)
> sd_Sij=lapply(f_Sij,sd)
> cv_Sij=mapply('/',sd_Sij,mean_Sij)
> prior=mapply(FUN='*',f_Sij,cv_Sij)
> posterior=lapply(prior,sum)
> g_Sij=mapply('/',prior,posterior)
> g_Sij_sum=lapply(g_Sij,sum)
>
> R=1-pgeom(unlist(lapply(g_Sij,sum)),unlist(lapply(p,mean)))
> ```

### (c) Running a Bayesian Regression
Some specific libraries are used. 
> ```ruby
> require(GGally)
> library(GGally)
> require(CCA)
> library(CCA)
> library(CCP)
> library(BAS)
> library(tidyverse)
> library(caret)
> ```
>
Consider the region of interest "Edema" (ED). We run a Bayesian regression and then compare the results with a GLM via some cross validation tools. 
> ```ruby
> reg1 = bas.lm(VOLUME_ED~., data=data_ED, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
> v11 = round(abs(predict(reg1)$fit))
> hist_ED = data_ED[,3:6]
> spat_ED = data_ED[,7:15]
> vol_ED = data_ED[,1]
> Vol_hist_ED = cbind(vol_ED,hist_ED)
> ggpairs(Vol_hist_ED)
> Vol_spat_ED = cbind(vol_ED,spat_ED)
> ggpairs(Vol_spat_ED)
> Hist_Spat_ED = cbind(hist_ED, spat_ED)
> ggpairs(Hist_Spat_ED)
> q <- length(spat_ED)
> cc3 = cc(spat_ED,vol_ED)
> cc4 = comput(spat_ED, vol_ED, cc3)
> rho2 = cc3$cor
> n2 = dim(vol_ED)[1]
> p2 = length(vol_ED)
> p.asym(rho2,n2,p2,q,tstat = "Wilks")
> cc5 = cc(hist_ED,vol_ED)
> cc6 = comput(hist_ED,vol_ED,cc5)
> rho3 = cc5$cor
> q3 = length(hist_ED)
> p.asym(rho3,n2,p2,q3,tstat="Wilks")
> ED.glm <- glm(log(VOLUME_ED) ~ ., data = data_ED)
> cv.err <- cv.glm(data_ED, ED.glm)$delta
> cv.err.6 <- cv.glm(data_ED, ED.glm, K = 6)$delta
> ```
>
The same process has been repeated for other regions of interests. 
Let us next consider the region of interest "GD-Enhancing Tumor" (ET).

> ```ruby
> reg2 = bas.lm(VOLUME_ET~., data=data_ET, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
> v21 = round(abs(predict(reg2)$fit))
> hist_ET = data_ET[,2:5]
> spat_ET = data_ET[,8:16]
> Hist_Spat_ET = cbind(hist_ET, spat_ET)
> ggpairs(Hist_Spat_ET)
> vol_ET = data_ET[,1]
> Vol_hist_ET = cbind(vol_ET,hist_ET)
> ggpairs(Vol_hist_ET)
> Vol_spat_ET = cbind(vol_ET,spat_ET)
> ggpairs(Vol_spat_ET)
> q <- length(spat_ET)
> cc3 = cc(spat_ET,vol_ET)
> cc4 = comput(spat_ET, vol_ET, cc3)
> rho2 = cc3$cor
> n2 = dim(vol_ET)[1]
> p2 = length(vol_ET)
> p.asym(rho2,n2,p2,q,tstat = "Wilks")
> cc5 = cc(hist_ET,vol_ET)
> cc6 = comput(hist_ET,vol_ET,cc5)
> rho3 = cc5$cor
> q3 = length(hist_ET)
> p.asym(rho3,n2,p2,q3,tstat="Wilks")
> ET.glm <- glm(log(VOLUME_ET) ~ ., data = data_ET)
> cv.err <- cv.glm(data_ET, ET.glm)$delta
> cv.err.6 <- cv.glm(data_ET, ET.glm, K = 6)$delta
> ```
>

Let us next consider the region of interest "Non-Enhancing Tumor" (NET).
> ```ruby
> reg3 = bas.lm(VOLUME_NET~., data=data_NET, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
> summary(reg3)
> v31 = round(abs(predict(reg3)$fit))
> hist_NET = data_NET[,2:5]
> spat_NET = data_NET[,8:16]
> Hist_Spat_NET = cbind(hist_NET, spat_NET)
> ggpairs(Hist_Spat_NET)
> vol_NET = data_NET[,1]
> Vol_hist_NET = cbind(vol_NET,hist_NET)
> ggpairs(Vol_hist_NET)
> Vol_spat_NET = cbind(vol_NET,spat_NET)
> ggpairs(Vol_spat_NET)
> q <- length(spat_NET)
> cc3 = cc(spat_NET,vol_NET)
> cc4 = comput(spat_NET, vol_NET, cc3)
> rho2 = cc3$cor
> n2 = dim(vol_NET)[1]
> p2 = length(vol_NET)
> p.asym(rho2,n2,p2,q,tstat = "Wilks")
> cc5 = cc(hist_NET,vol_NET)
> cc6 = comput(hist_NET,vol_NET,cc5)
> rho3 = cc5$cor
> q3 = length(hist_NET)
> p.asym(rho3,n2,p2,q3,tstat="Wilks")
> NET.glm <- glm(log(VOLUME_NET) ~ ., data = data_NET)
> summary(NET.glm)
> cv.err <- cv.glm(data_NET, NET.glm)$delta
> cv.err.6 <- cv.glm(data_NET, NET.glm, K = 6)$delta
> ```

Let us next consider the region of interest "Tumor Core" (TC).
> ```ruby
> reg4 = bas.lm(VOLUME_TC~., data=data_TC, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
> summary(reg4)
> v41 = round(abs(predict(reg4)$fit))
> spat_TC = data_TC[,3:9]
> vol_TC = data_TC[,1]
> Vol_spat_TC = cbind(vol_TC,spat_TC)
> ggpairs(Vol_spat_TC)
> cc1 = cc(spat_TC,vol_TC)
> cc2 = comput(spat_TC, vol_TC, cc1)
> rho = cc1$cor
> n = dim(vol_TC)[1]
> p = length(vol_TC)
> q = length(spat_TC)
> p.asym(rho,n,p,q,tstat = "Wilks")
> TC.glm <- glm(log(VOLUME_TC) ~ ., data = data_TC)
> summary(TC.glm)
> cv.err <- cv.glm(data_TC, TC.glm)$delta
> cv.err.6 <- cv.glm(data_TC, TC.glm, K = 6)$delta
> ```
>


Let us next consider the region of interest "Whole Tumor" (WT).
> ```ruby
> reg5 = bas.lm(VOLUME_WT~., data=data_WT, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
> summary(reg5)
> v51 = round(abs(predict(reg5)$fit))
> spat_WT = data_WT[,2:10]
> vol_WT = data_WT[,1]
> Vol_spat_WT = cbind(vol_WT,spat_WT)
> ggpairs(Vol_spat_WT)
> cc1 = cc(spat_WT,vol_WT)
> cc2 = comput(spat_WT, vol_WT, cc1)
> rho = cc1$cor
> n = dim(vol_WT)[1]
> p = length(vol_WT)
> q = length(spat_WT)
> p.asym(rho,n,p,q,tstat = "Wilks")
> WT.glm <- glm(log(VOLUME_WT) ~ ., data = data_WT)
> summary(WT.glm)
> cv.err <- cv.glm(data_WT, WT.glm)$delta
> cv.err.6 <- cv.glm(data_WT, WT.glm, K = 6)$delta
> ```

The detailed codes are attached in the files test_run.R and gbm_model.R. 

## Data Source: 
> Scarpace, L. The Cancer Imaging Archive http://doi.org/10.7937/K9/TCIA.2016.RNYFUYE9 (2016)


## Acknowledgement: 
> Preprocessing of images using parallel segmentation and masking of the tumor location is from https://github.com/lelechen63/MRI-tumor-segmentation-Brats.
  

