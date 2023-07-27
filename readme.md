## Abstract:
Predicting the eventual volume of tumor cells, that might proliferate from a given tumor, can help in cancer early detection and medical procedures planning to prevent their migration to other organs. In this work, a new statistical framework is proposed using Bayesian techniques for detecting the eventual volume of cells expected to proliferate from a Glioblastoma (GBM) tumor. Specifically, the tumor region is first extracted using a parallel image segmentation algorithm. Once the tumor region is determined, we are interested in the number of cells that can proliferate from this tumor until its survival time. For this, we construct the posterior distribution of the tumor cell numbers based on proposed likelihood function and a certain prior. We also determine the corresponding probability that no tumor cell goes undetected when we find the ultimate eventual volume. The model so developed gives excellent results on our dataset. It is expected that the model will work on any dataset where the changes are not measured with regards time. This is particularly crucial for predicting tumor volumes for a deadly tumor like glioblastoma. Furthermore, we extend the detection model and conduct a Bayesian regression analysis by incorporating radiomic features to discover %that their inclusion enhances the chances of those non-tumor cells remaining undetected. The main focus of the study was to develop a time-independent prediction model that can reliably predict the ultimate volume a malignant tumor of the fourth-grade severity can attain, and can also determine if incorporation of radiomic properties of the tumor enhances the chances of no malignant cells remaining undetected.

## Code Usage:
### (a) Simulation Run:

> Total number of voxels tested = 10
> 
> Frequency of Malignant cells detected at each voxel:
>
> ```ruby
> aj=c(61,9,34,83,14,21,5,11,8,5,4) ## malignant cells detected at each test
> ```
>
> 
> Frequency of Overall tumor cells detected at each voxel:
>
> ```ruby
> nj=c(226,107,119,307,74,21,33,225,268,125) ## overall tumor cells detected at each test
> ```
>
> 
> Choosing a sample of tests conducted until a malignant cell gets detected:
>
> ```ruby
> sub_aj = c(sample(aj[1],1),sample(aj[2],1),sample(aj[3],1),sample(aj[4],1),sample(aj[5],1),sample(aj[6],1),
>             sample(aj[7],1),sample(aj[8],1),sample(aj[9],1),sample(aj[10],1))
> ```
>
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
>
> Deriving the final distribution:
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
>
> 
