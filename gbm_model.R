library(readxl)
data = read_excel("F:\\Pre-operative_TCGA_GBM_NIfTI_and_Segmentations\\TCGA_GBM_Volume.xlsx")
data2 = read_excel("F:\\Pre-operative_TCGA_GBM_NIfTI_and_Segmentations\\TCGA_GBM_Intensity.xlsx")
data_ED = read_excel("F:\\Pre-operative_TCGA_GBM_NIfTI_and_Segmentations\\gbm.xlsx",sheet="ED")
data_ET = read_excel("F:\\Pre-operative_TCGA_GBM_NIfTI_and_Segmentations\\gbm.xlsx",sheet="ET")
data_NET = read_excel("F:\\Pre-operative_TCGA_GBM_NIfTI_and_Segmentations\\gbm.xlsx",sheet="NET")
data_TC = read_excel("F:\\Pre-operative_TCGA_GBM_NIfTI_and_Segmentations\\gbm.xlsx",sheet="TC")
data_WT = read_excel("F:\\Pre-operative_TCGA_GBM_NIfTI_and_Segmentations\\gbm.xlsx",sheet="WT")

summary(data)
summary(data2)

require(ggplot2)
library(ggplot2)
require(GGally)
library(GGally)
require(CCA)
library(CCA)
library(CCP)
library(BAS)
library(tidyverse)
library(caret)

aj = data[,5] ## VOLUME_ET
nj = data[,10] ## Volume of Whole Brain

##j=7
##aj=c(61,9,34,83,14,21,4)
##nj=c(226,107,119,307,74,21,33)
Aj=cumsum(aj)
A=sum(aj)
##d=1

## Considering 5 ROIs
Si1 = data[[5]] ## VOLUME_ET
Si2 = data[[6]] ## VOLUME_NET
Si3 = data[[7]] ## VOLUME_ED
Si4 = data[[8]] ## VOLUME_TC
Si5 = data[[9]] ## VOLUME_WT
Sij=list(Si1,Si2,Si3,Si4,Si5)

Sij_sum=lapply(Sij,sum)
p=mapply(FUN='/',Sij,Sij_sum)
sum(p[,1])
pij_sum=lapply(p,sum)

f_Si1=dbinom(Si1,sum(Si1),p[,1])
f_Si2=dbinom(Si2,sum(Si2),p[,2])
f_Si3=dbinom(Si3,sum(Si3),p[,3])
f_Si4=dbinom(Si4,sum(Si4),p[,4])
f_Si5=dbinom(Si5,sum(Si5),p[,5])

f_Sij=list(f_Si1,f_Si2,f_Si3,f_Si4,f_Si5)
f_Sij_sum=lapply(f_Sij,sum)
mean_Sij=lapply(f_Sij,mean)
sd_Sij=lapply(f_Sij,sd)
cv_Sij=mapply('/',sd_Sij,mean_Sij)
prior=mapply(FUN='*',f_Sij,cv_Sij)
posterior=lapply(prior,sum)
g_Sij=mapply('/',prior,posterior)
g_Sij_sum=lapply(g_Sij,sum)

R=1-pgeom(unlist(lapply(g_Sij,sum)),unlist(lapply(p,mean)))
V = unlist(Sij)
region = rep(c("ET", "NET", "ED", "TC", "WT"), each = 102)
##new_data = data.frame(R[1:102], R[103:204], R[205:306], R[307:408], R[409:510], Si1, Si2, Si3, Si4, Si5, V)
new_data = data.frame(R, V, region)
vol_region = rep(c("VOLUME_ET","VOLUME_NET","VOLUME_ED","VOLUME_TC","VOLUME_WT"), each = 102)
##new_data = data.frame(R[1:102], R[103:204], R[205:306], R[307:408], R[409:510], Si1, Si2, Si3, Si4, Si5, V)
new_data = data.frame(R, V, region)
data_fit = data.frame(R, V, region, vol_region)


mean_R = mean(R)
sd_R = sd(R)
mean_V = mean(V)
sd_V = sd(V)

## for simulation purpose
V1 = c(9525,68592,5899,31614,7338,17679,34935,70998,83517,117105,86271,37513)
R1 = c(0.9937805,0.9961458,0.9556447,0.9907852,0.8806554,0.9778719,0.9738203,
		0.9777224,0.9890068,0.9839206,0.9632160,0.9502363) 
lcl_V = V1 - 1.95996*(sd_V/1000)
ucl_V = V1 + 1.95996*(sd_V/1000)
lcl_R = R1 - 1.95996*(sd_R/100)
ucl_R = R1 + 1.95996*(sd_R/100)



## graphical presentation
theme_set(theme_minimal())
 fig = ggplot(data = new_data, aes(x = V, y = sort(R), color = region, group = region)) +
  geom_line(size = 1) + xlab("Eventual volume of undetected cancer cells in cubic mm.") + ylab("Probability that no cancer cells remain undetected ")
fig + theme(
  panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white"), 
axis.text = element_text(face="bold"),
axis.title.x = element_text(color="blue", size=14, face="bold"),
axis.title.y = element_text(color="#993333", size=14, face="bold"), 
legend.title = element_text(size=14, face="bold")
  )


##surv = rep(data$time,5)
##fit = lm(V~surv)
##summary(fit)


## Bayesian Regression & Cross Validation

## ED
fit1 = lm(VOLUME_ED~., data=data_ED)
summary(fit1)
v1 = round(abs(predict(fit1)))
reg1 = bas.lm(VOLUME_ED~., data=data_ED, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
v11 = round(abs(predict(reg1)$fit))
hist_ED = data_ED[,3:6]
##ggpairs(hist_ED)
spat_ED = data_ED[,7:15]
##ggpairs(spat_ED)
vol_ED = data_ED[,1]
Vol_hist_ED = cbind(vol_ED,hist_ED)
ggpairs(Vol_hist_ED)
Vol_spat_ED = cbind(vol_ED,spat_ED)
ggpairs(Vol_spat_ED)
Hist_Spat_ED = cbind(hist_ED, spat_ED)
ggpairs(Hist_Spat_ED)
##cc1 = cc(spat_ED, hist_ED)
##cc2 = comput(spat_ED, hist_ED, cc1)
##rho <- cc1$cor
##n <- dim(hist_ED)[1]
##p <- length(hist_ED)
q <- length(spat_ED)
##p.asym(rho, n, p, q, tstat = "Wilks")
cc3 = cc(spat_ED,vol_ED)
cc4 = comput(spat_ED, vol_ED, cc3)
rho2 = cc3$cor
n2 = dim(vol_ED)[1]
p2 = length(vol_ED)
p.asym(rho2,n2,p2,q,tstat = "Wilks")
cc5 = cc(hist_ED,vol_ED)
cc6 = comput(hist_ED,vol_ED,cc5)
rho3 = cc5$cor
q3 = length(hist_ED)
p.asym(rho3,n2,p2,q3,tstat="Wilks")
ED.glm <- glm(log(VOLUME_ED) ~ ., data = data_ED)
cv.err <- cv.glm(data_ED, ED.glm)$delta
cv.err.6 <- cv.glm(data_ED, ED.glm, K = 6)$delta


## ET
fit2 = lm(VOLUME_ET~., data=data_ET)
summary(fit2)
v2 = round(abs(predict(fit2)))
reg2 = bas.lm(VOLUME_ET~., data=data_ET, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
v21 = round(abs(predict(reg2)$fit))
hist_ET = data_ET[,2:5]
##ggpairs(hist_ET)
spat_ET = data_ET[,8:16]
##ggpairs(spat_ET)
Hist_Spat_ET = cbind(hist_ET, spat_ET)
ggpairs(Hist_Spat_ET)
vol_ET = data_ET[,1]
Vol_hist_ET = cbind(vol_ET,hist_ET)
ggpairs(Vol_hist_ET)
Vol_spat_ET = cbind(vol_ET,spat_ET)
ggpairs(Vol_spat_ET)
##cc1 = cc(spat_ET, hist_ET)
##cc2 = comput(spat_ET, hist_ET, cc1)
##rho <- cc1$cor
##n <- dim(hist_ET)[1]
##p <- length(hist_ET)
q <- length(spat_ET)
##p.asym(rho, n, p, q, tstat = "Wilks")
cc3 = cc(spat_ET,vol_ET)
cc4 = comput(spat_ET, vol_ET, cc3)
rho2 = cc3$cor
n2 = dim(vol_ET)[1]
p2 = length(vol_ET)
p.asym(rho2,n2,p2,q,tstat = "Wilks")
cc5 = cc(hist_ET,vol_ET)
cc6 = comput(hist_ET,vol_ET,cc5)
rho3 = cc5$cor
q3 = length(hist_ET)
p.asym(rho3,n2,p2,q3,tstat="Wilks")
ET.glm <- glm(log(VOLUME_ET) ~ ., data = data_ET)
cv.err <- cv.glm(data_ET, ET.glm)$delta
cv.err.6 <- cv.glm(data_ET, ET.glm, K = 6)$delta




## NET
fit3 = lm(VOLUME_NET~., data=data_NET)
summary(fit3)
v3 = round(abs(predict(fit3)))
reg3 = bas.lm(VOLUME_NET~., data=data_NET, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
summary(reg3)
v31 = round(abs(predict(reg3)$fit))
hist_NET = data_NET[,2:5]
##ggpairs(hist_NET)
spat_NET = data_NET[,8:16]
Hist_Spat_NET = cbind(hist_NET, spat_NET)
ggpairs(Hist_Spat_NET)
##ggpairs(spat_NET)
vol_NET = data_NET[,1]
Vol_hist_NET = cbind(vol_NET,hist_NET)
ggpairs(Vol_hist_NET)
Vol_spat_NET = cbind(vol_NET,spat_NET)
ggpairs(Vol_spat_NET)
##cc1 = cc(spat_NET, hist_NET)
##cc2 = comput(spat_NET, hist_NET, cc1)
##rho <- cc1$cor
##n <- dim(hist_NET)[1]
##p <- length(hist_NET)
q <- length(spat_NET)
##p.asym(rho, n, p, q, tstat = "Wilks")
cc3 = cc(spat_NET,vol_NET)
cc4 = comput(spat_NET, vol_NET, cc3)
rho2 = cc3$cor
n2 = dim(vol_NET)[1]
p2 = length(vol_NET)
p.asym(rho2,n2,p2,q,tstat = "Wilks")
cc5 = cc(hist_NET,vol_NET)
cc6 = comput(hist_NET,vol_NET,cc5)
rho3 = cc5$cor
q3 = length(hist_NET)
p.asym(rho3,n2,p2,q3,tstat="Wilks")
NET.glm <- glm(log(VOLUME_NET) ~ ., data = data_NET)
summary(NET.glm)
cv.err <- cv.glm(data_NET, NET.glm)$delta
cv.err.6 <- cv.glm(data_NET, NET.glm, K = 6)$delta




## TC
fit4 = lm(VOLUME_TC~., data=data_TC)
summary(fit4)
v4 = round(abs(predict(fit4)))
reg4 = bas.lm(VOLUME_TC~., data=data_TC, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
summary(reg4)
v41 = round(abs(predict(reg4)$fit))
spat_TC = data_TC[,3:9]
##ggpairs(spat_TC)
vol_TC = data_TC[,1]
Vol_spat_TC = cbind(vol_TC,spat_TC)
ggpairs(Vol_spat_TC)
cc1 = cc(spat_TC,vol_TC)
cc2 = comput(spat_TC, vol_TC, cc1)
rho = cc1$cor
n = dim(vol_TC)[1]
p = length(vol_TC)
q = length(spat_TC)
p.asym(rho,n,p,q,tstat = "Wilks")
TC.glm <- glm(log(VOLUME_TC) ~ ., data = data_TC)
summary(TC.glm)
cv.err <- cv.glm(data_TC, TC.glm)$delta
cv.err.6 <- cv.glm(data_TC, TC.glm, K = 6)$delta





## WT
fit5 = lm(VOLUME_WT~., data=data_WT)
summary(fit5)
v5 = round(abs(predict(fit5)))
reg5 = bas.lm(VOLUME_WT~., data=data_WT, prior = "BIC", modelprior = Bernoulli(0.5), include.always=~., n.models=1)
summary(reg5)
v51 = round(abs(predict(reg5)$fit))
spat_WT = data_WT[,2:10]
##ggpairs(spat_WT)
vol_WT = data_WT[,1]
Vol_spat_WT = cbind(vol_WT,spat_WT)
ggpairs(Vol_spat_WT)
cc1 = cc(spat_WT,vol_WT)
cc2 = comput(spat_WT, vol_WT, cc1)
rho = cc1$cor
n = dim(vol_WT)[1]
p = length(vol_WT)
q = length(spat_WT)
p.asym(rho,n,p,q,tstat = "Wilks")
WT.glm <- glm(log(VOLUME_WT) ~ ., data = data_WT)
summary(WT.glm)
cv.err <- cv.glm(data_WT, WT.glm)$delta
cv.err.6 <- cv.glm(data_WT, WT.glm, K = 6)$delta


## Correlation factor among Volume features
Vol = cbind(data_ED[,1],data_ET[,1],data_NET[,1],data_TC[,1],data_WT[,1])
ggpairs(Vol)



## Test run subjecting Eventual Volume obtained above as Response
Sij1=list(v1,v2,v3,v4,v5)
Sij1_sum=lapply(Sij1,sum)
p1=mapply(FUN='/',Sij1,Sij1_sum)
sum(p1[,1])
pij1_sum=lapply(p1,sum)

Sij2=list(v11,v21,v31,v41,v51)
Sij2_sum=lapply(Sij2,sum)
p2=mapply(FUN='/',Sij2,Sij2_sum)
sum(p2[,1])
pij2_sum=lapply(p2,sum)


f_Si11=dbinom(v1,sum(v1),p1[,1])
f_Si21=dbinom(v2,sum(v2),p1[,2])
f_Si31=dbinom(v3,sum(v3),p1[,3])
f_Si41=dbinom(v4,sum(v4),p1[,4])
f_Si51=dbinom(v5,sum(v5),p1[,5])

f_Si12=dbinom(v11,sum(v11),p2[,1])
f_Si22=dbinom(v21,sum(v21),p2[,2])
f_Si32=dbinom(v31,sum(v31),p2[,3])
f_Si42=dbinom(v41,sum(v41),p2[,4])
f_Si52=dbinom(v51,sum(v51),p2[,5])


f_Sij1=list(f_Si11,f_Si21,f_Si31,f_Si41,f_Si51)
f_Sij1_sum=lapply(f_Sij1,sum)
mean_Sij1=lapply(f_Sij1,mean)
sd_Sij1=lapply(f_Sij1,sd)
cv_Sij1=mapply('/',sd_Sij1,mean_Sij1)
prior1=mapply(FUN='*',f_Sij1,cv_Sij1)
posterior1=lapply(prior1,sum)
g_Sij1=mapply('/',prior1,posterior1)
g_Sij_sum1=lapply(g_Sij1,sum)


f_Sij2=list(f_Si12,f_Si22,f_Si32,f_Si42,f_Si52)
f_Sij2_sum=lapply(f_Sij2,sum)
mean_Sij2=lapply(f_Sij2,mean)
sd_Sij2=lapply(f_Sij2,sd)
cv_Sij2=mapply('/',sd_Sij2,mean_Sij2)
prior2=mapply(FUN='*',f_Sij2,cv_Sij2)
posterior2=lapply(prior2,sum)
g_Sij2=mapply('/',prior2,posterior2)
g_Sij_sum2=lapply(g_Sij2,sum)


R1=1-pgeom(unlist(lapply(g_Sij1,sum)),unlist(lapply(p1,mean)))
R2=1-pgeom(unlist(lapply(g_Sij2,sum)),unlist(lapply(p2,mean)))


V1 = unlist(Sij1)
new_data1 = data.frame(R1, V1, region)

V2 = unlist(Sij2)
new_data2 = data.frame(R2, V2, region)


theme_set(theme_minimal())
 fig1 = ggplot(data = new_data1, aes(x = V1, y = sort(R1), color = region, group = region)) +
  geom_line(size = 1) + xlab("Eventual volume of undetected cancer cells in cubic mm.") + ylab("Probability that no cancer cells remain undetected ")
fig1 + theme(
  panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white"), 
axis.text = element_text(face="bold"),
axis.title.x = element_text(color="blue", size=14, face="bold"),
axis.title.y = element_text(color="#993333", size=14, face="bold"), 
legend.title = element_text(size=14, face="bold")
  )



theme_set(theme_minimal())
fig2 = ggplot(data = new_data2, aes(x = V2, y = sort(R2), color = region, group = region)) +
  		geom_line(size = 1) + xlab("Eventual volume of undetected cancer cells in cubic mm.") + ylab("Probability that no cancer cells remain undetected ")
fig2 + theme(
  panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white"), 
axis.text = element_text(face="bold"),
axis.title.x = element_text(color="blue", size=14, face="bold"),
axis.title.y = element_text(color="#993333", size=14, face="bold"), 
legend.title = element_text(size=14, face="bold")
  )




#j=seq(1:510)
#plot(j,sort(R),type="l",xlab="Number of tests",ylab="Probability of tumor cells detected")


##sd_new=mapply('*',qnorm(mean(R)),sd_Sij)
##Sij_new=mapply('+',mean_Sij,sd_new)
##p_new=mapply('/',Sij_new,sum(Sij_new))
##cv_new=mapply('/',sd_new,mean_Sij)
##prior_new=mapply('*',f_Sij,cv_new)
##posterior_new=lapply(prior_new,sum)
##g_Sij_new=mapply('/',prior_new,posterior_new)
##g_Sij_sum_new=lapply(g_Sij_new,sum)
##R_new=1-pgeom(unlist(g_Sij_sum_new),unlist(lapply(p_new,mean)))
