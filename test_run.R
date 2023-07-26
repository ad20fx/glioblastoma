j=10 ## 10 tests to be done
aj=c(61,9,34,83,14,21,5,11,8,5,4) ## malignant cells detected at each test
nj=c(226,107,119,307,74,21,33,225,268,125) ## overall tumor cells detected at each test
Aj=cumsum(aj)
A=sum(aj)
Si1=c(rep(c(38,22),c(30,30)),1)
Si2=c(9,9,9,9,2,2,2,1,1)
Si3=c(rep(c(11,18,15),c(11,11,11)),1)
Si4=rep(c(74,9,1,1,2,1),c(22,21,20,15,3,2))
Si5=rep(c(2,2,2,1,1,1,1),each=2)
Si6=rep(c(2,1,1),each=7)
Si7=c(3,1,1,1,4,5)
Si8=c(5,4,1,2,3,1,4,6,1,2,1)
Si9=c(6,8,5,4,1)
Si10=c(2,1,3,1)
Sij=list(Si1,Si2,Si3,Si4,Si5,Si6,Si7,Si8, Si9, Si10)

Sij_sum=lapply(Sij,sum)
p=mapply(FUN='/',Sij,Sij_sum)
pij_sum=lapply(p,sum)

f_Si1=dbinom(Si1,sum(Si1),p[[1]])
f_Si2=dbinom(Si2,sum(Si2),p[[2]])
f_Si3=dbinom(Si3,sum(Si3),p[[3]])
f_Si4=dbinom(Si4,sum(Si4),p[[4]])
f_Si5=dbinom(Si5,sum(Si5),p[[5]])
f_Si6=dbinom(Si6,sum(Si6),p[[6]])
f_Si7=dbinom(Si7,sum(Si7),p[[7]])
f_Si8=dbinom(Si8,sum(Si8),p[[8]])
f_Si9=dbinom(Si9,sum(Si9),p[[9]])
f_Si10=dbinom(Si10,sum(Si10),p[[10]])

f_Sij=list(f_Si1,f_Si2,f_Si3,f_Si4,f_Si5,f_Si6,f_Si7, f_Si8, f_Si9, f_Si10)
f_Sij_sum=lapply(f_Sij,sum)
mean_Sij=lapply(f_Sij,mean)
sd_Sij=lapply(f_Sij,sd)
cv_Sij=mapply('/',sd_Sij,mean_Sij)
prior=mapply(FUN='*',f_Sij,cv_Sij)
posterior=lapply(prior,sum)
g_Sij=mapply('/',prior,posterior)
g_Sij_sum=lapply(g_Sij,sum)

R=1-pgeom(unlist(lapply(g_Sij,sum)),unlist(lapply(p,mean)))


j=c(1038,1346,3528,5008,6224,6559,8587,8990,9001,9101)
plot(j,sort(R),type="l",xlab="Eventual volume",ylab="Probability")

sd = c(0.03878859,0.04529278,0.04607087,0.04733609,0.05313459,0.08066568,
		0.08291579,0.08381236,0.09560712,0.141277)

lcl = j - 1.95996*sd
ucl = j + 1.95996*sd

sd_p = lapply(p,sd)
mean_p = lapply(p,mean)
lcl_p = R - 1.95996*unlist(sd_p)
ucl_p = R + 1.95996*unlist(sd_p)
sort(lcl_p)
sort(ucl_p)

d = data.frame(j,sort(R),rep("Estimate",10))
colnames(d) = c("PosteriorMean","Probability","Parameters")
d1 = data.frame(lcl,sort(lcl_p),rep("LCL",10))
colnames(d1) = c("PosteriorMean","Probability","Parameters")
d2 = data.frame(ucl,sort(ucl_p),rep("UCL",10))
colnames(d2) = c("PosteriorMean","Probability","Parameters")

est = rbind(d,d1,d2)
theme_set(theme_minimal())
 fig1 = ggplot(data = est, aes(x = PosteriorMean, y = Probability, color = Parameters, group = Parameters)) +
  geom_line(size = 1) + xlab("Posterior Mean") + ylab("Probability")
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


sd_new=mapply('*',qnorm(mean(R)),sd_Sij)
Sij_new=mapply('+',mean_Sij,sd_new)
p_new=mapply('/',Sij_new,sum(Sij_new))
cv_new=mapply('/',sd_new,mean_Sij)
prior_new=mapply('*',f_Sij,cv_new)
posterior_new=lapply(prior_new,sum)
g_Sij_new=mapply('/',prior_new,posterior_new)
g_Sij_sum_new=lapply(g_Sij_new,sum)
R_new=1-pgeom(unlist(g_Sij_sum_new),unlist(lapply(p_new,mean)))


theme_set(theme_minimal())
 fig2 = ggplot(data = d, aes(x = j, y = sort(R))) +
  geom_line(size = 1) + xlab("Posterior Mean") + ylab("Probability")
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
