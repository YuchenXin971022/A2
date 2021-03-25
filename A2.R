library(pastecs)
library(tidyverse)
library(data.table)
library(stringr)
library(mfx)
library(MASS)
library(bayesm)
library(rmarkdown)


## Exercise 1 Data Desription

library(bayesm)
data(margarine)
choiceprice=margarine$choicePrice
demos=margarine$demos
all=merge(choiceprice,demos,by="hhid")

avg=apply(as.matrix(choiceprice[,3:12]),2, mean)
print(avg)

std<-function(x) sd(x)/sqrt(length(x))
stdc=apply(as.matrix(choiceprice[,3:12]),2, std)
print(stdc)

#Average and dispersion in product characteristics
library(tidyr)
library(dplyr)
#Compute all summary statistics including both chosen products and non-chosen products
means <- t(all %>% summarise_at(3:12,mean))
mins <- t(all %>% summarise_at(3:12,min))
maxs <- t(all %>% summarise_at(3:12,max))
sds <- t(all %>% summarise_at(3:12,sd))
vars <- t(all %>% summarise_at(3:12,var))
des1=cbind(means,mins,maxs,sds,vars)
label1=c("mean","min","max","sd","var")
colnames(des1) <- label1
des1=round(des1,digits=3)
des1

#Market Share :

ms<-as.data.frame(choiceprice %>% count(choice))
ms$n<-ms$n/nrow(all)
ms

#Market Share by Characteristics:

#finding the price 
for (i in 1:nrow(all)){
  all$choiceprice[i]=all[i,all$choice[i]+2]
  all$choices[i]=colnames(all)[all$choice[i]+2]
} 

means <- as.data.frame((t(all %>% summarise_at(3:12,mean))))
means <- cbind(choices = rownames(means), means)
rownames(means) <- 1:nrow(means)
means
all=left_join(all, unique(means), by = c("choices"="choices"))

all$indict<-all$choiceprice>all$V1
sum(all$indict)

all%>%
  filter(all$indict==TRUE) %>%
  count(choice)%>%
  mutate(n=n/sum(all$indict))

all%>%
  filter(all$indict==TRUE) %>%
  count(choice)%>%
  mutate(n=n/(nrow(all)-sum(all$indict)))

#consider mapping on choices and attributes 

des3= all %>%
  group_by(choice) %>%
  summarize(
    famsize1_2=sum(Fs3_4 == 0 & Fs5.==0),
    famsize3_4=sum(Fs3_4 == 1 & Fs5.==0),
    famsize5.=sum(Fs3_4 == 0 & Fs5.==1),
    college=sum(college==1),
    whtcollar=sum(whtcollar==1),
    retired=sum(retired==1)
  )
notdes3= all %>%
  group_by(choice) %>%
  summarize(
    notcollege=sum(college==0),
    notwhtcollar=sum(whtcollar==0),
    notretired=sum(retired==0)
  )
des3=merge(des3,notdes3)
des3


## Exercise 2 First Model
#Likelihood and Optimization: 

library(mlogit)
library(stargazer)
library(texreg)
library(survival)
library(nnet)
price=margarine$choicePrice
colnames(price)[3:12]=paste0("price",1:10)
price<-price[-c(1)]
clogit0=mlogit.data(price,varying=2:11,shape="wide",sep="",choice="choice")
clogit1=mlogit(choice~price,data=clogit0)
summary(clogit1)

ni=nrow(price)
nj=ncol(price[,2:11])
Y=class.ind(all$choice)

clogit_ll<-function(beta){
  
  #Create the constant as instructed
  intercept=cbind(0,matrix(rep(beta[1:nj-1],each=ni),ni,nj-1))
  
  #Use the lecture definition of conditional logit to compute the likelihood
  XB=price[,2:11]*beta[nj]
  XB=intercept+XB
  eXB=exp(XB)
  teXB=rowSums(eXB)
  prob=eXB/teXB
  
  #Compute the neg log likelihood for each choice using the choice matrix
  llik=sum(Y*log(prob))
  return(-llik)
}

set.seed(0)
clogit <- optim(runif(10,-0.1,0.1),clogit_ll,method="BFGS")



clogit$par

## Exercise 3 Second Model

colnames(all)[3:12]=paste0("price",1:10)
mlogit0=mlogit.data(all,varying=3:12,shape="wide",sep="",choice="choice")
mlogit1=mlogit(choice~0 | Income,data=mlogit0)
summary(mlogit1)

mlogit_ll<-function(beta){
  
  #Create the constant as instructed
  intercept1=cbind(0,matrix(rep(beta[1:9],each=ni),ni,nj-1))
  intercept2=cbind(0,matrix(rep(beta[10:18],each=ni),ni,nj-1))
  #Use the lecture definition of conditional logit to compute the likelihood
  XB=intercept1+intercept2*cbind(all[,13],all[,13],all[,13],all[,13],all[,13],all[,13],all[,13],all[,13],all[,13],all[,13])
  eXB=exp(XB)
  teXB=rowSums(eXB)
  prob=eXB/teXB
  
  #Compute the neg log likelihood for each choice using the choice matrix
  ll=-sum(Y*log(prob))
  return(ll)
}

set.seed(0)
model2 <- optim(runif(18,-0.1,0.1),mlogit_ll,method="BFGS")
model2$par[10:18]



beta=matrix(0,nrow=1,ncol=11)
beta[1,1]=0
beta[1,2:11]=clogit$par
delta<-diag(1,10,10)

mf=matrix(0,nrow=10,ncol=10)
for (i in 1:nrow(all)){
  intercept=beta[1:10]
  
  #Use the lecture definition of conditional logit to compute the likelihood
  xibi=as.matrix(exp(price[i,2:11]*beta[11]+intercept))
  xib=rowSums(xibi)
  pij=xibi/xib
  dpij<-matrix(NA,10,10)
  delta<-diag(1,10,10)
  for (k in 1:10){
    delta[,k]=delta[,k]-pij
  }
  for (m in 1:10){
    dpij[m,]=pij*delta[m,]*beta[11]
  }
  mf=mf+dpij
}

mf/4470

z<-with(clogit0,data.frame(price=tapply(price,idx(clogit0,2),mean)))
effects(clogit1,covariate="price",data=z)


beta=matrix(0,nrow=1,ncol=20)
beta[1,1]=0
beta[1,2:10]=model2$par[1:9]
beta[1,11]=0
beta[1,12:20]=model2$par[10:18]


#Marginal Effects of MNL:

intercept1=beta[1:10]
intercept2=beta[11:20]
mf=matrix(0,nrow=1,ncol=10)
for (i in 1:nrow(all)){
  xibi=exp(intercept1+intercept2*cbind(all[i,13],all[i,13],all[i,13],all[i,13],all[i,13],all[i,13],all[i,13],all[i,13],all[i,13],all[i,13]))
  xib=rowSums(xibi)
  pij=xibi/xib
  betai=rowSums(pij*intercept2)
  dpixi=pij*(intercept2-betai)
  mf=mf+dpixi
}

#check with packages

z<-with(mlogit0,data.frame(Income=tapply(Income,idx(clogit0,2),mean)))
effects(mlogit1,covariate="Income",data=z)

mf/4470


## Exercise 5 Mixed Logit

mixed_ll<-function(para){
  beta1=matrix(0,nrow=1,ncol=11)
  gamma1=matrix(0,nrow=1,ncol=20)
  beta1[1,1]=0
  beta1[1,2:10]=para[1:9]
  beta1[1,11]=para[10]
  gamma1[1,1]=0
  gamma1[1,2:10]=para[1:9]
  gamma1[1,11]=0
  gamma1[1,12:20]=para[11:19]
  intercept=cbind(0,matrix(rep(beta1[2:10],each=ni),ni,nj-1))/2
  intercept1=cbind(0,matrix(rep(gamma1[2:10],each=ni),ni,nj-1))/2
  intercept2=cbind(0,matrix(rep(gamma1[12:20],each=ni),ni,nj-1))
  XB=intercept1+intercept2*cbind(all[,13],all[,13],all[,13],all[,13],all[,13],all[,13],all[,13],all[,13],all[,13],all[,13])
  WG=price[,2:11]*beta1[11]+intercept
  e=exp(XB+WG)
  se=rowSums(e)
  prob=e/se
  
  #Compute the neg log likelihood for each choice using the choice matrix
  ll=sum(Y*log(prob))
  return(-ll)
}

set.seed(0)
model3 <-  optim(runif(19,-0.1,0.1),mixed_ll,method="BFGS", hessian = FALSE)

model3$par

#consider alternative specification:

price[,12]=all$Income
colnames(price)[12] <- "Income"
logit0=mlogit.data(price,varying=2:11,shape="wide",sep="",choice="choice")
logit1=mlogit(choice~price|Income,data=logit0)
summary(logit1)

#removing one choice, removing choice10 data
pricer<-price %>%
  filter(price$choice!=1)
pricer<-pricer[-c(2)]

all1<-all %>%
  filter(all$choice!=1)
ni=nrow(pricer)
nj=ncol(pricer[,2:10])
Y=class.ind(all1$choice)

#computing ll
mixed_ll2<-function(para){
  beta1=matrix(0,nrow=1,ncol=10)
  gamma1=matrix(0,nrow=1,ncol=18)
  beta1[1,1]=0
  beta1[1,2:9]=para[1:8]
  beta1[1,10]=para[9]
  gamma1[1,1]=0
  gamma1[1,2:9]=para[1:8]
  gamma1[1,10]=0
  gamma1[1,11:18]=para[10:17]
  intercept=cbind(0,matrix(rep(beta1[2:9],each=ni),ni,nj-1))/2
  intercept1=cbind(0,matrix(rep(gamma1[2:9],each=ni),ni,nj-1))/2
  intercept2=cbind(0,matrix(rep(gamma1[11:18],each=ni),ni,nj-1))
  XB=intercept1+intercept2*cbind(all1[,13],all1[,13],all1[,13],all1[,13],all1[,13],all1[,13],all1[,13],all1[,13],all1[,13])
  WG=all1[,4:12]*para[9]+intercept
  e=exp(XB+WG)
  se=rowSums(e)
  prob=e/se
  ll=sum(Y*log(prob))
  return(-ll)
}

set.seed(0)
model4 <-optim(runif(17,-0.1,0.1),mixed_ll2,method="BFGS", hessian = FALSE)

model4$par

logit01=mlogit.data(pricer,varying=2:10,shape="wide",sep="",choice="choice")
logit11=mlogit(choice~price|Income,data=logit01)
summary(logit11)

ni=nrow(price)
nj=ncol(price[,2:11])
Y=matrix(0,ni,nj)
for (i in 1:nj){
  for (j in 2:ni){
    if (price$choice[j]==i){
      Y[j,i]=1
    }
  }
}
mixed_ll(model3$par)