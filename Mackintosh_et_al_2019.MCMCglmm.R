#MCMCGLMM PHYLOGENETIC MODEL

rm(list=ls())

#load packages

library(ape)
library(MCMCglmm)
library(phytools)

#import data

dat<-read.table(file.choose(),header=T)

#check data 

names(dat)

#transform predictors by z transformation
#this allows effect sizes to be compared (mean is 0, SD is 1)

dat<-cbind(dat, scale(dat$size))

names(dat)[names(dat)=="scale(dat$size)"]<-"scalesize"

#import species tree in newick format without node labels 

phylo<-read.newick(file.choose())

#reroot tree by 'mori' then drop this outgroup and other taxa

phylo<-reroot(phylo, interactive=TRUE)

mori<-cbind("mori")

phylo<-drop.tip(phylo,mori)

plot(phylo)

#calculate the inverse matrix of branch lengths 

inv.phylo<-inverseA(phylo,nodes='ALL', scale=FALSE)

#set prior for random effect of phylogeny and units

prior1<-list(G=list(G1=list(V=diag(2),nu=2, alpha.mu=c(0,0), 
alpha.V=diag(2)*1000)),R=list(V=diag(2),nu=0.002))

#run a bivariate MCMCglmm where each "pi" has been log transformed 
#error family is gaussian as "log_pi" is normally distributed

model_1<-MCMCglmm(cbind(log_4, log_0)~trait-1+trait:scalesize+trait:scalechromosome, rcov=~us(trait):units,random=~us(trait):phylo,
family=cbind("gaussian", "gaussian"),ginverse=list(phylo=inv.phylo$Ainv),prior=prior1,
data=dat,nitt=13000*10,burnin=3000*10,thin=10*5)

#Model summary

summary(model_1)

#Plot sampling of posterior, adjust "thin" if there is autocorrelation 

plot(model_1)


