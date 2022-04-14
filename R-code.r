setwd("../Supplementary Material")

#The first set of analyses
#the proportion of â€œthreatenedâ€? species between Temperate and Tropics

#Wilcoxon rank sum exact test
SR <- read.csv("p_two_regions.csv")
wilcox.test(SR$Temperate,SR$Tropics,alternative="two.sided")

#data:  SR$Temperate and SR$Tropics
#W = 13, p-value = 4.28e-05
#alternative hypothesis: true location shift is not equal to 0

# Figure 2
library(ggplot2)
data<- read.csv("proportion.csv")
wilcox<-ggplot(data,aes(group,Mean), fill=group)+
      stat_boxplot(geom = 'errorbar',width=0.2,cex=1)+
      geom_boxplot(width=0.8,cex=1,fill='grey',outlier.shape = NA)+xlab("") +
      ylab("The proportion of threatened species") +
      scale_y_continuous(limits = c(0,0.4))+
      scale_fill_manual(values=c("#377EB8","#E41A1C"))+
      theme(axis.text.x = element_text(size = 18, color = "black"))+
      theme(axis.text.y = element_text(size = 15, color = "black"))+
      theme(axis.title.y = element_text(size = 18, color = "black"))+
      geom_text(x=1.5,y=0.4,label="p<0.001", cex=8,colour = "black")
wilcox   
ggsave("Figure 2.pdf",width = 15,height = 15,units = "cm",dpi=600) 


#The second set of analyses
#To analyze the relationship of the proportion of threatened species with the six evolutionary dynamic predictors 
#for all 27 regions being considered as a whole and for each of the regions

#Packages
library(rstan)
library(ape)
library(MASS)
library(brms)

# Code for standardized
scale <- function(x){
  return( (x-mean(x))/sd(x))
}

#test multicollinearity of variables
library(car)
data <- read.csv("0_overall_global_biodiversity_hotspots.csv")

#Subsetting to vetted species
data<- subset(data,n.assess>0)

# Scaling continuous predictors
data$species_richness_scaled <- scale(log(data$Richness))
data$diversification_stem_scaled <- scale(log(data$logStem.DR+1))

with(data, cor(species_richness_scaled, diversification_stem_scaled))#0.863


#Bayesian binomial-logit univariate regression models
###Load sorted data files
data <- read.csv("0_overall_global_biodiversity_hotspots.csv")

## Subsetting to vetted species (for brm models)
data<- subset(data,n.assess>0)

#Import phylogenetic tree
tree<- read.tree("Ramirez_etal_pruned.tre")

# These families we have data for, but are not in tree
setdiff(data$FAMILY, tree$tip.label) 

# Identify families in tree with no data
setdiff(tree$tip.label, data$FAMILY)

# remove any families in tree with no data 
tree_vetted <- drop.tip(tree, setdiff(tree$tip.label, as.character(data$FAMILY)))

# phylogenetic correlation structure
phyloMat <- vcv(tree_vetted, corr=TRUE)

# creating another family level variable (for non-phylogenetic family level effects)
data$family_name <- data$FAMILY

####Hypothesis 1 - Species richness
# Scaling continuous predictor 1
data$species_richness_scaled <- scale(log(data$Richness))

#Model 1 - threat probability as a function of species richness 
Model_1<- brm(n.thrt|trials(n.assess) ~ species_richness_scaled +(1|FAMILY)+(1|family_name), 
           data =data, 
           family = binomial(), 
           iter=10000, thin=4, cov_ranef = list(FAMILY = phyloMat),
           control=list(adapt_delta=0.99,max_treedepth=12), cores=4)

summary(Model_1)
outcome <- summary(Model_1$fit)
outcome <- summary(Model_1$fit, prob=c(0.025,0.05,0.95,0.975))
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
to_print <- data.frame(to_print)
to_print

####Hypothesis 2 - familiy age
# Scaling continuous predictor 2
data$stem_age_scaled <- scale(log(data$StemAge))

#Model 2.1 - threat probability as a function of stem age
Model_2.1<- brm(n.thrt|trials(n.assess) ~ age_scaled +(1|FAMILY)+(1|family_name), 
           data =data, 
           family = binomial(), 
           iter=10000, thin=4, cov_ranef = list(FAMILY = phyloMat),
           control=list(adapt_delta=0.98,max_treedepth=12), cores=4)

summary(Model_2.1)
outcome <- summary(Model_2.1$fit)
outcome <- summary(Model_2.1$fit, prob=c(0.025,0.05,0.95,0.975))
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
to_print <- data.frame(to_print)
to_print

####Hypothesis 2 - familiy age
# Scaling continuous predictor 3
data<- subset(data,CrownAge>0)
data$crown_age_scaled <- scale(log(data$CrownAge))

#Model 2.2 - threat probability as a function of crown age
Model_2.2<- brm(n.thrt|trials(n.assess) ~ age_scaled +(1|FAMILY)+(1|family_name), 
           data =data, 
           family = binomial(), 
           iter=10000, thin=6, cov_ranef = list(FAMILY = phyloMat),
           control=list(adapt_delta=0.99,max_treedepth=12), cores=4)

summary(Model_2.2)
outcome <- summary(Model_2.2$fit)
outcome <- summary(Model_2.2$fit, prob=c(0.025,0.05,0.95,0.975))
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
to_print <- data.frame(to_print)
to_print

#####Hypothesis 3 - Diversification Rate
#Model 3.1
# Scaling continuous predictor 4
#data$logStem.DR <- with(data, log(Richness)/StemAge)
data$diversification_stem_scaled <- scale(log(data$logStem.DR+1))

#Model 3.1 - threat probability as a function of Diversification Rate(based on stem age)
Model_3.1<- brm(n.thrt| trials(n.assess) ~ diversification_stem_scaled +(1|FAMILY)+(1|family_name), 
           data =data, 
           family = binomial(), 
           iter=10000, thin=4, cov_ranef = list(FAMILY = phyloMat),
           control=list(adapt_delta=0.98,max_treedepth=12), cores=4)

summary(Model_3.1)
outcome <- summary(Model_3.1$fit)
outcome <- summary(Model_3.1$fit, prob=c(0.025,0.05,0.95,0.975))
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
to_print <- data.frame(to_print)
to_print

#e=0,e=0.9, stem age
###Load sorted data files
data <- read.csv("0_overall_global_biodiversity_hotspots.csv")

#Calculate rates assuming different extinction fractions (e)
e<-0
r0<-log(data$Richness*(1-e)+e)/data$StemAge

e<-0.9
r0.9<-log(data$Richness*(1-e)+e)/data$StemAge

data<-cbind(r0,r0.9)
data<-as.data.frame(data)
#Correlation test
pairs(data[,c("r0","r0.9")])
with(data, cor(r0, r0.9))

#Model 3.2
# Scaling continuous predictor 5
data$diversification_crown_scaled <- scale(log(data$logCrown.DR+1))

#Model 3.2 - threat probability as a function of Diversification Rate(based on crown age)
Model_3.3<- brm(n.thrt|trials(n.assess) ~ diversification_crown_scaled +(1|FAMILY)+(1|family_name), 
           data =data, 
           family = binomial(), 
           iter=10000, thin=4, cov_ranef = list(FAMILY = phyloMat),
           control=list(adapt_delta=0.98,max_treedepth=12), cores=4)

summary(Model_3.3)
outcome <- summary(Model_3.3$fit)
outcome <- summary(Model_3.3$fit, prob=c(0.025,0.05,0.95,0.975))
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
to_print <- data.frame(to_print)
to_print

#e=0, e=0.9, crown age
###Load sorted data files
data <- read.csv("0_overall_global_biodiversity_hotspots.csv")
data<- subset(data,CrownAge>0)

#Calculate rates assuming different extinction fractions (e)
data$div_simple <- with(data, log(Richness)/CrownAge)

e<-0
r0<-log(data$Richness*(1-e)+e)/data$CrownAge

e<-0.9
r0.9<-log(data$Richness*(1-e)+e)/data$CrownAge

data<-cbind(r0,  r0.9)
data<-as.data.frame(data)
#Correlation test
pairs(data[,c("r0","r0.9")])
with(data, cor(r0, r0.9))


####Hypothesis 4 Lineage Turnover
# Scaling continuous predictor 6
data$Turnover_scaled <- scale(log(data$Turnover))

#Model 4 - threat probability as a function of crown-stem ratio 
Model_4<- brm(n.thrt|trials(n.assess) ~ Turnover_scaled +(1|FAMILY)+(1|family_name), 
           data =data, 
           family = binomial(), 
           iter=10000, thin=4, cov_ranef = list(FAMILY = phyloMat),
           control=list(adapt_delta=0.99,max_treedepth=12), cores=4)
summary(Model_4)
outcome <- summary(Model_4$fit)
outcome <- summary(Model_4$fit, prob=c(0.025,0.05,0.95,0.975))
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
to_print <- data.frame(to_print)
to_print


#Bayesian binomial-logit multivariable regression models
# stem age + diversification rate based on stem age+ Turnover
###Load sorted data files
data <- read.csv("0_overall_global_biodiversity_hotspots.csv")

## Subsetting to vetted species (for brm models)
data<- subset(data,n.assess>0)

#Import phylogenetic tree
tree<- read.tree("Ramirez_etal_pruned.tre")

# These families we have data for, but are not in tree
setdiff(data$FAMILY, tree$tip.label) 

# Identify families in tree with no data
setdiff(tree$tip.label, data$FAMILY)

# remove any families in tree with no data 
tree_vetted <- drop.tip(tree, setdiff(tree$tip.label, as.character(data$FAMILY)))

# phylogenetic correlation structure
phyloMat <- vcv(tree_vetted, corr=TRUE)

# creating another family level variable (for non-phylogenetic family level effects)
data$family_name <- data$FAMILY

####
# Scaling continuous predictors
data<- subset(data,CrownAge>0)
data$stem_age_scaled <- scale(log(data$StemAge))
data$crown_age_scaled <- scale(log(data$CrownAge))
data$diversification_scaled <- scale(log(data$logCrown.DR+1))
data$diversification_scaled <- scale(log(data$logCrown.DR+1))
data$Turnover_scaled <- scale(log(data$Turnover))

#Model mix 1 - threat probability as a function of species richness 
Model_mix_1<- brm(n.thrt|trials(n.assess) ~ stem_age_scaled+ diversification_scaled+ Turnover_scaled+ (1|FAMILY)+ (1|family_name), 
                data =data, 
                family = binomial(), 
                iter=10000, thin=4, cov_ranef = list(FAMILY = phyloMat),
                control=list(adapt_delta=0.99,max_treedepth=12), cores=4)

summary(Model_mix_1)
outcome <- summary(Model_mix_1$fit)
outcome <- summary(Model_mix_1$fit, prob=c(0.025,0.05,0.95,0.975))
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
to_print <- data.frame(to_print)
to_print

# crown age + diversification rate based on crown age+ Turnover
#Model mix 2 - threat probability as a function of species richness 
Model_mix_2<- brm(n.thrt|trials(n.assess) ~ crown_age_scaled+ diversification_scaled+ Turnover_scaled+ (1|FAMILY)+ (1|family_name), 
                data =data, 
                family = binomial(), 
                iter=10000, thin=4, cov_ranef = list(FAMILY = phyloMat),
                control=list(adapt_delta=0.99,max_treedepth=12), cores=4)
summary(Model_mix_2)
outcome <- summary(Model_mix_2$fit)
outcome <- summary(Model_mix_2$fit, prob=c(0.025,0.05,0.95,0.975))
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
to_print <- data.frame(to_print)
to_print



##The third set of analyses- null models
library(dplyr)
library(tidyr)
a<-list.files(path="../Supplementary Material/2_data for brms and null models",pattern=".csv$",full.names =TRUE,include.dirs = TRUE)
a

#Set the matrix of calculation results for each area
final_result<-matrix(ncol=28, nrow=4)
dimnames(final_result)[[1]]<- c('obs.mean','null.mean','null.sd','SES')

#this is a example for Turnover Hypothesis, for other examples, change "Turnover" with Richness, age or diversification rate.   

for (i in 1:28)
{
  data<- read.csv(a[i])
  sum_rarinum<-sum(data$n.thrt)
  data1<-data[rep(1:nrow(data),data$Richness),]
  mat2<-matrix(ncol=1000,nrow=(dim(data)[1]))
  mat<-matrix(ncol=1,nrow=(dim(data)[1]))
  mat<-as.data.frame(mat)
  mat$FAMILY<-data$FAMILY
  for (j in 1:1000)
  {
    set.seed(j)
    rarinum<- sample_n(data1,sum_rarinum,replace=FALSE)
    b<-unique(rarinum$FAMILY)
    b<-as.data.frame(b)
    b$checked<-1
    dimnames(b)[[2]]<-c('FAMILY','checked')
    z<-merge(data,b,by='FAMILY',all.x =TRUE,sort=FALSE)
    mat1<-merge(mat,z,by='FAMILY',all.x=TRUE,sort=FALSE)
    mat2[,j]<-mat1$checked}
  mat3<-as.data.frame(mat2)
  mat3[is.na(mat3)] <- 0
  result<- data.frame(data, mat3)
  result <- result %>% drop_na(Turnover)
  result2 <-data.frame(sample=dimnames(result)[[2]][11:dim(result)[2]])
  result2$sum<-c(apply(result[,11:dim(result)[2]],2,sum))
  result2$mean<- as.numeric(apply(matrix(result$Turnover,nrow=1)* result[,11:dim(result)[2]], 2, sum)/result2$sum)
  null.mean<-mean(result2$mean)
  null.sd<-sd(result2$mean)
  data.obs<-as.data.frame(subset(data,data$n.thrt>0))
  data.obs <- data.obs %>% drop_na(Turnover)
  obs.mean<-mean(data.obs$Turnover)
  SES<-(obs.mean-null.mean)/null.sd
  final_result[,i]<-c(obs.mean,null.mean,null.sd,SES)
}

write.csv(final_result,"results.csv")
