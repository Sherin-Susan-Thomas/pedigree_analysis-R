setwd("C:/R")

group6_data <- read.csv("C:/R/group6_data.csv", header=FALSE)
View(group6_data)
Data <- as.data.frame(read.table(file = "group6_data.csv",header = TRUE, sep = ","))
names(Data)[1] <- "animal"
names(Data)[2] <- "age"
names(Data)[3] <- "sex"
names(Data)[4] <- "pca"
Data$animal <- as.factor(Data$animal)
Data$age <- as.factor(Data$age)
Data$sex <- as.factor(Data$sex)
Data$pca <- as.numeric(Data$pca)
head(Data)
Data

beagle_pedigree <- read.csv("C:/R/beagle_pedigree.csv")
View(beagle_pedigree)
Ped <- as.data.frame(read.table(file = "beagle_pedigree.csv", header = TRUE, sep = ","))
head (Ped)
library(MCMCglmm)
prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

model1.1 <- MCMCglmm(pca ~ 1, random = ~animal, pedigree = Ped,  data = Data, prior = prior1.1)
plot(model1.1$Sol)
plot(model1.1$VCV)
model1.1 <- MCMCglmm(pca ~ 1, random = ~animal, pedigree = Ped,data = Data, nitt = 65000, thin = 50, burnin = 15000, prior = prior1.1)

autocorr(model1.1$VCV)
posterior.mode(model1.1$VCV)
HPDinterval(model1.1$VCV)

#Estimating heritability
posterior.heritability1.1 <- model1.1$VCV[, "animal"]/(model1.1$VCV[, "animal"] + model1.1$VCV[, "units"])
HPDinterval(posterior.heritability1.1, 0.95)
posterior.mode(posterior.heritability1.1)
plot(posterior.heritability1.1)


#Adding fixed effects sex
model1.2 <- MCMCglmm(pca ~ sex, random = ~animal, pedigree = Ped,data = Data, prior = prior1.1, nitt = 65000, thin = 50, burnin = 15000, verbose = FALSE)
posterior.mode(model1.2$Sol[, "sex2"])

HPDinterval(model1.2$Sol[, "sex2"], 0.95)

posterior.mode(model1.2$VCV)
posterior.heritability1.2 <- model1.2$VCV[, "animal"]/(model1.2$VCV[,"animal"] + model1.2$VCV[, "units"])
posterior.mode(posterior.heritability1.2)
HPDinterval(posterior.heritability1.2, 0.95)


#adding random effects-age + fixed factor-sex

prior1.3 <- list(G = list(G1 = list(V = 1, n = 0.002), G2 = list(V = 1,n = 0.002)), R = list(V = 1, n = 0.002))

model1.3 <- MCMCglmm(pca ~ sex, random = ~animal + age, pedigree = Ped, data = Data, nitt = 65000, thin = 50, burnin = 15000, prior = prior1.3,verbose = FALSE)

posterior.mode(model1.3$VCV)
posterior.heritability1.3 <- model1.3$VCV[, "animal"]/(model1.3$VCV[,"animal"] + model1.3$VCV[, "age"] + model1.3$VCV[, "units"])


posterior.mode(posterior.heritability1.2)
posterior.mode(posterior.heritability1.3)
# adding age as random factor
prior1.4 <- list(G = list(G1 = list(V = 1, n = 0.002), G2 = list(V = 1,n = 0.002)), R = list(V = 1, n = 0.002))
model1.4 <- MCMCglmm(pca ~ 1, random = ~animal + age, pedigree = Ped, data = Data, nitt = 65000, thin = 50, burnin = 15000, prior = prior1.4,verbose = FALSE)

posterior.mode(model1.4$VCV)
posterior.heritability1.4 <- model1.4$VCV[, "animal"]/(model1.4$VCV[,"animal"] + model1.4$VCV[, "age"] + model1.3$VCV[, "units"])
posterior.mode(posterior.heritability1.4)
#Testing significance of variance components
model1.2$DIC
model1.3$DIC
model1.1$DIC
model1.4$DIC

