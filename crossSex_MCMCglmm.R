library(MCMCglmm)
library(tidyverse)
library(AtchleyMice)
library(sommer)
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}
if(!require(nadiv)){install.packages("nadiv"); library(nadiv)}
if(!require(tictoc)){install.packages("tictoc"); library(tictoc)}
if(!require(pedtools)){install.packages("pedtools"); library(pedtools)}

# Simulate some data

set.seed(123)
n <- 100
p <- 4
sex <- sample(c(0, 1), n, replace = TRUE)
age <- rnorm(n, mean = 50, sd = 10)
x <- rnorm(n, mean = 0, sd = 1)
a <- rnorm(n, mean = 0, sd = 1)

ped = randomPed(n, 10, seed = 2)
plot(ped)
ped = data.frame(ID = ped$ID, sire = ped$FIDX, dam = ped$MIDX)
ped = prepPed(ped)
A <- as.matrix(nadiv::makeA(ped))
#ped <- mice_pedigree
#F6_ID = mice_info$F6$ID
#A <- as.matrix(nadiv::makeA(ped))[F6_ID, F6_ID]
L_A = chol(A)

# 4x4 correlation matrix
corrG <- matrix(c(
  1.0,  0.8,  0.6,  0.3,  
  0.8,  1.0,  0.3,  0.6,  
  0.6,  0.3,  1.0,  0.8,  
  0.3,  0.6,  0.8,  1.0   
), nrow = 4, ncol = 4) + diag(0.1, 4) 
#corrG <- cov2cor(corrG) # Check the correlation matrix
# Ensure the matrix is positive definite
if (!all(eigen(corrG)$values > 0)) {
  stop("The correlation matrix is not positive definite.")
}

a <- t(L_A) %*% matrix(rnorm(n*p), n, p) %*% chol(corrG)

beta_sex = c(0, 0, 1, 1) * 1/2
beta_age = c(1, 1, 1, 1) * 0.2
beta_x = rep(rnorm(2, mean = 0, sd = 1), 2)

betas <- rbind(beta_age, beta_sex, beta_x)
X <- cbind(age, sex, x)

Ytotal <- X %*% betas + a + rnorm(4, mean = 0, sd = 1)
Y <- Ytotal
Y[sex == 1, 1:2] <- Y[sex == 1, 3:4]
Y <- Y[, 1:2]
Y

lm(Y ~ x + age + sex)
dat = data.frame(id = 1:n, Y, X)
names(dat)[2:3] <- c("y1", "y2")

# fixed effects model
m1 <- MCMCglmm(cbind(y1, y2) ~ trait + trait:x + trait:age + trait:sex, 
         rcov = ~ us(trait):units,
         data = dat, 
         family = rep("gaussian", 2))
summary(m1)

# Adding the random effects
m2 <- MCMCglmm(cbind(y1, y2) ~ trait + trait:x + trait:age + trait:sex, 
               random = ~ us(trait):id,  # Random effects for each trait by individual
               rcov = ~ us(trait):units, # Residual covariance structure
               data = dat, 
               family = rep("gaussian", 2),
               prior = list(
                 R = list(V = diag(2), nu = 2),  # Residual prior
                 G = list(G1 = list(V = diag(2), nu = 2))  # Random effect prior
               ),
               verbose = FALSE)
summary(m2)

# Adding the random effects covariance structure
# Invert the A matrix and convert it to sparse format
Ainv <- nadiv::makeAinv(ped)$Ainv

# Update the MCMCglmm model to use the inverted A matrix
m3 <- MCMCglmm(cbind(y1, y2) ~ trait + trait:x + trait:age + trait:sex, 
               random = ~ us(trait):id,  
               ginverse = list(id = Ainv),
               rcov = ~ us(trait):units, 
               data = dat, 
               family = rep("gaussian", 2),
               prior = list(
                 R = list(V = diag(2), nu = 2),  # Residual prior
                 G = list(G1 = list(V = diag(2), nu = 2))  # Random effect prior
               ), 
               verbose = FALSE)
summary(m3)

# extracting the posterior samples for the random effects
post_samples <- m3$VCV
id_col = grep("id", colnames(post_samples))
matrix(colMeans(post_samples[, id_col]), 2, 2) |> cov2cor()
corrG
