if(!require(pak)){install.packages("pak")}
if(!require(tidyverse)){pak::pkg_install("tidyverse"); library(tidyverse)}
if(!require(MCMCglmm)){pak::pkg_install("MCMCglmm"); library(MCMCglmm)}
if(!require(mvtnorm)){pak::pkg_install("mvtnorm"); library(mvtnorm)}
if(!require(nadiv)){pak::pkg_install("nadiv"); library(nadiv)}
if(!require(tictoc)){pak::pkg_install("tictoc"); library(tictoc)}
if(!require(pedtools)){pak::pkg_install("pedtools"); library(pedtools)}
if(!require(patchwork)){pak::pkg_install("patchwork"); library(patchwork)}
if(!require(cowplot)){pak::pkg_install("cowplot"); library(cowplot)}
if(!require(ggthemes)){pak::pkg_install("ggthemes"); library(ggthemes)}
if(!require(RColorBrewer)){pak::pkg_install("RColorBrewer"); library(RColorBrewer)}
if(!require(corrplot)){pak::pkg_install("corrplot"); library(corrplot)}

# Nice matrix plotting function
corrPlot = function(M, title = ""){
melted_cormat <- reshape2::melt(M)
    heatmap = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
        geom_tile() +
        scale_fill_gradientn(colours=brewer.pal(11, "RdBu"), 
                             limits = c(-1, 1), 
                              breaks=c(-1, -0.5, 0 , 0.5, 1))  + 
        labs(x = "Traits", y = "") +
        theme_tufte() + ggtitle(title) + 
        theme(legend.position = "bottom",
              legend.key.width= unit(1.5, 'cm'), 
              legend.title = element_blank(),
              plot.title = element_text(size = 14),
              axis.title = element_text(size = 12),
              axis.text.y = element_text(size = 6),
              axis.text.x = element_text(size = 6, angle = 90))
    heatmap
}

# Simulate some data

n = 100
ped = randomPed(n, 10, seed = 2)
#plot(ped)
ped = data.frame(ID = ped$ID, sire = ped$FIDX, dam = ped$MIDX)
ped = prepPed(ped)
A <- as.matrix(nadiv::makeA(ped))
Ainv <- nadiv::makeAinv(ped)$Ainv
# library(AtchleyMice)
# ped <- mice_pedigree
# F6_ID = mice_info$F6$ID
# A = as.matrix(nadiv::makeA(ped))
# Acols = match(F6_ID, colnames(A))
# A <- A[Acols, Acols]
# Ainv <- nadiv::makeAinv(ped)$Ainv[Acols, Acols]

L_A = chol(A)
n = dim(A)[1]

set.seed(123)
t <- 2
p <- t*2
sex <- sample(c(0, 1), n, replace = TRUE)
age <- rnorm(n, mean = 50, sd = 10)
x <- rnorm(n, mean = 0, sd = 1)

# 4x4 correlation matrix
corrG <- matrix(c(
  1.0,  0.8,  0.6,  0.3,  
  0.8,  1.0,  0.3,  0.6,  
  0.6,  0.3,  1.0,  0.8,  
  0.3,  0.6,  0.8,  1.0   
), nrow = p, ncol = p) + diag(0.1, p) 
#corrG <- cov2cor(corrG) # Check the correlation matrix
# Ensure the matrix is positive definite
if (!all(eigen(corrG)$values > 0)) {
  stop("The correlation matrix is not positive definite.")
}

a <- t(L_A) %*% matrix(rnorm(n*p), n, p) %*% chol(corrG)
rownames(a) <- 1:n

beta_sex = c(0, 0, 1, 1) 
beta_age = c(1, 1, 1, 1) * 0.2
beta_x = rep(rnorm(t, mean = 0, sd = 1), p/t)

betas <- rbind(beta_age, beta_sex, beta_x)
X <- cbind(age, sex, x)

Ytotal <- X %*% betas + a + rnorm(p, mean = 0, sd = 1)
Y <- Ytotal
# carefull here if you increase the number of traits
Y[sex == 1, 1:2] <- Y[sex == 1, 3:4] 
Y <- Y[, 1:2]
Y

sex %*% matrix(beta_sex, ncol  =4) 


dat = data.frame(id = 1:n, Y, X)
names(dat) <- c("id", "y1", "y2", "age", "sex", "x")
dat$age = dat$age - mean(dat$age)
dat$x = dat$x - mean(dat$x)

lm(Y ~ x + age + sex, data = dat)

# fixed effects model
m1 <- MCMCglmm(cbind(y1, y2) ~ trait + trait:x + trait:age + trait:sex - 1, 
         rcov = ~ us(trait):units,
         data = dat, 
         family = rep("gaussian", t),
         verbose = FALSE)
summary(m1)

# Adding the random effects
m2 <- MCMCglmm(cbind(y1, y2) ~ trait + trait:x + trait:age + trait:sex - 1, 
               random = ~ us(trait):id,  # Random effects for each trait by individual
               rcov = ~ us(trait):units, # Residual covariance structure
               data = dat, 
               family = rep("gaussian", t),
               prior = list(
                 R = list(V = diag(t), nu = t),  # Residual prior
                 G = list(G1 = list(V = diag(t), nu = t))  # Random effect prior
               ),
               verbose = FALSE)
summary(m2)

# Adding the random effects covariance structure
colnames(Ainv) <- rownames(Ainv) <- dat$id
# Update the MCMCglmm model to use the inverted A matrix
m3 <- MCMCglmm(cbind(y1, y2) ~ trait + trait:x + trait:age + trait:sex - 1, 
               random = ~ us(trait):id,  
               ginverse = list(id = Ainv),
               rcov = ~ us(trait):units, 
               data = dat, 
               family = rep("gaussian", t),
               prior = list(
                 R = list(V = diag(t), nu = t),  # Residual prior
                 G = list(G1 = list(V = diag(t), nu = t))  # Random effect prior
               ), 
               verbose = FALSE)
summary(m3) 

# extracting the posterior samples for the random effects
post_samples <- m3$VCV
id_col = grep("id", colnames(post_samples))
matrix(colMeans(post_samples[, id_col]), 2, 2) |> cov2cor()
corrG

# Now lets add the across sex correlations
# This is done by doubling the dimensionality of the G matrix, such that 
# each level of the sex variable has its own G matrix on the block diagonal,
# and the between-sex correlations are estimated outside the block diagonal.

# Because no individual is of both sexes, the residual structure is just the
# block diagonal elements, so there are two separate residual matrices, 
# one per sex. 

# This difference in the genetic and covariance strucutre is reflected in the
# prior specifications: one big G matrix and two small R matrices.
prior = list(R = list(R1 = list(V = diag(t), nu = 1.002),
                      R2 = list(V = diag(t), nu = 1.002)),
             G = list(G1 = list(V = diag(p) * 1.02, nu = p+1)))
dat$sex = factor(dat$sex)
m4 <- MCMCglmm(cbind(y1, y2) ~ trait + trait:x + trait:age + trait:sex - 1, 
               random = ~ us(trait:sex):id,  
               ginverse = list(id = Ainv),
               rcov = ~ us(trait:at.level(sex, "0")):units + 
                        us(trait:at.level(sex, "1")):units, 
               data = dat, 
               family = rep("gaussian", t),
               prior = prior, 
               verbose = FALSE)
summary(m4)
post_samples <- m4$VCV
id_col = grep("id", colnames(post_samples))
cG = matrix(colMeans(post_samples[, id_col]), p, p) |> cov2cor()

corrPlot(cov2cor(corrG),  "True Correlations") +
  corrPlot(cG, "Estimated Correlation")  + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom', plot.title = element_text(family = "serif"))

par(mfrow = c(1, 2))
corrplot.mixed(cov2cor(corrG), main = "True Correlation", upper = "ellipse")
corrplot.mixed(cG, main = "Estimated Correlation", upper = "ellipse")
