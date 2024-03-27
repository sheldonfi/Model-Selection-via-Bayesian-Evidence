## working with a Gamma - Poisson conjugate model
## We set our parameters lambda, a and b for Poisson(lambda) and Gamma(a,b)

########################################################
### Data simulations
set.seed(1)

sets <- 50   # the number of sets of data
n <- 50 # number of data samples
lambda <- 2 # lambda value used for generating data

y <- matrix(0,sets,n) # y is a matrix to hold sets of data 

for (i in 1:sets) 
{ y[i,] <- rpois(n,lambda)
  }

ybar <- apply(y,1,mean) # ybar is a vector which holds the means of each data set

#######################################################
### Model Specifications

sets=50   # the number of sets of data
n = 50    # number of data samples in each set
N=10000  # size of posterior and prior samples
a = 0.1 # Gamma prior shape parameter
b = 0.1 # Gamma prior scale parameter

########################################################
#### FUNCTIONS
### We create the functions here that we need for our methods

logmarglikelihood <- function(data,ybar,n,a,b) {a*log(b)+log(gamma(n*ybar+a))-sum(log(factorial(data)))-log(gamma(a))-(n*ybar+a)*log(n+b)}

likelihood <- function(x,data) {prod(dpois(data,x, log = FALSE))}
loglikelihood <- function(x,data) {sum(dpois(data,x, log = TRUE))}

prior <- function(x) {dgamma(x,a,b, log = FALSE)}
logprior <- function(x) {dgamma(x,a,b, log = TRUE)}

########################################################
#### True model evidence values 

results <- numeric(sets) # results is a vector which holds our true model evidence values
for (i in 1:sets) { results[i] <- logmarglikelihood(y[i,],ybar[i],n,a,b) }

########################################################
#### Sampling from the posterior 
set.seed(1)

postsims <- t(apply(y,1,function(x) {rgamma(N, a + sum(x), b+n)})) # postsims is a matrix which holds the posterior simulations

########################################################
#### Sampling from the prior
set.seed(1)

priorsims <- matrix(rgamma(sets*N,a,b),sets,N) # priorsims is a matrix which holds the prior simulations

#priorsims[priorsims==0] <- .Machine$double.xmin

########################################################
### Prior vs Posterior density plots
library(ggplot2)

x.plot <- seq(0,5,0.001)
y.prior <- dgamma(x.plot,a,b)
y.post <- dgamma(x.plot,a+sum(y[1,]), b + n)

data <- data.frame(x = x.plot, y = y.prior)
ggplot(data = data, aes(x = x.plot, y = y.prior)) +
  geom_line(color="#69b3a2", size = 1) +
  labs(x = "Values", y = "Density") +
  scale_y_continuous(limits = c(0, 2.1))

data <- data.frame(x = x.plot, y = y.post)
ggplot(data = data, aes(x = x.plot, y = y.post)) +
  geom_line(color="#69b3a2", size = 1) +
  labs(x = "Values", y = "Density") +
  scale_y_continuous(limits = c(0, 2.1))

data <- data.frame(x = x.plot, y1 = y.prior, y2 = y.post)
ggplot(data = data, aes(x = x.plot)) +
  geom_line(aes(y=y1,colour="Prior"), size = 1) +
  geom_line(aes(y=y2,colour="Posterior"), size = 1) + 
  scale_colour_manual("", breaks = c("Prior","Posterior"),
                      values = c("#69b3a2", "red")) +
  labs(x = "Values", y = "Density") +
  scale_y_continuous(limits = c(0, 2.1)) +
  scale_x_continuous(limits = c(0, 4))

############################################
#### METHODS
############################################
### MONTE CARLO ESTIMATOR
set.seed(1)

MClikelihoods <- matrix(0,sets,N)

for(i in 1:sets){
MClikelihoods[i,] <- sapply(priorsims[i,],likelihood,data = y[i,])
}

MCestimates <- apply(MClikelihoods,1, function(x) {log(1/N*sum(x))})

MCerrors <- abs(MCestimates-results)
MCmean <- mean(MCerrors)
MCsd <- sqrt(var(MCerrors))
boxplot(MCerrors)

############################################
#### HARMONIC MEAN ESTIMATOR
set.seed(1)

HMlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  HMlikelihoods[i,] <- sapply(postsims[i,],likelihood,data = y[i,])
}

HMestimates <- apply(HMlikelihoods,1,function(x) {-log(1/N*sum(1/x))})

HMerrors <- abs(HMestimates-results)
HMmean <- mean(HMerrors)
HMsd <- sqrt(var(HMerrors))
boxplot(HMerrors)

############################################
#### GENERALISED HARMONIC MEAN ESTIMATOR
set.seed(1)
GHMtiming <- system.time({
GHMproposaldensities <- matrix(0,sets,N)
for (i in 1:sets){
  mean = mean(postsims[i,])
  var = var(postsims[i,])
  GHMproposaldensities[i,] <- sapply(postsims[i,],dnorm,mean = mean, sd = sqrt(var))
}

GHMpriors <- t(apply(postsims,1,dgamma, shape = a, rate = b))

matrix <- GHMproposaldensities/(HMlikelihoods*GHMpriors)

GHMestimates <- apply(matrix,1,function(x) {-log(1/N*sum(x))})
})
print(GHMtiming)

GHMerrors <- abs(GHMestimates-results)
GHMmean <- mean(GHMerrors)
GHMsd <- sqrt(var(GHMerrors))
boxplot(GHMerrors)

############################################################################################################################
## NEWTON-RAFTERY ESTIMATOR
set.seed(1)

delta = 0.2
NRiterations = 1000

indices <- matrix(0,sets,delta*N)
for (i in 1:sets) {
  indices[i,] <- sample(N,delta*N)
}

NRsims <- postsims

for (i in 1:sets){
  for(j in indices[i,]){
    NRsims[i,j] <- priorsims[i,j]
  }}


NRlikelihoods <- matrix(0,sets,N)
for (i in 1:sets){
  NRlikelihoods[i,] <- sapply(NRsims[i,],likelihood,data = y[i,])
}

NRestimates <- rep(1,50)

sequence <- matrix(0,sets,NRiterations)
for(i in 1:sets){
  sequence[i,1] <- NRestimates[i]
}

for (i in 1:sets){
  for (j in 1:NRiterations) {
    NRestimates[i] <- (sum(NRlikelihoods[i,]*(delta*NRestimates[i]+(1-delta)*NRlikelihoods[i,])^-1))/(sum((delta*NRestimates[i]+(1-delta)*NRlikelihoods[i,])^-1))
  sequence[i,j] <- NRestimates[i]
    }
}

NRestimates <- log(NRestimates)

NRerrors <- abs(NRestimates-results)
NRmean <- mean(NRerrors)
NRsd <- sqrt(var(NRerrors))
boxplot(NRerrors)

par(mar = c(4, 2.5, 2, 1.5))

plot(10:NRiterations,sequence[1,10:NRiterations], xlab = 'iterations', ylab = '')
############################################
#### LAPLACE-METROPOLIS ESTIMATOR
set.seed(1)

logposterior <- matrix(0,sets,N)
for ( i in 1:sets) {
  for( j in 1:N) {
  logposterior[i,j] <- loglikelihood(postsims[i,j],y[i,])+logprior(postsims[i,j])
}}

maxlambda <- numeric(sets)
for ( i in 1:sets ) {maxlambda[i] <-  which.max(logposterior[i,])}

S <- numeric(sets)
for (i in 1:sets) {S[i] <- var(postsims[i,])}

d <- 1

LMestimates <- numeric(sets)
for (i in 1:sets) { LMestimates[i] <- (d/2)*log(2*pi)+0.5*log(S[i])+logposterior[i,][maxlambda[i]]}

LMestimates

LMerrors <- abs(LMestimates-results)
LMmean <- mean(LMerrors)
LMsd <- sqrt(var(LMerrors))

boxplot(LMerrors)

############################################
#### BRIDGE SAMPLING ESTIMATOR
set.seed(1)
BStiming <- system.time({
BSiterations=1

N1=N
N2=N

s1=N1/(N1+N2)
s2=N2/(N1+N2)

BSproposalsims <- matrix(0,sets,N2)
for (i in 1:sets){
  mean = mean(postsims[i,])
  var = var(postsims[i,])
  BSproposalsims[i,] <- rnorm(N2,mean = mean, sd = sqrt(var))
}

likelihoods2 <- matrix(0,sets,N2)

for (i in 1:sets){
  likelihoods2[i,] <- sapply(BSproposalsims[i,],likelihood,data = y[i,])
}

priors2 <- t(apply(BSproposalsims,1,dgamma, shape = a, rate = b))

proposaldensities2 <- matrix(0,sets,N2)
for (i in 1:sets){
  mean = mean(postsims[i,])
  var = var(postsims[i,])
  proposaldensities2[i,] <- sapply(BSproposalsims[i,],dnorm,mean = mean, sd = sqrt(var))
}

likelihoods1 <- HMlikelihoods

priors1 <- GHMpriors

proposaldensities1 <- GHMproposaldensities

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(1,N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(1,N1)
}

sequence <- matrix(0,BSiterations,sets)
for(i in 1:BSiterations){
matrix2 <- (likelihoods2*priors2)/(s1*likelihoods2*priors2+s2*x2*proposaldensities2)
matrix1 <- proposaldensities1/(s1*likelihoods1*priors1+s2*x1*proposaldensities1)
numerator <- apply(matrix2,1,function(x) {1/N2*sum(x)})
denominator <- apply(matrix1,1,function(x) {1/N1*sum(x)})
BSestimates <- numerator/denominator

sequence[i,] <- BSestimates

x2 <- matrix(0,sets,N2)
for (i in 1:sets){
  x2[i,] <- rep(exp(BSestimates[i]),N2)
}

x1 <- matrix(0,sets,N1)
for (i in 1:sets){
  x1[i,] <- rep(exp(BSestimates[i]),N1)
}

}
})
print(BStiming)
BSestimates <- log(BSestimates)
BSerrors <- abs(BSestimates-results)
BSmean <- mean(BSerrors)
BSsd <- sqrt(var(BSerrors))
boxplot(BSerrors)

plot(1:BSiterations,sequence[1:BSiterations,1], xlab = 'iterations', ylab = '')
############################################
#### FOURIER INTEGRAL ESTIMATOR
set.seed(1)
FItiming <- system.time({
N <- N
R <- 20

postmeans <- apply(postsims,1,mean)

FIestimates <- numeric(sets)
postestimates <- numeric(sets)

for(i in 1:sets){
  postestimates[i] <- 1/(N*pi)*sum(sin(R*(postmeans[i]-postsims[i,]))/(postmeans[i]-postsims[i,]))
}

for( i in 1:sets){
FIestimates[i] <- loglikelihood(postmeans[i],y[i,])+logprior(postmeans[i])-log(postestimates[i])
}
})
print(FItiming)

FIerrors <- abs(FIestimates-results)
boxplot(abs(FIestimates-results))
FImean <- mean(FIerrors)
FIsd <- sqrt(var(FIerrors))

FImatrix <- matrix(0,100,50)
for (R in 1:100){
  FIestimates <- numeric(sets)
  postestimates <- numeric(sets)
  
  for(i in 1:sets){
    postestimates[i] <- 1/(N*pi)*sum(sin(R*(postmeans[i]-postsims[i,]))/(postmeans[i]-postsims[i,]))
  }
  
  for( i in 1:sets){
    FIestimates[i] <- loglikelihood(postmeans[i],y[i,])+logprior(postmeans[i])-log(postestimates[i])
  }
  FImatrix[R,] <- FIestimates
}

x.plot <- 1:100
FIdata <- data.frame(x = x.plot, y1 = FImatrix[,1], y2 = FImatrix[,2],y3 = FImatrix[,3], y4 = FImatrix[,4], y5 = FImatrix[,5])
ggplot(data = FIdata, aes(x = x)) +
  geom_line(aes(y=y1,colour="estimate 1"), size = 1) +
  geom_line(aes(y=y2,colour="estimate 2"), size = 1) + 
  geom_line(aes(y=y3,colour="estimate 3"), size = 1) +
  geom_line(aes(y=y4,colour="estimate 4"), size = 1) +
  geom_line(aes(y=y5,colour="estimate 5"), size = 1) +
  scale_colour_manual("", breaks = c("estimate 1","estimate 2","estimate 3","estimate 4", "estimate 5"),
                      values = c("#69b3a2", "lightblue","pink","orange","red")) +
  labs(x = "R", y = "density values") 

############################################################################################################################
## box plot of errors
par(mar = c(2.2,4.3,2,0.4))

boxplot(MCerrors,HMerrors,GHMerrors,NRerrors,LMerrors,BSerrors,FIerrors, xlab = '', ylab = 'Absolute Error', names = c('Monte Carlo','H. Mean', 'G. H. Mean', 'N-Raftery','L-Metropolis','B. Sampling', 'F. Integral'), outine.col = 'red', pch = 20 , 
        boxfill = "#69b3a2", whiskcol = 'black', staplecol = 'black', border = 'black', cex.lab = 1.4, cex.main = 1.3)

## box plot of errors for best estimators
boxplot(GHMerrors,LMerrors,BSerrors,FIerrors, xlab = '', ylab = 'Absolute Error', names = c('Gen. Harm. Mean','Laplace-Metropolis','Bridge S.', 'Fourier I.'), outine.col = 'red', pch = 20 , 
        boxfill = "#69b3a2", whiskcol = 'black', staplecol = 'black', border = 'black', cex.lab = 1.4, cex.main = 1.3)

############################################################################################################################
## Summaries

df <- data.frame(
  'Methods' = c("Mean", "Std. Dev."),
  'MC' = c(MCmean,MCsd),
  'HM' = c(HMmean,HMsd),
  'GHM' = c(GHMmean,GHMsd),
  'NR' = c(NRmean,NRsd),
  'LM' = c(LMmean,LMsd),
  'BS' = c(BSmean,BSsd),
  'FI' = c(FImean,FIsd)
)

knitr::kable(df)
