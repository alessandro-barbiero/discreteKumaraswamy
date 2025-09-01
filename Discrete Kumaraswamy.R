#library(Rsolnp) # for optimization routines
#library(nloptr) # for optimization routines
library(bbmle)  # for MLE
library(CUB)    # for dissimilarity measure
#######################################
# a dissimilarity measure between     #
# two vectors of probabilities        #
#######################################
L2 <- function(f,p)
{
  k <- length(p)
  1/(1+1/k*sum((f/p-1)^2))
}
#######################################
# continuous Kumaraswamy distribution #
#######################################
library(extraDistr)
# plots of the pdf
k <- 5
# 5 x 5
op<-par()
par(mfrow=c(5,5), mar=c(1,1,1,1), mai=c(0.25,0.25,0.1,0.2), cex=.6)
alpha <- c(0.5, 0.75, 1, 1.5, 2.5) # values of the first shape par.
for(i in 1:length(alpha))
{
  for(j in 1:length(alpha))
  {
    cat("i=",i,"j=",j,"\n")
    plot(function(x) dkumar(x,alpha[i],alpha[j]),ylim=c(0,2))
    text.a <- as.expression(bquote(a == .(alpha[i])))
    text.b <- as.expression(bquote(b  == .(alpha[j])))
    text(x=.5,y=1.8, labels = text.a)
    text(x=.5,y=1.5, labels = text.b)
  }
}
par(op)
##################################
#### Kumaraswamy distribution ####
##################################
############## cdf ###############
##################################
F.K <- Vectorize(function(x, a, b)
{
  if(x<0 | x>1) return(0)
  else return(1-(1-x^a)^b)
}
)
F.K(0.5,2,2)
############################
###### integer moments #####
############################
M.K <- function(alpha,beta,r) beta*beta(1+r/alpha,beta)
# skewness
sk.K <- function(alpha,beta)
{
  E <- M.K(alpha,beta,1)
  V <- M.K(alpha,beta,2) - E^2
  E3 <- M.K(alpha,beta,3)
  (E3-3*E*V-E^3)/V^1.5
}
sk.K(3,3)

######################################
## discrete Kumaraswamy distribution #
#### with shape parameters a and b ###
# and parameter k (n. of categories) #
######################################
################ pmf #################
######################################
ddiscreteK <- function(k, a, b)
{
i <- 1:k
return(F.K((i)/k,a,b)-F.K((i-1)/k,a,b))
}
ddiscreteK(5,2,3)
######################################
################# cdf ################
######################################
pdiscreteK <- function(i,k,a,b)
{
  1-(1-(i/k)^a)^b
}
######################################
######## quantile function ###########
######################################
qdiscreteK <- function(u,k,a,b)
{
  ceiling(k*(1-(1-u)^(1/b))^(1/a))
}
######################################
######### random generation ##########
######################################
rdiscreteK <- function(n, k, a, b)
{
  u <- runif(n)
  qdiscreteK(u,k,a,b)
}
#
# graphs PMF
#
# 5 x 5
k <- 7
op<-par()
par(mfrow=c(5,5), mar=c(1,1,1,1), mai=c(0.25,0.25,0.1,0.2), cex=.6)
alpha<-c(0.5, 0.75, 1, 1.5, 2.5)
for(i in 1:length(alpha))
{
  for(j in 1:length(alpha))
  {
    cat("i=",i,"j=",j,"\n")
    plot(ddiscreteK(k,alpha[i],alpha[j]),type="h",ylim=c(0,0.95))
    points(1:k,ddiscreteK(k,alpha[i],alpha[j]),pch=19)
    text.a <- as.expression(bquote(a == .(alpha[i])))
    text.b <- as.expression(bquote(b  == .(alpha[j])))
    text(x=4,y=.85, labels = text.a)
    text(x=4,y=.7, labels = text.b)
  }
}
par(op)
#
###########################################
################# MOMENTS #################
###########################################
E.K <- Vectorize(function(a,b,k)
{
  i <- 0:(k-1)
  sum(1-F.K(i/k,a,b))
}
)
E.K(1,1,5)
#
E.K2 <- Vectorize(function(a,b,k)
{
  i <- 1:(k-1)
  1 + sum((2*i+1)*(1-F.K(i/k,a,b)))
}
)
E.K2(1,1,5)
#
V.K <- function(a,b,k) E.K2(a,b,k)-(E.K(a,b,k))^2
V.K(1,1,5)
Sk.K <- Vectorize(function(a,b,k)
{
  i <- 1:k
  sum((i-E.K(a,b,k))^3*ddiscreteK(k,a,b))/(E.K2(a,b,k)-(E.K(a,b,k))^2)^1.5
}
)
# graph of Skewness
# Create a sequence for a and b in the range (0, 4)
a <- seq(0.1, 4, by = 0.1)
b <- seq(0.1, 4, by = 0.1)
# Create a grid of points
zS5 <- outer(a, b, Sk.K, k=5)
zS7 <- outer(a, b, Sk.K, k=7)
zS9 <- outer(a, b, Sk.K, k=9)
# example of asymmetric discrete K with zero skewness
Sk.K.fixed <- function(b) Sk.K(a=2,b,k=7)
b <- uniroot(function(x) Sk.K.fixed(x), lower=1,upper=4)$root
ddiscreteK(k=7, a=2, b=b)

# Plot the level curves of E, V, and DI=V/E
k <- 5
f.E <- function(a, b, k=k) {
  E.K(a, b, k)
}
f.V <- function(a, b, k=k) {
  E.K2(a, b, k) - (E.K(a, b, k))^2
}
# Create a sequence for a and b in the range (0, 4)
a <- seq(0.01, 4, by = 0.01)
b <- seq(0.01, 4, by = 0.01)
# Create a grid of points
zE <- outer(a, b, f.E, k=5)
zV <- outer(a, b, f.V, k=5)
# Plot the level curves of E
op<-par()
par(mai=c(0.6,0.6,0.25,0.1), mgp=c(2.5,1,0), mfrow=c(1,3),cex=0.75)
contour(a, b, zE, xlab = "a", ylab = "b", main = "Level curves of expected value, k=5")
contour(a, b, zV, xlab = "a", ylab = "b", main = "Level curves of variance, k=5")
contour(a, b, zV/zE, xlab = "a", ylab = "b", main = "Level curves of DI, k=5")
par(op)

#################################################
############## E S T I M A T I O N ##############
#################################################
library(bbmle)
est.dK <- function(x, method="ML")
{
#k <- max(x)
if(method=="MM") # Method of moments
{
f.mom <- function(par,k)
{
  a <- par[1]
  b <- par[2]
  (mean(x)-E.K(a,b,k))^2 + (mean(x^2)-E.K2(a,b,k))^2
}
res <- optim(par=c(1,1), k=k, fn=f.mom, control=list(abstol=1e-10))
return(res$par)
}
else if(method=="MP") # Method of Proportions
{
f.prop <- function(b)
{
  1-(1-p1)^(1/b) - (1-pk^(1/b))^(-log(k)/(log(k-1)-log(k)))
}
p1<-mean(x==1)
pk<-mean(x==k)
if(p1==0 | pk==0) return(c(NA,NA))
else
{
eps<-.05 # or another arbitrary small value
b.P <- uniroot(f.prop, lower=eps, upper=10)$root
b.P
a.P <- log(1-(1-p1)^(1/b.P))/(-log(k))
a.P
return(c(a.P,b.P))
}
}
else if (method=="MCS") # Method of minimum Chi-squared
{
chi.square <- function(arg,k)
{
  a <- arg[1]
  b <- arg[2]
  sum((ddiscreteK(k, a, b) - sapply(1:k, function(val) mean(x == val)))^2/ddiscreteK(k, a, b))
}
res <- optim(par=c(1,1), fn=chi.square, k=k)
return(res$par)
}
else # Maximum Likelihood Method
{
log.lik.discreteK <- function(a,b,k)
{
  elle <- -sum(log(ddiscreteK(k,a,b)[x]))
  ifelse(elle==Inf, 999, elle)
}
if (length(unique(x)) == 2) return(list(c(NA,NA),matrix(NA,2,2)))
else
{
res<-mle2(minuslogl=log.lik.discreteK, fixed=list(k=k), start=list(a=1,b=1),method="L-BFGS-B", lower=c(0.01,0.01), upper=c(10,10))
CI <- confint(res)
CI[which(is.na(CI))]<-0 # !!!
return(list(res@coef,CI))
}
}
}

########################
####  DISCRETE BETA ####
# Sciandra et al. 2024 #
########################
# pmf
ddiscreteB <- function(k, a, b)
{
  i   <- 1:k
  pbeta(i/k,a,b) - pbeta((i-1)/k,a,b)
}
# Expectation
E.B <- function(a,b,k)
{
  j   <- 1:(k-1)
  k - sum(pbeta(j/k,a,b))
}
# Variance
V.B <- function(a,b,k)
{
  j   <- 1:(k-1)
  k^2 - sum((1+2*j)*pbeta(j/k,a,b)) - (E.B(a,b,k))^2
}
# log-likelihood function discrete Beta
log.lik.discreteB <- function(a,b,k)
{
  -sum(log(ddiscreteB(k,a,b)[x]))
}
# log-likelihood function discrete K.
log.lik.discreteK <- function(a,b,k)
{
  elle <- -sum(log(ddiscreteK(k,a,b)[x]))
  ifelse(elle==Inf, 99999, elle)
}

#########################
# Example of estimation #
####   EX SYMM K=5   ####
#########################
k<-5
x <- rep(1:5,c(10,15,50,15,10))
res<-mle2(minuslogl=log.lik.discreteK, fixed=list(k=k), start=list(a=1,b=1))
summary(res)
p.hat <- ddiscreteK(k,a=res@coef[1], b=res@coef[2])
p.hat
barplot(table(x)/sum(table(x)))
points(1:5+seq(-0.3,0.5,0.2), p.hat)
res<-mle2(minuslogl=log.lik.discreteB, fixed=list(k=k), start=list(a=1,b=1))
summary(res)

#########################
####   EX SYMM K=7   ####
#########################

k <- 7
x <- rep(1:7,c(5,10,20,30,20,10,5))
res <- mle2(minuslogl=log.lik.discreteK, fixed=list(k=k), start=list(a=1,b=1))
summary(res)
p.hat <- ddiscreteK(k,a=res@coef[1], b=res@coef[2])
p.hat
barplot(table(x)/sum(table(x)))
points(1:7+seq(-0.3,0.9,0.2), p.hat)
res <- mle2(minuslogl=log.lik.discreteB, fixed=list(k=k), start=list(a=1,b=1))
summary(res)


########################################################
############ REPARAMETRIZATION & REGRESSION ############
########################################################

# CDF (F.K)
F.K <- Vectorize(function(x, a, b) {
  if (x < 0 | x > 1) return(0)
  return(1 - (1 - x^a)^b)
})

# Compute a and b from eta and gamma
compute_ab <- function(eta, gamma) {
  a <- exp(gamma) * exp(eta) / (1 + exp(eta))
  b <- exp(gamma) * 1 / (1 + exp(eta))
  return(list(a = a, b = b))
}
# Compute eta and gamma from a and b
compute_eg <- function(a, b) {
  eta <- log(a) - log(b)
  gamma <- log(a+b)
  return(list(eta = eta, gamma = gamma))
}
# PMF
ddiscreteK_pmf <- Vectorize(function(x, k, eta, gamma,log=FALSE) {
  params <- compute_ab(eta, gamma)
  a <- params$a
  b <- params$b
  p <- F.K(x / k, a, b) - F.K((x - 1) / k, a, b)
  return(ifelse(log==FALSE,p,log(p)))
})

# Example of simulation for the regression model
# with discrete k response variable
# Define the covariates and the data
set.seed(12345) # for reproducibility
n <- 100
h <- 3          # N. of covariates (excluding the intercept)

# Randomly generate the covariates
X <- matrix(rnorm(n * h), nrow = n, ncol = h)
# Parameters and Observations
eta.v   <- -0.5 - 0.25*X[,1] + X[,2] - 0.1*X[,3]
gamma.v <- 0.5  + 0.25*X[,2] + X[,2] - 0.1*X[,3]
par     <- compute_ab(eta.v, gamma.v)
y       <- rdiscreteK(n, k=5, par$a, par$b)
mle_fit <- mle2(y~ddiscreteK_pmf(eta,gamma,k=5),     ## specify the distribution of the response y
                data=data.frame(y,X),                ## need to specify as data frame
                parameters=list(eta~X,gamma~X),      ## linear model for eta and gamma
                start=list(eta=0,gamma=0))

summary(mle_fit)
####################################################
########## R E G R E S S I O N   M O D E L #########
# S T A R T   O F    S I M U L A T I O N   P L A N #
####################################################
n <- 250
S <- 1000
estimates <- matrix(0,S,6)
A <- matrix(0, S, n)
B <- matrix(0, S, n)
set.seed(12345)
for(s in 1:S)
{
  # sampling x1 from a standard normal
  x1 <- rnorm(n)
  # sampling x2 from a uniform in (-1,1)
  x2 <- runif(n,-1,1)
  # regression equations for the transformed parameters eta and gamma
  eta   <- 0      + 0.5*x1 + 0.75*x2
  gamma <- log(2) - 0.25*x1 + 0.5*x2
  # computing the corresponding parameters a and b
  param <- compute_ab(eta,gamma)
  A[s,]<-param$a
  B[s,]<-param$b
  # sampling from the discrete Kumaraswamy
  y <- rdiscreteK(n,param$a,param$b,k=5)
  x <- data.frame(x1,x2)
  # computing MLEs
  mle_fit <- mle2(y~ddiscreteK_pmf(eta,gamma,k=max(y)),
                  data=data.frame(y,x),
                  parameters=list(eta~unlist(x1)+unlist(x2),gamma~unlist(x1)+unlist(x2)),           ## linear model for loglambda
                  start=list(eta=0,gamma=0))
  estimates[s,] <- mle_fit@coef
  print(s)
}
apply(estimates,2,mean)
save(estimates,A,B,file="est.regrnew.Rdata")
names <- c(expression(lambda[0]), expression(lambda[1]), expression(lambda[2]),expression(omega[0]),
           expression(omega[1]), expression(omega[2]))
boxplot(estimates,names=names)
# regression coefficients used
par.reg <- c(0,0.5,0.75,log(2),-0.25,0.5)
for(i in 1:6)
{
  segments(x0=i-0.525,x1=i+0.525,y0=par.reg[i],col="blue",lty=2)
}
AA <- as.vector(A)
BB <- as.vector(B)
summary(AA)
summary(BB)
mean(AA>1 & BB>1)
mean(AA<1 & BB<1)
mean(AA>1 & BB<1)
mean(AA<1 & BB>1)
cor(AA,BB)
####################################################
########## R E G R E S S I O N   M O D E L #########
### E N D   O F    S I M U L A T I O N   P L A N ###
####################################################

#############################################
########### E X A M P L E S   O F ###########
######### D A T A   A N A L Y S I S #########
#############################################

#############################################
################## EUROBAR ##################
#############################################
# Eurobarometer survey data can be found at
# https://data.europa.eu/data/datasets/s703_72_1_ebs322?locale=en
# here I directly construct the full EU data set
# related to the question
# `How serious a problem do you think climate change is at this moment?
# Please use a scale from 1 to 10
# by providing
# absolute frequencies ni
ni <- c(408,357,826,1216,3095,3219,4691,4998,2295,4757)
# of the 10 ordered categories xi
xi <- 1:10
x  <- rep(xi, ni)
k  <- length(xi)
# ML estimation of discrete Kumaraswamy
resK <- mle2(minuslogl=log.lik.discreteK, fixed=list(k=k), start=list(a=1,b=1))
summary(resK)
AIC(resK)
pe <- table(x)/sum(table(x))
barplot(table(x)/sum(table(x)),ylim=c(0,0.225))
pK<-ddiscreteK(k, resK@coef[1],resK@coef[2])
points(1:10+seq(-0.3,1.5,0.2),pK,pch=19,col="black",type="b")
# ML estimation of discrete Beta
resBeta<-mle2(minuslogl=log.lik.discreteB, fixed=list(k=k), start=list(a=1,b=1))
summary(resBeta)
AIC(resBeta)
pB<-ddiscreteB(k, resBeta@coef[1],resBeta@coef[2])
points(1:10+seq(-0.3,1.5,0.2),pB,pch=12,col="red",type="b")
legend("topleft", c("dK","dBeta"),bty="n", pch=c(19,12),col=c("black","red"))
mtext(expression(f[i]),side=2,las=1,line=1.5, at=0.125)
legend("topright","EU",bty="n")
box()
#
# the data set for ITALIANS
ni <- c(9,12,27,40,83,173,186,204,95,168)
xi <- 1:10
x <- rep(xi, ni)
k <- length(xi)
resK<-mle2(minuslogl=log.lik.discreteK, fixed=list(k=k), start=list(a=1,b=1))
summary(resK)
AIC(resK)
pe <- table(x)/sum(table(x))
barplot(table(x)/sum(table(x)),ylim=c(0,0.225))
pK<-ddiscreteK(k, resK@coef[1],resK@coef[2])
points(1:10+seq(-0.3,1.5,0.2),pK,pch=19,type="b")
resBeta<-mle2(minuslogl=log.lik.discreteB, fixed=list(k=k), start=list(a=1,b=1))
summary(resBeta)
AIC(resBeta)
pB<-ddiscreteB(k, resBeta@coef[1],resBeta@coef[2])
points(1:10+seq(-0.3,1.5,0.2),pB,pch=12,col="red",type="b")
legend("topleft", c("dK","dBeta"),bty="n", pch=c(19,12),col=c("black","red"))
legend("topright","Italy",bty="n")
mtext(expression(f[i]),side=2,las=1,line=1.5, at=0.125)
box()

#############################################
################## PREFMOD ##################
#############################################
library(prefmod)
data(issp2000)
m <- issp2000#[issp2000$CNTRY==2,]
#
CAR <- m$CAR
table(CAR)
x <- CAR
k <- 5
resK <- mle2(minuslogl=log.lik.discreteK, fixed=list(k=k), start=list(a=1,b=1))
summary(resK)
op<-par()
par(mai=c(0.5,0.5,0.1,0.1), las=1)
pe <- table(x)/sum(table(x))
barplot(table(x)/sum(table(x)),ylim=c(0,0.475))
pK<-ddiscreteK(k, resK@coef[1],resK@coef[2])
points(1:5+seq(-0.3,0.5,0.2),pK,pch=19,type="b")
mtext(expression(f[i]),side=2,las=1,line=1.5)
#
resBeta<-mle2(minuslogl=log.lik.discreteB, fixed=list(k=k), start=list(a=1,b=1))
summary(resBeta)
pB<-ddiscreteB(5, resBeta@coef[1],resBeta@coef[2])
points(1:5+seq(-0.3,0.5,0.2),pB,pch=12,col="red",type="b")
legend("topright", c("dK","dBeta"),bty="n", pch=c(19,12),col=c("black","red"))
box()
#
dissim(pe,pK)
dissim(pe,pB)
# estimates, SE, and max log.lik.
c(resK@coef[1],sqrt(diag(resK@vcov)[1]),resK@coef[2],sqrt(diag(resK@vcov)[2]),-resK@min)
c(resBeta@coef[1],sqrt(diag(resBeta@vcov)[1]),resBeta@coef[2],sqrt(diag(resBeta@vcov)[2]),-resBeta@min)
#
CAR <- issp2000$CAR
SEX <- factor(issp2000$SEX)
URB <- factor(issp2000$URB)
AGE <- factor(issp2000$AGE,ordered=TRUE)
CNTRY <- factor(issp2000$CNTRY)
EDU <- factor(issp2000$EDU, ordered=TRUE)
contrasts(EDU) <- contr.treatment(2)
contrasts(AGE) <- contr.treatment(3)
# regression model including all covariates
X <- data.frame(SEX,URB,AGE,CNTRY,EDU)
mle_fit <- mle2(CAR~ddiscreteK_pmf(eta,gamma,k=max(CAR)),     ## use log link/exp inverse-link
                data=data.frame(CAR,X),                       ## need to specify as data frame
                parameters=list(eta~unlist(SEX)+unlist(URB)+unlist(AGE)+unlist(CNTRY)+unlist(EDU),gamma~unlist(SEX)+unlist(URB)+unlist(AGE)+unlist(CNTRY)+unlist(EDU)),
                start=list(eta=0,gamma=0))
summary(mle_fit)
# ADDED - godness-of-fit check
eta.coef   <- mle_fit@fullcoef[1:8]
gamma.coef <- mle_fit@fullcoef[9:16]
A<-model.matrix(~SEX+URB+AGE+CNTRY+EDU,data=X)
A%*%eta.coef
A%*%gamma.coef
parK.est <- compute_ab(A%*%eta.coef,A%*%gamma.coef)
set.seed(538)
n <- length(x)
y.sim <- rdiscreteK(n=n, k=5, parK.est[[1]], parK.est[[2]])
table(y.sim)
table(x)

# regression model with country and sex as covariates
X <- data.frame(CNTRY,SEX)
mle_fit2 <- mle2(CAR~ddiscreteK_pmf(eta,gamma,k=max(CAR)),      ## use log link/exp inverse-link
                data=data.frame(CAR,CNTRY),                     ## need to specify as data frame
                parameters=list(eta~unlist(CNTRY)+unlist(SEX),gamma~unlist(CNTRY)+unlist(SEX)),
                start=list(eta=0,gamma=0))
summary(mle_fit2)
AIC(mle_fit)
AIC(mle_fit2)
#

###########################################
#################### UTAH  ################
###########################################

# Data from the Utah Air Quality Risk and Behavioral Action Survey
# are freely available in the ICPSR repository
# https://www.openicpsr.org/openicpsr/project/117904/version/V1/view?path=/openicpsr/117904/fcr:versions/V1&type=project
# we focus on item Q18_04
# and directly provide the distribution of responses
n <- c(369,291,206,172,65,33,10,6)
h <- 1:8
x <- rep(h, n)
k <- length(h)
# fitting the discrete Kumaraswamy distribution
resK<-mle2(minuslogl=log.lik.discreteK, fixed=list(k=k), start=list(a=1,b=1))
AICK <- 2*resK@min+2*2
AICK
# fitting the discrete Beta distribution
resB<-mle2(minuslogl=log.lik.discreteB, fixed=list(k=k), start=list(a=1,b=1))
AICB <- 2*resB@min+2*2
AICB

c(resK@coef[1],sqrt(diag(resK@vcov)[1]),resK@coef[2],sqrt(diag(resK@vcov)[2]),-resK@min)
c(resB@coef[1],sqrt(diag(resB@vcov)[1]),resB@coef[2],sqrt(diag(resB@vcov)[2]),-resB@min)

# barplot with points
f <- n/sum(n)
par(mai=c(0.5,0.6,0.1,0.1), las=1)
barplot(f, names.arg = 1:8, ylim=c(0,max(f)*1.1))
points((1:k)+seq(-0.3,1.1,0.2),ddiscreteK(k, resK@coef[1], resK@coef[2]),pch=19,type="o")
points((1:k)+seq(-0.3,1.1,0.2),ddiscreteB(k, resB@coef[1], resB@coef[2]),pch=12,type="o",col="red")
mtext(expression(f[i]),side=2,las=1,line=1.5,at=0.225)
legend("topright", c("dK","dBeta"), pch=c(19,12),col=c("black","red"))
box()
#almost the same
pK <- ddiscreteK(k, resK@coef[1], resK@coef[2])
pB <- ddiscreteB(k, resB@coef[1], resB@coef[2])
dissim(f,pK);L2(f,pK)
dissim(f,pB);L2(f,pB)

library(readxl)
# reading from the Excel file in my directory
dataset <- read_excel("UTAH.xlsx")
y <- dataset[[7]] # MAJOR INDUSTRY
index.na <- which(y=="MD")
# delete NA in y
dataset <- dataset[-index.na,]
dim(dataset)
y <- dataset[[7]]
y <- as.numeric(y)
levels <-c("Less than $25,000","$25,000 to $34,999","$35,000 to $49,999","$50,000 to $74,999","$75,000 to $99,999","$100,000 to $149,999","Greater than $150,000","Prefer not to answer")
X1 <- factor(dataset[[1]],ordered=FALSE,levels=levels)
X2 <- factor(dataset[[2]],levels=c("Democrat","Independent","Libertarian","Republican","Other","No preference"))

X <- data.frame(X1,X2)
mle_fit <- mle2(y~ddiscreteK_pmf(eta,gamma,k=max(y)),     ## specify the response distribution
                data=data.frame(y,X),                     ## need to specify data as data frame
                parameters=list(eta~unlist(X1)+unlist(X2),## linear model for transformed par. eta and gamma
                gamma~unlist(X1)+unlist(X2)),
                start=list(eta=0,gamma=0))

summary(mle_fit)
AIC(mle_fit) # AIC=3780.976

# ADDED - godness-of-fit check
n <- length(y)
eta.coef   <- mle_fit@fullcoef[1:13]
gamma.coef <- mle_fit@fullcoef[14:26]
A<-model.matrix(~X1+X2,data=X)
A%*%eta.coef
A%*%gamma.coef
parK.est <- compute_ab(A%*%eta.coef,A%*%gamma.coef)
set.seed(538)
y.sim <- rdiscreteK(n=n, k=8, parK.est[[1]], parK.est[[2]])
table(y.sim)
table(y)

#######################################################
##################  S T A R T   O F  ##################
############ S I M U L A T I O N   P L A N ############
#######################################################
set.seed(12345)
nS <- 10000 # n. of MC runs
A <- matrix(0, nS, 4)
B <- matrix(0, nS, 4)
CI.A <- matrix(0,nS,2)
CI.B <- matrix(0,nS,2)
n <- 200 # sample size
X <- matrix(0, nS, n)
a <- 1.5   # (0.5, 0.75, 1, 1.5, 2)
b <- 0.75 # (0.5, 0.75, 1, 1.5, 2)
k <- 7   # 5 or 7
for(i in 1:nS)
{
  cat("Simulazione n. ",i,"\n")
  x <- rdiscreteK(n, k, a, b)
  X[i,] <- x
  res.P   <- est.dK(x, "MP")
  res.M   <- est.dK(x, "MM")
  res.MCS <- est.dK(x, "MCS")
  res.ML  <- est.dK(x, "ML")
  A[i, 1] <- res.P[1]
  B[i, 1] <- res.P[2]
  A[i, 2] <- res.M[1]
  B[i, 2] <- res.M[2]
  A[i, 3] <- res.MCS[1]
  B[i, 3] <- res.MCS[2]
  A[i, 4] <- res.ML[[1]][1]
  B[i, 4] <- res.ML[[1]][2]
  CI.A[i,] <- res.ML[[2]][1,]
  CI.B[i,] <- res.ML[[2]][2,]
}
missing.rate <- mean(is.na(A[,1]))
lab <- c("P","M","MCS","ML")
colnames(A) <- lab
colnames(B) <- lab

m.A <- apply(A,2,mean, na.rm=TRUE)
m.B <- apply(B,2,mean, na.rm=TRUE)
# MC bias
bias.A <- apply(A-a,2,mean, na.rm=TRUE)
bias.B <- apply(B-b,2,mean, na.rm=TRUE)
# MC variance
var.A  <- apply(A^2, 2, mean, na.rm=TRUE) - apply(A, 2, mean, na.rm=TRUE)^2
var.B <- apply(B^2, 2, mean, na.rm=TRUE) - apply(B, 2, mean, na.rm=TRUE)^2
# MC root mean squared error
rmse.A <- sqrt(apply((A-a)^2,2,mean, na.rm=TRUE))
rmse.B <- sqrt(apply((B-b)^2,2,mean, na.rm=TRUE))
# coverage
cov.A <- mean(CI.A[,1] <= a & a <= CI.A[,2], na.rm=TRUE)
cov.B <- mean(CI.B[,1] <= b & b <= CI.B[,2], na.rm=TRUE)
# length
len.A <- mean(CI.A[,2] - CI.A[,1], na.rm=TRUE)
len.B <- mean(CI.B[,2] - CI.B[,1], na.rm=TRUE)
#
filename <- paste("discreteK a",a,"b",b,"k",k,"n",n)
save(file=paste(filename,".Rdata",sep=""),list=ls())
# boxplots
boxplot(A, main=paste("estimators of parameter a"))
abline(h=a, lty=3)
boxplot(B, main=paste("estimators of parameter b"))
abline(h=b, lty=3)
# scatter plots
op<-par()
par(mfrow=c(2,2), mai=c(0.65,0.7,0.2,0.1), mgp=c(2.5,1,0))
for(i in 1:4)
{
  plot(A[,i],B[,i],main=lab[i],xlab=expression(hat(a)),ylab=expression(hat(b)))
  lab.rho <- bquote(r==.(round(cor(A[,i],B[,i],use="complete.obs"),3)))
  text(x=min(A[,i],na.rm=TRUE), y=.98*max(B[,i],na.rm=TRUE), labels=lab.rho,pos=4)
}
par(op)
# printing on a txt file the main results
filetxt <- paste(filename,"txt",sep=".")
sink(file=filetxt)
cat("a=",a,"\n")
cat("b=",b,"\n")
cat("k=",k,"\n")
cat("n=",n,"\n")
cat("\n")
cat("MC mean of a and b:","\n")
m.A
m.B
cat("\n")
cat("MC BIAS of a and b:","\n")
bias.A
bias.B
cat("\n")
cat("MC Variance of a and b:","\n")
var.A
var.B
cat("\n")
cat("MC RMSE of a and b:","\n")
rmse.A
rmse.B
cat("\n")
cat("Coverage of CIs for a and b:","\n")
cov.A
cov.B
cat("\n")
cat("Length of CIs for a and b:","\n")
len.A
len.B
cat("\n")
cat("Rate of missing data for the method of proportion:","\n")
missing.rate
sink()
#######################################################
##################### E N D   O F #####################
############ S I M U L A T I O N   P L A N ############
#######################################################
