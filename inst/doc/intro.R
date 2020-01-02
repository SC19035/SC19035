## ----eval=FALSE----------------------------------------------------------
#  mytable <-function(x,y){
#    x=as.factor(x)
#    y=as.factor(y)
#    xi=as.integer(x)
#    xm=as.integer(max(xi))
#    yi=(as.integer(y)-1L)*xm
#    ym=as.integer(max(yi)+xm)
#    matrix(.Internal(tabulate(xi+yi,ym)),xm)
#  }

## ----eval=FALSE----------------------------------------------------------
#  cdf <- function(x,p,q) {
#    n<-length(x)
#    res<-numeric(n)
#    for(i in 1:n){
#      res[i]<-integrate(f, lower=-Inf, upper=x[i],
#                        rel.tol=.Machine$double.eps^0.25,
#                        p=p,q=q)$value
#    }
#    return(res)
#  }

## ------------------------------------------------------------------------
paste("R","for Beginning")
paste("Today is",date())

## ------------------------------------------------------------------------
A<-array(1:9,dim=(c(3,3)))
B<-array(9:1,dim=(c(3,3)))
C<-A%*%B;C

## ------------------------------------------------------------------------
X1<-c(35,42,46,37,38,42,41)
X2<-c(60,72,64,62,71,69,75)
plot(X1,X2)
hist(X1)

## ------------------------------------------------------------------------
rrayleigh<-function(n,sigma){
  sigma<-1:1000
  u<-runif(n)
  x<-sqrt(-2*sigma**{2}*log(1-u))
  return (x)
}

n1<-1000
sigma1<-4
rr=head(rrayleigh(n1,sigma1),10)
hist(rrayleigh(n1,sigma1),prob=TRUE)

## ------------------------------------------------------------------------
n<-1000
u<-runif(n)
sigma<-10
x<-sqrt(-2*sigma**{2}*log(1-u))
hist(x,prob=TRUE,main = expression(f(x)==x/(sigma**2)*exp(-((x/sigma)**2)/2)))

## ------------------------------------------------------------------------
n<-1000
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
u<-runif(n)
p<-as.integer(u>0.25)
x<-p*x1+(1-p)*x2 ##the mixture

## ------------------------------------------------------------------------
hist(x, prob = TRUE)

## ------------------------------------------------------------------------
q1<-as.integer(u>0.5)
y1<-q1*x1+(1-q1)*x2
hist(y1,prob=TRUE)
q2<-as.integer(u>0.42)
y2<-q2*x1+(1-q2)*x2
hist(y2,prob=TRUE)

## ------------------------------------------------------------------------
require(graphics)
runif.wis<-function(sigma,n){
  a<-ncol(sigma)
  b<-chol(sigma)
  X<-matrix(ncol=a,chol=a)
  for(i in 1:a){
    for(j in 1:a){
      X[i,j]=rnorm(2,4)
    }
  }
 X[lower.tri(X)]=0
 C<-numeric(a)
 for(i in 1:a){
 C[i]=rchisq(1,n-i+1)
  X[i,i]=C[i]
 }
  d<-b%*%X
  e<-t(X)%*%t(b)
  f<-d%*%e
 return(f)
}

## ------------------------------------------------------------------------
m <- 10000
x <- runif(m)
theta.hat <- mean(sin(pi/3*x)) * (pi/3)
print(theta.hat)
print(1-cos(pi/3))

## ------------------------------------------------------------------------
## the definition of function
MC.Phi <- function(x, R = 10000, antithetic = TRUE) {
u <- runif(R/2)
if (!antithetic) v <- runif(R/2) else
v <- 1 - u
u <- c(u, v)
cdf <- numeric(length(x))
for (i in 0:length(x)) {
g <-mean(exp(-u * x)/(1+(u*x)^2))
}
g
}
## Use Monte Carlo integration with antithetic variables to estimate
set.seed(200)
MC1<-MC.Phi(1,anti=FALSE)
MC2<-MC.Phi(1)
print(MC1)
print(MC2)

## find the approximate reduction in variance
m <- 1000
MC1 <- MC2<-numeric(m)
x<-1
for (i in 1:m) {
MC1[i] <- MC.Phi(x,R = 1000, anti = FALSE)
MC2[i] <- MC.Phi(x,R = 1000)
}
 print(sd(MC1))
 print(sd(MC2))
 print((var(MC1)-var(MC2))/var(MC1))

## ------------------------------------------------------------------------
m <- 10000
theta.hat <- se <- numeric(5)
g <- function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
# choose the importance function : f3=exp（-x）/(1+exp（-1）)
u <- runif(m) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg<- function(x){ g(x) / (exp(-x) / (1 - exp(-1)))
}
#stratify
M <- 10000 #number of replicates
k <- 5 #number of strata
r <- M / k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
estimates <- matrix(0, N, 2)
# estimate
for (i in 1:N) {
estimates[i, 1] <- mean(fg(runif(M)))
for (j in 1:k)
T2[j] <- mean(fg(runif(M/k, (j-1)/k, j/k)))
estimates[i, 2] <- mean(T2)
}
apply(estimates, 2, mean)
apply(estimates, 2, var)

## ------------------------------------------------------------------------

alpha <- .5 
n <- 20 
m <- 1000 
beta<- seq(1, 40, 1)
N <- length(beta) 
pwr <- numeric(N) 
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))#critical value for the skewness test 


sk <- function(x) { #computes the sample skewness coeff. 
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) 
  } 

for (j in 1:N) { #for each beta 
  b<- beta[j] 
  sktests <- numeric(m) 
  for (i in 1:m) { #for each replicate 
    x <- rbeta(n,b,b)
    sktests[i] <- as.integer(abs(sk(x)) >= cv) 
    } 
  pwr[j] <- mean(sktests) } #plot power vs beta 


plot(beta, pwr, type = "b", xlab = bquote(beta), ylim = c(0.3,0.5))
abline(h = .1, lty = 3) 
se <- sqrt(pwr * (1-pwr) / m) #add standard errors 
lines(beta, pwr+se, lty = 3) 
lines(beta, pwr-se, lty = 3)



## ------------------------------------------------------------------------
n <- 20 
alpha <- .01 
mu0 <- 1 
m <- 10000 #number of replicates 
p <- numeric(m) #storage for p-values 
for (j in 1:m) { 
  x <- rchisq(n,1)
  ttest <- t.test(x, alternative = "greater", mu = mu0) 
  p[j] <- ttest$p.value 
  }
p.hat <- mean(p < alpha) 
se.hat <- sqrt(p.hat * (1 - p.hat) / m) 
print(c(p.hat, se.hat))

## ------------------------------------------------------------------------
n <- 20 
alpha <- .01 
mu0 <- 1 
m <- 10000 #number of replicates 
p <- numeric(m) #storage for p-values 
for (j in 1:m) { 
  x <- runif(n,0,2)
  ttest <- t.test(x, alternative = "greater", mu = mu0) 
  p[j] <- ttest$p.value 
  }
p.hat <- mean(p < alpha) 
se.hat <- sqrt(p.hat * (1 - p.hat) / m) 
print(c(p.hat, se.hat))

## ------------------------------------------------------------------------
n <- 20 
alpha <- .01 
mu0 <- 1 
m <- 10000 #number of replicates 
p <- numeric(m) #storage for p-values 
for (j in 1:m) { 
  x <- rexp(n,1)
  ttest <- t.test(x, alternative = "greater", mu = mu0) 
  p[j] <- ttest$p.value 
  }
p.hat <- mean(p < alpha) 
se.hat <- sqrt(p.hat * (1 - p.hat) / m) 
print(c(p.hat, se.hat))

## ------------------------------------------------------------------------
 library(bootstrap)
data(scor)
pairs(scor[1:5], main = "the scatter diagram ",bg = c("red", "green", "blue","yellow","black"))

## ------------------------------------------------------------------------
library('GGally')
set.seed(1836)
sc <- scor[,c(1:5)]
ggpairs(sc[sample.int(nrow(sc),20),])

## ------------------------------------------------------------------------
# Estimate standard error
#set up the bootstrap 
B <- 100 #number of replicates 
n <- nrow(scor) #sample size 
R12 <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R 
for (b in 1:B) { 
        i <- sample(1:n, size = n, replace = TRUE) 
        MEC <- scor$mec[i] #i is a vector of indices 
        VEC <- scor$vec[i] 
        R12[b] <- cor(MEC, VEC) 
 } 
print(se.R12 <- sd(R12))

## ------------------------------------------------------------------------
B <- 100 
n <- nrow(scor)  
R34 <- numeric(B)  
for (b in 1:B) { 
        i <- sample(1:n, size = n, replace = TRUE) 
        ALG <- scor$alg[i] 
        ANA <- scor$ana[i] 
        R34[b] <- cor(ALG, ANA) 
 } 
print(se.R34 <- sd(R34))

## ------------------------------------------------------------------------
B <- 100 
n <- nrow(scor)  
R35 <- numeric(B)  
for (b in 1:B) { 
        i <- sample(1:n, size = n, replace = TRUE) 
        ALG <- scor$alg[i] 
        STA <- scor$sta[i] 
        R35[b] <- cor(ALG, STA) 
 } 
print(se.R35 <- sd(R35))

## ------------------------------------------------------------------------
B <- 100 
n <- nrow(scor)  
R45 <- numeric(B)  
for (b in 1:B) { 
        i <- sample(1:n, size = n, replace = TRUE) 
        ANA <- scor$ana[i] 
        STA <- scor$sta[i] 
        R45[b] <- cor(ANA, STA) 
 } 
print(se.R45 <- sd(R45))

## ----exercise 7.B--------------------------------------------------------

set.seed(1111)
library(boot)
n <- 100
x <- rnorm(n, 0, 1)

sk <- function(x,i) {
xbar <- mean(x[i])
m3 <- mean((x[i] - xbar)^3)
m2 <- mean((x[i] - xbar)^2)
return( m3 / m2^1.5 )
}

m<-50     
# sample size
norm1_cil<-norm1_cir<-numeric(n)
basic1_cil<-basic1_cir<-numeric(n)
perc1_cil<-perc1_cir<-numeric(n)
# Monte Carlo and bootstrap
for (i in 1:n) {
  x<-rnorm(m)
  boot.obj<-boot(x,statistic=sk,R=2000)
  ci<-boot.ci(boot.obj,type=c('norm','basic','perc'))
  norm1_cil[i]<-ci$normal[2]
  norm1_cir[i]<-ci$normal[3]
  basic1_cil[i]<-ci$basic[4]
  basic1_cir[i]<-ci$basic[5]
  perc1_cil[i]<-ci$percent[4]
  perc1_cir[i]<-ci$percent[5]
}
# compute coverage rates 
print(c(mean((norm1_cil<0)*(norm1_cir>0)),mean((basic1_cil<0)*(basic1_cir>0)),
        mean((perc1_cil<0)*(perc1_cir>0))))


## ------------------------------------------------------------------------
set.seed(1111)
library(boot)
n <- 100
x <- rgamma(n, 0, 5)

sk <- function(x,i) {
xbar <- mean(x[i])
m3 <- mean((x[i] - xbar)^3)
m2 <- mean((x[i] - xbar)^2)
return( m3 / m2^1.5 )
}

m<-50     
# sample size
norm2_cir<-numeric(n)
basic2_cir<-numeric(n)
perc2_cir<-numeric(n)
# Monte Carlo and bootstrap
for (i in 1:n) {
  x <- rgamma(n,5)
  boot.obj<-boot(x,statistic=sk,R=2000)
  ci<-boot.ci(boot.obj,type=c('norm','basic','perc'))
  norm2_cir[i]<-ci$normal[3]
  basic2_cir[i]<-ci$basic[3]
  perc2_cir[i]<-ci$percent[3]
}
# compute coverage rates 
print(c(mean(norm2_cir>0),mean(basic2_cir>0),
        mean(perc2_cir>0)))

## ------------------------------------------------------------------------

library(bootstrap)
data(scor)
scor_pca <- princomp(scor, cor = F)
print(scor_pca )
theta.hat <- 0.619115#the proportion of variance explained by the first principal component is 0.619115

#jackknife
n <- nrow(scor)
theta.hat.j <- numeric(n)
for (i in 1:n) {
scor_j <- scor[-i, ]
scor_j_cov <- cov(scor_j)
theta.hat.j[i] <-eigen(scor_j_cov)$value[1] / sum(eigen(scor_j_cov)$value)
}

#bias
bias <- (n - 1) * (mean(theta.hat.j) - theta.hat)
print(bias)

#standard error
se <-sqrt((n - 1) * mean((theta.hat.j - mean(theta.hat.j)) ^ 2))
print(se)



## ------------------------------------------------------------------------

library(DAAG)
attach(ironslag)

L1 <-  lm(magnetic ~ chemical)
L2 <-  lm(magnetic ~ chemical + I(chemical ^ 2))
L3 <-  lm(log(magnetic) ~ chemical)
L4 <-  lm(magnetic ~ chemical +I(chemical^2)+ I(chemical^3)) # replacing the Log-Log model with a cubic polynomial model

n <- length(ironslag$magnetic) #in DAAG ironslag 
e1 <- e2 <- e3 <- e4 <- numeric(n)
for (k in 1:n) { 
  y <- magnetic[-k] 
  x <- chemical[-k]
  
  J1 <- lm(y ~ x) 
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k] 
  e1[k] <- magnetic[k] - yhat1

  J2 <- lm(y ~ x + I(x^2)) 
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2 
  e2[k] <- magnetic[k] - yhat2

  J3 <- lm(log(y) ~ x) 
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k] 
  yhat3 <- exp(logyhat3) 
  e3[k] <- magnetic[k] - yhat3

  J4 <- lm(y ~ x +I(x^2)+ I(x^3)) 
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k]+ J4$coef[3]*chemical[k]^2  +J4$coef[4]*chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
} 
print(c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)))

## ------------------------------------------------------------------------
print(J1)#Y=1.7881+0.8895X
print(summary(J1))#Multiple R-squared is 0.5411

print(J2)#Y=25.98525-1.52948*X+0.05707*X^2
print(summary(J2))#Multiple R-squared is 0.607


print(J3)#LOG(Y)=2.1274+0.04057*X
print(summary(J3))#Multiple R-squared is 0.5307

print(J4)#Y=-2.7637+2.9084*X+(-0.1609)*X^2+0.0034*X^3
print(summary(J4))#Multiple R-squared is 0.6176



## ------------------------------------------------------------------------

n1 <- 10
n2 <- 20
mu1 <- mu2 <- 0 
sigma1 <- sigma2 <- 1 
m <- 1000


##Estimate the m

alpha<-0.06
mx<-log(alpha/2)/log(n1/(n1+n2))
my<-log(alpha/2)/log(n2/(n1+n2))
print(c(mx,my))#choose the bigger my 


# the critical value for both Cx,cy is m=8
counttest <- function(x, y) { 
  X <- x - mean(x) 
  Y <- y - mean(y) 
  outx <- sum(X > max(Y)) + sum(X < min(Y)) 
  outy <- sum(Y > max(X)) + sum(Y < min(X)) # return 1 (reject) or 0 (do not reject H0) 
  return(as.integer(max(c(outx, outy)) > 8)) 
  }


tests <- replicate(m, expr = {
  x <- rnorm(n1, mu1, sigma1) 
  y <- rnorm(n2, mu2, sigma2) 
  x <- x - mean(x) #centered by sample mean 
  y <- y - mean(y) 
  counttest(x, y) })
alphahat <- mean(tests) 
print(alphahat) 


## ------------------------------------------------------------------------
library(boot)
library(Ball)
library(MASS)
library(ggplot2)
alpha <- 0.05
n<-c(seq(1,100,5))
k<-50
#denefite the statistic of cov,ndcov2

dCov <- function(x, y) { 
  x <- as.matrix(x) 
  y <- as.matrix(y) 
  n <- nrow(x) 
  m <- nrow(y) 
  if (n != m || n < 2) stop("Sample sizes must agree") 
  if (! (all(is.finite(c(x, y))))) stop("Data contains missing or infinite values")
  
  Akl <- function(x) {
    d <- as.matrix(dist(x)) 
    m <- rowMeans(d) 
    M <- mean(d) 
    a <- sweep(d, 1, m) 
    b <- sweep(a, 2, m) 
    return(b + M) 
    } 
  A <- Akl(x) 
  B <- Akl(y) 
  dCov <- sqrt(mean(A * B)) 
  dCov
}


ndCov2 <- function(z, ix, dims) { #dims contains dimensions of x and y 
  p <- dims[1] 
  q <- dims[2] 
  d <- p + q 
  x <- z[ , 1:p] #leave x as is 
  y <- z[ix, -(1:p)] #permute rows of y 
  return(nrow(z) * dCov(x, y)^2) 
  } 


set.seed(1001) 


pow_dCor_Model1<-pow_ball_Model1<-pow_dCor_Model2<-pow_ball_Model2<-numeric(length(n))
dcor1<-dcor2<-p_ball1<-p_ball2<-numeric(k)

for (j in 1:length(n)) {
  dcor1 <-p_ball1< -numeric(k)
  for (i in 1:k) {
    set.seed(i)
    # the function "rmvnorm" is used to 
    # generate the multidimensional normal data
    X<-mvrnorm(n[j],rep(0,2),diag(2))
    e<-mvrnorm(n[j],rep(0,2),diag(1,2))
    Y1<-(X/4)+e
    Z1<-cbind(X,Y1)
   
     t1<-bcov.test(X,Y1,R=99)
    p_ball1[i]<-t1$p.value
    boot.obj1<-boot(data=Z1,statistic=ndCov2,R=99,sim="permutation",dims=c(2, 2))
    temp1<-c(boot.obj1$t0, boot.obj1$t)
    dcor1[i]<-mean(temp1>=temp1[1])
    }
  pow_dCor_Model1[j]<-mean(dcor1<alpha)
  pow_ball_Model1[j]<-mean(p_ball1<alpha)
}
dat1<-data.frame(pow_dCor_Model1,pow_ball_Model1)

ggplot(dat1,aes(n))+geom_point(y=pow_dCor_Model1,fill="white")+geom_line(y=pow_dCor_Model1,colour="red")+geom_point(y=pow_ball_Model1,fill="white")+geom_line(y=pow_ball_Model1,colour="blue")


for (j in 1:length(n)) {
dcor2<- p_ball2<-numeric(k)
for (i in 1:k) {
    set.seed(i)
    # the function "rmvnorm" is used to 
    # generate the multidimensional normal data
    X<-mvrnorm(n[j],rep(0,2),diag(2))
    e<-mvrnorm(n[j],rep(0,2),diag(2))
    Y2<-(X/4)*e
    Z2<-cbind(X,Y2)
   
 t2<-bcov.test(X,Y2,R=99)
    p_ball2[i]<-t2$p.value
    boot.obj2<-boot(data=Z2,statistic=ndCov2,R=200,sim="permutation",dims=c(2, 2))
    temp2<-c(boot.obj2$t0, boot.obj2$t)
    dcor2[i]<-mean(temp2>=temp2[1])
}
  pow_dCor_Model2[j]<-mean(dcor2<alpha)
  pow_ball_Model2[j]<-mean(p_ball2<alpha)
  }
dat2<-data.frame(pow_dCor_Model2,pow_ball_Model2)


ggplot(dat2,aes(n))+geom_point(y=pow_dCor_Model1,fill="white")+geom_line(y=pow_dCor_Model1,colour="red")+geom_point(y=pow_ball_Model1,fill="white")+geom_line(y=pow_ball_Model1,colour="blue")

## ------------------------------------------------------------------------
x<-c(10,20,50,100,1000)
log(exp(x))==exp(log(x))
exp(log(x))==x
log(exp(x))==x

## ------------------------------------------------------------------------
x<-c(10,20,50,100,1000)
isTRUE(all.equal(log(exp(x)), x,exp(log(x)))) 
all.equal(log(exp(x)),x,exp(log(x))) 

## ------------------------------------------------------------------------
#11.5

f = function(a,k){
 
 ck = sqrt(a^2*k/(k+1-a^2))
 l1 = integrate(function(u){(1+u^2/k)^(-(k+1)/2)},0,ck)$value
 l2 = 2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
 l1*l2
}

solve1 = function(k){
  
  output = uniroot(function(a){f(a,k)-f(a,k-1)},lower=1,upper=2)
  output$root
}

k = c(4,25,100,500,1000)
root1 = matrix(0,2,length(k))

for (i in 1:length(k)){
  root1[2,i]=round(solve1(k[i]),4)
}

root1[1,] = k
rownames(root1) = c('k','root')
root1

#11.4

S = function(a,k){
 
 ck = sqrt(a^2*k/(k+1-a^2))
 pt(ck,df=k,lower.tail=FALSE)
}

solve2 = function(k){
  output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
  output$root
}

root2 = matrix(0,2,length(k))

for (i in 1:length(k)){
  root2[2,i]=round(solve2(k[i]),4)
}

root2[1,] = k
rownames(root2) = c('k','A(k)')
root2

## ------------------------------------------------------------------------

nA. <- 28; nB. <- 24; nOO <- 41; nAB <- 70
p <- q <- r <- numeric(100)
p[1] <- 0.2
q[1] <- 0.2
r[1] <- (1- p[1]- q[1])

f <- function(a,b) {
  return((nB.*b/(2-b-2*a)+nB.+nAB)/(nA.*a/(2-a-2*b)+nA.+nAB))
}
g <- function(a,b) {
((1-a/(2-a-2*b))*nA.+(1-b/(2-b-2*a))*nB.+2*nOO)/((nB.*b/(2-b-2*a)+nB.+nAB))
}
threshold <- 1e-5

for (k in 2:100) {
   p[k] <- 1/(1+f(p[k-1],q[k-1])*(1+g(p[k-1],q[k-1])))
   q[k] <- f(p[k-1],q[k-1])/(1+f(p[k-1],q[k-1])*(1+g(p[k-1],q[k-1])))
   r[k] <- 1- p[k] - q[k]
   if((p[k]-p[k-1] <= threshold) & (q[k]-q[k-1] <= threshold) &(r[k]-r[k-1] <= threshold))
   {print(c(k, p[k], q[k],r[k]))
       break
    }
}
x <- seq(1,k,1)
plot(x, p[1:k], "b", col = "green",ylim=c(0,0.7), main = "The log-maximum likelihood values in M-steps" , xlab = "The number of iteration", ylab = "The value of iteration")
lines(x, q[1:k], "b", col = "yellow")
lines(x, r[1:k], "b", col = "red")
legend("topright", legend = c("p", "q", "r"),lty = 1, col = c("green", "yellow", "red"))


## ------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
# lapply
la1 <- lapply(formulas, lm, data = mtcars)
la2 <- lapply(formulas, function(x) lm(formula = x, data = mtcars))
# for loop
lf1 <- vector("list", length(formulas))
for (i in seq_along(formulas)){
lf1[[i]] <- lm(formulas[[i]], data = mtcars)
}


## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
# lapply without anonymous function
la <- lapply(bootstraps, lm, formula = mpg ~ disp)
# for loop
lf <- vector("list", length(bootstraps))
for (i in seq_along(bootstraps)){
lf[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}


## ------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
sapply(la1, rsq)
sapply(la2, rsq)
sapply(lf1, rsq)

## ------------------------------------------------------------------------
sapply(la, rsq)
sapply(lf, rsq)

## ------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
# anonymous function:
sapply(trials, function(x) x[["p.value"]])
# without anonymous function:
sapply(trials, "[[", "p.value")

## ------------------------------------------------------------------------
library(parallel)
mcsapply<-function(x){
  unlist(mclapply(x, sqrt, mc.cores = 3))
}

