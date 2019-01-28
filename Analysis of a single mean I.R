rm(list = ls())  # clear the memory
###########################

X = as.matrix(read.table("datasets/T5-2.DAT"))

n = nrow(X)
p = ncol(X)

colnames(X) = c("SS&H", "Vb", "Sc")
xbar = colMeans(X)
S = cov(X)
S

########################################
# Bonferroni intervals on the three marginal means; 95%
alpha = .05
m = 3
alpham = alpha/m
int = matrix(nrow = m, ncol = 2, data = NA)
for(i in 1:m) {
psihat = xbar[i]
sd_psihat = sqrt(S[i,i]/n)
qT = qt(alpham/2, n-1, lower.tail = 0)
int[i,] = c(psihat - qT*sd_psihat, psihat + qT*sd_psihat)
}
int_Bon = cbind(int[,1], xbar, int[,2])
colnames(int_Bon) = c("lower", "xbar", "upper")

cat("Bonferroni intervals on the three marginal means","\n")
int_Bon

########################################
# Simultaneous intervals on these AND ALL OTHERS
alpha = .05
m = 3
int = matrix(nrow = m, ncol = 2, data = NA)
for(i in 1:m) {
I = diag(m)
a = I[,i]  # Any other vectors a could be used without lowering the overall confidence
psihat = t(a)%*%xbar
sd_psihat = sqrt(t(a)%*%S%*%a/n)
qF = sqrt(((n-1)*p/(n-p))*qf(alpha, p, n-p, lower.tail = 0))

int[i,] = c(psihat - qF*sd_psihat, psihat + qF*sd_psihat)
}
int_Scheffe = cbind(int[,1], xbar, int[,2])
colnames(int_Scheffe) = c("lower", "xbar", "upper")
int_Scheffe

# Values in text not very accurate - lots of rounding

########################################

# 95% confidence ellipsoid on means of VB and Sc:
X23 = X[, 2:3]
n = nrow(X23)
p = ncol(X23)
alpha = .05
xbar23 = colMeans(X23)
S23 = cov(X23)
U = chol(S23) # S23 = U'U
c = sqrt(((n-1)*p/(n-p))*qf(alpha, p, n-p, lower.tail = 0))

phi = 2*pi*seq(from = 0, to = 1, length = 201) # 2*pi*(0,1/200, 2/200, ... 199/200, 1)

z = rbind(cos(phi), sin(phi))  # 201 columns, each of norm 1
mu = xbar23 + (c/sqrt(n))*t(U)%*%z
dev.new()
plot(mu[1,], mu[2,], xlab = "mean Vb", ylab = "mean Sc", type = 'l')

########################################

# Simultaneous intervals on these two means, and all other comparisons of these two means only; 
# should agree with the extremes of the confidence ellipsoid:
alpha = .05
m = 2
int = matrix(nrow = m, ncol = 2, data = NA)
for(i in 1:m) {
I = diag(m)
a = I[,i]
psihat = t(a)%*%xbar23
sd_psihat = sqrt(t(a)%*%S23%*%a/n)
qF = sqrt(((n-1)*p/(n-p))*qf(alpha, p, n-p, lower.tail = 0))

int[i,] = c(psihat - qF*sd_psihat, psihat + qF*sd_psihat)
}
int2_Scheffe = cbind(int[,1], xbar23, int[,2])
colnames(int2_Scheffe) = c("lower", "xbar", "upper")
int2_Scheffe
abline(v=int[1,]) # vertical lines at the endpoints of the interval on mean(Vb)
abline(h=int[2,]) # horizontal lines at the endpoints of the interval on mean(Sc)

extremes_of_ellipsoid = rbind(c(min(mu[1,]), max(mu[1,])),c(min(mu[2,]), max(mu[2,])))
extremes_of_ellipsoid


########################################

# 95% prediction region for a new (VB,Sc):
c = sqrt(((n-1)*p/(n-p))*qf(alpha, p, n-p, lower.tail = 0))
c = c*sqrt((n+1)/n)

phi = 2*pi*seq(from = 0, to = 1, length = 201) # 2*pi*(0,1/200, 2/200, ... 199/200, 1)

z = rbind(cos(phi), sin(phi))  # 201 columns, each of norm 1
mu = xbar23 + (c/sqrt(n))*t(U)%*%z
dev.new()
plot(mu[1,], mu[2,], xlab = "Vb", ylab = "Sc", type = 'l')


## Superimpose the previous confidence ellipsoid on the means:
c = sqrt(((n-1)*p/(n-p))*qf(alpha, p, n-p, lower.tail = 0))
phi = 2*pi*seq(from = 0, to = 1, length = 201) # 2*pi*(0,1/200, 2/200, ... 199/200, 1)
mu = xbar23 + (c/sqrt(n))*t(U)%*%z
#dev.new()
lines(mu[1,], mu[2,], type = 'l')


########################################
## Test mu = c(500, 50, 25) using Tsqd and the large sample test
p = 3
mu0 = c(500, 50, 25)
Tsqd = n*t(xbar-mu0)%*%solve(S,xbar-mu0)
Fcalc = (n-p)*Tsqd/(p*(n-1))
pval = pf(Fcalc, p, n-p, lower.tail = 0) # prob(F > Fcalc)

Lambda = (1+Tsqd/(n-1))^(-n/2)
chisqd.calc = -2*log(Lambda)
pval.approx = pchisq(chisqd.calc, p, lower.tail = 0)

pval
pval.approx