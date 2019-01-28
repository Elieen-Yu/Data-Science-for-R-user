rm(list = ls())  # clear the memory
###########################

########################################
########################################

## Matched pairs:

X = as.matrix(read.table("datasets/T6-1.dat"))

# Columns represent lab 1 (c1, c2) and lab 2(c3, c4)
diff = cbind(X[,1] - X[,3], X[,2] - X[,4])
colnames(diff) = c("BOD", "SS") # biochemical oxyden demand, suspended solids
dbar = colMeans(diff)
Sd = cov(diff)

dbar
Sd

n = nrow(diff)
p = ncol(diff)

########################################
# Bonferroni intervals on the two marginal means; 95%
alpha = .05
m = 2
alpham = alpha/m
int = matrix(nrow = m, ncol = 2, data = NA)
for(i in 1:m) {
psihat = dbar[i]
sd_psihat = sqrt(Sd[i,i]/n)
qT = qt(alpham/2, n-1, lower.tail = 0)
int[i,] = c(psihat - qT*sd_psihat, psihat + qT*sd_psihat)
}
int_Bon = cbind(int[,1], dbar, int[,2])
colnames(int_Bon) = c("lower", "dbar", "upper")

cat("Bonferroni intervals on the two marginal means","\n")
int_Bon

########################################
# Simultaneous intervals on these AND ALL OTHERS
alpha = .05
m = 2
int = matrix(nrow = m, ncol = 2, data = NA)
for(i in 1:m) {
I = diag(m)
a = I[,i]  # Any other vectors a could be used without lowering the overall confidence
psihat = t(a)%*%dbar
sd_psihat = sqrt(t(a)%*%Sd%*%a/n)
qF = sqrt(((n-1)*p/(n-p))*qf(alpha, p, n-p, lower.tail = 0))

int[i,] = c(psihat - qF*sd_psihat, psihat + qF*sd_psihat)
}
int_Scheffe = cbind(int[,1], dbar, int[,2])
colnames(int_Scheffe) = c("lower", "dbar", "upper")
int_Scheffe


########################################

# 95% confidence ellipsoid:
U = chol(Sd) # S = U'U
c = sqrt(((n-1)*p/(n-p))*qf(alpha, p, n-p, lower.tail = 0))

phi = 2*pi*seq(from = 0, to = 1, length = 301)

z = -rbind(cos(phi), sin(phi))  # 200 columns, each of norm 1
mu = dbar + (c/sqrt(n))*t(U)%*%z
plot(mu[1,], mu[2,], xlab = "mean BOD", ylab = "mean SS")
points(0,0, pch=22)
abline(h=0)
abline(v=0)
########################################

extremes_of_ellipsoid = rbind(c(min(mu[1,]), max(mu[1,])),c(min(mu[2,]), max(mu[2,])))
extremes_of_ellipsoid


########################################
########################################

## Repeated measures:


X = as.matrix(read.table("datasets/T6-2.dat",header=F))	
C = rbind(c(-1, -1, 1, 1), c(1, -1, 1, -1), c(1, -1, -1, 1))
Y = X%*%t(C)

colnames(Y) = c("H effect", "C effect", "interaction") 
ybar = colMeans(Y)
Sy = cov(Y)

ybar
Sy

n = nrow(Y)
q = ncol(Y)

########################################
# Test that all three contrasts vanish, i.e. mean(y) = c(0,0,0):

Tsqd = n*t(ybar)%*%solve(Sy,ybar)
pval = pf(((n-1)*q/(n-q))*Tsqd, q, n-q, lower.tail = 0)
pval

# Simultaneous intervals on these AND ALL OTHERS
alpha = .05
m = 3
int = matrix(nrow = m, ncol = 2, data = NA)
for(i in 1:m) {
I = diag(m)
a = I[,i]  # Any other vectors a could be used without lowering the overall confidence
psihat = t(a)%*%ybar
sd_psihat = sqrt(t(a)%*%Sy%*%a/n)
qF = sqrt(((n-1)*q/(n-q))*qf(alpha, q, n-q, lower.tail = 0))

int[i,] = c(psihat - qF*sd_psihat, psihat + qF*sd_psihat)
}
int_Scheffe = cbind(int[,1], ybar, int[,2])
colnames(int_Scheffe) = c("lower", "ybar", "upper")
int_Scheffe



