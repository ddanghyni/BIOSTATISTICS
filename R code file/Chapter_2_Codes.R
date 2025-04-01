
# packages ----------------------------------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("multtest")
#install.packages("nortest")
#install.packages("car")
#install.packages("outliers")
#install.packages("ape")
#install.packages("SNPassoc")
#install.packages("genetics")

library(multtest)
library(nortest)
library(car)
library(outliers)
library(ape)
library(SNPassoc)
library(genetics)


# Z-test ------------------------------------------------------------------

library(multtest)

data(golub, package = "multtest")

golubFactor <- factor(golub.cl, levels=0:1, labels=c("ALL","AML"))

x <- golub[2058, golubFactor=="ALL"]

n <- length(x)

sigma <- 0.25; mu0 <- 0
z.value <- sqrt(n)*(mean(x) - mu0)/sigma

2*pnorm(-abs(z.value)) # p-value


# Confidence Interval -----------------------------------------------------

mean(x) + qnorm(0.025) * sigma/sqrt(n)
mean(x) + qnorm(0.975) * sigma/sqrt(n)
mean(x) + c(-1, 1) * qnorm(0.975) * sigma / sqrt(n)


# One Sample T-test -------------------------------------------------------

#'모집단의 분산(표준편차)을 모르고, 표본이 정규성을 만족한다고 가정할 때 사용한다.

#'이때 검정 통계량은 t 분포를 따르고 t분포는 degree of freedom에 따라 모양이 달라진다.

#'모집단의 분산 대신 표본의 구한 표본분산을 이용하기 때문에 정규분포가 아닌 t 분포를 따르게 된다.(자세한 증명은 수통 때)


t.test(x, mu=0) # default two side test

#' mu를 기준으로 p-v는 오른쪽 부분 싸악
t.test(x, mu=0, alternative = "greater") # one side test

#' mu를 기준으로 p-v는 왼쪽 부분 싸악 
t.test(x, mu=0, alternative = "less") # one side test

# gred()으로 CCND3 gene의 index 찾기
ccnd3 <- grep("CCND3", golub.gnames[,2], ignore.case=TRUE) 
ccnd3

ALL <- golubFactor=="ALL" 

#' CCND3 gene에서 ALL로 검정을 하자..
#' H_0 = mu에 대한 가설검정 
t.test(golub[ccnd3, ALL], mu=0) 
t.test(golub[ccnd3, !ALL], mu=0)


# Two Sample T-test with Unequal Variances --------------------------------

#두 그룹의 평균을 비교할 때 사용한다.
#Two Sample T test를 진행하기 위해선
#1. 각 그룹 샘플에 대한 정규성 검정 실시
#2.Equal variance에 대한 검정 실시
#3.분산이 같다면 pooled t-test를 실시
#4.같지 않다면 welch's t-test를 실시(general version)
#5.그 다음 상황에 맞는 t test를 진행한다.

# t.test(data ~ group(factor), var.equal = ...)
t.test(golub[ccnd3,] ~ golubFactor, var.equal=FALSE) # result = > 두 집단은 다르다~



# Two Sample T-test with Equal Variances ----------------------------------

t.test(golub[ccnd3,] ~ golubFactor, var.equal=TRUE)



# Test for Equal Variances ------------------------------------------------

var.test(golub[ccnd3, ] ~ golubFactor) # 두 그룹의 분산이 같다고 할 수 있다->pooled 이용
var.test(golub[zyxin, ] ~ golubFactor) # 두 그룹의 부산이 같다고 할 수 없다->pooled 못 쓴다.. welch's 이용

bartlett.test(golub[ccnd3,] ~ golubFactor) 
bartlett.test(golub[zyxin,] ~ golubFactor)

fligner.test(golub[ccnd3,], golubFactor) 
fligner.test(golub[zyxin,], golubFactor)

leveneTest(golub[ccnd3,], golubFactor) 
leveneTest(golub[zyxin,], golubFactor)


# Q-Q plot ----------------------------------------------------------------

par(mfrow=c(1,2)) 
qqnorm(golub[ccnd3, golubFactor=="ALL"], pch=19, cex.lab=1.5, col="red", main=NULL) 
qqline(golub[ccnd3, golubFactor=="ALL"]) 

qqnorm(golub[zyxin, golubFactor=="ALL"], pch=19, cex.lab=1.5,col="red", main=NULL) 
qqline(golub[zyxin, golubFactor=="ALL"])



# Normality Tests ---------------------------------------------------------

# H_0 : \text{정규성을 따른다..!}
# Normality Tests는 각 집단 별로 다해야한다!!! -> 베리 중요

shapiro.test(golub[ccnd3, golubFactor=="ALL"]) 
shapiro.test(golub[ccnd3, golubFactor=="AML"]) 
shapiro.test(golub[zyxin, golubFactor=="ALL"]) # reject H_0
shapiro.test(golub[zyxin, golubFactor=="AML"]) 


# Outliers Test -----------------------------------------------------------

# H_0 : \text{Does not contain an outlier}

library(outliers)
grubbs.test(golub[ccnd3, golubFactor=="ALL"])
grubbs.test(golub[zyxin, golubFactor=="ALL"])


# Binomial Test -----------------------------------------------------------

# 이항분포 기반 가설검정은 특정 사건이 일어날 확률 $p$가 이론적으로 가정된 값 $p_0$와 같은지를 통계적으로 검정

# 예 : 어떤 microRNA 서열에서 퓨린이 나타날 확률을 $p_0$라고 가정하고, 실제로 관측된 퓨린 비율이 $p_0$와 다른지 판단.

binom.test(18, 22, p=0.7, alternative="greater", conf.level=0.95)


# R에서 분포와 관련된 함수 정리 -------------------------------------------------------

#' d : 확률밀도함수
#' 특정 분포의 확률밀도값 또는 확률질량값을 계산
dnorm(0, mean = 0, sd = 1)

#' p : 누적분포함수
#' 특정 분포에서 어떤 값 q 이하가 나올 확률을 계산
pnorm(0, mean = 0, sd = 1)

#' q : 분위수함수
#' 누적분포함수의 반대 개념으로, 누적확률 p에 해당하는 값을 계산
qnorm(0.975, mean = 0, sd = 1)
qnorm(0.5, mean = 0, sd = 1)

#' r : 난수생성
#' 특정 분포를 따르는 난수를 n개 생성하여 벡터로 반환
rnorm(10, mean = 0, sd = 1)



# Chi-squared Test --------------------------------------------------------

pi <- c(0.75,0.25) 
x <- c(5474, 1850) 
chisq.test(x, p=pi)



# Testing Indepnedence ----------------------------------------------------

#1. Cutoff 설정: 유전자 발현값이나 어떤 점수 등에 임의의 임계값을 정해, 이를 기준으로 예측

#2. Confusion matrix 구성

#3. 독립성 검정 수행

# H_0 : \text{cutoff 쓰레기, 실제 레이블과 예측 레이블이 독립이다.}

data(golub, package = "multtest") 
golubFactor <- factor(golub.cl, levels=0:1, labels=c("ALL","AML")) 
gdf5 <- grep("Gdf5", golub.gnames[,2], ignore.case=TRUE)

x <- golub[gdf5, ] 

# cutoff 설정
cutoff <- 0.1 

# cutoff에 따른 예측
pred <- ifelse(x < cutoff, "ALL", "AML") 
data.frame(predicted=pred, true=golubFactor) 
table(predicted=pred, true=golubFactor)

chisq.test(table(pred, golubFactor))

################################################################################
cutoff <- sort(x) 
pval <- 0

for (i in 1:length(cutoff)) { 
  pred <- ifelse(x < cutoff[i], "ALL", "AML") 
  pval[i] <- chisq.test(table(pred, golubFactor))$p.val 
}

plot(cutoff, pval, type="p", pch=20, ylab="pvalues") 
cutoff[pval < 0.05]


# Fisher’s Exact Test -----------------------------------------------------

fisher.test(data)


# Hardy Weinberg Equilibrium ----------------------------------------------

table(Geno)[1] # aa 
table(Geno)[2] # Aa

MAF = (table(Geno)[1] * 2 + table(Geno)[2] * 1) / (n * 2) # P_a
round(MAF, 4)

#' Chi-squ test gogo
table(Geno)

# expect count
Exp = c(n * MAF^2, 2 * n * MAF * (1 - MAF), n * (1 -MAF)^2) # Aa => Aa, aA -> 2배

# Chi-sqr test
sum((table(Geno) -  Exp)^2 / Exp )
round(sum((table(Geno) -  Exp)^2 / Exp ), 4)

t = c(MAF^2, 2*MAF*(1-MAF), (1-MAF)^2)
chisq.test(table(Geno), p = t)

pi = c(F^2, 2*F*(1-F), (1-F)^2)
chisq.test(table(Geno), p = pi)

table(Geno)


################################################################################

library(genetics)

# Geno 벡터를 유전자형 형태로 바꿔주기
# 0 = aa, 1 = Aa, 2 = AA
geno_vec = rep(NA, length(Geno))
geno_vec[Geno == 0] = "aa"
geno_vec[Geno == 1] = "Aa"
geno_vec[Geno == 2] = "AA"

geno_factor = genotype(geno_vec, sep = "")

# Fisher's exact HWE test (genetics 패키지)
HWE.exact(geno_factor)
HWE.chisq(geno_factor)




# Permutation Test --------------------------------------------------------

library(multtest) 
data(golub)

golubFactor <- factor(golub.cl, levels=0:1, labels=c("ALL","AML")) 
golubFactor

sample(golubFactor) 
sample(as.numeric(golubFactor))  # resampling...
sample(as.numeric(golubFactor))

set.seed(12345) 
sample(as.numeric(golubFactor))


y <- as.numeric(golubFactor)  # origin values
K <- 10000 # 10000번 rsampling 하자!
mat.y <- matrix(y, length(y), K) # 행은 y 즉 기존 데이터 수 (n+m), 열은 반복하는 만큼 ㄱㄱ
mat.y[,1:5]

per.y <- apply(mat.y, 2, sample) 
per.y[,1:5]

ccnd3 <- grep("CCND3", golub.gnames[,2], ignore.case=TRUE)
t.test(golub[ccnd3,] ~ per.y[,1], var.equal=FALSE) 
# per.y[,1] -> 첫 번째 resample로 했을때의 검정 통계량임..!
t.test(golub[ccnd3,] ~ per.y[,1], var.equal=FALSE)$stat

tobs <- t.test(golub[ccnd3,] ~ golubFactor, var.equal=FALSE)$stat
tobs


# 이제 검정 통계량 10000개 구해보자..!
fun <- function(t) t.test(golub[ccnd3,]~t, var.equal=FALSE)$stat 
T <- apply(per.y, 2, fun) 
T[1:100]

mean(abs(T) > abs(tobs))
(sum(abs(T) > abs(tobs))+1)/(K+1) # p-value
1/(K+1)


gdf5 <- grep("Gdf5", golub.gnames[,2], ignore.case=TRUE) 
fun2 <- function(t) t.test(golub[gdf5,]~t, var.equal=FALSE)$stat 
T2 <- apply(per.y, 2, fun2) 
tobs2 <- t.test(golub[gdf5,] ~ golubFactor, var.equal=FALSE)$stat 
(sum(abs(T2) > abs(tobs2))+1)/(K+1)



# Bootstrap Confidence Interval -------------------------------------------


K <- 10000
mat <- matrix(golub[ccnd3,], length(golub[ccnd3,]), K) 
mat[1:5, 1:5]
fun3 <- function(t) sample(t, replace=TRUE) 
boot <- apply(mat, 2, fun3) 
bmean <- apply(boot, 2, mean) 
quantile(bmean, c(0.025, 0.975)) # 95% 신뢰구간! 



# Wilcoxon Signed Rank Test -----------------------------------------------

# one sample t-test의 비모수 버전

x <- c(6003, 6304, 6478, 6245, 6134, 6204, 6150) 
wilcox.test(x, mu=6000)

nkr <- grep("Nkr", golub.gnames[, 2], ignore.case=TRUE) 
shapiro.test(golub[nkr, golubFactor=="ALL"]) # 정규성 안따르노 ㅋㅎ
shapiro.test(golub[nkr, golubFactor=="AML"]) # 정규성 따른다.



# Bootstrap Confidence Interval in nonparametric test. --------------------

set.seed(12345) 
K <- 100000 
nkr.ALL <- golub[nkr, golubFactor=="ALL"] 
nkr.ALL 
sample(nkr.ALL, replace=TRUE) 
mean(sample(nkr.ALL, replace=TRUE))

mat <- matrix(nkr.ALL, length(nkr.ALL), K) 
fun3 <- function(t) sample(t, replace=TRUE) 
boot <- apply(mat, 2, fun3) 
bmean <- apply(boot, 2, mean) # bootstrap 표본평균 계산
quantile(bmean, c(0.025, 0.975))


# Wilcoxon Rank Sum Test --------------------------------------------------

# two sample t test의 비모수 검정 방법.

igf <- grep("IGFBP5",golub.gnames[,2], ignore.case = TRUE) 
shapiro.test(golub[igf, golubFactor=="ALL"]) 
shapiro.test(golub[igf, golubFactor=="AML"]) 
wilcox.test(golub[igf,] ~ golubFactor)


# Applications to Multiple Genes ------------------------------------------

# 정규성 통과한 유전자 비율 예제

dim(golub) 
fun <- function(t) shapiro.test(t)$p.value

all <- apply(golub[, golubFactor=="ALL"], 1, fun)
aml <- apply(golub[, golubFactor=="AML"], 1, fun)

sum(all > 0.05)/nrow(golub) 
sum(aml > 0.05)/nrow(golub)

pval.t <- function(t) t.test(t ~ golubFactor)$p.value 
pval.w <- function(t) wilcox.test(t ~ golubFactor)$p.value

pt <- apply(golub, 1, pval.t) 
pw <- apply(golub, 1, pval.w)

pval <- data.frame(cbind(pt, pw)) 

# t-test는 유의(p<0.05)인데 wilcox는 유의하지 않은 유전자 수
nrow(pval[pt < 0.05 & pw >= 0.05, ]) # 서로 반대의 결과를 보여주는..

# wilcox는 유의(p<0.05)인데 t-test는 유의하지 않은 유전자 수
nrow(pval[pw < 0.05 & pt >= 0.05, ]) 

# 각각의 방법에서 p-value가 가장 작은 유전자의 인덱스 찾기
apply(pval, 2, function(x) which(x==min(x))) # 2124이 가장 signal 쌔다..?

# 뭐 나중에는 ranking 순으로 gene별 중요도 확인도 한대.






































