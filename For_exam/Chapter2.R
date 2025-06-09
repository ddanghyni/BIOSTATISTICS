library(multtest)
library(nortest)
library(car)
library(outliers)
library(ape)
library(SNPassoc)
library(genetics)

# The Z-test --------------------------------------------------------------

data(golub, package = "multtest")

golubFactor = factor(golub.cl, levels = c(0, 1), labels = c("ALL", "AML"))

x = golub[2058, golubFactor == 'ALL']

n = length(x)

sigma = 0.25; mu0 = 0

z.value = sqrt(n) * (mean(x) - mu0) / sigma
2*pnorm(-abs(z.value))


# One Sample T-test -------------------------------------------------------

Gdf5 = grep('Gdf5', golub.gnames[,2], ignore.case = TRUE)

x = golub(Gdf5, golubFactor = 'ALL')
mu0 = 0
n = 27
t.value = sqrt(n) * (mean(x) - mu0) / sd(x)

t.test(x, mu = 0)
t.test(x, mu = 0, alternative = 'greater')
t.test(x, mu = 0, alternative = 'less')


# One Sample T-test CCND3 gene expression ---------------------------------

ccnd3 = grep('ccnd3', golub.gnames[,2], ignore.case = TRUE)

par(mfrow=c(1,2)) 
stripchart(golub[ccnd3,] ~ golubFactor, method="jitter", cex.lab=1.5, ylab="", vertical = TRUE, col=c("red", "darkgreen"), xlab="Leukemia subtype")

boxplot(golub[ccnd3,] ~ golubFactor, cex.lab=1.5, main=NULL, xlab="Leukemia subtype", col=c("purple","green"), ylab="CCND3 (Cyclin D3) Expression")

ALL = golubFactor == 'ALL'

# ALL
t.test(golub[ccnd3, ALL], mu = 0, alternative = 'greater')

# AML
t.test(golub[ccnd3, !ALL], mu = 0, alternative = 'greater')


# Two Samplpe T-test with Unequal Variances -------------------------------

# welch's t-test

ccnd3 = grep('ccnd3', golub.gnames[,2], ignore.case = TRUE)

golubFactor = factor(golub.cl, levels = c(0, 1), labels = c('ALL', 'AML'))

t.test(golub[ccnd3, ] ~ golubFactor, var.equal = FALSE)

t.test(golub[ccnd3, ] ~ golubFactor, var.equal = FALSE, alternative = 'greater')

t.test(golub[ccnd3, ] ~ golubFactor, var.equal = FALSE, alternative = 'less')


# Two Sample T-test with Equal Variances ----------------------------------

# pooled t-test

golubFactor = factor(golub.cl, levels = c(0, 1), labels = c("ALL", "AML"))
ccnd3 = grep(ccnd3, golub.gnames[, 2], ignore.case = TRUE)
t.test(golub[ccnd3, ] ~ golubFactor, var.equal = TRUE)


# Test for equal variance -------------------------------------------------

# F-test

ccnd3 = grep('ccnd3', golub.gnames[, 2], ignore.case = TRUE)
zyxin = grep('zyxin', golub.gnames[, 2], ignore.case = TRUE)

tapply(golub[ccnd3, ], golubFactor, var) # -> 이놈은 비슷한데..
tapply(golub[zyxin, ], golubFactor, var) # -> 이놈은 확실히 다르쥬,,? 

var.test(golub[ccnd3, ] ~ golubFactor) # H_0 기각 못함.
var.test(golub[zyxin, ] ~ golubFactor) # H_0 기각.

bartlett.test(golub[ccnd3,] ~ golubFactor) 
bartlett.test(golub[zyxin,] ~ golubFactor)

fligner.test(golub[ccnd3,], golubFactor) 
fligner.test(golub[zyxin,], golubFactor)

library(car) 
leveneTest(golub[ccnd3,], golubFactor) 
leveneTest(golub[zyxin,], golubFactor)


# HIstogram and Q-Q plot --------------------------------------------------

par(mfrow=c(1,2)) 
hist(golub[ccnd3, golubFactor=="ALL"], cex.lab=1.5, col="orange",nclass=20, main=NULL, xlab="CCND3 (Cyclin D3) Expression") 
hist(golub[zyxin, golubFactor=="ALL"], cex.lab=1.5, col="orange",nclass=20, main=NULL, xlab="Zyxin Expression")

# q-q plot
par(mfrow=c(1,2)) 
qqnorm(golub[ccnd3, golubFactor=="ALL"], pch=19, cex.lab=1.5, col="red", main=NULL) 
qqline(golub[ccnd3, golubFactor=="ALL"]) 

qqnorm(golub[zyxin, golubFactor=="ALL"], pch=19, cex.lab=1.5,col="red", main=NULL) 
qqline(golub[zyxin, golubFactor=="ALL"])

#'이런 대각선에서의 벗어남은 정규분포와의 **비정규성(non-normality)**을 나타낸다.
#' 예시로 여기선 **양쪽 꼬리(tail)**가 정규분포 꼬리보다 **더 가늘다(skinny)**고 설명.
#' 
#' 꼬리 쪽 점이 수직적 -> 정규분포보다 더 skinny한 tail
#' 꼬리 쪽 점이 수평적 -> 정규분포보다 더 fat한 tail


# Normality test ----------------------------------------------------------

#' Shapiro Wilk test - > based on Q-Q plot
#' Anderson-Darling test -> based on dist of data
#' 각 집단 별로 다 해야한다!!

shapiro.test(golub[ccnd3, golubFactor == 'ALL'])
shapiro.test(golub[ccnd3, golubFactor == 'AML'])

shapiro.test(golub[zyxin, golubFactor == 'ALL']) # -> H_0 기각! -> 정규성 x
shapiro.test(golub[zyxin, golubFactor == 'AML'])

#' Anderson-Darling test

library(nortest)

ad.test(golub[ccnd3, golubFactor == 'ALL'])
ad.test(golub[ccnd3, golubFactor == 'AML'])

ad.test(golub[zyxin, golubFactor == 'ALL']) # -> H_0 기각! -> 정규성 x
ad.test(golub[zyxin, golubFactor == 'AML'])


# Outliers Test -----------------------------------------------------------

library(outliers)

#' 자동으로 설명해주는 alter hypo 잘 보기!

grubbs.test(golub[ccnd3, golubFactor == 'ALL'])
grubbs.test(golub[ccnd3, golubFactor == 'AML'])

grubbs.test(golub[zyxin, golubFactor == 'ALL']) # 이상치가 아니다!!
grubbs.test(golub[zyxin, golubFactor == 'AML'])


golub[ccnd3, golubFactor == 'ALL']


# Binomial Test -----------------------------------------------------------

#' Example : A microRNA of length 22 contains 18 purines.
#' Null hypo - > H_0 : p = 0.7

binom.test(18, 22, p = 0.7, alternative = "greater", conf.level = 0.95)


# Chi-squared test --------------------------------------------------------

library(ape) 
obs <- table(read.GenBank(c("X94991.1"),as.character=TRUE)) 
obs
  
e = rep(sum(obs) / 4, 4) # expected count
test = sum((obs-e)^2/e)

test

1-pchisq(test,3)

qchisq(0.95,3)

chisq.test(obs)
# observed count를 넣으면 됨!! table 자료
# pi를 넣지 않으면 default가 모든 그룹의 확률이 같다!

pi = c(0.75, 0.25)
x = c(5474, 1850)
chisq.test(x, p = pi)


# Testing independence ----------------------------------------------------

data(golub, package = 'multtest')
golubFactor = factor(golub.cl, levels = c(0, 1), labels = c("ALL", "AML"))
gdf5 = grep('gdf5', golub.gnames[, 2], ignore.case = TRUE)

x = golub[gdf5, ]
cutoff = 0.1 # -> 임의의 값
pred = ifelse(x < cutoff, 'ALL', 'AML')

data.frame(pred = pred, true = golubFactor)
table(pred = pred, true = golubFactor)

chisq.test(table(pred = pred, true = golubFactor))
#' p-v = 0.7 -> not reject H_0
#' 이 gene은 좋은 역할을 하지 못한다..!



# 모든 각 express 값을 cutoff로 해보자! --------------------------------------------

cutoff = sort(x)
pval = NULL

for(i in 1:length(cutoff)){
  pred = ifelse(x < cutoff[i], 'ALL', 'AML')
  pval[i] = chisq.test(table(pred, golubFactor))$p.val
}

plot(cutoff, pval, type = "p", pch = 20, ylab = 'pvalues')
cutoff[pval < 0.05]
which(x == cutoff[pval < 0.05])


ccnd3 <- grep("CCND3", golub.gnames[,2], ignore.case=TRUE)
x2 <- golub[ccnd3, ]
cutoff <- sort(x2)
pval2 <- 0
for (i in 1:length(cutoff)) {
  pred2 <- ifelse(x2 < cutoff[i], "ALL", "AML")
  pval2[i] <- chisq.test(table(pred2, golubFactor))$p.val
}
plot(cutoff, pval2, type="p", pch=20, ylab="pvalues")
cutoff[pval2 < 0.05]
boxplot(x2 ~ golubFactor, cex.lab=1.5, main=NULL,
        xlab="Leukemia subtype", col=c("lightblue", "orange"))
abline(h=cutoff[pval2 < 0.05], lty=2, col="gray")


# Fisher's Exact Test -----------------------------------------------------

data = matrix(c(300, 500, 3000, 7000), 2, byrow = TRUE)
fisher.test(data)


# Hardy-Weinberg Equilibruim ----------------------------------------------

library(SNPassoc)
data(asthma, package = 'SNPassoc')

dim(asthma)
lapply(asthma[, -c(1:6)], table)

snp1 = asthma$rs2274276
table(snp1)

# 관측빈도
ObsCount = table(snp1)
Nobs = sum(ObsCount)

# P_G 추정
FreqG = (2 * ObsCount[3] + ObsCount[2]) / (Nobs * 2)

# 기대빈도
ExpCount = c(Nobs*(1-FreqG)^2, 2*Nobs*FreqG*(1-FreqG), Nobs*FreqG^2)

rbind(ObsCount, ExpCount)

ChiSqStat = sum((ObsCount - ExpCount)^2/ExpCount)
pchisq(ChiSqStat, df=2, lower.tail=FALSE) 
# Lower.tail = FALSE = > 오른쪽 꼬리 확률 계산



# HWE.chisq(), HWE.exact() ------------------------------------------------

library(genetics)

# Sample 많
Snp1 = genotype(snp1, sep = "")
table(Snp1)
summary(Snp1)
HWE.chisq(Snp1)

# Sample 적
country = asthma$country
table(country)

Snp1bg = genotype(snp1[country == 'Belgium'], sep = "")
summary(Snp1bg)
HWE.chisq(Snp1bg)
HWE.exact(Snp1bg) # => sample이 적을때는 exact test !!!


# Permutation Test : Two-sample Test --------------------------------------

set.seed(123)
x = c(1, 3, 4)
y = c(2, 2, 1, 2)
sample(c(x, y))


# Permutation Test --------------------------------------------------------

library(multtest)
data(golub)

golubFactor = factor(golub.cl, levels = c(0, 1), labels = c('ALL', 'AML'))
golubFactor

sample(golubFactor)
sample(as.numeric(golubFactor))

set.seed(12345)
sample(as.numeric(golubFactor))

# original data
y = as.numeric(golubFactor)

# 반복수
k = 10000

# 실험 결과를 저장할 메트릭스
# 각 열이 각 permutation이다.
mat.y = matrix(y, length(y), k)
mat.y[, 1:4]

# 각 열에 대해서 sampling 
per.y = apply(mat.y, 2, sample)
per.y[, 1:4]

# 우리가 실험할 gene 선택
ccnd3 = grep('ccnd3', golub.gnames[,2], ignore.case = TRUE)

# two sample t test 진행 -> 검정통계량 값만 보면 된다잉~
t.test(golub[ccnd3, ] ~ per.y[,1], var.equal=FALSE)
t.test(golub[ccnd3, ] ~ per.y[,1], var.equal=FALSE)$stat

# 기존값
tobs = t.test(golub[ccnd3, ] ~ per.y[, 1], var.equal = FALSE)$stat
tobs # -> 차이가 거의 안나는걸 알 수 있다.

# 이제 10000번 반복해보자..
fun = function(t) t.test(golub[ccnd3, ] ~ t, var.equal = FALSE)$stat
T = apply(per.y, 2, fun)

mean(abs(T) > abs(tobs))
#'→ 실제보다 더 극단적인 t-값이 순열에서 얼마나 나오는지를 비율로 계산
#'→ 즉, Permutation 기반의 p-value다.
(sum(abs(T) > abs(tobs))+1)/(k+1)

#' → 순열된 그룹 중 46.2%는 실제보다 더 극단적인 t값을 만들었단 얘기
#' → 그러면 실제 관측된 tobs는 전혀 특별하지 않음
#' → 그러니까 "H₀ 하에서도 이런 값은 흔하다"
#' ⇒ H₀를 기각할 근거가 없다

hist(T, breaks=100, col="orange", main="", xlim=c(-7, 7), xlab="Null Distribution of Test Statistic")
x0 <- seq(-7, 7, len=1000) 
y0 <- dnorm(seq(-7, 7, len=1000)) 
lines(x0, y0*1000, col=2, lwd=3) 
abline(v=-abs(tobs), col=4, lty=2, lwd=2) 
abline(v=abs(tobs), col=4, lty=2, lwd=2) 
text(abs(tobs)-1, 350, col=4, paste("|T| = ", round(abs(tobs), 4), sep=""))


# example2 ----------------------------------------------------------------

gdf5 <- grep("Gdf5", golub.gnames[,2], ignore.case=TRUE) 
fun2 <- function(t) t.test(golub[gdf5,]~t, var.equal=FALSE)$stat 
T2 <- apply(per.y, 2, fun2) 
tobs2 <- t.test(golub[gdf5,] ~ golubFactor, var.equal=FALSE)$stat 
(sum(abs(T2) > abs(tobs2))+1)/(k+1)

hist(T2, breaks=100, col="orange", main="", xlab="Null Distribution of Test Statistic") 
x0 <- seq(min(T2), max(T2), len=1000) 
y0 <- dnorm(seq(min(T2), max(T2), len=1000)) 
lines(x0, y0*1000, col=2, lwd=3) 
abline(v=-abs(tobs2), col=4, lty=2, lwd=2) 
abline(v=abs(tobs2), col=4, lty=2, lwd=2) 
text(abs(tobs2)+1, 350, col=4, paste("|T| = ", round(abs(tobs2), 4), sep=""))



# Confidence Interval -----------------------------------------------------



































































