
# -------------------------------------------------------------------------

set.seed(111222) 
n <- 20; k <- 4; x <- NULL
m <- runif(k, -1, 1)
for (i in 1:k) x <- c(x, rnorm(n, m[i], 2))
fac <- gl(k, n)

#1. SSW
g = lm(x ~ fac)
anova(g)$Sum[1] #SSB
anova(g)$Sum[2] #SSW  F = SSB / SSW

#2. 
shapiro.test(residuals(g))

#3.
fac2 = relevel(fac, ref = 3)
g2 = lm(x~fac2)
summary(g2)$coef[4, 3]

#4.
g3 = TukeyHSD(aov(x~fac))
g4 = pairwise.t.test(x, fac, p.adjust.method = "bonferroni")
g5 = pairwise.t.test(x, fac, p.adjust.method = "holm")
mean(c(g3$fac[6, 4], g4$p.value[3,3], g5$p.value[3,3]))



# -------------------------------------------------------------------------

data(golub, package="multtest") 
zyxin <- grep("Zyxin", golub.gnames[,2], ignore.case=TRUE) 
cmyb <- grep("c-myb", golub.gnames[,2], ignore.case=TRUE) 
x <- golub[zyxin, ] 
y <- golub[cmyb, ] 

g <- lm(y ~ x) 
g$coef
coef(g)
summary(g)

# confidence interval -----------------------------------------------------
confint(g, level=0.90)

# anova, ssr, sst ---------------------------------------------------------
n <- length(x)
p <- 1
df1 <- n - p - 1
df2 <- g$df
df3 <- df.residual(g)
c(df1, df2, df3)
summary(g)$sigma^2 # MSE, summary(g)$sigma -> sqrt(mse)
sum(g$res^2)/g$df
deviance(g)/df.residual(g) # deviance -> SSE
anova(g)
anova(g)$Sum # SSR, SSE 각각 
sum(anova(g)$Sum) # SST
names(summary(g))


# R^2 ---------------------------------------------------------------------
summary(g)$r.squared
summary(g)$adj.r.squared # adj-r^2

# F-stat ------------------------------------------------------------------
names(summary(g))
summary(g)$fstat


# -------------------------------------------------------------------------

x1 <- golub[grep("Gdf5", golub.gnames[,2], ignore.case=TRUE),]
x2 <- golub[grep("CCND3", golub.gnames[,2], ignore.case=TRUE),]
x3 <- golub[grep("Nkr", golub.gnames[, 2], ignore.case=TRUE),]
x4 <- golub[grep("IGFBP5",golub.gnames[,2], ignore.case = TRUE),]
data <- data.frame(Gdf5=x1, CCND3=x2, Nkr=x3, IGFBP5=x4, cmyb=y)
data

g2 <- lm(cmyb ~ ., data) 
summary(g2)
summary(g2)$coef
confint(g2) # confidence interval of coefficients

n <- length(x)
p <- 4
df1 <- n - p - 1
df2 <- g2$df
df3 <- df.residual(g2) # degree of freedom of residuals
c(df1, df2, df3)

summary(g)$sigma
summary(g)$sigma^2 # MSE
# -------------------------------------------------------------------------

g0 <- lm(cmyb ~ 1, data) # intercept model (udner H_0)
summary(g0)

anova(g0, g2) # compare intercept model and full model
# > if H_0 is not rejected => intercept model is goood

sum(anova(g2)$F, na.rm = TRUE)/4

cor(data)
SLR <- vector("list", length=3)
names(SLR) <- c("bhat", "rsq", "sigma")
for (k in 1:(ncol(data)-1)) {
  g0 <- lm(data$cmyb ~ data[,k])
  SLR$bhat[k] <- coef(g0)[2] # \hat\beta_1
  SLR$rsq[k] <- summary(g0)$r.squared
  SLR$sigma[k] <- summary(g0)$sigma
}

SLR
data.frame(MLR=coef(g2)[-1], SLR=SLR$bhat)
data.frame(MLR=summary(g2)$r.squared, SLR=SLR$rsq)
data.frame(MLR=summary(g2)$sigma, SLR=SLR$sigma)

x0 <- data.frame(Gdf5=median(x1), CCND3=median(x2),
                 Nkr=median(x3), IGFBP5=median(x4))
x0
predict(g2, x0)

# 각 simple model에서의 예측값..
pred <- NULL
for (k in 1:length(x0)) {
  g0 <- lm(data$cmyb ~ data[,k])
  pred[k] <- sum(coef(g0)[1] + x0[k]*coef(g0)[2])
}
pred


# -------------------------------------------------------------------------

y <- c(2, 3, 1, 2, 8, 7, 9, 8, 11, 12, 13, 12)
factor <- gl(3, 4) # 4개 씩 3 그룹?
factor

data.frame(x = factor, y = y)
model.matrix(y ~ factor - 1) # no intercept term


# -------------------------------------------------------------------------

set.seed(123)
y <- c(rnorm(10, -1, 1), rnorm(10, 0, 1), rnorm(10, 1, 1))
y
factor <- gl(3, 10)
factor

model.matrix(y ~ factor - 1)
model.matrix(y ~ factor)


# ref ---------------------------------------------------------------------
factor <- relevel(factor, ref=2)
model.matrix(y ~ factor)
g3 <- lm(y ~ factor)




# All pairwise test -------------------------------------------------------

pairwise.t.test(y, factor, p.adjust.method="bonferroni")
pairwise.t.test(y, factor, p.adjust.method="holm")
TukeyHSD(aov(y ~ factor))
plot(TukeyHSD(aov(y ~ factor)))

set.seed(1234)
y <- rnorm(30, 1.9, 1)
y
factor <- gl(3, 10)
g <- lm(y ~ factor)
anova(g)
pairwise.t.test(y, factor, p.adjust.method="bonferroni")
TukeyHSD(aov(g))
plot(TukeyHSD(aov(g)))


# Multiple test in anova --------------------------------------------------

dim(ex[, B1B2B3])
fun <- function(t) anova(lm(t ~ factor))$Pr[1]
anova.pValues <- apply(ex[, B1B2B3], 1, fun)

pBF <- p.adjust(anova.pValues, method="bonferroni")
pHO <- p.adjust(anova.pValues, method="holm")
pBH <- p.adjust(anova.pValues, method="BH")

alpha <- 0.05
c(sum(pBF < alpha), sum(pHO < alpha), sum(pBH < alpha))


# Checking assumptions ----------------------------------------------------

#.1 정규성 검정
library(ALL)
data(ALL)

B1B2B3 <- ALL$BT %in% c("B1","B2","B3")
ex <- exprs(ALL)
y <- as.numeric(ex[row.names(ex)=="1866_g_at", B1B2B3])
factor <- factor(ALL$BT[B1B2B3], labels=c("B1","B2","B3"))

# 정규성 검정
# H_0 is normality okay
res <- residuals(lm(y ~ factor))
shapiro.test(res) # pv < 0.05 이므로 아노바 ㄴㄴ

#2. 등분산 검정
library(lmtest)
bptest(lm(y ~ factor), studentize=FALSE)
bartlett.test(y ~ factor)
fligner.test(y, factor)
library(car)
leveneTest(y, factor)


# Robust test -------------------------------------------------------------

#,잔차가 가정을 따르지 않을 때...

#, When only homoscedasticity is violated, Welch's ANOVA -> 정규성은 만족해야함!!!

#, When normality is violated, the Kruskal-Wallis rank sum test -> 정규성 x 일 때!!!

#1. 등분산 x 일때!
# Welch's ANOVA -> Games Howell test
oneway.test(y ~ factor, var.equal = FALSE)

#2. 정규성 x 일 때!
kruskal.test(y~factor)

f2 <- function(t) kruskal.test(t ~ factor)$p.value
kruskal.pValues <- apply(ex[, B1B2B3], 1, f2)
pBF <- p.adjust(kruskal.pValues, method="bonferroni")
pHO <- p.adjust(kruskal.pValues, method="holm")
pBH <- p.adjust(kruskal.pValues, method="BH")
alpha <- 0.05
c(sum(pBF < alpha), sum(pHO < alpha), sum(pBH < alpha))



# Two way ANOVA -----------------------------------------------------------

library(ALL)
data(ALL)
ex <- exprs(ALL)
dim(ex)
table(ALL$BT)
table(ALL$mol.biol)

w1 <- ALL$BT %in% c("B", "B1", "B2", "B3", "B4")
w2 <- ALL$mol.biol %in% c("BCR/ABL", "NEG")
ex12 <- ex["32069_at", w1 & w2]
length(ex12)
facB <- ceiling(as.integer(ALL$BT[w1 & w2])/3)
facB
data.frame(B.type=ALL$BT[w1 & w2], facB=facB)
fac1 <- factor(facB, levels=1:2, labels=c("B012", "B34"))
fac2 <- factor(ALL$mol.biol[w1 & w2])
table(fac1)
table(fac2)
tapply(ex12, fac1, mean)
tapply(ex12, fac2, mean)


anova(lm(ex12 ~ fac1 * fac2))
summary(lm(ex12 ~ fac1 * fac2))

data.frame(fac1, fac2)
nfac <- rep(1, length(fac1))
nfac[fac1=="B012" & fac2=="NEG"] <- 2
nfac[fac1=="B34" & fac2=="BCR/ABL"] <- 3
nfac[fac1=="B34" & fac2=="NEG"] <- 4
nfac <- factor(nfac, labels=c("f11", "f12", "f21", "f22"))
table(nfac)
table(fac1, fac2)
res <- residuals(lm(ex12 ~ nfac))
shapiro.test(res)
library(lmtest)
bptest(lm(ex12 ~ nfac), studentize=FALSE)
anova(lm(ex12 ~ nfac))
TukeyHSD(aov(lm(ex12 ~ nfac)))
plot(TukeyHSD(aov(lm(ex12 ~ nfac))))












































































































