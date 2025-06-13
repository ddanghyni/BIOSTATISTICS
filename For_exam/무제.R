library(multtest)
library(ALL)
library(hgu95av2.db)
library(ROCR)
library(rpart)
library(rpart.plot)
library(randomForest)
library(e1071)


# Classification ----------------------------------------------------------

data(golub, package = "multtest")
Labels <- factor(golub.cl, levels=0:1, 
                 labels=c("ALL", "notALL"))
ccnd3 <- grep("CCND3", golub.gnames[,2], ignore.case=TRUE)
sort(golub[ccnd3, ])

decision <- golub[ccnd3, ] >= 1.27
c(sum(decision), ncol(golub))

Pred <- factor(decision, levels=c("TRUE","FALSE"),
               labels=c("ALL", "notALL"))
data.frame(predict=Pred, label=Labels)
table(Pred, Labels)

perf <- function(pred, label) {
  tab <- table(pred, label)
  sensitivity <- tab[1, 1]/sum(tab[ ,1])
  specificity <- tab[2, 2]/sum(tab[ ,2])
  PV.positive <- tab[1, 1]/sum(tab[1, ])
  PV.negative <- tab[2, 2]/sum(tab[2, ])
  c(sensitivity, specificity, PV.positive, PV.negative)
}
perf(Pred, Labels)



# 모든 cut off 이용 -----------------------------------------------------------

cutoff <- sort(golub[ccnd3,])
cutoff <- c(cutoff, Inf)
cutoff

res <- matrix(0, length(cutoff), 4)
rownames(res) <- round(cutoff, 3)
colnames(res) <- c("sensitivity", "specificity", "PV.positive",
                   "PV.negative")
for (i in 1:length(cutoff)) {
  decision <- golub[ccnd3, ] >= cutoff[i]
  Pred <- factor(decision, levels=c("TRUE","FALSE"),
                 labels=c("ALL", "notALL"))
  tab <- table(Labels, Pred)
  res[i, ] <- perf(Pred, Labels)
}
res

data.frame(FPR=1-res[,2], TPR=res[,1])
plot(1-res[,2], res[,1], type="l", xlab="False postive rate",
     ylab="True postive rate", col="blue", lwd=3)




# ROC CURVE ---------------------------------------------------------------

library(ROCR)

true = factor(golub.cl,
              level = c(0, 1), 
              label = c("TRUE", "FALSE"))
# label 무조건 TRUE, FALSE로 작성!!

predccnd3 <- prediction(golub[ccnd3, ], true)
perfccnd3 <- performance(predccnd3, "tpr", "fpr" )
plot(perfccnd3, lwd=3, col="blue")

# auc
performance(predccnd3, "auc")@y.values

gdf5 <- grep("GDF5", golub.gnames[ ,2], ignore.case=TRUE)
predgdf5 <- prediction(golub[gdf5, ], true)
perfgdf5 <- performance(predgdf5, "tpr", "fpr" )
plot(perfgdf5, lwd=4, col="magenta")
performance(predgdf5, "auc")@y.values



# Logistic Regression -----------------------------------------------------

data(golub, package = "multtest")
Factor <- factor(golub.cl, levels=0:1, labels=c("ALL", "AML"))
ccnd3 <- grep("CCND3", golub.gnames[ ,2], ignore.case=TRUE)
data.frame(x=golub[ccnd3,], y=Factor)

# 모델 만들기
g <- glm(Factor ~ golub[ccnd3,], family="binomial")

#glm()에서 이진 종속변수를 factor로 넣으면, 
#predict(..., type = "response") 는 
#“두 번째(level == 1)인 수준” 이 일어날 확률을 돌려줍니다.

# 추정
tmp <- predict(g, type="response") # AML일 확률!! 
# level 1 인 놈 기준으로 반환
pred <- predict(g, type="response") < 0.5
est <- factor(pred, levels=c(TRUE, FALSE),
              labels=c("ALL", "not ALL"))
data.frame(pred=est, Label=Factor)
table(est, Factor)
perf(est, Factor) # 앞에 있는 함수 이용


library(ROCR)
true <- factor(golub.cl, levels=0:1, labels=c("FALSE", "TRUE"))
yprob <- predict(g, type="response")
pred <- prediction(yprob, true)
perf <- performance(pred, "tpr", "fpr" )
plot(perf, lwd=2, col="blue")
performance(pred, "auc")@y.values



# Test Errors -------------------------------------------------------------

set.seed(1234)
train <- sample(1:ncol(golub), 19)
test <- setdiff(1:ncol(golub), train)
Data <- data.frame(x=golub[ccnd3, ], y=Factor)
g <- glm(y ~ x, data=Data, family=binomial, subset=train)
p <- predict(g, Data, type="response")
pred.train <- p[train] > 0.5
pred.test <- p[test] > 0.5
t1 <- table(pred.train, Factor[train])
t1
1-sum(diag(t1))/sum(t1)
t2 <- table(pred.test, Factor[test])
t2
1-sum(diag(t2))/sum(t2)


set.seed(12345)
K <- 100
miss <- matrix(0, K, 2)
colnames(miss) <- c("training", "test")
for (k in 1:K) {
  train <- sample(1:ncol(golub), 19)
  test <- setdiff(1:ncol(golub), train)
  g <- glm(y ~ x, data=Data, family=binomial, subset=train)
  p <- predict(g, Data, type="response")
  pred.train <- p[train] > 0.5
  pred.test <- p[test] > 0.5
  t1 <- table(pred.train, Factor[train])
  t2 <- table(pred.test, Factor[test])
  miss[k,1] <- 1-sum(diag(t1))/sum(t1)
  miss[k,2] <- 1-sum(diag(t2))/sum(t2)
}
apply(miss, 2, summary)



# Classification Trees ----------------------------------------------------

gene <- golub[ccnd3,]
tree <- rpart(golubFactor ~ gene, method="class")
prp(tree, branch.lwd=4, branch.col="blue", extra=101)
tree
predict(tree, type="prob")
predict(tree, type="class")
pred <- predict(tree, type="class")
data.frame(pred=pred, golub=golubFactor)
table(pred, golubFactor)



# anova -------------------------------------------------------------------
library(ALL)
data(ALL)
ALLB123 <- ALL[ ,ALL$BT %in% c("B1","B2","B3")]
ALLB123$BT
table(ALLB123$BT)
table(ALL$BT)
names <- featureNames(ALL)
names
library(hgu95av2.db)
symb <- mget(names, env=hgu95av2SYMBOL)
unlist(symb)
unlist(symb)[1:100]
ALLBTnames <- ALLB123[names, ]
dim(ALLBTnames)
dim(ALL)
probeData <- as.matrix(exprs(ALLBTnames))
row.names(probeData) <- unlist(symb)
probeData[1:20, 1:5]

fun <- function(x) anova(lm(x ~ ALLB123$BT))$Pr[1]
anova.pValue <- apply(exprs(ALLB123), 1, fun)
ww <- anova.pValue < 0.00001
sum(ww)
diagnosed <- factor(ALLBTnames$BT)
diagnosed
Data <- data.frame(t(probeData[ww, ]))
dim(Data)
fit <- rpart(diagnosed ~ ., data=Data)
prp(fit, branch.lwd=4, branch.col="blue", extra=101)
fit

pred <- predict(fit, type="class")
table(pred, diagnosed)
prob <- predict(fit, type="prob")
out <- data.frame(prob, predicted=pred, diagnosis=diagnosed)
out
out[pred!=diagnosed, ]


# Classification Trees train test -----------------------------------------

set.seed(123)
train <- sample(1:78, 39, replace=FALSE)
test <- setdiff(1:78, train)
table(diagnosed[train])
table(diagnosed[test])
fit.tr <- rpart(diagnosed ~ ., data=Data, subset=train)
pred.tr <- predict(fit.tr, Data[train, ], type="class")
table(pred.tr, diagnosed[train])
mean(pred.tr != diagnosed[train])
pred.te <- predict(fit.tr, Data[test, ], type="class")
table(pred.te, diagnosed[test])
mean(pred.te != diagnosed[test])



# Random Forest -----------------------------------------------------------

library(randomForest)
xtr <- Data[train, ]
xte <- Data[test, ]
ytr <- diagnosed[train]
yte <- diagnosed[test]
rf1 <- randomForest(x=xtr, y=ytr, xtest=xte, ytest=yte,
                    ntree=1000, mtry=1)
rf1$test$confusion
rf1.conf <- rf1$test$confusion[1:3, 1:3]
1- sum(diag(rf1.conf))/sum(rf1.conf)


rf2 <- randomForest(x=xtr, y=ytr, xtest=xte, ytest=yte,
                    ntree=1000, mtry=3)
rf2.conf <- rf2$test$confusion[1:3, 1:3]
1- sum(diag(rf2.conf))/sum(rf2.conf)
rf3 <- randomForest(x=xtr, y=ytr, xtest=xte, ytest=yte,
                    ntree=1000, mtry=5)
rf3.conf <- rf3$test$confusion[1:3, 1:3]
1- sum(diag(rf3.conf))/sum(rf3.conf)
rf4 <- randomForest(x=xtr, y=ytr, xtest=xte, ytest=yte,
                    ntree=1000, mtry=10)
rf4.conf <- rf4$test$confusion[1:3, 1:3]
1- sum(diag(rf4.conf))/sum(rf4.conf)



# -------------------------------------------------------------------------
set.seed(111222) 
library(rpart)
data(golub, package = "multtest")
y <- factor(golub.cl, levels=0:1, labels=c("ALL", "notALL"))
x <- t(golub[sample(1:nrow(golub), 3), ])

g <- glm(y ~ x, family="binomial")
pred <- predict(g, type="response") < 0.3
est <- factor(pred, levels=c(TRUE, FALSE), labels=c("ALL", "not ALL"))
T <- table(est, y)

est

perf <- function(pred, label) {
  tab <- table(pred, label)
  sensitivity <- tab[1, 1]/sum(tab[ ,1])
  specificity <- tab[2, 2]/sum(tab[ ,2])
  PV.positive <- tab[1, 1]/sum(tab[1, ])
  PV.negative <- tab[2, 2]/sum(tab[2, ])
  round(c(sensitivity, specificity, PV.positive, PV.negative),4)
}
perf(est, y)

round(c(T[1,1]/sum(T[,1]), T[2,2]/sum(T[,2])), 4)


g2 <- rpart(y ~ x, method="class")
pred2 <- predict(g2, type="prob")[,2] < 0.3 # AML 기준
est2 <- factor(pred2, levels=c(TRUE, FALSE), labels=c("ALL", "not ALL"))
T2 <- table(est2, y)
perf(est2, y)
round(c(T2[1,1]/sum(T2[,1]), T2[2,2]/sum(T2[,2])), 4)



















