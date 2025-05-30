
# Computing on a data matrix ----------------------------------------------

gene1 = c(1.00, 1.50, 1.25)
gene1

gene2 = c(1.35, 1.55, 1.00)
gene3 = c(-1.10, -1.50, -1.25)
gene4 = c(-1.20, -1.30, -1.00)

# 한 줄로 다들어 간다..!
genes = c(gene1, gene2, gene3, gene4)
genes

# byrow = TRUE -> 데이터가 row 방향으로 쑥쑥 들어간다.
# gene을 행으로 만들고 sample을 col로 만들고 싶어!
geneData = matrix(genes, nrow = 4, ncol = 3, byrow = TRUE)
geneData


# gene을 열로 sample을 Col로
testData_ = matrix(genes, nrow = 3, ncol = 4)
testData_
t(geneData)

# 이상하쥬..?
testData = matrix(genes, nrow = 4, ncol = 3, byrow = FALSE)
testData

rownames(geneData) = c("gene1", "gene2", "gene3", "gene4")
colnames(geneData) = c("Eric", "Peter", "Anna")
geneData

class(geneData)


# indexing matrix ---------------------------------------------------------

geneData[3,]
geneData[c(1, 2), ]

# 이름으로도 indexing이 가능하다..!
geneData[c('gene1', 'gene2'),]

#' 행 열을 가진 데이터에서 단순 정수로 indexing을하면
#' row 방향 즉 위에서 아래로 내려가면서 순서대로
#' 4행 3열짜리 행렬을 백터처럼 다뤄서, 5번째 원소를 pick!
geneData[5]


# apply -------------------------------------------------------------------

apply(geneData, 1, mean) # 각 행을 가지고 놀자!
apply(geneData, 2, mean) # 각 열을 가지고 놀자!

#' 이게 R에서는 Matrix로 반환할때
#' 열로 쌓이게 되어있음 ㅇㅇ
apply(geneData, 1, function(x) x ** 2)

apply(apply(geneData, 1, function(x) x ** 2), 2, mean)
apply(geneData ** 2, 1, mean)


# lapply ------------------------------------------------------------------

# 리스틑 또는 벡터의 각 요소에 함수를 적용!
# 결과를 리스트로 반환
my_list = list(a = 1:5, b = 6:10)
my_list

#apply(my_list, 2, mean)-> error!

lapply(my_list, mean)

# my_list ** 2 -> error!
lapply(my_list, function(x) x ** 2)
lapply(lapply(my_list, function(x) x ** 2), mean)

# sapply ------------------------------------------------------------------

#lapply와 동일하지만 결과가 리스트가 아니라 벡터/행렬로 나옴

sapply(my_list, mean)
sapply(sapply(my_list, function(x) x ** 2, simplify = FALSE), mean)

# tapply ------------------------------------------------------------------

# 벡터의 값을 그룹 별로 나눠서 함수 적용
scores <- c(80, 90, 70, 95, 88, 92)
group <- factor(c("A", "A", "B", "B", "C", "C"))
tapply(scores, group, mean)



# ++ aggregate ------------------------------------------------------------

geneData
geneType = factor(c("ALL",'ALL', 'AML', 'AML'))
df = data.frame(geneData, geneType)
aggregate(geneData ~ geneType, data = df, mean)


# df indexing -------------------------------------------------------------

patients.df = data.frame(
  patientID = c("101", "102", "103", "104"),
  treatment = c("drug", "placebo", "drug", "placebo"),
  age = c(20, 30, 24, 22) 
)

patients.df

patients.df[patients.df$treatment == "drug", ]

# 행렬에선 rownames지만 df에선 row.names
row.names(patients.df) = patients.df$patientID

# 하나의 정수만 넣으면 해당 열!!!을 뽑는다.
list(patients.df[2],
      patients.df[[2]])



# Generating factor -------------------------------------------------------

# gl( # of group, # of element of each group)
gl(3, 5)
gl(2, 3, labels = c('ALL', 'AML'))

vector = c(1, 1, 0, 1, 1, 0)
factor(vector, levels = c(1, 0), labels = c("ALL", 'AML'))
factor(vector, levels = c(0, 1), labels = c("ALL", 'AML'))


# Functions ---------------------------------------------------------------

set.seed(123)

A = matrix(sample(1:100, 12), nrow = 3, ncol = 4)
A

apply(A, 2, function(x) sum((x - mean(x))^2) / (length(x) - 1))  

# 추가적인 arguments를 넣고 싶을 때
fun1 = function(x, a) sum(x - a)^2/ (length(x) - 1)
apply(A, 2, fun1, a = 40)

data = matrix(1:6, 2, 3)
data

lapply(data, function(x) x^2) # 굳이..?



# ALL data ----------------------------------------------------------------

library(ALL)
data(ALL)
str(ALL)

sampleNames(ALL) # 환자 이름
featureNames(ALL) # 각 환자가 가지고 있는 gene이름

ALL$BT
table(ALL$BT)

BT = ALL$BT %in% c('T', 'T1', 'T2', 'T3', 'T4')
as.numeric(BT)
y = factor(as.numeric(BT), levels = c(0, 1), labels = c('B', 'T'))
#factor(as.numeric(BT), labels = c('B', 'T'))
table(y)


# golub data --------------------------------------------------------------

library(multtest)

data(golub)

#'각 sample의 tumor class -> ALL, AML에 대한 정보
#' ALL -> 0, AML -> 1
#' 나중에 당연히 Factor로 변환해야겠지?
golub.cl

#' gene name 보관
#' gene index, ID, Name
dim(golub.gnames)
golub.gnames[1:10,]

#' gene expression이 담긴 데이터
#' 행이 Gene 열이 Sample
dim(golub)
t(golub)[1:5, 1:10]


# M92287_at ---------------------------------------------------------------

golub.gnames[1042, ]

# 각 sample의 CCND3 gene에 대한 gene expression
golub[1042, ]

# ALL, AML Factor화
golubFactor = factor(golub.cl, levels = c(0, 1), labels = c("ALL", "AML"))
golubFactor == 'ALL'

#'1042번째 gene(CCND3)에서 이제 암 상태가 ALL인 sample만 뽑!
golub[1042, golubFactor == 'ALL'] 

list(
ALL = mean(golub[1042, golubFactor == 'ALL']),
AML = mean(golub[1042, golubFactor == 'AML']))

#' 저 위에 놈을 모든 gene에 대해서
#' 각 gene 별 ALL과 AML의 평균
meanALL = apply(golub[, golubFactor == 'ALL'], 1, mean)
meanAML = apply(golub[, golubFactor == 'AML'], 1, mean)

# 각 sample 별
apply(golub[, golubFactor == 'ALL'], 2, mean)

############################################
# 데이터 구조를 잘 파악하는게 일단 중요!!! #
############################################


# grep --------------------------------------------------------------------

# index 반환...
cd33 = grep("CD33", golub.gnames[,2], ignore.case = TRUE)
grep("CCND3", golub.gnames[,2], ignore.case = TRUE)

golub.gnames[cd33, ]
golub[cd33, ]

