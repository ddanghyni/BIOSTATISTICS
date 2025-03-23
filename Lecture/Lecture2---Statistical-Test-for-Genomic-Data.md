Lecture2 : Statistical Test for Genomic Data
================
2025-03-13

- [Required Packages](#required-packages)
- [Statistical Hypothesis Testing](#statistical-hypothesis-testing)
- [The Z-test](#the-z-test)
- [One sample T-test](#one-sample-t-test)
- [Two Sample T-test with Unequal Variances(welch’s
  t-test)](#two-sample-t-test-with-unequal-varianceswelchs-t-test)
- [Two Sample T-test with Equal
  Variances](#two-sample-t-test-with-equal-variances)
- [Test for Equal Variances](#test-for-equal-variances)
- [Histogram and Q-Q plot(정규성
  검정)](#histogram-and-q-q-plot정규성-검정)
- [Normality Tests](#normality-tests)
- [Outliers Test](#outliers-test)
- [Binomial Test](#binomial-test)
- [Chi-squared Test](#chi-squared-test)
- [Fisher’s Exact Test](#fishers-exact-test)

### Required Packages

``` r
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("multtest")
#install.packages("nortest")
#install.packages("car")
#install.packages("outliers")
#install.packages("ape")
#install.packages("SNPassoc")
#install.packages("genetics")
```

### Statistical Hypothesis Testing

- Sample data 만든 후 그 data에서 **sample statistic** 계산 조진다.
- 그 다음 우리는 **Null huypothesis**가 **true**라 가정한다잉!
- 그 다음 $H_0$의 distribution에서 sample statistic의 위치(?)를 본다.
- 우리는 지금 $H_0$가 true라고 가정을 했기 때문에 sample data에서 나온
  ST가 $H_0$의 distribution에 잘 맞아야한다.
- 여기서 잘 맞아야한다라는 뜻은 그 ST가 터무니 없는 확률로 나온다면
  안된다는 뜻이다.
- 이때 그 확률을 우리는 **p-value**라고 한다.
- 그리고 우리는 $H_0$가 true일떄 $H_0$를 reject할 허용 오차를
  significance level $\alpha$라고 한다.

> p-value가 0.05일 때, 만약 $H_0$가 True일 때 관측된 결과(ST)가 우연히
> 나타날 확률이 50%이다! 라는 느낌! 즉 $\alpha = 0.05$는 내가 관측한
> 값이 우연히 나타날 확률이 0.05 보단 작아선 안된다! 라는 threshold를
> setting 한 느낌!

- 핵심은 sample statistic의 distribution을 알아야 한다!!!

### The Z-test

- 현실에선 쓸 일 없음 ㅋㅋㅋㅋㅋ

- 왜? 모표준편차를 어케아노..

- Under the Null

$$
 H_0: \mu = 0 \quad \text{vs} \quad H_1: \mu \neq 0
$$

``` r
library(multtest)
```

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

``` r
library(nortest)
library(car)
```

    ## Loading required package: carData

``` r
library(outliers)
library(ape)
library(SNPassoc)
```

    ## Registered S3 method overwritten by 'SNPassoc':
    ##   method            from       
    ##   summary.haplo.glm haplo.stats

    ## 
    ## Attaching package: 'SNPassoc'

    ## The following object is masked from 'package:ape':
    ## 
    ##     additive

``` r
library(genetics)
```

    ## Loading required package: combinat

    ## 
    ## Attaching package: 'combinat'

    ## The following object is masked from 'package:utils':
    ## 
    ##     combn

    ## Loading required package: gdata

    ## 
    ## Attaching package: 'gdata'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

    ## The following object is masked from 'package:stats':
    ## 
    ##     nobs

    ## The following object is masked from 'package:utils':
    ## 
    ##     object.size

    ## The following object is masked from 'package:base':
    ## 
    ##     startsWith

    ## Loading required package: gtools

    ## 
    ## Attaching package: 'gtools'

    ## The following object is masked from 'package:car':
    ## 
    ##     logit

    ## Loading required package: MASS

    ## Loading required package: mvtnorm

    ## Registered S3 methods overwritten by 'genetics':
    ##   method                      from    
    ##   print.LD                    SNPassoc
    ##   print.LD.data.frame         SNPassoc
    ##   print.summary.LD.data.frame SNPassoc
    ##   summary.LD.data.frame       SNPassoc

    ## 

    ## NOTE: THIS PACKAGE IS NOW OBSOLETE.

    ## 

    ##   The R-Genetics project has developed an set of enhanced genetics

    ##   packages to replace 'genetics'. Please visit the project homepage

    ##   at http://rgenetics.org for informtion.

    ## 

    ## 
    ## Attaching package: 'genetics'

    ## The following objects are masked from 'package:SNPassoc':
    ## 
    ##     LD, LDplot, LDtable

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     %in%, order

    ## The following objects are masked from 'package:base':
    ## 
    ##     %in%, as.factor, order

``` r
library(multtest)

data(golub, package = "multtest")

golubFactor <- factor(golub.cl, levels=0:1, labels=c("ALL","AML"))
  
x <- golub[2058, golubFactor=="ALL"]

n <- length(x)

sigma <- 0.25; mu0 <- 0
z.value <- sqrt(n)*(mean(x) - mu0)/sigma

2*pnorm(-abs(z.value)) # p-value
```

    ## [1] 0.9991094

``` r
f <- function(x) dnorm(x,0,1)
alpha <- 0.05
qz <- qnorm(1-alpha/2)
x1 <- seq(-4, -qz, 0.01); y1 <- dnorm(x1, 0, 1)
x2 <- seq(-qz, qz, 0.01); y2 <- dnorm(x2, 0, 1)
x3 <- seq(qz, 4, 0.01); y3 <- dnorm(x3, 0, 1)
plot(f, -4, 4, cex.lab=1.5, xlab="x", ylab="f(x)",
main="Normal probability density function f(x)")
polygon(c(-4, x1, -qz), c(0, y1, 0), col="red")
polygon(c(-qz, x2, qz), c(0, y2, 0), col="lightblue")
polygon(c(qz, x3, 4), c(0, y3, 0), col="red")
arrows(-3, 0.15, -3, 0.03)
text(-3, 0.23, "Rejection")
text(-3, 0.20, "Region")
text(-3, 0.17, expression(alpha/2))
arrows(3, 0.15, 3, 0.03)
text(3, 0.23, "Rejection")
text(3, 0.20, "Region")
text(3, 0.17, expression(alpha/2))
text(0, 0.23, "Acceptance")
text(0, 0.20, "Region")
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

#### Confidence Interval

- If we were to repeat the procedure thousands of times, we are 95%
  certain that the true mean falls in the confidence interval -\> 이게
  핵심!

``` r
mean(x) + qnorm(0.025) * sigma/sqrt(n)
```

    ## [1] -0.09424511

``` r
mean(x) + qnorm(0.975) * sigma/sqrt(n)
```

    ## [1] 0.09435251

``` r
mean(x) + c(-1, 1) * qnorm(0.975) * sigma / sqrt(n)
```

    ## [1] -0.09424511  0.09435251

### One sample T-test

- 모집단의 분산(표준편차)을 모르고, 표본이 **정규성을 만족한다고
  가정**할 때 사용한다.
- 이때 검정 통계량은 t 분포를 따르고 t분포는 degree of freedom에 따라
  모양이 달라진다.
- 모집단의 분산 대신 표본의 구한 표본분산을 이용하기 때문에 정규분포가
  아닌 t 분포를 따르게 된다.(자세한 증명은 수통 때)

``` r
f <- function(x) dt(x, 5)
alpha <- 0.05
qt5 <- qt(1-alpha/2, 5)
x1 <- seq(-4, -qt5, 0.01)
y1 <- f(x1)
x2 <- seq(-qt5, qt5, 0.01)
y2 <- f(x2)
x3 <- seq(qt5, 4, 0.01)
y3 <- f(x3)
plot(f, -4, 4, xlab="x", ylab="f(x)",
main="T-distribution probability density function f(x)")
polygon(c(-4, x1, -qt5), c(0, y1 , 0), col="red")
polygon(c(-qt5, x2, qt5), c(0, y2, 0), col="lightblue")
polygon(c(qt5, x3, 4), c(0, y3, 0), col = "red")

arrows(-3, 0.15, -3, 0.03)
text(-3, 0.23, "Rejection")
text(-3, 0.20, "Region")
text(-3, 0.17, expression(alpha/2))
arrows(3, 0.15, 3, 0.03)
text(3, 0.23, "Rejection")
text(3, 0.20, "Region")
text(3, 0.17, expression(alpha/2))
text(0, 0.23, "Acceptance")
text(0, 0.20, "Region")
mtext(expression(t[0.025]), side=1, at=-qt5, col="red")
mtext(expression(t[0.975]), side=1, at=qt5, col="red")
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

#### R에서 분포와 관련된 함수 정리

``` r
#' d : 확률밀도함수
#' 특정 분포의 확률밀도값 또는 확률질량값을 계산
dnorm(0, mean = 0, sd = 1)
```

    ## [1] 0.3989423

``` r
#' p : 누적분포함수
#' 특정 분포에서 어떤 값 q 이하가 나올 확률을 계산
pnorm(0, mean = 0, sd = 1)
```

    ## [1] 0.5

``` r
#' q : 분위수함수
#' 누적분포함수의 반대 개념으로, 누적확률 p에 해당하는 값을 계산
qnorm(0.975, mean = 0, sd = 1)
```

    ## [1] 1.959964

``` r
qnorm(0.5, mean = 0, sd = 1)
```

    ## [1] 0

``` r
#' r : 난수생성
#' 특정 분포를 따르는 난수를 n개 생성하여 벡터로 반환
rnorm(10, mean = 0, sd = 1)
```

    ##  [1] -0.1999557 -0.2893906  1.8621327  1.3250889 -0.8567678  0.6671946
    ##  [7] -1.1923083 -0.7168538  0.4241996 -0.6618517

``` r
x <- golub[2058, golubFactor=="ALL"]

mu0 <- 0

n <- 27

t.value <- sqrt(n)*(mean(x) - mu0)/sd(x)

t.value
```

    ## [1] 0.001076867

``` r
qt(0.975, 26)
```

    ## [1] 2.055529

``` r
2 * pt(-abs(t.value), 26) # p-value
```

    ## [1] 0.999149

``` r
alpha <- 0.05

mean(x) + c(qt(alpha/2, n-1), qt(1-alpha/2, n-1)) * sd(x)/sqrt(n)
```

    ## [1] -0.1024562  0.1025636

#### One Sample T-test를 위한 R function : `t.test(data= , mu= , alternative = )`

``` r
t.test(x, mu=0) # default two side test
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  x
    ## t = 0.0010769, df = 26, p-value = 0.9991
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1024562  0.1025636
    ## sample estimates:
    ##   mean of x 
    ## 5.37037e-05

$$
H_1: \mu > 0
$$

``` r
#' mu를 기준으로 p-v는 오른쪽 부분 싸악 
t.test(x, mu=0, alternative = "greater") # one side test
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  x
    ## t = 0.0010769, df = 26, p-value = 0.4996
    ## alternative hypothesis: true mean is greater than 0
    ## 95 percent confidence interval:
    ##  -0.085006       Inf
    ## sample estimates:
    ##   mean of x 
    ## 5.37037e-05

$$
H_1: \mu < 0
$$

``` r
#' mu를 기준으로 p-v는 왼쪽 부분 싸악 
t.test(x, mu=0, alternative = "less") # one side test
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  x
    ## t = 0.0010769, df = 26, p-value = 0.5004
    ## alternative hypothesis: true mean is less than 0
    ## 95 percent confidence interval:
    ##        -Inf 0.08511341
    ## sample estimates:
    ##   mean of x 
    ## 5.37037e-05

#### Example of One sample T-test

- Using a box-and-whiskers plot and a side-by-side plot, we can
  investigate the range of the ALL and AML gene expression values of
  CCND3 (CyclinD3), which are collected in row 1042 of the data matrix
  golub.

> 그니깐 CCND3 gene에서 ALL vs AML을 하자.

``` r
# gred()으로 CCND3 gene의 index 찾기
ccnd3 <- grep("CCND3", golub.gnames[,2], ignore.case=TRUE) 
ccnd3
```

    ## [1] 1042

``` r
par(mfrow=c(1,2)) 
stripchart(golub[ccnd3,] ~ golubFactor, method="jitter", cex.lab=1.5, ylab="", vertical = TRUE, col=c("red", "darkgreen"), xlab="Leukemia subtype")

boxplot(golub[ccnd3,] ~ golubFactor, cex.lab=1.5, main=NULL, xlab="Leukemia subtype", col=c("purple","green"), ylab="CCND3 (Cyclin D3) Expression")
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
ALL <- golubFactor=="ALL" 

#' CCND3 gene에서 ALL로 검정을 하자..
#' H_0 = mu에 대한 가설검정 
t.test(golub[ccnd3, ALL], mu=0) 
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  golub[ccnd3, ALL]
    ## t = 20.06, df = 26, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  1.699817 2.087948
    ## sample estimates:
    ## mean of x 
    ##  1.893883

- t = 20.06 존나 큰 값 -\> 무조건 기각 즉! ALL의 평균 값이 0이
  아니다라고 할 수 있다.

``` r
t.test(golub[ccnd3, !ALL], mu=0)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  golub[ccnd3, !ALL]
    ## t = 3.6249, df = 10, p-value = 0.004651
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  0.2449118 1.0262701
    ## sample estimates:
    ## mean of x 
    ## 0.6355909

- 0.004이라서 기각이긴한데 살짝 부족함.

``` r
t.test(golub[ccnd3, ALL], mu=0, alternative="greater")
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  golub[ccnd3, ALL]
    ## t = 20.06, df = 26, p-value < 2.2e-16
    ## alternative hypothesis: true mean is greater than 0
    ## 95 percent confidence interval:
    ##  1.732853      Inf
    ## sample estimates:
    ## mean of x 
    ##  1.893883

``` r
t.test(golub[ccnd3, !ALL], mu=0, alternative="greater") # two side랑 p-v만 다르다.
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  golub[ccnd3, !ALL]
    ## t = 3.6249, df = 10, p-value = 0.002326
    ## alternative hypothesis: true mean is greater than 0
    ## 95 percent confidence interval:
    ##  0.3177962       Inf
    ## sample estimates:
    ## mean of x 
    ## 0.6355909

``` r
t.test(golub[ccnd3, !ALL], mu=0, alternative="less")
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  golub[ccnd3, !ALL]
    ## t = 3.6249, df = 10, p-value = 0.9977
    ## alternative hypothesis: true mean is less than 0
    ## 95 percent confidence interval:
    ##       -Inf 0.9533856
    ## sample estimates:
    ## mean of x 
    ## 0.6355909

``` r
t.test(golub[ccnd3, ALL], mu=0, alternative="less")
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  golub[ccnd3, ALL]
    ## t = 20.06, df = 26, p-value = 1
    ## alternative hypothesis: true mean is less than 0
    ## 95 percent confidence interval:
    ##      -Inf 2.054912
    ## sample estimates:
    ## mean of x 
    ##  1.893883

### Two Sample T-test with Unequal Variances(welch’s t-test)

- 두 그룹의 평균을 비교할 때 사용한다.
- Two Sample T test를 진행하기 위해선
  1.  각 그룹 샘플에 대한 정규성 검정 실시
  2.  Equal variance에 대한 검정 실시
  3.  분산이 같다면 pooled t-test를 실시
  4.  같지 않다면 welch’s t-test를 실시(general version)
  5.  그 다음 상황에 맞는 t test를 진행한다.

``` r
ccnd3 <- grep("CCND3", golub.gnames[,2], ignore.case=TRUE) 
ccnd3
```

    ## [1] 1042

``` r
par(mfrow=c(1,2)) 
stripchart(golub[ccnd3,] ~ golubFactor, method="jitter", cex.lab=1.5, ylab="", vertical = TRUE, col=c("red", "darkgreen"), xlab="Leukemia subtype")

boxplot(golub[ccnd3,] ~ golubFactor, cex.lab=1.5, main=NULL, xlab="Leukemia subtype", col=c("purple","green"), ylab="CCND3 (Cyclin D3) Expression")
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
data.frame(golub[ccnd3,], y =golubFactor)
```

    ##    golub.ccnd3...   y
    ## 1         2.10892 ALL
    ## 2         1.52405 ALL
    ## 3         1.96403 ALL
    ## 4         2.33597 ALL
    ## 5         1.85111 ALL
    ## 6         1.99391 ALL
    ## 7         2.06597 ALL
    ## 8         1.81649 ALL
    ## 9         2.17622 ALL
    ## 10        1.80861 ALL
    ## 11        2.44562 ALL
    ## 12        1.90496 ALL
    ## 13        2.76610 ALL
    ## 14        1.32551 ALL
    ## 15        2.59385 ALL
    ## 16        1.92776 ALL
    ## 17        1.10546 ALL
    ## 18        1.27645 ALL
    ## 19        1.83051 ALL
    ## 20        1.78352 ALL
    ## 21        0.45827 ALL
    ## 22        2.18119 ALL
    ## 23        2.31428 ALL
    ## 24        1.99927 ALL
    ## 25        1.36844 ALL
    ## 26        2.37351 ALL
    ## 27        1.83485 ALL
    ## 28        0.88941 AML
    ## 29        1.45014 AML
    ## 30        0.42904 AML
    ## 31        0.82667 AML
    ## 32        0.63637 AML
    ## 33        1.02250 AML
    ## 34        0.12758 AML
    ## 35       -0.74333 AML
    ## 36        0.73784 AML
    ## 37        0.49470 AML
    ## 38        1.12058 AML

#### R Function for Unequal Variances in two sample t-test : t.test(data ~ groupfactor, var.equal = FALSE)

``` r
t.test(golub[ccnd3,] ~ golubFactor, var.equal=FALSE) # result = > 두 집단은 다르다~
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  golub[ccnd3, ] by golubFactor
    ## t = 6.3186, df = 16.118, p-value = 9.871e-06
    ## alternative hypothesis: true difference in means between group ALL and group AML is not equal to 0
    ## 95 percent confidence interval:
    ##  0.8363826 1.6802008
    ## sample estimates:
    ## mean in group ALL mean in group AML 
    ##         1.8938826         0.6355909

- From the result, we can say that the data provide strong evidence that
  the population means do differ.

``` r
t.test(golub[ccnd3,] ~ golubFactor, var.equal=FALSE, alternative="greater")
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  golub[ccnd3, ] by golubFactor
    ## t = 6.3186, df = 16.118, p-value = 4.935e-06
    ## alternative hypothesis: true difference in means between group ALL and group AML is greater than 0
    ## 95 percent confidence interval:
    ##  0.9107706       Inf
    ## sample estimates:
    ## mean in group ALL mean in group AML 
    ##         1.8938826         0.6355909

``` r
t.test(golub[ccnd3,] ~ golubFactor, var.equal=FALSE, alternative="less")
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  golub[ccnd3, ] by golubFactor
    ## t = 6.3186, df = 16.118, p-value = 1
    ## alternative hypothesis: true difference in means between group ALL and group AML is less than 0
    ## 95 percent confidence interval:
    ##      -Inf 1.605813
    ## sample estimates:
    ## mean in group ALL mean in group AML 
    ##         1.8938826         0.6355909

### Two Sample T-test with Equal Variances

- 두 집단의 분산에 대한 검정 결과로 $H_0$를 기각하지 못하면 즉, 두
  집단의 분산이 같으면 **pooled sample variance**를 이용한 **pooled
  t-test**를 실시.

#### \`t.test(data ~ groupFactor, var.equal=TRUE)

``` r
t.test(golub[ccnd3,] ~ golubFactor, var.equal=TRUE)
```

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  golub[ccnd3, ] by golubFactor
    ## t = 6.7983, df = 36, p-value = 6.046e-08
    ## alternative hypothesis: true difference in means between group ALL and group AML is not equal to 0
    ## 95 percent confidence interval:
    ##  0.8829143 1.6336690
    ## sample estimates:
    ## mean in group ALL mean in group AML 
    ##         1.8938826         0.6355909

### Test for Equal Variances

- F 분포로 조지는 놈!

``` r
ccnd3 <- grep("CCND3", golub.gnames[,2], ignore.case=TRUE) 
zyxin <- grep("Zyxin", golub.gnames[,2], ignore.case=TRUE)
```

``` r
# 각 그룹별 분산의 차이 
tapply(golub[ccnd3,], golubFactor, var) 
```

    ##       ALL       AML 
    ## 0.2406642 0.3381806

``` r
tapply(golub[zyxin,], golubFactor, var)
```

    ##       ALL       AML 
    ## 0.5224983 0.1351442

#### R function for variance test : `var.test(data~groupFactor)`

``` r
var.test(golub[ccnd3, ] ~ golubFactor) # 두 그룹의 분산이 같다고 할 수 있다->pooled 이용
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  golub[ccnd3, ] by golubFactor
    ## F = 0.71164, num df = 26, denom df = 10, p-value = 0.4652
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.2127735 1.8428387
    ## sample estimates:
    ## ratio of variances 
    ##          0.7116441

``` r
var.test(golub[zyxin, ] ~ golubFactor) # 두 그룹의 부산이 같다고 할 수 없다->pooled 못 쓴다.. welch's 이용
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  golub[zyxin, ] by golubFactor
    ## F = 3.8662, num df = 26, denom df = 10, p-value = 0.02968
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##   1.155958 10.011795
    ## sample estimates:
    ## ratio of variances 
    ##           3.866228

- 분산에 대한 검정은 F검정 말고도 여러 종류가 있다! 참고!!

``` r
bartlett.test(golub[ccnd3,] ~ golubFactor) 
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  golub[ccnd3, ] by golubFactor
    ## Bartlett's K-squared = 0.42236, df = 1, p-value = 0.5158

``` r
bartlett.test(golub[zyxin,] ~ golubFactor)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  golub[zyxin, ] by golubFactor
    ## Bartlett's K-squared = 5.036, df = 1, p-value = 0.02483

``` r
fligner.test(golub[ccnd3,], golubFactor) 
```

    ## 
    ##  Fligner-Killeen test of homogeneity of variances
    ## 
    ## data:  golub[ccnd3, ] and golubFactor
    ## Fligner-Killeen:med chi-squared = 0.14115, df = 1, p-value = 0.7071

``` r
fligner.test(golub[zyxin,], golubFactor)
```

    ## 
    ##  Fligner-Killeen test of homogeneity of variances
    ## 
    ## data:  golub[zyxin, ] and golubFactor
    ## Fligner-Killeen:med chi-squared = 7.2671, df = 1, p-value = 0.007023

``` r
library(car) 
leveneTest(golub[ccnd3,], golubFactor) 
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1336 0.7169
    ##       36

``` r
leveneTest(golub[zyxin,], golubFactor)
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1  9.2119 0.004448 **
    ##       36                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

- 여기까지 복습 완료

### Histogram and Q-Q plot(정규성 검정)

``` r
# histogram
par(mfrow=c(1,2)) 
hist(golub[ccnd3, golubFactor=="ALL"], cex.lab=1.5, col="orange",nclass=20, main=NULL, xlab="CCND3 (Cyclin D3) Expression") 
hist(golub[zyxin, golubFactor=="ALL"], cex.lab=1.5, col="orange",nclass=20, main=NULL, xlab="Zyxin Expression")
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
# q-q plot
par(mfrow=c(1,2)) 
qqnorm(golub[ccnd3, golubFactor=="ALL"], pch=19, cex.lab=1.5, col="red", main=NULL) 
qqline(golub[ccnd3, golubFactor=="ALL"]) 

qqnorm(golub[zyxin, golubFactor=="ALL"], pch=19, cex.lab=1.5,col="red", main=NULL) 
qqline(golub[zyxin, golubFactor=="ALL"])
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

### Normality Tests

``` r
shapiro.test(golub[ccnd3, golubFactor=="ALL"]) 
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  golub[ccnd3, golubFactor == "ALL"]
    ## W = 0.94663, p-value = 0.1774

``` r
shapiro.test(golub[ccnd3, golubFactor=="AML"]) 
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  golub[ccnd3, golubFactor == "AML"]
    ## W = 0.91756, p-value = 0.2989

``` r
shapiro.test(golub[zyxin, golubFactor=="ALL"]) # reject H_0
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  golub[zyxin, golubFactor == "ALL"]
    ## W = 0.89198, p-value = 0.00881

``` r
shapiro.test(golub[zyxin, golubFactor=="AML"]) 
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  golub[zyxin, golubFactor == "AML"]
    ## W = 0.96378, p-value = 0.8178

``` r
library(nortest) 
ad.test(golub[ccnd3, golubFactor=="ALL"]) 
```

    ## 
    ##  Anderson-Darling normality test
    ## 
    ## data:  golub[ccnd3, golubFactor == "ALL"]
    ## A = 0.52154, p-value = 0.1683

``` r
ad.test(golub[zyxin, golubFactor=="ALL"])
```

    ## 
    ##  Anderson-Darling normality test
    ## 
    ## data:  golub[zyxin, golubFactor == "ALL"]
    ## A = 1.013, p-value = 0.009583

### Outliers Test

``` r
library(outliers)
grubbs.test(golub[ccnd3, golubFactor=="ALL"])
```

    ## 
    ##  Grubbs test for one outlier
    ## 
    ## data:  golub[ccnd3, golubFactor == "ALL"]
    ## G = 2.92639, U = 0.65796, p-value = 0.0183
    ## alternative hypothesis: lowest value 0.45827 is an outlier

``` r
grubbs.test(golub[zyxin, golubFactor=="ALL"])
```

    ## 
    ##  Grubbs test for one outlier
    ## 
    ## data:  golub[zyxin, golubFactor == "ALL"]
    ## G = 1.63480, U = 0.89326, p-value = 1
    ## alternative hypothesis: lowest value -1.47649 is an outlier

### Binomial Test

``` r
sum(dbinom(18:22, 22, 0.7)) # d => pdf value
```

    ## [1] 0.1645488

``` r
1 - pbinom(17, 22, 0.7) # p => cdf
```

    ## [1] 0.1645488

``` r
binom.test(18, 22, p=0.7, alternative="greater", conf.level=0.95)
```

    ## 
    ##  Exact binomial test
    ## 
    ## data:  18 and 22
    ## number of successes = 18, number of trials = 22, p-value = 0.1645
    ## alternative hypothesis: true probability of success is greater than 0.7
    ## 95 percent confidence interval:
    ##  0.6309089 1.0000000
    ## sample estimates:
    ## probability of success 
    ##              0.8181818

### Chi-squared Test

``` r
library(ape) 
obs <- table(read.GenBank(c("X94991.1"),as.character=TRUE)) 
obs
```

    ## X94991.1
    ##   a   c   g   t 
    ## 410 789 573 394

``` r
e <- rep(sum(obs)/4, 4) 
e
```

    ## [1] 541.5 541.5 541.5 541.5

``` r
test <- sum((obs-e)^2/e) 
test
```

    ## [1] 187.0674

``` r
1-pchisq(test, 3) 
```

    ## [1] 0

``` r
qchisq(0.95, 3)
```

    ## [1] 7.814728

``` r
chisq.test(obs)
```

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  obs
    ## X-squared = 187.07, df = 3, p-value < 2.2e-16

``` r
pi <- c(0.75,0.25) 
x <- c(5474, 1850) 
chisq.test(x, p=pi)
```

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  x
    ## X-squared = 0.26288, df = 1, p-value = 0.6081

``` r
data(golub, package = "multtest") 
golubFactor <- factor(golub.cl, levels=0:1, labels=c("ALL","AML")) 
gdf5 <- grep("Gdf5", golub.gnames[,2], ignore.case=TRUE)
```

``` r
x <- golub[gdf5, ] 
cutoff <- 0.1 
pred <- ifelse(x < cutoff, "ALL", "AML") 
data.frame(predicted=pred, true=golubFactor) 
```

    ##    predicted true
    ## 1        AML  ALL
    ## 2        ALL  ALL
    ## 3        AML  ALL
    ## 4        ALL  ALL
    ## 5        ALL  ALL
    ## 6        ALL  ALL
    ## 7        AML  ALL
    ## 8        AML  ALL
    ## 9        ALL  ALL
    ## 10       ALL  ALL
    ## 11       ALL  ALL
    ## 12       ALL  ALL
    ## 13       ALL  ALL
    ## 14       ALL  ALL
    ## 15       ALL  ALL
    ## 16       AML  ALL
    ## 17       ALL  ALL
    ## 18       AML  ALL
    ## 19       ALL  ALL
    ## 20       ALL  ALL
    ## 21       ALL  ALL
    ## 22       ALL  ALL
    ## 23       AML  ALL
    ## 24       ALL  ALL
    ## 25       AML  ALL
    ## 26       ALL  ALL
    ## 27       AML  ALL
    ## 28       ALL  AML
    ## 29       ALL  AML
    ## 30       AML  AML
    ## 31       ALL  AML
    ## 32       ALL  AML
    ## 33       AML  AML
    ## 34       ALL  AML
    ## 35       AML  AML
    ## 36       ALL  AML
    ## 37       AML  AML
    ## 38       AML  AML

``` r
table(predicted=pred, true=golubFactor)
```

    ##          true
    ## predicted ALL AML
    ##       ALL  18   6
    ##       AML   9   5

``` r
chisq.test(table(pred, golubFactor))
```

    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect

    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  table(pred, golubFactor)
    ## X-squared = 0.11005, df = 1, p-value = 0.7401

``` r
cutoff <- sort(x) 
pval <- 0
```

``` r
for (i in 1:length(cutoff)) { 
  pred <- ifelse(x < cutoff[i], "ALL", "AML") 
  pval[i] <- chisq.test(table(pred, golubFactor))$p.val 
  }
```

    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred, golubFactor)): Chi-squared approximation may
    ## be incorrect

``` r
plot(cutoff, pval, type="p", pch=20, ylab="pvalues") 
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

``` r
cutoff[pval < 0.05]
```

    ## [1] -0.6839

``` r
pred <- ifelse(x < cutoff[1], "ALL", "AML") 
table(predicted=pred, true=golubFactor)
```

    ##          true
    ## predicted ALL AML
    ##       AML  27  11

``` r
boxplot(x ~ golubFactor, cex.lab=1.5, main=NULL, 
        xlab="Leukemia subtype", col=c("lightblue", "orange"))
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

<hr>

``` r
ccnd3 <- grep("CCND3", golub.gnames[,2], ignore.case=TRUE) 
x2 <- golub[ccnd3, ]
```

``` r
cutoff <- sort(x2) 
pval2 <- 0 
for (i in 1:length(cutoff)) {
  pred2 <- ifelse(x2 < cutoff[i], "ALL", "AML")
  pval2[i] <- chisq.test(table(pred2, golubFactor))$p.val 
}
```

    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect
    ## Warning in chisq.test(table(pred2, golubFactor)): Chi-squared approximation may
    ## be incorrect

``` r
plot(cutoff, pval2, type="p", pch=20, ylab="pvalues") 
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

``` r
cutoff[pval2 < 0.05]
```

    ##  [1] -0.74333  0.45827  0.63637  0.73784  0.82667  0.88941  1.02250  1.10546
    ##  [9]  1.12058  1.27645  1.32551  1.36844  1.45014  1.52405  1.78352  1.80861
    ## [17]  1.81649  1.83051  1.83485  1.85111  1.90496  1.92776  1.96403  1.99391
    ## [25]  1.99927

``` r
boxplot(x2 ~ golubFactor, cex.lab=1.5, main=NULL, xlab="Leukemia subtype", 
        col=c("lightblue", "orange")) 
abline(h=cutoff[pval2 < 0.05], lty=2, col="gray")
```

![](Lecture2---Statistical-Test-for-Genomic-Data_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

### Fisher’s Exact Test

``` r
data <- matrix(c(300,500,3000,7000), 2, byrow=TRUE) 
fisher.test(data)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  data
    ## p-value = 1.336e-05
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  1.201492 1.629240
    ## sample estimates:
    ## odds ratio 
    ##   1.399977
