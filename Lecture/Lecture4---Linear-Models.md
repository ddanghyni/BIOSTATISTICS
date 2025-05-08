Lecture4 : Linear Models
================
KIM SANG HYUN(202211545)
2025-05-06

- [01. Introduction to linear Models](#01-introduction-to-linear-models)
  - [Definition of Linear Models](#definition-of-linear-models)
  - [Least Square Regression Fit](#least-square-regression-fit)
  - [Residual Sum of Squares (RSS,
    SSE)](#residual-sum-of-squares-rss-sse)
  - [Least Square Estimation](#least-square-estimation)
  - [Hypothesis Testing](#hypothesis-testing)
  - [Confidence Intervals](#confidence-intervals)
  - [Sum of Squares](#sum-of-squares)
  - [ANOVA Table](#anova-table)
  - [Least Square Estimation](#least-square-estimation-1)
  - [R Squared](#r-squared)
  - [F test Statistic](#f-test-statistic)
  - [Multiple linear Regression](#multiple-linear-regression)
  - [Overall F-test in ANOVA table](#overall-f-test-in-anova-table)
- [2. Linear model with a factor](#2-linear-model-with-a-factor)
- [3. One-way ANOVA](#3-one-way-anova)
  - [All Pairwise Test](#all-pairwise-test)
  - [Example of One-way ANOVA](#example-of-one-way-anova)
- [4. Multiple testing for ANOVA](#4-multiple-testing-for-anova)
- [5. Checking assumptions](#5-checking-assumptions)
- [6. Robust tests](#6-robust-tests)
- [7. Two-way ANOVA](#7-two-way-anova)

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
library(ALL)
library(lmtest)
```

    ## Loading required package: zoo

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

## 01. Introduction to linear Models

- We learned how statistical tests can be used to discover genes with
  different means with respect to **two groups**.

- We’ll learn how to perform similar tests between **three or more
  groups**.

- A technique is called **analysis of variance (ANOVA)**

> 다중집단 비교해보자!! by ANOVA

- ANOVA is based on the assumption that gene expression values are

  - Normally distributed

  - Have Equal variance (homogeneity)

    > 집단 간의 분산이 같아야 한다..

- across the groups of patients, samples, or experiments.

### Definition of Linear Models

- Given a continuous value of $y_i$, a basic form of the linear model is
  :

$$
y_i = \beta_0 + x_i\beta_1 + e_i, \text{for i = 1,...,n}
$$

$$
e \text{~} N(0,\sigma^2)
$$

- It’s commonly assumed that the error variables $e_1,...,e_n$ are
  **independent and normally distributed** with zero mean!!

$$
y_i \text{~}N(\beta_0 + x_i\beta_1,\sigma^2)
$$

### Least Square Regression Fit

- The LS regression line provides the best fit of the data points.

> RSS, 즉 residual의 합을 가\~~장 작게 만들어주는 그 선의 기울기와 절편!

![](images/clipboard-1768613021.png)

### Residual Sum of Squares (RSS, SSE)

- Let $\hat{y_i} = \hat{\beta_0} + \hat{\beta_1x_i}$ be the prediction
  for $Y$ based on the $i$th value of $X$. Then,

$$
e_i = y_i - \hat{y_i}
$$

- represents the $i$th **residual**.

- We define the RSS as

$$
\text{RSS} = \sum_{i=1}^{n} e_i^2 = \sum_{i=1}^{n}(y_i - \hat{y_i})^2
$$

### Least Square Estimation

- 위에서 언급한 RSS를 가장 작게 만들어주는 $\beta_1, \beta_2$를
  구하자..!

$$
(\hat{\beta_0}, \hat{\beta^1}) = \text{argmin}_{\beta_0,\beta_1} RSS
$$

### Hypothesis Testing

- If one wish to ony estimate $\beta_0, \beta_1$ , it is \*\*not
  essential to assume any distirbutional form\*\* for the errors!

> 회귀계수만 추정할때는 에러가 정규성 가성 따를 필요없음!! 걍 미분
> 박는거니깐..

- For confidence intervals or hypothesis tests, one need to assume that
  the errors are \*\*normally distributed\*\*.

> 하지만 회귀계수의 검정이나 신뢰구간을 구할 때는 T 분포 등을 이용하므로
> 정규성 검정을 꼭해야한다..!

- Mathematically, this corresponds to testing

$$
H_0 : \beta_1 = 0  \quad \text{VS}\quad H_1 : \beta_1 \neq 0
$$

- since if $\beta_1 = 0$ (귀무가설 기각 못함.. t-stat이 작음 -\> pv가
  0.05보다 큼) then the model reduces to

$$
Y = \beta_0 + e
$$

- and $X$ is not associated with $Y$.

- To test the null hypo, we compute a t-statistic.

![](images/clipboard-1979367006.png)

> MSE = $\hat{\sigma^2}$ !!

### Confidence Intervals

- 이거 모르면 자퇴해야겠지?

![](images/clipboard-3006532652.png)

### Sum of Squares

- Sum of squares total (SST) : $y_i$가 평균을 중심으로 얼마나 펴져
  있는가? (**Total variation**)

$$
\text{SST} = \sum_{i=1}^n (y_i - \bar{y_i})^2
$$

> 이 Total variation 즉, 전체적인 변동을 두 가지 regression과 error로
> 설명이 가능하다.

- Sum of squares regression (SSR) : fitted된 model이 평균을 중심으로
  얼마나 펴져 있는가? (Regression)

$$
\text{RSS} = \sum_{i = 1}^{n}(\hat{y_i}-\bar{y_i})^2
$$

- SST = SSR + SSE 가 된다.

> 간단하게 생각해보면 SSE가 크면 SSR이 작겠지? -\> fitted model이
> 별로다..!

### ANOVA Table

- $F$ test for $H_0 : \beta_1 \quad\text{VS}\quad H_1:\beta_1 \neq  0$

> t test와 다르게 f test 즉 아노바는 모든 애들이 동시에 0으로 가는지
> 체크
>
> 즉, F-test는 Regression Model이 의미가 있는 지 없는 지 testing!!

$$
F = {\text{MSR}\over \text{MSE}}
$$

> 식을 보면 일단 기각 시킬려면 MSR이 MSE보다 커야함!
>
> MSR이 크다!?! -\> Regression Model is gooood~ -\> Reject H_0

![](images/clipboard-2449914760.png)

### Least Square Estimation

``` r
data(golub, package="multtest") 
zyxin <- grep("Zyxin", golub.gnames[,2], ignore.case=TRUE) 
cmyb <- grep("c-myb", golub.gnames[,2], ignore.case=TRUE) 
x <- golub[zyxin, ] 
y <- golub[cmyb, ] 
```

``` r
g <- lm(y ~ x) 
g$coef
```

    ## (Intercept)           x 
    ##   1.5552103  -0.5881673

``` r
coef(g)
```

    ## (Intercept)           x 
    ##   1.5552103  -0.5881673

``` r
summary(g)
```

    ## 
    ## Call:
    ## lm(formula = y ~ x)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.51086 -0.40011 -0.04658  0.44662  1.03868 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.5552     0.1100  14.132 2.92e-16 ***
    ## x            -0.5882     0.1012  -5.814 1.23e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6603 on 36 degrees of freedom
    ## Multiple R-squared:  0.4843, Adjusted R-squared:   0.47 
    ## F-statistic: 33.81 on 1 and 36 DF,  p-value: 1.23e-06

``` r
plot(x, y, pch=19, xlab="Relative Zyxin gene expression",
ylab="Relative c-MYB gene expression", cex.lab=1.5, col="blue")
abline(g$coef, lwd=3, lty=2, col="red")
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
sum(g$res^2) # SSE (RSS)
```

    ## [1] 15.6939

``` r
sum(residuals(g)^2)
```

    ## [1] 15.6939

``` r
deviance(g)
```

    ## [1] 15.6939

``` r
plot(fitted(g), residuals(g), pch=19, xlab="Fitted values",
ylab="Residuals")
abline(h=0, lty=2, col="grey")
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
summary(g)$coef
```

    ##               Estimate Std. Error   t value     Pr(>|t|)
    ## (Intercept)  1.5552103  0.1100494 14.131930 2.919560e-16
    ## x           -0.5881673  0.1011560 -5.814457 1.230053e-06

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
g$coef
```

    ## (Intercept)           x 
    ##   1.5552103  -0.5881673

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
g$coef - qt(0.975, g$df)* summary(g)$coef[,2]
```

    ## (Intercept)           x 
    ##   1.3320198  -0.7933212

``` r
g$coef + qt(0.975, g$df)* summary(g)$coef[,2]
```

    ## (Intercept)           x 
    ##   1.7784008  -0.3830134

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
confint(g)
```

    ##                  2.5 %     97.5 %
    ## (Intercept)  1.3320198  1.7784008
    ## x           -0.7933212 -0.3830134

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
confint(g, level=0.90) # confidence intervals
```

    ##                    5 %       95 %
    ## (Intercept)  1.3694142  1.7410064
    ## x           -0.7589488 -0.4173859

``` r
n <- length(x)
p <- 1
df1 <- n - p - 1
df2 <- g$df
df3 <- df.residual(g)
c(df1, df2, df3)
```

    ## [1] 36 36 36

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
summary(g)$sigma^2 # MSE, summary(g)$sigma -> sqrt(mse)
```

    ## [1] 0.4359417

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
sum(g$res^2)/g$df
```

    ## [1] 0.4359417

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
deviance(g)/df.residual(g) # deviance -> SSE
```

    ## [1] 0.4359417

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
anova(g)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y
    ##           Df Sum Sq Mean Sq F value   Pr(>F)    
    ## x          1 14.738 14.7383  33.808 1.23e-06 ***
    ## Residuals 36 15.694  0.4359                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
anova(g)$Sum # SSR, SSE 각각 
```

    ## [1] 14.73828 15.69390

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
sum(anova(g)$Sum) # SST
```

    ## [1] 30.43218

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
print("====================================================")
```

    ## [1] "===================================================="

``` r
names(summary(g))
```

    ##  [1] "call"          "terms"         "residuals"     "coefficients" 
    ##  [5] "aliased"       "sigma"         "df"            "r.squared"    
    ##  [9] "adj.r.squared" "fstatistic"    "cov.unscaled"

### R Squared

- The coefficient of determination $R^2$ is number between 0 and 1 that
  indicates

- How well the data points fit the statistical model.

> 그니깐 나의 모델이 주어진 데이터에 얼마나 적합한지…를 나타내는 측도..

$$
R^2 = 1 - {\text{SSE}\over\text{SST}}
$$

> 1에 가까울수록 나의 모델이 데이터에 잘 맞다는 뜻… -\> 즉 뺴지는 텀의
> 분자가 작을 수록
>
> 즉, SSE가 작을수록 -\> 회귀 모델로 설명하지 못하는 부분이 작을수록..

- In gneneral, as we keep adding predictor variables to our model,

- the $R^2$ will impove

> feature가 늘어나면 늘어날 수록 R 값은 무조건 1에 가까워짐!
>
> 이걸 방지하기 위해서 feature의 수가 늘어나면 패널티를 주는 adjusted
> $R^2$을 쓰자..!

``` r
summary(g)$r.squared
```

    ## [1] 0.4842992

``` r
summary(g)$adj.r.squared # adj-r^2
```

    ## [1] 0.4699741

### F test Statistic

``` r
names(summary(g))
```

    ##  [1] "call"          "terms"         "residuals"     "coefficients" 
    ##  [5] "aliased"       "sigma"         "df"            "r.squared"    
    ##  [9] "adj.r.squared" "fstatistic"    "cov.unscaled"

``` r
summary(g)$fstat
```

    ##    value    numdf    dendf 
    ## 33.80792  1.00000 36.00000

### Multiple linear Regression

![](images/clipboard-1567527149.png)

- We interpret $\beta_j$ as the average effect on $y$ of one unit
  increase in $x_j$ , **holding all other predictors fixed**.

- 여기서 feature들 간의 상관관계가 없는게 가장 이상적이다..

- 상관관계가 있으면 회귀계수를 추정할때 상관이 있는 변수의 유무에 따라
  계수의 값이 변할 수 있기 때문이다.

- 변수들간의 상관관계가 있으면..

  - The variance of all coefficients tends to increase..

  - Interpertations become hazardous..

![](images/clipboard-1367441817.png)

- RSS는 simple linear regression에서 단순히 추정할 회귀계수가 늘어날
  뿐이다..

- The unbiased estiamte of $\sigma^2$ is

$$
\hat{\sigma^2} = {\text{RSS}\over{n-p-1}} = \text{MSE}
$$

- 이때 p는 number of parameter이다.

### Overall F-test in ANOVA table

- Hypothesis

  - $H_0 : \beta_1 = ... = \beta_p = 0$

  - $H_1:$ At least one $\beta_j$ is not equal to 0.

> 당연이 귀무가설 기각하지 못하면 $y= \beta_0 + e$

![](images/clipboard-870843951.png)

``` r
x1 <- golub[grep("Gdf5", golub.gnames[,2], ignore.case=TRUE),]
x2 <- golub[grep("CCND3", golub.gnames[,2], ignore.case=TRUE),]
x3 <- golub[grep("Nkr", golub.gnames[, 2], ignore.case=TRUE),]
x4 <- golub[grep("IGFBP5",golub.gnames[,2], ignore.case = TRUE),]
data <- data.frame(Gdf5=x1, CCND3=x2, Nkr=x3, IGFBP5=x4, cmyb=y)
data
```

    ##        Gdf5    CCND3      Nkr   IGFBP5     cmyb
    ## 1   0.24735  2.10892 -1.45769 -1.45769  1.71135
    ## 2  -0.26002  1.52405 -1.39420 -1.39420  0.74335
    ## 3   0.23335  1.96403 -1.46227 -1.46227  1.89287
    ## 4   0.01026  2.33597 -1.40715 -1.40715  2.34702
    ## 5   0.04565  1.85111 -1.42668 -1.42668  1.85889
    ## 6  -0.14913  1.99391 -1.21719 -1.21719  2.03111
    ## 7   0.16400  2.06597 -1.37386 -1.37386  2.05606
    ## 8   0.54271  1.81649 -1.36832 -1.36832  2.14942
    ## 9   0.01893  2.17622 -1.47649 -1.07665  1.67872
    ## 10 -0.10892  1.80861 -1.21583 -1.21583  1.36812
    ## 11 -0.30807  2.44562 -1.28137  0.38217  1.75256
    ## 12  0.06702  1.90496 -1.03209 -1.03209  0.62576
    ## 13 -0.04994  2.76610 -1.36149 -1.36149  2.05013
    ## 14 -0.20145  1.32551  2.07770 -1.39979  1.12207
    ## 15 -0.12313  2.59385 -1.39503 -1.39503  2.10825
    ## 16  0.16253  1.92776 -1.40095 -1.40095  2.17166
    ## 17 -0.04562  1.10546 -1.56783 -1.56783  2.04341
    ## 18  0.23475  1.27645 -1.20466 -1.20466  1.16980
    ## 19 -0.32004  1.83051 -1.24482 -1.24482  2.13521
    ## 20 -0.18194  1.78352 -1.60767 -1.60767  1.83307
    ## 21 -0.68390  0.45827 -1.06221 -1.06221  3.08642
    ## 22  0.08487  2.18119 -1.12665 -1.12665  2.55785
    ## 23  0.12644  2.31428 -1.20963 -1.20963  2.11331
    ## 24 -0.13465  1.99927 -1.48332 -1.48332  2.32874
    ## 25  0.10843  1.36844 -1.25268 -1.25268  1.46383
    ## 26 -0.03124  2.37351 -1.27619 -1.27619  1.94300
    ## 27  0.55321  1.83485 -1.23051 -1.23051  1.93723
    ## 28  0.06908  0.88941 -1.43337 -1.43337  0.45261
    ## 29 -0.44772  1.45014 -1.08902 -1.08902  1.30993
    ## 30  0.29207  0.42904 -1.29865 -0.36744  0.00467
    ## 31 -0.07202  0.82667 -1.26183 -1.26183 -0.00368
    ## 32 -0.03390  0.63637 -1.44434 -1.44434  0.65624
    ## 33  0.19627  1.02250 -1.47218 -1.47218 -0.43408
    ## 34  0.04091  0.12758 -1.34158 -1.34158  0.87417
    ## 35  0.27934 -0.74333 -1.22961 -1.22961  0.88220
    ## 36  0.08552  0.73784 -1.39906 -1.39906 -0.10793
    ## 37  0.26222  0.49470 -1.34579 -1.34579 -0.81267
    ## 38  0.21272  1.12058 -1.32403 -1.32403  0.41329

- `lm(cmyb ~ . , data)`

  - data라는 df 안에서 cmyb는 y(target)으로 ~ . (df 안에 있는 나머지)는
    feature로..

``` r
g2 <- lm(cmyb ~ ., data) 
print("=====================================")
```

    ## [1] "====================================="

``` r
summary(g2)
```

    ## 
    ## Call:
    ## lm(formula = cmyb ~ ., data = data)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.41378 -0.32811 -0.00134  0.35798  1.83861 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.009281   0.583479  -0.016 0.987405    
    ## Gdf5        -0.889039   0.506074  -1.757 0.088239 .  
    ## CCND3        0.681214   0.157580   4.323 0.000133 ***
    ## Nkr         -0.069515   0.216890  -0.321 0.750602    
    ## IGFBP5      -0.247653   0.357863  -0.692 0.493757    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7337 on 33 degrees of freedom
    ## Multiple R-squared:  0.4162, Adjusted R-squared:  0.3455 
    ## F-statistic: 5.883 on 4 and 33 DF,  p-value: 0.001093

``` r
print("=====================================")
```

    ## [1] "====================================="

``` r
summary(g2)$coef
```

    ##                 Estimate Std. Error     t value    Pr(>|t|)
    ## (Intercept) -0.009280741  0.5834787 -0.01590588 0.987405260
    ## Gdf5        -0.889038700  0.5060739 -1.75673695 0.088239293
    ## CCND3        0.681214185  0.1575804  4.32296159 0.000133379
    ## Nkr         -0.069515384  0.2168898 -0.32051017 0.750602433
    ## IGFBP5      -0.247653487  0.3578633 -0.69203382 0.493756938

``` r
print("=====================================")
```

    ## [1] "====================================="

``` r
confint(g2) # confidence interval of coefficients
```

    ##                  2.5 %    97.5 %
    ## (Intercept) -1.1963772 1.1778157
    ## Gdf5        -1.9186538 0.1405764
    ## CCND3        0.3606144 1.0018140
    ## Nkr         -0.5107810 0.3717502
    ## IGFBP5      -0.9757318 0.4804248

``` r
n <- length(x)
p <- 4
df1 <- n - p - 1
df2 <- g2$df
df3 <- df.residual(g2) # degree of freedom of residuals
c(df1, df2, df3)
```

    ## [1] 33 33 33

``` r
summary(g)$sigma
```

    ## [1] 0.6602588

``` r
summary(g)$sigma^2
```

    ## [1] 0.4359417

``` r
g0 <- lm(cmyb ~ 1, data) # intercept model (udner H_0)
summary(g0)
```

    ## 
    ## Call:
    ## lm(formula = cmyb ~ 1, data = data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2209 -0.6322  0.3237  0.6463  1.6782 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.4083     0.1471   9.572 1.49e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9069 on 37 degrees of freedom

``` r
anova(g0, g2) # compare intercept model and full model
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: cmyb ~ 1
    ## Model 2: cmyb ~ Gdf5 + CCND3 + Nkr + IGFBP5
    ##   Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
    ## 1     37 30.432                                
    ## 2     33 17.765  4    12.667 5.8827 0.001093 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# > if H_0 is not rejected => intercept model is goood
```

``` r
summary(g2)$fstat
```

    ##     value     numdf     dendf 
    ##  5.882707  4.000000 33.000000

``` r
anova(g2)
```

    ## Analysis of Variance Table
    ## 
    ## Response: cmyb
    ##           Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## Gdf5       1  2.2097  2.2097  4.1047 0.0509121 .  
    ## CCND3      1 10.1448 10.1448 18.8450 0.0001266 ***
    ## Nkr        1  0.0550  0.0550  0.1022 0.7511839    
    ## IGFBP5     1  0.2578  0.2578  0.4789 0.4937569    
    ## Residuals 33 17.7649  0.5383                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
sum(anova(g2)$F, na.rm = TRUE)/4
```

    ## [1] 5.882707

``` r
cor(data)
```

    ##               Gdf5       CCND3         Nkr      IGFBP5        cmyb
    ## Gdf5    1.00000000 -0.08778002 -0.19143305 -0.16513812 -0.26946109
    ## CCND3  -0.08778002  1.00000000 -0.06126599  0.04850305  0.59879652
    ## Nkr    -0.19143305 -0.06126599  1.00000000  0.02785847 -0.03527205
    ## IGFBP5 -0.16513812  0.04850305  0.02785847  1.00000000 -0.02646746
    ## cmyb   -0.26946109  0.59879652 -0.03527205 -0.02646746  1.00000000

``` r
SLR <- vector("list", length=3)
names(SLR) <- c("bhat", "rsq", "sigma")
for (k in 1:(ncol(data)-1)) {
  g0 <- lm(data$cmyb ~ data[,k])
  SLR$bhat[k] <- coef(g0)[2] # \hat\beta_1
  SLR$rsq[k] <- summary(g0)$r.squared
  SLR$sigma[k] <- summary(g0)$sigma
}
SLR
```

    ## $bhat
    ## [1] -0.98843495  0.70404156 -0.05627507 -0.07019522
    ## 
    ## $rsq
    ## [1] 0.0726092764 0.3585572699 0.0012441173 0.0007005266
    ## 
    ## $sigma
    ## [1] 0.8854145 0.7363668 0.9188507 0.9191007

``` r
SLR
```

    ## $bhat
    ## [1] -0.98843495  0.70404156 -0.05627507 -0.07019522
    ## 
    ## $rsq
    ## [1] 0.0726092764 0.3585572699 0.0012441173 0.0007005266
    ## 
    ## $sigma
    ## [1] 0.8854145 0.7363668 0.9188507 0.9191007

``` r
data.frame(MLR=coef(g2)[-1], SLR=SLR$bhat)
```

    ##                MLR         SLR
    ## Gdf5   -0.88903870 -0.98843495
    ## CCND3   0.68121418  0.70404156
    ## Nkr    -0.06951538 -0.05627507
    ## IGFBP5 -0.24765349 -0.07019522

``` r
data.frame(MLR=summary(g2)$r.squared, SLR=SLR$rsq)
```

    ##         MLR          SLR
    ## 1 0.4162477 0.0726092764
    ## 2 0.4162477 0.3585572699
    ## 3 0.4162477 0.0012441173
    ## 4 0.4162477 0.0007005266

``` r
data.frame(MLR=summary(g2)$sigma, SLR=SLR$sigma)
```

    ##        MLR       SLR
    ## 1 0.733709 0.8854145
    ## 2 0.733709 0.7363668
    ## 3 0.733709 0.9188507
    ## 4 0.733709 0.9191007

``` r
x0 <- data.frame(Gdf5=median(x1), CCND3=median(x2),
Nkr=median(x3), IGFBP5=median(x4))
x0
```

    ##      Gdf5   CCND3       Nkr    IGFBP5
    ## 1 0.04328 1.81255 -1.343685 -1.343685

``` r
predict(g2, x0)
```

    ##        1 
    ## 1.613151

``` r
# 각 simple model에서의 예측값..
pred <- NULL
for (k in 1:length(x0)) {
  g0 <- lm(data$cmyb ~ data[,k])
  pred[k] <- sum(coef(g0)[1] + x0[k]*coef(g0)[2])
}
pred
```

    ## [1] 1.388527 1.607442 1.414129 1.414688

## 2. Linear model with a factor

> $y$는 연속형이지만 , $X$가 범주형일 때.. -\> factor 처리 ㄱㄱㅎ!

- Example : Modeling 3 groups of genes..

  - Suppose we have the following artificial gene expression values : 2,
    3, 1, 2 in group 1

  - 8, 7, 9, 8 in Group 2, and 11, 12, 13, 12 in group 3

``` r
y <- c(2, 3, 1, 2, 8, 7, 9, 8, 11, 12, 13, 12)
factor <- gl(3, 4) # 4개 씩 3 그룹?
factor
```

    ##  [1] 1 1 1 1 2 2 2 2 3 3 3 3
    ## Levels: 1 2 3

``` r
data.frame(x = factor, y = y)
```

    ##    x  y
    ## 1  1  2
    ## 2  1  3
    ## 3  1  1
    ## 4  1  2
    ## 5  2  8
    ## 6  2  7
    ## 7  2  9
    ## 8  2  8
    ## 9  3 11
    ## 10 3 12
    ## 11 3 13
    ## 12 3 12

``` r
model.matrix(y ~ factor - 1) # no intercept term
```

    ##    factor1 factor2 factor3
    ## 1        1       0       0
    ## 2        1       0       0
    ## 3        1       0       0
    ## 4        1       0       0
    ## 5        0       1       0
    ## 6        0       1       0
    ## 7        0       1       0
    ## 8        0       1       0
    ## 9        0       0       1
    ## 10       0       0       1
    ## 11       0       0       1
    ## 12       0       0       1
    ## attr(,"assign")
    ## [1] 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$factor
    ## [1] "contr.treatment"

``` r
lm(y ~ factor - 1)$coef # 평균...
```

    ## factor1 factor2 factor3 
    ##       2       8      12

``` r
summary(lm(y ~ factor - 1))$coef
```

    ##         Estimate Std. Error   t value     Pr(>|t|)
    ## factor1        2  0.4082483  4.898979 8.488414e-04
    ## factor2        8  0.4082483 19.595918 1.086960e-08
    ## factor3       12  0.4082483 29.393877 2.979521e-10

$$
\hat{y} = 2x_1 + 8x_2 + 12x_3
$$

- 여기에 (1, 0, 0) 대입해보자… -\> group의 평균이 나온다.

## 3. One-way ANOVA

- A frequent problem is that of testing the null hypo that three or more
  population means are equal to each other.

- ANOVA can perform this test by comparing **within group variance** to
  **between group variances**.

- Null hypothesis

$$
H_0 : \mu_1 = \mu_2 = ... = \mu_k
$$

> 여기서 주의!! = 0 인 term이 없다!! 즉 F 검정과 다르다..!!
>
> 0인지 아닌지는 신경 안쓰고 단순히 그룹 간의 평균이 같냐 아니냐만
> 판단….

- The three sample means per patient group $\bar{y_j}$ and the overall
  mean $\bar{y}$ can be calculated by :

$$
\bar{y_j} = {1\over n}\sum_{i=1}^{n}y_{ij} \quad \text{for}\quad j = 1, 2, 3 
$$

> 각 집단 별 평균

$$
\bar{y} = {1\over n}(\sum_{j=1}^3\sum_{i=1}^ny_{ij})
$$

> 전체 평균

- Sum of squares within (SSW) : 그룹 안에서의 변동

$$
\text{SSW} = \sum_{j=1}^{g}\sum_{i=1}^{n}(y_{ij}-\bar{y_j})^2
$$

> g : is number of groups

- Sum of squares between (SSB) : 그룹 간의 변동

$$
\text{SSB} = \sum_{j=1}^{g}\sum_{i=1}^n(\bar{y_j}-\bar{y})^2
$$

- F-statististic

$$
f = {\text{SSB}/(g-1)\over\text{SSW} / (N-g)}
$$

> 자자 식을 보면 SSB 즉 그룹 간의 변동이 작으면 F 값이 작아져서
> 귀무가설을 기각 못하고 그럼 그룹들의 평균은 같다라는 결론..

![](images/clipboard-3169972116.png)

<img src="images/clipboard-3610410082.png" style="width:50.0%" />

> 아노바는 반드시 동분산 가정이 성립 되어야지만 사용 할 수 있다…!

![](images/clipboard-2472942022.png)

``` r
set.seed(123)
y <- c(rnorm(10, -1, 1), rnorm(10, 0, 1), rnorm(10, 1, 1))
y
```

    ##  [1] -1.56047565 -1.23017749  0.55870831 -0.92949161 -0.87071226  0.71506499
    ##  [7] -0.53908379 -2.26506123 -1.68685285 -1.44566197  1.22408180  0.35981383
    ## [13]  0.40077145  0.11068272 -0.55584113  1.78691314  0.49785048 -1.96661716
    ## [19]  0.70135590 -0.47279141 -0.06782371  0.78202509 -0.02600445  0.27110877
    ## [25]  0.37496073 -0.68669331  1.83778704  1.15337312 -0.13813694  2.25381492

``` r
factor <- gl(3, 10)
factor
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3
    ## Levels: 1 2 3

``` r
model.matrix(y ~ factor - 1)
```

    ##    factor1 factor2 factor3
    ## 1        1       0       0
    ## 2        1       0       0
    ## 3        1       0       0
    ## 4        1       0       0
    ## 5        1       0       0
    ## 6        1       0       0
    ## 7        1       0       0
    ## 8        1       0       0
    ## 9        1       0       0
    ## 10       1       0       0
    ## 11       0       1       0
    ## 12       0       1       0
    ## 13       0       1       0
    ## 14       0       1       0
    ## 15       0       1       0
    ## 16       0       1       0
    ## 17       0       1       0
    ## 18       0       1       0
    ## 19       0       1       0
    ## 20       0       1       0
    ## 21       0       0       1
    ## 22       0       0       1
    ## 23       0       0       1
    ## 24       0       0       1
    ## 25       0       0       1
    ## 26       0       0       1
    ## 27       0       0       1
    ## 28       0       0       1
    ## 29       0       0       1
    ## 30       0       0       1
    ## attr(,"assign")
    ## [1] 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$factor
    ## [1] "contr.treatment"

``` r
model.matrix(y ~ factor)
```

    ##    (Intercept) factor2 factor3
    ## 1            1       0       0
    ## 2            1       0       0
    ## 3            1       0       0
    ## 4            1       0       0
    ## 5            1       0       0
    ## 6            1       0       0
    ## 7            1       0       0
    ## 8            1       0       0
    ## 9            1       0       0
    ## 10           1       0       0
    ## 11           1       1       0
    ## 12           1       1       0
    ## 13           1       1       0
    ## 14           1       1       0
    ## 15           1       1       0
    ## 16           1       1       0
    ## 17           1       1       0
    ## 18           1       1       0
    ## 19           1       1       0
    ## 20           1       1       0
    ## 21           1       0       1
    ## 22           1       0       1
    ## 23           1       0       1
    ## 24           1       0       1
    ## 25           1       0       1
    ## 26           1       0       1
    ## 27           1       0       1
    ## 28           1       0       1
    ## 29           1       0       1
    ## 30           1       0       1
    ## attr(,"assign")
    ## [1] 0 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$factor
    ## [1] "contr.treatment"

``` r
g1 <- lm(y ~ factor - 1)
g2 <- lm(y ~ factor)
```

``` r
summary(g1)$coef
```

    ##           Estimate Std. Error    t value    Pr(>|t|)
    ## factor1 -0.9253744   0.308421 -3.0003606 0.005740635
    ## factor2  0.2086220   0.308421  0.6764194 0.504528252
    ## factor3  0.5754411   0.308421  1.8657648 0.072974278

``` r
print("======================================================")
```

    ## [1] "======================================================"

``` r
summary(g2)$coef
```

    ##               Estimate Std. Error   t value    Pr(>|t|)
    ## (Intercept) -0.9253744  0.3084210 -3.000361 0.005740635
    ## factor2      1.1339963  0.4361732  2.599876 0.014937889
    ## factor3      1.5008155  0.4361732  3.440870 0.001901063

``` r
g1$coef[2] - g1$coef[1]
```

    ##  factor2 
    ## 1.133996

``` r
g1$coef[3] - g1$coef[1]
```

    ##  factor3 
    ## 1.500815

``` r
anova(g1)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## factor     3 12.310  4.1032  4.3136 0.01311 *
    ## Residuals 27 25.683  0.9512                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
print("======================================================")
```

    ## [1] "======================================================"

``` r
anova(g2)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## factor     2 12.243  6.1216  6.4354 0.005183 **
    ## Residuals 27 25.683  0.9512                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

![](images/clipboard-3558157797.png)

``` r
c(summary(g1)$sigma, summary(g2)$sigma)
```

    ## [1] 0.975313 0.975313

``` r
cbind(residuals(g1), residuals(g2))
```

    ##            [,1]         [,2]
    ## 1  -0.635101291 -0.635101291
    ## 2  -0.304803134 -0.304803134
    ## 3   1.484082670  1.484082670
    ## 4  -0.004117253 -0.004117253
    ## 5   0.054662091  0.054662091
    ## 6   1.640439343  1.640439343
    ## 7   0.386290562  0.386290562
    ## 8  -1.339686879 -1.339686879
    ## 9  -0.761478496 -0.761478496
    ## 10 -0.520287614 -0.520287614
    ## 11  1.015459837  1.015459837
    ## 12  0.151191866  0.151191866
    ## 13  0.192149490  0.192149490
    ## 14 -0.097939245 -0.097939245
    ## 15 -0.764463096 -0.764463096
    ## 16  1.578291176  1.578291176
    ## 17  0.289228517  0.289228517
    ## 18 -2.175239117 -2.175239117
    ## 19  0.492733941  0.492733941
    ## 20 -0.681413369 -0.681413369
    ## 21 -0.643264833 -0.643264833
    ## 22  0.206583958  0.206583958
    ## 23 -0.601445575 -0.601445575
    ## 24 -0.304332356 -0.304332356
    ## 25 -0.200480395 -0.200480395
    ## 26 -1.262134438 -1.262134438
    ## 27  1.262345918  1.262345918
    ## 28  0.577931991  0.577931991
    ## 29 -0.713578064 -0.713578064
    ## 30  1.678373794  1.678373794

> 그래 이제 mu_2 - mu1 = beta_1인건 알겠어.. 그럼 다른건??? 예를 들어
> mu_3 - mu_2 같은 거….

``` r
factor <- relevel(factor, ref=2)
model.matrix(y ~ factor)
```

    ##    (Intercept) factor1 factor3
    ## 1            1       1       0
    ## 2            1       1       0
    ## 3            1       1       0
    ## 4            1       1       0
    ## 5            1       1       0
    ## 6            1       1       0
    ## 7            1       1       0
    ## 8            1       1       0
    ## 9            1       1       0
    ## 10           1       1       0
    ## 11           1       0       0
    ## 12           1       0       0
    ## 13           1       0       0
    ## 14           1       0       0
    ## 15           1       0       0
    ## 16           1       0       0
    ## 17           1       0       0
    ## 18           1       0       0
    ## 19           1       0       0
    ## 20           1       0       0
    ## 21           1       0       1
    ## 22           1       0       1
    ## 23           1       0       1
    ## 24           1       0       1
    ## 25           1       0       1
    ## 26           1       0       1
    ## 27           1       0       1
    ## 28           1       0       1
    ## 29           1       0       1
    ## 30           1       0       1
    ## attr(,"assign")
    ## [1] 0 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$factor
    ## [1] "contr.treatment"

``` r
g3 <- lm(y ~ factor)
```

``` r
summary(g3)$coef
```

    ##               Estimate Std. Error    t value   Pr(>|t|)
    ## (Intercept)  0.2086220  0.3084210  0.6764194 0.50452825
    ## factor1     -1.1339963  0.4361732 -2.5998761 0.01493789
    ## factor3      0.3668192  0.4361732  0.8409942 0.40773784

``` r
factor <- relevel(factor, ref=3)
g_tmp <- lm(y ~ factor)
summary(g_tmp)$coef
```

    ##               Estimate Std. Error    t value    Pr(>|t|)
    ## (Intercept)  0.5754411  0.3084210  1.8657648 0.072974278
    ## factor2     -0.3668192  0.4361732 -0.8409942 0.407737841
    ## factor1     -1.5008155  0.4361732 -3.4408703 0.001901063

``` r
summary(g2)$coef
```

    ##               Estimate Std. Error   t value    Pr(>|t|)
    ## (Intercept) -0.9253744  0.3084210 -3.000361 0.005740635
    ## factor2      1.1339963  0.4361732  2.599876 0.014937889
    ## factor3      1.5008155  0.4361732  3.440870 0.001901063

``` r
summary(g3)$coef
```

    ##               Estimate Std. Error    t value   Pr(>|t|)
    ## (Intercept)  0.2086220  0.3084210  0.6764194 0.50452825
    ## factor1     -1.1339963  0.4361732 -2.5998761 0.01493789
    ## factor3      0.3668192  0.4361732  0.8409942 0.40773784

``` r
anova(g2)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## factor     2 12.243  6.1216  6.4354 0.005183 **
    ## Residuals 27 25.683  0.9512                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(g3)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## factor     2 12.243  6.1216  6.4354 0.005183 **
    ## Residuals 27 25.683  0.9512                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

> F-value는 똑같다,,

``` r
factor = gl(3, 10)
summary(lm(y ~ factor))
```

    ## 
    ## Call:
    ## lm(formula = y ~ factor)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.17524 -0.64122 -0.05103  0.46612  1.67837 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  -0.9254     0.3084  -3.000  0.00574 **
    ## factor2       1.1340     0.4362   2.600  0.01494 * 
    ## factor3       1.5008     0.4362   3.441  0.00190 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9753 on 27 degrees of freedom
    ## Multiple R-squared:  0.3228, Adjusted R-squared:  0.2727 
    ## F-statistic: 6.435 on 2 and 27 DF,  p-value: 0.005183

``` r
anova(lm(y ~ factor))
```

    ## Analysis of Variance Table
    ## 
    ## Response: y
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## factor     2 12.243  6.1216  6.4354 0.005183 **
    ## Residuals 27 25.683  0.9512                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### All Pairwise Test

- 모든 가능한 그룹 쌍에 대해 가설을 각각 검정

$$
H_0 : \mu_i = \mu_j
$$

``` r
pairwise.t.test(y, factor, p.adjust.method="bonferroni")
```

    ## 
    ##  Pairwise comparisons using t tests with pooled SD 
    ## 
    ## data:  y and factor 
    ## 
    ##   1      2     
    ## 2 0.0448 -     
    ## 3 0.0057 1.0000
    ## 
    ## P value adjustment method: bonferroni

``` r
pairwise.t.test(y, factor, p.adjust.method="holm")
```

    ## 
    ##  Pairwise comparisons using t tests with pooled SD 
    ## 
    ## data:  y and factor 
    ## 
    ##   1      2     
    ## 2 0.0299 -     
    ## 3 0.0057 0.4077
    ## 
    ## P value adjustment method: holm

``` r
TukeyHSD(aov(y ~ factor))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = y ~ factor)
    ## 
    ## $factor
    ##          diff         lwr      upr     p adj
    ## 2-1 1.1339963  0.05254072 2.215452 0.0384525
    ## 3-1 1.5008155  0.41935988 2.582271 0.0052283
    ## 3-2 0.3668192 -0.71463643 1.448275 0.6812485

``` r
plot(TukeyHSD(aov(y ~ factor)))
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

- Example: No difference between the means

``` r
set.seed(1234)
y <- rnorm(30, 1.9, 1)
y
```

    ##  [1]  0.6929343  2.1774292  2.9844412 -0.4456977  2.3291247  2.4060559
    ##  [7]  1.3252600  1.3533681  1.3355480  1.0099622  1.4228073  0.9016136
    ## [13]  1.1237461  1.9644588  2.8594941  1.7897145  1.3889905  0.9888046
    ## [19]  1.0628283  4.3158352  2.0340882  1.4093141  1.4594521  2.3595894
    ## [25]  1.2062798  0.4517951  2.4747557  0.8763443  1.8848617  0.9640514

``` r
factor <- gl(3, 10)
g <- lm(y ~ factor)
anova(g)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y
    ##           Df  Sum Sq Mean Sq F value Pr(>F)
    ## factor     2  0.4767 0.23837  0.2778 0.7596
    ## Residuals 27 23.1692 0.85812

``` r
pairwise.t.test(y, factor, p.adjust.method="bonferroni")
```

    ## 
    ##  Pairwise comparisons using t tests with pooled SD 
    ## 
    ## data:  y and factor 
    ## 
    ##   1 2
    ## 2 1 -
    ## 3 1 1
    ## 
    ## P value adjustment method: bonferroni

``` r
TukeyHSD(aov(g))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = g)
    ## 
    ## $factor
    ##             diff        lwr       upr     p adj
    ## 2-1  0.264986701 -0.7621747 1.2921481 0.7997113
    ## 3-1 -0.004789407 -1.0319509 1.0223720 0.9999263
    ## 3-2 -0.269776109 -1.2969376 0.7573853 0.7932794

``` r
plot(TukeyHSD(aov(g)))
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

### Example of One-way ANOVA

``` r
library(ALL)
data(ALL)
dim(ALL)
```

    ## Features  Samples 
    ##    12625      128

``` r
ALL[1:10, 1:5]
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 10 features, 5 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: 01005 01010 ... 04007 (5 total)
    ##   varLabels: cod diagnosis ... date last seen (21 total)
    ##   varMetadata: labelDescription
    ## featureData: none
    ## experimentData: use 'experimentData(object)'
    ##   pubMedIds: 14684422 16243790 
    ## Annotation: hgu95av2

``` r
str(ALL)
```

    ## Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
    ##   ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
    ##   .. .. ..@ name             : chr "Chiaretti et al."
    ##   .. .. ..@ lab              : chr "Department of Medical Oncology, Dana-Farber Cancer Institute, Department of Medicine, Brigham and Women's Hospi"| __truncated__
    ##   .. .. ..@ contact          : chr ""
    ##   .. .. ..@ title            : chr "Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with"| __truncated__
    ##   .. .. ..@ abstract         : chr "Gene expression profiles were examined in 33 adult patients with T-cell acute lymphocytic leukemia (T-ALL). Non"| __truncated__
    ##   .. .. ..@ url              : chr ""
    ##   .. .. ..@ pubMedIds        : chr [1:2] "14684422" "16243790"
    ##   .. .. ..@ samples          : list()
    ##   .. .. ..@ hybridizations   : list()
    ##   .. .. ..@ normControls     : list()
    ##   .. .. ..@ preprocessing    : list()
    ##   .. .. ..@ other            : list()
    ##   .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
    ##   .. .. .. .. ..@ .Data:List of 1
    ##   .. .. .. .. .. ..$ : int [1:3] 1 0 0
    ##   .. .. .. .. ..$ names: chr "MIAME"
    ##   ..@ assayData        :<environment: 0x112cffc88> 
    ##   ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
    ##   .. .. ..@ varMetadata      :'data.frame':  21 obs. of  1 variable:
    ##   .. .. .. ..$ labelDescription: chr [1:21] " Patient ID" " Date of diagnosis" " Gender of the patient" " Age of the patient at entry" ...
    ##   .. .. ..@ data             :'data.frame':  128 obs. of  21 variables:
    ##   .. .. .. ..$ cod           : chr [1:128] "1005" "1010" "3002" "4006" ...
    ##   .. .. .. ..$ diagnosis     : chr [1:128] "5/21/1997" "3/29/2000" "6/24/1998" "7/17/1997" ...
    ##   .. .. .. ..$ sex           : Factor w/ 2 levels "F","M": 2 2 1 2 2 2 1 2 2 2 ...
    ##   .. .. .. ..$ age           : int [1:128] 53 19 52 38 57 17 18 16 15 40 ...
    ##   .. .. .. ..$ BT            : Factor w/ 10 levels "B","B1","B2",..: 3 3 5 2 3 2 2 2 3 3 ...
    ##   .. .. .. ..$ remission     : Factor w/ 2 levels "CR","REF": 1 1 1 1 1 1 1 1 1 1 ...
    ##   .. .. .. ..$ CR            : chr [1:128] "CR" "CR" "CR" "CR" ...
    ##   .. .. .. ..$ date.cr       : chr [1:128] "8/6/1997" "6/27/2000" "8/17/1998" "9/8/1997" ...
    ##   .. .. .. ..$ t(4;11)       : logi [1:128] FALSE FALSE NA TRUE FALSE FALSE ...
    ##   .. .. .. ..$ t(9;22)       : logi [1:128] TRUE FALSE NA FALSE FALSE FALSE ...
    ##   .. .. .. ..$ cyto.normal   : logi [1:128] FALSE FALSE NA FALSE FALSE FALSE ...
    ##   .. .. .. ..$ citog         : chr [1:128] "t(9;22)" "simple alt." NA "t(4;11)" ...
    ##   .. .. .. ..$ mol.biol      : Factor w/ 6 levels "ALL1/AF4","BCR/ABL",..: 2 4 2 1 4 4 4 4 4 2 ...
    ##   .. .. .. ..$ fusion protein: Factor w/ 3 levels "p190","p190/p210",..: 3 NA 1 NA NA NA NA NA NA 1 ...
    ##   .. .. .. ..$ mdr           : Factor w/ 2 levels "NEG","POS": 1 2 1 1 1 1 2 1 1 1 ...
    ##   .. .. .. ..$ kinet         : Factor w/ 2 levels "dyploid","hyperd.": 1 1 1 1 1 2 2 1 1 NA ...
    ##   .. .. .. ..$ ccr           : logi [1:128] FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##   .. .. .. ..$ relapse       : logi [1:128] FALSE TRUE TRUE TRUE TRUE TRUE ...
    ##   .. .. .. ..$ transplant    : logi [1:128] TRUE FALSE FALSE FALSE FALSE FALSE ...
    ##   .. .. .. ..$ f.u           : chr [1:128] "BMT / DEATH IN CR" "REL" "REL" "REL" ...
    ##   .. .. .. ..$ date last seen: chr [1:128] NA "8/28/2000" "10/15/1999" "1/23/1998" ...
    ##   .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
    ##   .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
    ##   .. .. .. .. ..@ .Data:List of 1
    ##   .. .. .. .. .. ..$ : int [1:3] 1 1 0
    ##   .. .. .. .. ..$ names: chr "AnnotatedDataFrame"
    ##   ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
    ##   .. .. ..@ varMetadata      :'data.frame':  0 obs. of  1 variable:
    ##   .. .. .. ..$ labelDescription: logi(0) 
    ##   .. .. ..@ data             :'data.frame':  12625 obs. of  0 variables
    ##   .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
    ##   .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
    ##   .. .. .. .. ..@ .Data:List of 1
    ##   .. .. .. .. .. ..$ : int [1:3] 1 1 0
    ##   .. .. .. .. ..$ names: chr "AnnotatedDataFrame"
    ##   ..@ annotation       : chr "hgu95av2"
    ##   ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
    ##   .. .. ..@ varMetadata      :'data.frame':  0 obs. of  1 variable:
    ##   .. .. .. ..$ labelDescription: chr(0) 
    ##   .. .. ..@ data             :'data.frame':  128 obs. of  0 variables
    ##   .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
    ##   .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
    ##   .. .. .. .. ..@ .Data:List of 1
    ##   .. .. .. .. .. ..$ : int [1:3] 1 1 0
    ##   .. .. .. .. ..$ names: chr "AnnotatedDataFrame"
    ##   ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
    ##   .. .. ..@ .Data:List of 4
    ##   .. .. .. ..$ : int [1:3] 2 10 0
    ##   .. .. .. ..$ : int [1:3] 2 5 5
    ##   .. .. .. ..$ : int [1:3] 1 3 0
    ##   .. .. .. ..$ : int [1:3] 1 0 0
    ##   .. .. ..$ names: chr [1:4] "R" "Biobase" "eSet" "ExpressionSet"

``` r
exprs(ALL)[1:10, 1:5]
```

    ##              01005     01010    03002    04006    04007
    ## 1000_at   7.597323  7.479445 7.567593 7.384684 7.905312
    ## 1001_at   5.046194  4.932537 4.799294 4.922627 4.844565
    ## 1002_f_at 3.900466  4.208155 3.886169 4.206798 3.416923
    ## 1003_s_at 5.903856  6.169024 5.860459 6.116890 5.687997
    ## 1004_at   5.925260  5.912780 5.893209 6.170245 5.615210
    ## 1005_at   8.570990 10.428299 9.616713 9.937155 9.983809
    ## 1006_at   3.656143  3.853979 3.646808 3.874289 3.547361
    ## 1007_s_at 7.623562  7.543604 7.916954 6.816397 7.516981
    ## 1008_f_at 8.903547  9.903953 8.494499 9.533983 8.871669
    ## 1009_at   9.371888  9.322177 9.304982 9.135370 9.627175

``` r
dim(exprs(ALL))
```

    ## [1] 12625   128

``` r
table(ALL$BT)
```

    ## 
    ##  B B1 B2 B3 B4  T T1 T2 T3 T4 
    ##  5 19 36 23 12  5  1 15 10  2

``` r
B1B2B3 <- ALL$BT %in% c("B1","B2","B3")
ex <- exprs(ALL)
y <- as.numeric(ex[row.names(ex)=="1866_g_at", B1B2B3])
factor <- factor(ALL$BT[B1B2B3], labels=c("B1","B2","B3"))
```

``` r
col <- c("orange","darkgreen","blue")
xlab <- "B-cell ALL stage"
ylab <- "SKI-like oncogene expression"
par(mfrow=c(1, 2))
stripchart(y ~ factor, method="jitter", cex.lab=1.5, xlab=xlab,
pch=19, col=col, vertical=TRUE, ylab=ylab)
boxplot(y ~ factor, cex.lab=1.5, main=NULL, boxwex=0.3, col=col,
xlab=xlab, ylab=ylab)
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

``` r
g <- lm(y ~ factor)
summary(g)
```

    ## 
    ## Call:
    ## lm(formula = y ~ factor)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.59070 -0.25251 -0.07008  0.16395  1.29635 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  4.58222    0.08506  53.873  < 2e-16 ***
    ## factorB2    -0.43689    0.10513  -4.156 8.52e-05 ***
    ## factorB3    -0.72193    0.11494  -6.281 2.00e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3707 on 75 degrees of freedom
    ## Multiple R-squared:  0.3461, Adjusted R-squared:  0.3287 
    ## F-statistic: 19.85 on 2 and 75 DF,  p-value: 1.207e-07

``` r
anova(g)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y
    ##           Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## factor     2  5.4563 2.72813  19.848 1.207e-07 ***
    ## Residuals 75 10.3091 0.13745                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(g))
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## factor       2  5.456  2.7281   19.85 1.21e-07 ***
    ## Residuals   75 10.309  0.1375                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
pairwise.t.test(y, factor, p.adjust.method="bonferroni")
```

    ## 
    ##  Pairwise comparisons using t tests with pooled SD 
    ## 
    ## data:  y and factor 
    ## 
    ##    B1      B2     
    ## B2 0.00026 -      
    ## B3 6e-08   0.01553
    ## 
    ## P value adjustment method: bonferroni

``` r
TukeyHSD(aov(g))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = g)
    ## 
    ## $factor
    ##             diff        lwr         upr     p adj
    ## B2-B1 -0.4368885 -0.6882697 -0.18550726 0.0002489
    ## B3-B1 -0.7219306 -0.9967602 -0.44710107 0.0000001
    ## B3-B2 -0.2850422 -0.5216834 -0.04840097 0.0141733

``` r
plot(TukeyHSD(aov(g)))
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

``` r
y2 <- as.numeric(ex[row.names(ex)=="1242_at", B1B2B3])
ylab <- "Ets2 expression"

par(mfrow=c(1, 2))
stripchart(y2 ~ factor, method="jitter", cex.lab=1.5, pch=19,
vertical=TRUE, xlab=xlab, ylab=ylab, col=col)
boxplot(y2 ~ factor, cex.lab=1.5, main=NULL, xlab=xlab,
boxwex=0.3, ylab=ylab, col=col)
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

``` r
g <- lm(y2 ~ factor)
summary(g)
```

    ## 
    ## Call:
    ## lm(formula = y2 ~ factor)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.48408 -0.13461 -0.01297  0.12232  0.67138 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  6.55083    0.05673 115.483   <2e-16 ***
    ## factorB2     0.03331    0.07011   0.475    0.636    
    ## factorB3    -0.04675    0.07665  -0.610    0.544    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2473 on 75 degrees of freedom
    ## Multiple R-squared:  0.01925,    Adjusted R-squared:  -0.006898 
    ## F-statistic: 0.7362 on 2 and 75 DF,  p-value: 0.4823

``` r
anova(g)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y2
    ##           Df Sum Sq  Mean Sq F value Pr(>F)
    ## factor     2 0.0900 0.045012  0.7362 0.4823
    ## Residuals 75 4.5854 0.061138

``` r
pairwise.t.test(y2, factor, p.adjust.method="bonferroni")
```

    ## 
    ##  Pairwise comparisons using t tests with pooled SD 
    ## 
    ## data:  y2 and factor 
    ## 
    ##    B1   B2  
    ## B2 1.00 -   
    ## B3 1.00 0.69
    ## 
    ## P value adjustment method: bonferroni

``` r
TukeyHSD(aov(g))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = g)
    ## 
    ## $factor
    ##              diff        lwr        upr     p adj
    ## B2-B1  0.03330952 -0.1343427 0.20096171 0.8832341
    ## B3-B1 -0.04675441 -0.2300449 0.13653604 0.8151498
    ## B3-B2 -0.08006393 -0.2378856 0.07775778 0.4492574

``` r
plot(TukeyHSD(aov(g)))
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

``` r
y3 <- as.numeric(ex[row.names(ex)=="1242_at", ])
g <- lm(y3 ~ ALL$BT)
anova(g)
```

    ## Analysis of Variance Table
    ## 
    ## Response: y3
    ##            Df Sum Sq  Mean Sq F value Pr(>F)
    ## ALL$BT      9 0.2929 0.032543  0.5573 0.8294
    ## Residuals 118 6.8904 0.058393

``` r
xlab <- "B-cell and T-cell stage"
par(mfrow=c(2, 1))
stripchart(y3 ~ ALL$BT, method="jitter", cex.lab=1.5, pch=19,
vertical=TRUE, xlab=xlab, ylab=ylab, col=1:10)
boxplot(y3 ~ ALL$BT, cex.lab=1.5, main=NULL, xlab=xlab,
boxwex=0.3, ylab=ylab, col=1:10)
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

## 4. Multiple testing for ANOVA

``` r
dim(ex[, B1B2B3])
```

    ## [1] 12625    78

``` r
fun <- function(t) anova(lm(t ~ factor))$Pr[1]
anova.pValues <- apply(ex[, B1B2B3], 1, fun)
```

``` r
pBF <- p.adjust(anova.pValues, method="bonferroni")
pHO <- p.adjust(anova.pValues, method="holm")
pBH <- p.adjust(anova.pValues, method="BH")
```

``` r
alpha <- 0.05
c(sum(pBF < alpha), sum(pHO < alpha), sum(pBH < alpha))
```

    ## [1]  53  54 815

## 5. Checking assumptions

- 정규성 검정

``` r
library(ALL)
data(ALL)
```

``` r
B1B2B3 <- ALL$BT %in% c("B1","B2","B3")
ex <- exprs(ALL)
y <- as.numeric(ex[row.names(ex)=="1866_g_at", B1B2B3])
factor <- factor(ALL$BT[B1B2B3], labels=c("B1","B2","B3"))
```

``` r
# 정규성 검정
# H_0 is normality okay
res <- residuals(lm(y ~ factor))
shapiro.test(res) # pv < 0.05 이므로 아노바 ㄴㄴ
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  res
    ## W = 0.93462, p-value = 0.0005989

``` r
par(mfrow=c(1, 2))
qqnorm(res, pch=19, cex.lab=1.5, col="red", main=NULL)
qqline(res)
hist(res, nclass=20, col="orange", xlab="residuals", main="")
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

- 등분산 검정

- H_0 : 등분산이노

``` r
library(lmtest)
bptest(lm(y ~ factor), studentize=FALSE)
```

    ## 
    ##  Breusch-Pagan test
    ## 
    ## data:  lm(y ~ factor)
    ## BP = 8.7311, df = 2, p-value = 0.01271

``` r
bartlett.test(y ~ factor)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  y by factor
    ## Bartlett's K-squared = 7.3117, df = 2, p-value = 0.02584

``` r
fligner.test(y, factor)
```

    ## 
    ##  Fligner-Killeen test of homogeneity of variances
    ## 
    ## data:  y and factor
    ## Fligner-Killeen:med chi-squared = 5.7718, df = 2, p-value = 0.05581

``` r
library(car)
```

    ## Loading required package: carData

``` r
leveneTest(y, factor)
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.0666 0.05246 .
    ##       75                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
g <- lm(y ~ factor)
plot(fitted(g), res, pch=19, col="red",
xlab="Fitted values", ylab = "Residuals")
abline(h=c(max(res), 0, min(res)), lty=2, col="grey")
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

## 6. Robust tests

- 잔차가 가정을 따르지 않을 때…

- When only homoscedasticity is violated, **Welch’s ANOVA** -\> 정규성은
  만족해야함!!!

- When normality is violated, the **Kruskal-Wallis rank sum test** -\>
  정규성 x 일 때!!!

#### Example : Unequal variances

``` r
# Welch's ANOVA -> Games Howell test
oneway.test(y ~ factor, var.equal = FALSE)
```

    ## 
    ##  One-way analysis of means (not assuming equal variances)
    ## 
    ## data:  y and factor
    ## F = 14.157, num df = 2.000, denom df = 36.998, p-value = 2.717e-05

#### Example : Non-normallay distributed residuals

``` r
kruskal.test(y~factor)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  y by factor
    ## Kruskal-Wallis chi-squared = 30.667, df = 2, p-value = 2.192e-07

``` r
f2 <- function(t) kruskal.test(t ~ factor)$p.value
kruskal.pValues <- apply(ex[, B1B2B3], 1, f2)
```

``` r
pBF <- p.adjust(kruskal.pValues, method="bonferroni")
pHO <- p.adjust(kruskal.pValues, method="holm")
pBH <- p.adjust(kruskal.pValues, method="BH")
alpha <- 0.05
c(sum(pBF < alpha), sum(pHO < alpha), sum(pBH < alpha))
```

    ## [1]  15  15 696

## 7. Two-way ANOVA

``` r
library(ALL)
data(ALL)
```

``` r
ex <- exprs(ALL)
dim(ex)
```

    ## [1] 12625   128

``` r
table(ALL$BT)
```

    ## 
    ##  B B1 B2 B3 B4  T T1 T2 T3 T4 
    ##  5 19 36 23 12  5  1 15 10  2

``` r
table(ALL$mol.biol)
```

    ## 
    ## ALL1/AF4  BCR/ABL E2A/PBX1      NEG   NUP-98  p15/p16 
    ##       10       37        5       74        1        1

``` r
w1 <- ALL$BT %in% c("B", "B1", "B2", "B3", "B4")
w2 <- ALL$mol.biol %in% c("BCR/ABL", "NEG")
ex12 <- ex["32069_at", w1 & w2]
length(ex12)
```

    ## [1] 79

``` r
facB <- ceiling(as.integer(ALL$BT[w1 & w2])/3)
facB
```

    ##  [1] 1 1 2 1 1 1 1 1 1 2 2 1 2 1 1 2 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1
    ## [39] 1 2 1 1 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 1 1 2 2 2 2 1 1 1
    ## [77] 1 1 1

``` r
data.frame(B.type=ALL$BT[w1 & w2], facB=facB)
```

    ##    B.type facB
    ## 1      B2    1
    ## 2      B2    1
    ## 3      B4    2
    ## 4      B2    1
    ## 5      B1    1
    ## 6      B1    1
    ## 7      B1    1
    ## 8      B2    1
    ## 9      B2    1
    ## 10     B3    2
    ## 11     B3    2
    ## 12     B2    1
    ## 13     B3    2
    ## 14      B    1
    ## 15     B2    1
    ## 16     B3    2
    ## 17     B2    1
    ## 18     B3    2
    ## 19     B2    1
    ## 20     B2    1
    ## 21     B2    1
    ## 22     B1    1
    ## 23     B2    1
    ## 24     B2    1
    ## 25     B2    1
    ## 26      B    1
    ## 27      B    1
    ## 28     B2    1
    ## 29     B2    1
    ## 30     B2    1
    ## 31     B2    1
    ## 32     B2    1
    ## 33     B2    1
    ## 34     B2    1
    ## 35     B2    1
    ## 36     B4    2
    ## 37     B2    1
    ## 38     B2    1
    ## 39     B2    1
    ## 40     B4    2
    ## 41     B2    1
    ## 42     B2    1
    ## 43     B3    2
    ## 44     B3    2
    ## 45     B3    2
    ## 46     B3    2
    ## 47     B4    2
    ## 48     B3    2
    ## 49     B3    2
    ## 50     B1    1
    ## 51     B1    1
    ## 52     B3    2
    ## 53     B3    2
    ## 54     B3    2
    ## 55     B3    2
    ## 56     B3    2
    ## 57     B3    2
    ## 58     B3    2
    ## 59     B3    2
    ## 60     B3    2
    ## 61     B1    1
    ## 62     B2    1
    ## 63     B2    1
    ## 64     B1    1
    ## 65     B3    2
    ## 66     B4    2
    ## 67     B4    2
    ## 68     B2    1
    ## 69     B2    1
    ## 70     B3    2
    ## 71     B4    2
    ## 72     B4    2
    ## 73     B4    2
    ## 74     B2    1
    ## 75     B2    1
    ## 76     B2    1
    ## 77     B1    1
    ## 78     B2    1
    ## 79      B    1

``` r
fac1 <- factor(facB, levels=1:2, labels=c("B012", "B34"))
fac2 <- factor(ALL$mol.biol[w1 & w2])
table(fac1)
```

    ## fac1
    ## B012  B34 
    ##   48   31

``` r
table(fac2)
```

    ## fac2
    ## BCR/ABL     NEG 
    ##      37      42

``` r
tapply(ex12, fac1, mean)
```

    ##     B012      B34 
    ## 6.438820 6.190026

``` r
tapply(ex12, fac2, mean)
```

    ##  BCR/ABL      NEG 
    ## 6.552861 6.154722

``` r
col <- c("orange","blue")
xlab1 <- "B-cell ALL stage"
xlab2 <- "Molecular biology"
ylab <- "NEDD4 expression"
par(mfrow=c(1,2))
stripchart(ex12 ~ fac1, method="jitter", cex.lab=1.5, pch=19,
vertical=TRUE, xlab=xlab1, ylab=ylab, col=col)
stripchart(ex12 ~ fac2, method="jitter", cex.lab=1.5, pch=19,
vertical=TRUE, xlab=xlab2, ylab=ylab, col=col)
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

``` r
boxplot(ex12 ~ fac1, cex.lab=1.5, main=NULL, boxwex=0.3,
ylab=ylab, col=col, xlab=xlab1)
boxplot(ex12 ~ fac2, cex.lab=1.5, main=NULL, xlab=xlab2, ylab=ylab, col=col,boxwex=0.3)
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-71-2.png)<!-- -->

``` r
interaction.plot(fac1, fac2, ex12, type="b", col=col, xlab=xlab1,
pch=c(16, 17), lty=1, lwd=2, legend=F, ylab=ylab)
legend("topright", c("BCR/ABL","NEG"), bty="n", lty=1, lwd=2,
pch=c(16,17), col=col, inset=0.02)
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

``` r
interaction.plot(fac2, fac1, ex12, type="b", col=col, xlab=xlab2,
pch=c(16, 17), lty=1, lwd=2, legend=F, ylab=ylab)
legend("topright", c("B012","B34"), bty="n", lty=1, lwd=2,
pch=c(16,17), col=col, inset=0.02)
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-72-2.png)<!-- -->

``` r
tapply(ex12[fac2=="BCR/ABL"], fac1[fac2=="BCR/ABL"], mean)
```

    ##     B012      B34 
    ## 6.764918 6.241844

``` r
tapply(ex12[fac2=="NEG"], fac1[fac2=="NEG"], mean)
```

    ##     B012      B34 
    ## 6.162891 6.141447

``` r
tapply(ex12[fac1=="B012"], fac2[fac1=="B012"], mean)
```

    ##  BCR/ABL      NEG 
    ## 6.764918 6.162891

``` r
tapply(ex12[fac1=="B34"], fac2[fac1=="B34"], mean)
```

    ##  BCR/ABL      NEG 
    ## 6.241844 6.141447

``` r
anova(lm(ex12 ~ fac1 * fac2))
```

    ## Analysis of Variance Table
    ## 
    ## Response: ex12
    ##           Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## fac1       1  1.1659  1.1659  4.5999 0.0352127 *  
    ## fac2       1  3.2162  3.2162 12.6891 0.0006433 ***
    ## fac1:fac2  1  1.1809  1.1809  4.6592 0.0340869 *  
    ## Residuals 75 19.0094  0.2535                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(lm(ex12 ~ fac1 * fac2))
```

    ## 
    ## Call:
    ## lm(formula = ex12 ~ fac1 * fac2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.21905 -0.31589 -0.05137  0.24404  1.52015 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       6.7649     0.1073  63.026  < 2e-16 ***
    ## fac1B34          -0.5231     0.1686  -3.103   0.0027 ** 
    ## fac2NEG          -0.6020     0.1458  -4.128 9.39e-05 ***
    ## fac1B34:fac2NEG   0.5016     0.2324   2.159   0.0341 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5034 on 75 degrees of freedom
    ## Multiple R-squared:  0.2264, Adjusted R-squared:  0.1954 
    ## F-statistic: 7.316 on 3 and 75 DF,  p-value: 0.0002285

``` r
data.frame(fac1, fac2)
```

    ##    fac1    fac2
    ## 1  B012 BCR/ABL
    ## 2  B012     NEG
    ## 3   B34 BCR/ABL
    ## 4  B012     NEG
    ## 5  B012     NEG
    ## 6  B012     NEG
    ## 7  B012     NEG
    ## 8  B012     NEG
    ## 9  B012 BCR/ABL
    ## 10  B34 BCR/ABL
    ## 11  B34     NEG
    ## 12 B012     NEG
    ## 13  B34 BCR/ABL
    ## 14 B012     NEG
    ## 15 B012 BCR/ABL
    ## 16  B34 BCR/ABL
    ## 17 B012 BCR/ABL
    ## 18  B34 BCR/ABL
    ## 19 B012     NEG
    ## 20 B012 BCR/ABL
    ## 21 B012 BCR/ABL
    ## 22 B012     NEG
    ## 23 B012 BCR/ABL
    ## 24 B012     NEG
    ## 25 B012 BCR/ABL
    ## 26 B012     NEG
    ## 27 B012 BCR/ABL
    ## 28 B012     NEG
    ## 29 B012 BCR/ABL
    ## 30 B012 BCR/ABL
    ## 31 B012     NEG
    ## 32 B012 BCR/ABL
    ## 33 B012 BCR/ABL
    ## 34 B012 BCR/ABL
    ## 35 B012     NEG
    ## 36  B34 BCR/ABL
    ## 37 B012     NEG
    ## 38 B012     NEG
    ## 39 B012     NEG
    ## 40  B34 BCR/ABL
    ## 41 B012 BCR/ABL
    ## 42 B012 BCR/ABL
    ## 43  B34     NEG
    ## 44  B34     NEG
    ## 45  B34     NEG
    ## 46  B34     NEG
    ## 47  B34 BCR/ABL
    ## 48  B34 BCR/ABL
    ## 49  B34     NEG
    ## 50 B012     NEG
    ## 51 B012     NEG
    ## 52  B34     NEG
    ## 53  B34 BCR/ABL
    ## 54  B34     NEG
    ## 55  B34     NEG
    ## 56  B34     NEG
    ## 57  B34     NEG
    ## 58  B34     NEG
    ## 59  B34 BCR/ABL
    ## 60  B34 BCR/ABL
    ## 61 B012     NEG
    ## 62 B012     NEG
    ## 63 B012 BCR/ABL
    ## 64 B012 BCR/ABL
    ## 65  B34     NEG
    ## 66  B34     NEG
    ## 67  B34     NEG
    ## 68 B012     NEG
    ## 69 B012 BCR/ABL
    ## 70  B34     NEG
    ## 71  B34 BCR/ABL
    ## 72  B34 BCR/ABL
    ## 73  B34 BCR/ABL
    ## 74 B012     NEG
    ## 75 B012     NEG
    ## 76 B012 BCR/ABL
    ## 77 B012     NEG
    ## 78 B012 BCR/ABL
    ## 79 B012 BCR/ABL

``` r
nfac <- rep(1, length(fac1))
nfac[fac1=="B012" & fac2=="NEG"] <- 2
nfac[fac1=="B34" & fac2=="BCR/ABL"] <- 3
nfac[fac1=="B34" & fac2=="NEG"] <- 4
nfac <- factor(nfac, labels=c("f11", "f12", "f21", "f22"))
table(nfac)
```

    ## nfac
    ## f11 f12 f21 f22 
    ##  22  26  15  16

``` r
table(fac1, fac2)
```

    ##       fac2
    ## fac1   BCR/ABL NEG
    ##   B012      22  26
    ##   B34       15  16

``` r
res <- residuals(lm(ex12 ~ nfac))
shapiro.test(res)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  res
    ## W = 0.98762, p-value = 0.648

``` r
library(lmtest)
bptest(lm(ex12 ~ nfac), studentize=FALSE)
```

    ## 
    ##  Breusch-Pagan test
    ## 
    ## data:  lm(ex12 ~ nfac)
    ## BP = 5.9348, df = 3, p-value = 0.1148

``` r
anova(lm(ex12 ~ nfac))
```

    ## Analysis of Variance Table
    ## 
    ## Response: ex12
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## nfac       3  5.563 1.85432  7.3161 0.0002285 ***
    ## Residuals 75 19.009 0.25346                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov(lm(ex12 ~ nfac)))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = lm(ex12 ~ nfac))
    ## 
    ## $nfac
    ##                diff        lwr         upr     p adj
    ## f12-f11 -0.60202738 -0.9852329 -0.21882189 0.0005352
    ## f21-f11 -0.52307392 -0.9660222 -0.08012565 0.0140737
    ## f22-f11 -0.62347091 -1.0581108 -0.18883102 0.0018080
    ## f21-f12  0.07895345 -0.3499591  0.50786599 0.9624938
    ## f22-f12 -0.02144353 -0.4417703  0.39888328 0.9991323
    ## f22-f21 -0.10039699 -0.5758245  0.37503056 0.9449440

``` r
plot(TukeyHSD(aov(lm(ex12 ~ nfac))))
```

![](Lecture4---Linear-Models_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->
