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
summary(g)$sigma^2 # MSE
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
print("====================================================")
```

    ## [1] "===================================================="

``` r
names(summary(g))
```

    ##  [1] "call"          "terms"         "residuals"     "coefficients" 
    ##  [5] "aliased"       "sigma"         "df"            "r.squared"    
    ##  [9] "adj.r.squared" "fstatistic"    "cov.unscaled"
