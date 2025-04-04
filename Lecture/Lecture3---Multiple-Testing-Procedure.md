Lecture3 : Multiple Testing Procedure
================
KIM SANG HYUN(202211545)
2025-04-01

- [Introduction to Multiple Testing](#introduction-to-multiple-testing)
- [Type l and Type ll Errors](#type-l-and-type-ll-errors)
- [Multiple Testing Problems](#multiple-testing-problems)

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
library(ISLR2)
library(qqman)
```

    ## 

    ## For example usage please run: vignette('qqman')

    ## 

    ## Citation appreciated but not required:

    ## Turner, (2018). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. Journal of Open Source Software, 3(25), 731, https://doi.org/10.21105/joss.00731.

    ## 

``` r
library(RColorBrewer)
```

## Introduction to Multiple Testing

- 기존의 검정은 single null hypothesis를 testing을 하였다.

- 이번에는 huge amounts of data에서 **great many null hypotheses**
  test를 진행한다.

- We might want to test $m$ null hypotheses.

$$
H_{01}, H_{02},H_{03},H_{04},...,H_{0m}
$$

- $H_{01}$ -\> 1 번째 gene에 대해서 ALL, AML 두 집단의 평균이 같을
  것이다…

- 이런 검정을 모든 Gene에 대해서 다 진행을 한다.

- 여기서 문제를 우리가 이런 multiple testing에 대한 해석을 어떻게 할
  것인가 이다! (즉, error control를 어떻게 할까?)

- **False Discovery Rate (FDR)**로 조지자!

- FDR를 공부하기전에 빌드업 조지자…

## Type l and Type ll Errors

In single test…

- True 값은 우리가 알 수 없다.

|                         | **Truth(우리가 알 수 없는 값)** |                   |
|-------------------------|---------------------------------|-------------------|
| **Decision(판단 결과)** | **H₀ (True Null)**              | **H₁ (True Alt)** |
| Accept H₀               | ✅ (Correct)                    | ❌ Type II Error  |
| Reject H₀               | ❌ Type I Error                 | ✅ (Correct)      |

- 오류가 발생하는 Case

  - $H_0$가 참인데 우리의 데이터를 기반으론 P-value가 0.05보다 낮은
    상황이 나와서 $H_0$를 기각을 한다면 우리는 그것을 **Type I error**
    라고 한다.

    > $H_0$ is in fact true -\> reject $H_0$ -\> Type I error.

  - $H_0$가 거짓이고 이때는 무조건 기각을 해야하는데 못했어 엉엉 ㅠㅠ
    -\> **Type II error**

  - Ideally we would like both the Type I and Type II error rates to be
    small..

    > 하지만 이놈들은 Trade-off 관계쥬,,

![](images/clipboard-753359818.png)

- α = 0.05는 **1종 오류의 허용 수준**.

- 다시 말해, H₀가 참일 때에도 **5% 확률로 H₀를 잘못 기각할 수 있음**을
  감수하겠다는 의미.

- **p-value = 0.03**은 **“귀무가설이 참일 때 이 정도 극단적인 데이터가
  나올 확률이 3%”** 라는 뜻이지,

- **“내가 틀릴 확률이 3%”는 아니다!**

- α는 null hypothesis가 참이라는 가정 하에서, 검정 통계량이 우연히
  극단적인 값으로 나올 확률의 기준선(상한)“이다.

## Multiple Testing Problems
