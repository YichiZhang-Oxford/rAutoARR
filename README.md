# rAutoARR

**Auto Adaptive Robust Regression R Package**

## Description

This R package implements the Alternating Gradient Descent, Alternating Gradient Descent with Barzilai-Borwein Method and Alternating Gradient Descent with Backtracking Method. It also includes the Huber Mean Estimation, Huber Covariance Matrix Estimation, Huber Regression and Adaptive Huber Regression from `R` library [FarmTest](https://CRAN.R-project.org/package=FarmTest), written by [Xiaoou Pan](https://www.math.ucsd.edu/~xip024/).

## Installation

Install `rAutoARR` from GitHub:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("YichiZhang-Oxford/rAutoARR")
library(rAutoARR)
```

## Functions

There are seven functions in this package:

-   `agd`: Alternating Gradient Descent
-   `agdBB`: Alternating Gradient Descent with Barzilai-Borwein Method
-   `agdBacktracking`: Alternating Gradient Descent with Backtracking Method
-   `huberMean`: Huber Mean Estimation
-   `huberCov`: Huber Covariance Matrix Estimation
-   `huberReg`: Huber Regression
-   `adaHuberReg`: Adaptive Huber Regression

## Examples 

First, we present an example of mean estimation about Huber and Alternating Gradient Descent related methods. We generate data from a log-normal distribution, which is asymmetric and heavy-tailed.

```r
library(rAutoARR)
n = 1000
X = rlnorm(n, 0, 1.5) - exp(1.5^2 / 2)
huberMean = huberMean(X)
agd = agd(X)
agdBB = agdBB(X)
```

Second, for each setting, we generate an independent sample of size *n = 100* and compute four mean estimators: the Sample Mean, the Huber estimator, the Alternating Gradient Descent estimator, and the Alternating Gradient Descent with Barzilai-Borwein Method. **Figure 1** displays the &alpha;-quantile of the estimation error, with &alpha; ranging from *0.5* to *1* based on 2000 simulations.

The four mean estimators perform almost identically for the normal data. For the heavy-tailed skewed distributions, the deviation of the sample mean from the population mean grows rapidly with the confidence level, in striking contrast to the DA-Huber estimator, the Alternating Gradient Descent estimator, and the Alternating Gradient Descent with Barzilai-Borwein Method.

![figure 1](https://github.com/YichiZhang-Oxford/rAutoARR/blob/main/example/figure_1.png)

**Figure 1**: Estimation error versus confidence level for the sample mean, the DA-Huber, and the Alternating Gradient Descent estimator, and the Alternating Gradient Descent with Barzilai-Borwein estimator based on *2000* simulations

Finally, in **Figure 2**, we examine the *99%*-quantile of the estimation error versus a distribution parameter measuring the tail behavior and the skewness. That is, for normal data we let &sigma; vary between *1* and *4*; for skewed generalized *t* distributions, we increase the shape parameter *q* from *2.5* to *4*; for the lognormal and Pareto distributions, the shape parameters &sigma; and &alpha; vary from *0.25* to *2* and *1.5* to *3*, respectively.

The DA-Huber, the Alternating Gradient Descent estimator, and the Alternating Gradient Descent with Barzilai-Borwein estimator show substantial improvement in the deviations from the population mean because the distribution tends to have heavier tails and becomes more skewed.

![figure 2](https://github.com/YichiZhang-Oxford/rAutoARR/blob/main/example/figure_2.png)

**Figure 2**: Empirical *99%*-quantile of the estimation error versus a parameter measuring
the tails and skewness for the sample mean, the DA-Huber, and the Alternating Gradient Descent estimator, and the Alternating Gradient Descent with Barzilai-Borwein estimator

## License

MIT

## Author(s)

Yichi Zhang <yichi.zhang@worc.ox.ac.uk>, Qiang Sun <qiang.sun@utoronto.ca>

## References

Sun, Q. (2021). Do we need to estimate the variance in robust mean estimation? [Paper](https://arxiv.org/pdf/2107.00118.pdf) 

Bose, K., Fan, J., Ke, Y., Pan, X. and Zhou, W.-X. (2020). FarmTest: An R package for factor-adjusted robust multiple testing. *R. J.* **12** 372-387. [Paper](https://journal.r-project.org/archive/2021/RJ-2021-023/index.html)

Fan, J., Ke, Y., Sun, Q. and Zhou, W.-X. (2019). FarmTest: Factor-adjusted robust multiple testing with approximate false discovery control. *J. Amer. Statist. Assoc.* **114** 1880-1893. [Paper](https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527700) 

Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. *J. Amer. Stat. Assoc.* **115** 254-265. [Paper](https://doi.org/10.1080/01621459.2018.1543124)

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2020). A new principle for tuning-free Huber regression. *Stat. Sinica* to appear. [Paper](https://www.math.ucsd.edu/~wez243/tfHuber.pdf)