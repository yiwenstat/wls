
## WLS

### Introduction

The weighted leverage score (WLS) method is a model-free variable screening method for high dimensional data regression. It was developed by utilizing both the left and right singular matrices to evaluate predictors’ importance under the semi-parametric model assumption. We call the importance score of each predictor the WLS. Theoretically, the WLS is consistent in ranking and selecting predictors. That is, the weighted leverage scores of non-redundant predictors have higher rankings than redundant predictors even in the "large p, small n" case, and the screening procedure using WLS can consistently exclude all redundant predictors. Compared with other feature screening methods, the WLS screening approach provides great ﬂexibility in model speciﬁcation and achieves high computational efﬁciency. Thus, it can be applied to a wide range of large-scale and high dimensional data collected from many ﬁelds of science.

The package includes the following file

- F_wls.R: the function for calculating the weighted leverage score.


### Example

#### An example dataset
Let the sample size n=1000 and the number of predictor p=1500. We generated the example dataset using the following code.

```
library(MASS)
set.seed(1234)
n=1000
p=1500
x=mvrnorm(n, rep(0,p), diag(1,p))
y=x[,1] + x[,2] + x[,3] + x[,4] + x[,5] + x[,6] + rnorm(n)
```

#### Calculate the WLS
To calculate the WLS of each predictor, the ```dr``` package is required.
```
source("F_wls.R")
calcwls=wlsFunc(x, y)
```

The selected predictors are
```
calcwls$select

```
```

[1]    3    5    6    4    2    1 1478 1373   31  813  380  553   54 1198  703

[16] 1357  525  210  225  121  113  162  749  412  522 1029  619 1166  587  825

[31] 1483  795  581 1132 1454 1498  706  693  421  592   17  352 1288  933 1355

[46] 1472  287 1204  689   96  361  957  215  901   39 1434  701 1450  943  608

[61] 1143
```

### Reference

Zhong, W., Liu, Y., and Zeng, P. (2019). A model-free variable screening method based on leverage score.

