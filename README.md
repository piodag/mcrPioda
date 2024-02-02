mcrPioda is a fork of the mcr package with additional functionalities:

 - MDeming regression
 - MMDeming regression
 - plotBoxEllipses for Mahalanobis distance hypothesis testing

The M-Deming and MM-Deming iterative procedures are still written in R and not compiled with C.

There is an urgent need for M-Deming since **all Passing Bablok regression are biased** when
the low precision data set (2 and 3 significant digits only).

Reference: [https://arxiv.org/pdf/2105.04628.pdf](https://arxiv.org/pdf/2105.04628.pdf)
See also: [https://ssmtstatistica.wordpress.com/2024/01/27/equivariant-passing-bablok-regression-27-01-2024/](https://ssmtstatistica.wordpress.com/2024/01/27/equivariant-passing-bablok-regression-27-01-2024/)
