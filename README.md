mcrPioda is a fork of the mcr package with additional functionalities:

 - MDeming regression - With bootstrap CI and jackknife CI + SE 
 - MMDeming regression - With bootstrap CI and jackknife CI + SE
 - plotBoxEllipses for Mahalanobis distance hypothesis testing

The M-Deming and MM-Deming iterative procedures are still written in R and not compiled with C.

There is an urgent need for M-Deming since **all Passing Bablok regression are biased** when
the low precision data set (2 and 3 significant digits only).

Power testing for the jackknife related methods are ongoing. Jackknife could be an
attractive alternative to bootstrap.

Worth mentioning that M-Deming relies on the very same recursive method proposed by Linnet
the for WDeming, just with an M- algorithm for robust weights.

MM-Deming algorithm is more complex. It is also a recursive method but relies on MM- methods. The most
experimental part of this algorithm is the determination of the starting values.
In this proposal (tested in 2021) the starting values for slope and intercept are calculated
as average of the two possible covariances (flipping X and Y) calculated with Rocke's method.
The subsequent iterative process is MM- like.

Reference: [https://arxiv.org/pdf/2105.04628.pdf](https://arxiv.org/pdf/2105.04628.pdf)
See also: [https://ssmtstatistica.wordpress.com/2024/01/27/equivariant-passing-bablok-regression-27-01-2024/](https://ssmtstatistica.wordpress.com/2024/01/27/equivariant-passing-bablok-regression-27-01-2024/)
