mcrPioda is a fork of the mcr package with additional functionalities:

 - MDeming regression - With bootstrap CI and jackknife CI + SE 
 - MMDeming regression - With bootstrap CI and jackknife CI + SE
 - NgMMDeming regression - With bootstrap CI and jackknife CI + SE
 - PiMMDeming regression - With bootstrap CI and jackknife CI + SE
 - plotBoxEllipses for Mahalanobis distance hypothesis testing
 
All regression functions are written in C out of the MM-Deming which is kept
for reproducibility and should be considered deprecated.

There is an urgent need for M-Deming since **all Passing Bablok regression are biased** with
low precision data sets (especially with 2 and 3 significant digits only).

Power testing for the jackknife and bootstrap Cis are ongoing. Jackknife could be an
attractive alternative to bootstrap for small sample size.

For the smallest samples the Bayesian Deming regression can be a better option.
The same is true with heteroscedastic data sets. Check package *rstanbdp*.

Worth mentioning that M-Deming relies on the very same recursive method proposed by Linnet
the for WDeming, just with an M- algorithm for robust weight.

MM-Deming algorithm is more complex. It is also a recursive method but relies on MM- methods and
uses bi-square re-descending weights. The new algorithms NgMM- and PiMM- refresh the mad() dispersion of the residuals at each iterations, making the end results much less sensible to the starting values.

Reference: [http://arxiv.org/pdf/2105.04628.pdf](http://arxiv.org/pdf/2105.04628)
See also: [https://ssmtstatistica.wordpress.com/2024/01/27/equivariant-passing-bablok-regression-27-01-2024/](https://ssmtstatistica.wordpress.com/2024/01/27/equivariant-passing-bablok-regression-27-01-2024/)
