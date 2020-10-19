ROC
========

`ROC` is a user-friendly R package for displaying and modeling ROC curves. 

Author
-------
Jeffrey D. Blume  
Vanderbilt University  
Professor of Biostatistics, Biomedical Informatics and Biochemistry  
Vice Chair for Education, Biostatistics  
Director of Graduate Education, Data Science Institute  
<i class="fas fa-envelope"></i>  j.blume@vanderbilt.edu  
https://www.statisticalevidence.com/

News
----
Version 0.1.0

Installation
------------

``` r
# install.packages("devtools")
devtools::install_github("statevidence/ROC")
```

Help
----

Use the `help()` function to see the full list of functions and documentation. 

``` r 
help(package='ROC')
```

Example
-------

The packages contains a set of functions for graphing, modeling and analyzing ROC curves and the area under them - the AUC.

``` r
library(ROC)

data(paca)

pcm.1 <- set.roc(log(paca$ca199), paca$case)
pcm.2 <- set.roc(log(paca$ca125), paca$case)


## Historgams of Data with smoothing
hist(pcm.1, smooth(pcm.1, adj=0.8), breaks=seq(0, 12, length.out = 20))
hist(pcm.2, smooth(pcm.2, adj=1), breaks=seq(0, 12, length.out = 20))


## Smoothed ROC
pcm1.smooth <- smooth(pcm.1, adj=0.8)
pcm2.smooth <- smooth(pcm.2, adj=1)


## Plot empirical and smoothed ROC curve
plot(pcm.1, points.col='dodgerblue3', flip.xaxis=TRUE)
plot(pcm.2, points.col='firebrick', show.plot=FALSE)

legend("right", c("Pa-Ca Marker CA-199","Pa-Ca Marker CA-125"),
				col=c('dodgerblue3','firebrick'),
				lty=1, lwd=2, pch=16, bty="n")

lines(pcm1.smooth, col="dodgerblue3", lwd=2)
lines(pcm2.smooth, col="firebrick", lwd=2)

```

explain

``` r
## Empirical and Smoothed AUCs
auc(pcm.1, level=0.95)
auc(pcm1.smooth, level=0.95)

auc(pcm.2, level=0.95)
auc(pcm2.smooth, level=0.95)
```

explain 

``` r 
## Binormal fit ROCs
pcm1.binorm <- fit(pcm.1)
summary(pcm1.binorm)
auc(pcm1.binorm)

pcm2.binorm <- fit(pcm.2)
summary(pcm2.binorm)
auc(pcm2.binorm)

hold <- rbind(cbind(a=pcm1.binorm$a$coef, b=pcm1.binorm$b$coef, 
		t0=pcm1.binorm$t0$coef, t1=pcm1.binorm$t1$coef, 
		auc=auc(pcm1.binorm)$auc),
		cbind(a=pcm2.binorm$a$coef, b=pcm2.binorm$b$coef, 
		t0=pcm2.binorm$t0$coef, t1=pcm2.binorm$t1$coef, 
		auc=auc(pcm2.binorm)$auc))
hold <- rbind(hold,(hold[1,]-hold[2,]))								
rownames(hold) <- c("M1", "M2", "Diff")
hold

lines(curve(pcm1.binorm), col="dodgerblue3", lwd=2, lty=2)
lines(curve(pcm2.binorm), col="firebrick", lwd=2, lty=2)
```
 explain
 
 ```r
 cutplot(pcm.1,pcm1.binorm)
 cutplot(pcm.2,pcm2.binorm)
 ```
 
 explain
 
 ```r
 ## Using ranks

 pcm.1r <- set.roc(rank(log(paca$ca199)), paca$case)
 pcm.2r <- set.roc(rank(log(paca$ca125)), paca$case)

 pcm1.binormr <- fit(pcm.1r)
 auc(pcm1.binormr)

 pcm2.binormr <- fit(pcm.2r)
 auc(pcm2.binormr)

 lines(curve(pcm1.binormr), col="black", lwd=2, lty=2)
 lines(curve(pcm2.binormr), col="black", lwd=2, lty=2)
 ```


References
----------

There should be some references...eventually.  


