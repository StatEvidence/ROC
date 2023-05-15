ROC
========

`ROC` is a user-friendly R package for displaying and modeling ROC curves. 
The package contains tools for graphing, modeling and analyzing ROC data. 

Author
-------
Jeffrey D. Blume  
School of Data Science
University of Virginia
<i class="fas fa-envelope"></i>  blume@virginia.edu 
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

Here is an example using cancer serum biomarker data from 141 individuals. 
51 control patients and 90 cases with pancreatic cancer have data on two 
biomarkers (CA-125 and CA-19-9). See Wieand et al. (1989).

We read the data in and set-up two roc-objects.

``` r
library(ROC)

data(paca)

pcm.1 <- set.roc(log(paca$ca199), paca$case)
pcm.2 <- set.roc(log(paca$ca125), paca$case)
```

Historgams that display the roc data along with smoothed empirical densities.

``` r
hist(pcm.1, smooth(pcm.1, adj=0.8), breaks=seq(0, 12, length.out = 20))
hist(pcm.2, smooth(pcm.2, adj=1), breaks=seq(0, 12, length.out = 20))
```

Compute the (smoothed) ROC curve that corresponde to the smoothed empirical densities.

``` r
pcm1.smooth <- smooth(pcm.1, adj=0.8)
pcm2.smooth <- smooth(pcm.2, adj=1)
```

Plot empirical and smoothed ROC curves.

``` r
plot(pcm.1, points.col='dodgerblue3', flip.xaxis=TRUE)
plot(pcm.2, points.col='firebrick', show.plot=FALSE)

legend("right", c("Pa-Ca Marker CA-199","Pa-Ca Marker CA-125"),
				col=c('dodgerblue3','firebrick'),
				lty=1, lwd=2, pch=16, bty="n")

lines(pcm1.smooth, col="dodgerblue3", lwd=2)
lines(pcm2.smooth, col="firebrick", lwd=2)
```

Compute the empirical and smoothed area under the ROC curve. 

``` r
auc(pcm.1, level=0.95)
auc(pcm1.smooth, level=0.95)

auc(pcm.2, level=0.95)
auc(pcm2.smooth, level=0.95)
```

Fit the Binormal model to the ROC data. Add their curves to the plot.
 
``` r 
pcm1.binorm <- fit(pcm.1)
summary(pcm1.binorm)
auc(pcm1.binorm)

pcm2.binorm <- fit(pcm.2)
summary(pcm2.binorm)
auc(pcm2.binorm)

lines(curve(pcm1.binorm), col="dodgerblue3", lwd=2, lty=2)
lines(curve(pcm2.binorm), col="firebrick", lwd=2, lty=2)
```

The Binormal fits are not as good as the smoothed fit. This is because our fitting routine jointly 
fits the distribution of the biomarkers scores - not just the ROC curve - and the normal 
assumption is more consequential because the empirical distributions are skewed.

Collect the ROC parameter estimates for comparions across biomarkers. 

``` r 
hold <- rbind(cbind(a=pcm1.binorm$a$coef, b=pcm1.binorm$b$coef, 
					t0=pcm1.binorm$t0$coef, t1=pcm1.binorm$t1$coef, 
					auc=auc(pcm1.binorm)$auc),
		      cbind(a=pcm2.binorm$a$coef, b=pcm2.binorm$b$coef, 
				  	t0=pcm2.binorm$t0$coef, t1=pcm2.binorm$t1$coef, 
					auc=auc(pcm2.binorm)$auc))
hold <- rbind(hold,(hold[2,]-hold[1,]))								
rownames(hold) <- c("M1", "M2", "Diff")
hold

lines(curve(pcm1.binorm), col="dodgerblue3", lwd=2, lty=2)
lines(curve(pcm2.binorm), col="firebrick", lwd=2, lty=2)
```

But rather than fitting two seperate ROC models, we can combine them using regression ideas.
To do this, we combine the data into one dataframe using an indicator variables.
 
```r
y=data.frame(score =c( log(paca$ca199), log(paca$ca125)), 
  			 status=c( paca$case,paca$case), 
  			 marker=c( rep("one",length(paca$ca199)), 
			 		rep("two",length(paca$ca125))) )

pcc   <- set.roc(y$score, y$status, data.frame("marker"=y$marker))
```

Then fit a fully specified model with 8 parameters.

```r
pcc.binorm <- fit(pcc, f.1 = a ~ marker, f.2 = b ~ marker,
					f.3 = t0 ~ marker, f.4 = t1 ~ marker)
summary(pcc.binorm)
hold # for comparison
auc(pcc.binorm)
```

To get the AUC for the second biomarker, we can use the predict.roc function.

```r
at = data.frame('marker'=factor(c('one', 'two'), levels=c('one','two')))
keep <- predict(pcc.binorm, at, level=0.95)

keep$auc
auc(pcm1.binorm)
auc(pcm2.binorm)

## Draw predicted curves from model
lines(curve(keep, caserow=1), col="purple", lwd=2)
lines(curve(keep, caserow=2), col="green", lwd=2)
```

To reduce the effect of modeling assumptions, consider modeling the ranks instead as a type of 
non-parametric approach.  

```r
pcm.1r <- set.roc(rank(log(paca$ca199)), paca$case)
pcm.2r <- set.roc(rank(log(paca$ca125)), paca$case)

pcm1.binormr <- fit(pcm.1r)
auc(pcm1.binormr)

pcm2.binormr <- fit(pcm.2r)
auc(pcm2.binormr)

lines(curve(pcm1.binormr), col="black", lwd=2, lty=2)
lines(curve(pcm2.binormr), col="black", lwd=2, lty=2)
```
 
Use a "cutplot" to show the relationship between Sensitivity and Specificity and biomarker score.

```r
cutplot(pcm.1, pcm1.binorm)
cutplot(pcm.2, pcm2.binorm)

cutplot(pcc, pcc.binorm, 
at = data.frame('marker'=factor('two', levels=c('one','two'))))

cutplot(pcc, pcc.binorm, 
at = data.frame('marker'=factor('one', levels=c('one','two'))))
par(mfrow=c(1,1))

x=cutplot(pcm.1, pcm1.binorm)
head(x$graph)
```

Or switch it up and try a bilogistic model.

```r
pcc.bilog <- fit(pcc, f.1 = a ~ marker, f.2 = b ~ marker,
					f.3 = t0 ~ marker, 
					f.4 = t1 ~ marker, model="bilogistic")
summary(pcc.binorm)

at = data.frame('marker'=factor(c('one', 'two'), levels=c('one','two')))
keep <- predict(pcc.bilog, at, level=0.95)
					
plot(pcm.1, points.col='dodgerblue3', flip.xaxis=TRUE)
plot(pcm.2, points.col='firebrick', show.plot=FALSE)

legend("right", c("Pa-Ca Marker CA-199","Pa-Ca Marker CA-125"), bty="n",
				col=c('dodgerblue3','firebrick'),
				lty=1, lwd=2, pch=16)

lines(curve(keep, caserow=1), col="black", lwd=2)
lines(curve(keep, caserow=2), col="black", lwd=2)	

lines(curve(pcm1.binorm), col="dodgerblue3", lwd=2, lty=2)
lines(curve(pcm2.binorm), col="firebrick", lwd=2, lty=2)				
```

References
----------

There should be some references...really, there should.  


