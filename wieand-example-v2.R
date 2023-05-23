########################################################
##
## Wieand (1989) Pancreatic Cancer Marker Data
##
########################################################

##### libraries
library(devtools)

#### Load package from Github
devtools::install_github("statevidence/ROC")
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


## Empirical
plot(pcm.1, points.col='dodgerblue3', flip.xaxis=TRUE)
plot(pcm.2, points.col='firebrick', show.plot=FALSE)

legend("right", c("Pa-Ca Marker CA-199","Pa-Ca Marker CA-125"), bty="n",
				col=c('dodgerblue3','firebrick'),
				lty=1, lwd=2, pch=16)

lines(pcm1.smooth, col="dodgerblue3", lwd=2)
lines(pcm2.smooth, col="firebrick", lwd=2)


## Empirical and Smoothed AUCs
auc(pcm.1, level=0.95)
auc(pcm1.smooth, level=0.95)

auc(pcm.2, level=0.95)
auc(pcm2.smooth, level=0.95)

cutplot(pcm.1,smooth=TRUE)
cutplot(pcm.1,pcm1.binorm) ## need to run code below first for binorm fit.

## Binormal fit ROCs
pcm1.binorm <- fit(pcm.1)
summary(pcm1.binorm)
auc(pcm1.binorm)

pcm2.binorm <- fit(pcm.2)
summary(pcm2.binorm)
auc(pcm2.binorm)

hold <- rbind(cbind(a=pcm1.binorm$a$coef, b=pcm1.binorm$b$coef, 
					t0=pcm1.binorm$t0$coef, t1=pcm1.binorm$t1$coef, auc=auc(pcm1.binorm)$auc),
		      cbind(a=pcm2.binorm$a$coef, b=pcm2.binorm$b$coef, 
				  	t0=pcm2.binorm$t0$coef, t1=pcm2.binorm$t1$coef, auc=auc(pcm2.binorm)$auc))
hold <- rbind(hold,(hold[2,]-hold[1,]))								
rownames(hold) <- c("M1", "M2", "Diff")
hold

lines(curve(pcm1.binorm), col="dodgerblue3", lwd=2, lty=2)
lines(curve(pcm2.binorm), col="firebrick", lwd=2, lty=2)


## Using ranks

pcm.1r <- set.roc(rank(log(paca$ca199)), paca$case)
pcm.2r <- set.roc(rank(log(paca$ca125)), paca$case)

pcm1.binormr <- fit(pcm.1r)
auc(pcm1.binormr)

pcm2.binormr <- fit(pcm.2r)
auc(pcm2.binormr)

lines(curve(pcm1.binormr), col="black", lwd=2, lty=2)
lines(curve(pcm2.binormr), col="black", lwd=2, lty=2)

## 

x= cutplot(pcm.1,pcm1.binorm, plot=FALSE)

head(x$graph)

cutplot(pcm.1,pcm1.binorm)
cutplot(pcm.2,pcm2.binorm)

par(mfrow=c(1,1))

## Combined ROC model for both markers

y=data.frame(score =c( log(paca$ca199), log(paca$ca125)), 
			 status=c( paca$case,paca$case), 
			 marker=factor(c( rep("one",length(paca$ca199)), rep("two",length(paca$ca125)))) )

pcc   <- set.roc(y$score, y$status, data.frame("marker"=y$marker))
pcc.r <- set.roc(rank(y$score), y$status, data.frame("marker"=y$marker))

## Pool AUC and plot

hist(pcc, breaks=seq(0, 12, length.out = 20))
hist(pcc, smooth(pcc, adj=0.8), breaks=seq(0, 12, length.out = 20))

auc(pcc)
auc(pcc.r)

pcc.fit <- fit(pcc)
pcc.fitr <- fit(pcc.r)

plot(pcc)

lines(curve(pcc.fit), col="darkorange", lwd=2)
lines(curve(pcc.fitr), col="purple", lwd=2)

legend("right", c("Pooled Data", "Ranked", "Original"), bty="n",
				col=c('Black','purple','darkorange'),
				lty=c(NA,1,1), lwd=2, pch=c(16,NA,NA))


## Model Pooled ROC with indicators for Marker

pcc.binorm <- fit(pcc, f.1 = a ~ marker, f.2 = b ~ marker,
					f.3 = t0 ~ marker, f.4 = t1 ~ marker)

pcc.binorm
fit(pcc, f.1 = a ~ marker, f.2 = b ~ marker,
					f.3 = t0 ~ marker, f.4 = t1 ~ marker, method="Nelder-Mead", maxit=20000)
					
fit(pcc, f.1 = a ~ marker, f.2 = b ~ marker,
					f.3 = t0 ~ marker, f.4 = t1 ~ marker, method="SANN", maxit=20000)
pcc.fit

auc(pcc.binorm)
auc(pcc.fit)

hold.all <- rbind(	cbind(a=pcc.binorm$a$coef['(Intercept)'], b=pcc.binorm$b$coef['(Intercept)'], 
						t0=pcc.binorm$t0$coef['(Intercept)'], t1=pcc.binorm$t1$coef['(Intercept)'], auc=NA),
					cbind(a=sum(pcc.binorm$a$coef), b=sum(pcc.binorm$b$coef), 
						t0=sum(pcc.binorm$t0$coef), t1=sum(pcc.binorm$t1$coef), auc=NA),
			  	  	cbind(a=pcc.binorm$a$coef['markertwo'], b=pcc.binorm$b$coef['markertwo'], 
				  		t0=pcc.binorm$t0$coef['markertwo'], t1=pcc.binorm$t1$coef['markertwo'], auc=NA))
rownames(hold.all) <- c("M1", "M2", "Diff")
hold
hold.all


## 

at = data.frame('marker'=factor(c('one', 'two'), levels=c('one','two')))
keep <- predict(pcc.binorm, at, level=0.95)

keep$auc
auc(pcm1.binorm)
auc(pcm2.binorm)

plot(pcm.1, points.col='dodgerblue3', flip.xaxis=TRUE)
plot(pcm.2, points.col='firebrick', show.plot=FALSE)

legend("right", c("Pa-Ca Marker CA-199","Pa-Ca Marker CA-125"), bty="n",
				col=c('dodgerblue3','firebrick'),
				lty=1, lwd=2, pch=16)

lines(pcm1.smooth, col="dodgerblue3", lwd=2)
lines(pcm2.smooth, col="firebrick", lwd=2)

lines(curve(pcm1.binorm), col="black", lwd=2, lty=2)
lines(curve(pcm2.binorm), col="black", lwd=2, lty=2)

lines(curve(keep, caserow=1), col="purple", lwd=2)
lines(curve(keep, caserow=2), col="green", lwd=2)

cutplot(pcc, pcc.binorm)

cutplot(pcc, pcc.binorm, 
	at = data.frame('marker'=factor('two', levels=c('one','two'))))
	
cutplot(pcc, pcc.binorm, 
	at = data.frame('marker'=factor('one', levels=c('one','two'))))

## Model Pooled ROC with indicators for Marker
pcc.bilog <- fit(pcc, f.1 = a ~ marker, f.2 = b ~ marker,
					f.3 = t0 ~ marker, f.4 = t1 ~ marker, model="bilogistic")

pcc.bilog

at = data.frame('marker'=factor(c('one', 'two'), levels=c('one','two')))
keep <- predict(pcc.bilog, at, level=0.95)

par(mfrow=c(1,1))					
plot(pcm.1, points.col='dodgerblue3', flip.xaxis=TRUE)
plot(pcm.2, points.col='firebrick', show.plot=FALSE)

legend("right", c("Pa-Ca Marker CA-199","Pa-Ca Marker CA-125"), bty="n",
				col=c('dodgerblue3','firebrick'),
				lty=1, lwd=2, pch=16)

lines(curve(keep, caserow=1), col="purple", lwd=2)
lines(curve(keep, caserow=2), col="green", lwd=2)

cutplot(pcc, pcc.bilog)

####

####					
###
##
#