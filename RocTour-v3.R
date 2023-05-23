################################################
## Tour of the ROC Package
################################################
##### libraries
library(devtools)

#### Load package from Github
devtools::install_github("statevidence/ROC")
library(ROC)

#################################
#### Generate fake ROC data

### Sample size
n 	  <- 200
pie.0 <- 0.5

n.0 <- ceiling(pie.0*n)
n.1 <- n-n.0

## Group Indicator
d 	 <- c(rep(0,n.0), rep(1,n.1))

## Covariates (age centered, age dichotomized at mean, factor)
age  <- rnorm(n, mean=35, sd=8)
age  <- age-mean(age)
over <- c(rep(0, ceiling(n/2)), rep(1, n-ceiling(n/2)))	
size <- sample(c("small", "medium", "big"), n, replace=TRUE)

## Model
mu 	<- ( 0 + 0.05*age + 0.25*(size=="medium") + 0.5*(size=="big") + 
			d*(1 + 0.05*age + 0.5*(size=="medium") + 1*(size=="big") ) )
sig <- 0.9 + 0.1*d + 0.1*over + 0.1*d*over

## Score
y	 <- round( rnorm(n, mean=mu, sd=sig) , 1)

## Dataframe for covariates
x <- data.frame(age, size, over) 


#################################
#### Set ROC object

myroc <- set.roc(y, d, x)
myroc

with(myroc, cbind(cutpoint, fp, fp.ci, tp, tp.ci))

## Fraction of group 1 scores at or above the overall median score
cut.med <- median(myroc$data.all$score)
with(myroc, cbind(cutpoint, tp, tp.ci))[myroc$cutpoint>cut.med,][1,]


#################################
#### Plotting the ROC curve

## Show empirical CI region ; flip x-axis
plot(myroc, flip.xaxis=TRUE, ci.region=TRUE)

## Add smoothed curve (from kernal smoothed densities)
myroc.smooth <- smooth(myroc, adj=0.8)
lines(myroc.smooth, col="purple", lwd=2)

## Add fitted curves 
binorm.fit <- fit(myroc)
binorm.crv <- curve(binorm.fit)

## Empirical
plot(myroc)
lines(binorm.crv, col="blue", lwd=1.5)
with(binorm.crv, matlines(x, cbind(ci.lo,ci.hi), col="firebrick", lty=1, lwd=1.5))


#################################
### Obtain the Area under the ROC curve

auc(myroc, level=0.95)
auc(binorm.fit, level=0.95)
auc(smooth(myroc))


#################################
### Histogram of data

hist(myroc, smooth(myroc, adj=1), breaks=seq(-30, 30, length.out = 81))

hist(myroc, myroc.smooth, binorm.fit, breaks=seq(-30, 30, length.out = 81))


#################################
### Fit an ROC curve

summary(binorm.fit)


## Fitting problems... (Add over to see)
mymod = fit.roc(roc.obj = myroc, f.1 = a ~ age + size, 
				f.2 = t0 ~ age + size, f.3 = b ~ 1, 
				f.4 = t1 ~ 1)
				
summary(mymod)
auc(mymod)


#################################
### Predict for covarite profiles 

spacing <- seq(floor(range(age)[1]), ceiling(range(age)[2]), 1)
at <- NULL

for (i in 1:length(spacing)) {

at <- rbind( at,
		data.frame(age=spacing[i],  size='small', over=0)	
			)
	}


keep <- predict(mymod, at, level=0.95)
keep

########################################################
## plotting

plot(0.5, 0.5, type="n",
		ylim=c(0,1),xlim=c(floor(min(age)),ceiling(max(age))), 
		ylab="AUC", xlab="Age")
lines(keep$auc$age, keep$auc$est, lwd=2, col="blue")
lines(keep$auc$age, keep$auc$lo, lwd=2, col="red")
lines(keep$auc$age, keep$auc$hi, lwd=2, col="red")
abline(h=0.5, lty=2, lwd=1.5, col="gray57")

legend( 'bottomright', c("AUC(age)", "95% CI"), col=c("blue","red"),lwd=2, lty=1, bty="n" )


#### For observed data set...
at = mymod$data.all[,-(1:2)]
keep <- predict(mymod, at, level=0.95)

plot(0.5, 0.5, type="n",
		ylim=c(0,1),xlim=c(floor(min(age)),ceiling(max(age))), 
		ylab="AUC", xlab="Age")
abline(h=0.5, lty=2, lwd=1.5, col="gray57")

with(keep$auc[over==1,],
{points(age[size=="big"], est[size=="big"], pch=16, col="blue",cex=1.2)
points(age[size=="medium"], est[size=="medium"], pch=16, col="forestgreen",cex=1.2)
points(age[size=="small"], est[size=="small"], pch=16, col="firebrick",cex=1.2)})

with(keep$auc[over==0,],
{points(age[size=="big"], est[size=="big"], pch=17, col="blue",cex=1.2)
points(age[size=="medium"], est[size=="medium"], pch=17, col="forestgreen",cex=1.2)
points(age[size=="small"], est[size=="small"], pch=17, col="firebrick",cex=1.2)})

legend( 'bottomright', c("Big", "Medium", "Small","Over=1","Over=0"), 
		col=c("blue", "forestgreen", "firebrick", "black", "black"),
		lwd=c(2,2,2,NA,NA), lty=c(1,1,1,NA,NA), pch=c(NA,NA,NA,17,16),bty="n" )


########################################################
## Showing the influence of cutpoints

cutplot(myroc)

try = fit.roc(roc.obj = myroc)
cutplot(myroc, try)

cutplot(myroc, mymod, at=data.frame(age=15, size= 'large', over= 1))


########################################################
## Plotting specific predicted ROC curves


plot(myroc)
zz <- curve(keep, caserow=1)
xx <- curve(keep, caserow=200)

lines(zz, col="purple", lwd=1.5)
lines(xx, col="green", lwd=1.5)


######
# sapply(keep,'[[', 1)

###
##
#





