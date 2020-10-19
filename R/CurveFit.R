#' Curve fitobj
#' 
#' @export
#'
curve.fitobj <- function(fit.obj, level=0.95, acc=0.01, ...) {

################################################
## Function:	curve.roc
## Version:		1.0		(old version = ROCfn.R)
## Purpose:		Compute Sensitivity for fitted
##				at 'fp' based on (a.hat, b.hat)
##				and SE.hat
##
## Model:		Sens = G(a + b * F^(-1)(1-spec))
##				G,F ~ CDFs grp.1, grp.0 scores
##
## Author: 		JD Blume
## Date: 		April 2020 
################################################

################################################
## fit.obj		Fitted ROC object
##
## acc			accuracy/density along x-axis
##
################################################

#### Errors / Warnings

	if (is(fit.obj,"fitobj")!=TRUE) {stop("Object not fitobj class.")}

#### Set-up
	
	dist.0 <- fit.obj$model$dist.0
	dpar.0 <- fit.obj$model$dpar.0
	
	dist.1 <- fit.obj$model$dist.1
	dpar.1 <- fit.obj$model$dpar.1

	if (fit.obj$model$model=="binormal") {
		
		dist.0 <- "normal"
		dpar.0 <- c(0,1)
		dist.1 <- "normal"
		dpar.1 <- c(0,1)
		
	}

	## WARNING intercept case only
	
	a <- fit.obj$a$coef[1]
	b <- fit.obj$b$coef[1]
	
	vcov <- fit.obj$vcov[1:2, 1:2]  

#### Normal Quantile for CIs
 
	z.alpha <-  -qnorm((1-level)/2)

#### Compute quantile for fp point (x-axis)

	fp <- seq(0.0001, 0.9999, acc)

	if (dist.0=="normal") 	{fp.dev <- qnorm(fp, mean=dpar.0[1], sd=dpar.0[2]) }
	if (dist.0=="logistic") {fp.dev <- qlogis(fp, location=dpar.0[1], scale=dpar.0[2]) }
	if (dist.0=="gamma") 	{fp.dev <- qgamma(fp, shape=dpar.0[1], scale=dpar.0[2]) } 

#### Compute ROC line and pointwise CIs (for sensitivity at given fp) 
#### 		in quantile-quantile space
		
	roc.line <- a + b*fp.dev
	line.sd  <- sqrt(vcov[1,1]+(fp.dev^2)*vcov[2,2]+2*fp.dev*vcov[1,2])

	roc.lo <- roc.line - z.alpha*line.sd
	roc.hi <- roc.line + z.alpha*line.sd

#### Transform line back into probabiltiy space
	
	if (dist.1=="normal") {
		tp.fit <- pnorm(roc.line, mean=dpar.1[1], sd=dpar.1[2])
		tp.lo  <- pnorm(roc.lo  , mean=dpar.1[1], sd=dpar.1[2])
		tp.hi  <- pnorm(roc.hi  , mean=dpar.1[1], sd=dpar.1[2])
							}
						
	if (dist.1=="logistic") {
		tp.fit <- plogis(roc.line, location=dpar.1[1], scale=dpar.1[2])
		tp.lo  <- plogis(roc.lo  , location=dpar.1[1], scale=dpar.1[2])
		tp.hi  <- plogis(roc.hi  , location=dpar.1[1], scale=dpar.1[2])
							}

	if (dist.1=="gamma") {
		tp.fit <- pgamma(roc.line, shape=dpar.1[1], scale=dpar.1[2])
		tp.lo  <- pgamma(roc.lo  , shape=dpar.1[1], scale=dpar.1[2])
		tp.hi  <- pgamma(roc.hi  , shape=dpar.1[1], scale=dpar.1[2])
							}		

### Prepare output

	fp <- c(0,fp,1)
	tp.fit <- c(0, tp.fit, 1)
	tp.lo  <- c(0, tp.lo, 1)
	tp.hi  <- c(0, tp.hi, 1)
	
	outlist <- list(x=fp, y=tp.fit, ci.lo=tp.lo, ci.hi=tp.hi, level=level)
	
	class(outlist) <- c("curve", class(outlist))	# added 10/12 (problems?)
	return(outlist)
}

###
##
#