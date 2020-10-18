################################################
## Function:	fitted.roc
## Version:		1.0
## Purpose:		Compute fitted values for ROC model
##	 			Parameters (a, b, theta.0, theta.1)
##
## Model:		Sens = G(a + b * F^(-1)(1-spec))
##				G,F ~ CDFs grp.1, grp.0 scores
##				H(scores.0; theta) ~ F  
##
## Current:		H() = theta.0 + theta.1*scores
##
## Author: 		JD Blume
## Date: 		May 2020 
################################################
#'
#' Compute fitted values for an ROC model
#'
#' @description This function will compute the fitted values for an ROC model.
#'
#' @param parms values for parameters \code{c(a, b, theta.0, theta.1)}
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param mm.list List of Model Matrix and numer of parameters for each formula. Obtained from \code{get.modelmatrix}.
#'
#' @details Computes fitted values, fitted variance, and MSE for ROC model.
#'
#' @return list of fittedvalues (NEED MORE HERE)
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{www.statisticalevidence.com} \code{\link{smooth.roc}}
#' @keywords ROC curve, AUC
#' 
#' @export
#'
#' @examples
#'

fitted.roc <- function(parms, roc.obj, mm.list, ...) {

#### Errors / Warnings

	if (is(roc.obj,"roc")!=TRUE) {stop("Error in fitted.roc: object not ROC class.")}
	
#### Set data for fitting
		
	data.all 	<- data.frame(roc.obj$data.all)	

#### Organize model matrix

	mm.a  <- mm.list$mm.a
	mm.b  <- mm.list$mm.b
	mm.t0 <- mm.list$mm.t0
	mm.t1 <- mm.list$mm.t1

	p.a   <- mm.list$p.a
	p.b   <- mm.list$p.b
	p.t0  <- mm.list$p.t0
	p.t1  <- mm.list$p.t1

#### Compute ROC Parameters 

	data.all$a <- 0
	data.all$b <- 1
	
	data.all[data.all$status==1,]$a	 <- mm.a  %*% parms[1:p.a]
	data.all[data.all$status==1,]$b	 <- mm.b  %*% parms[(p.a+1):(p.a+p.b)]

#### Compute transformation parameters 
	
	data.all$theta.0 <- mm.t0 %*% parms[(p.a+p.b+1):(p.a+p.b+p.t0)]
	data.all$theta.1 <- mm.t1 %*% parms[(p.a+p.b+p.t0+1):(p.a+p.b+p.t0+p.t1)]

#### Model Kernal/Transformation
#### Currently assumes WLOG that mu=0 and sigma=1

	data.all$yhat <- with(data.all, a/(b*theta.1) - theta.0/theta.1 )
	data.all$vhat <- with(data.all, 1/((b*theta.1)^2) )

	data.all$resid <- data.all$score-data.all$yhat

	mse <- mean(data.all$resid^2)
	rse <- sqrt(sum(data.all$resid^2)/(length(data.all$resid)-(p.a+p.b+p.t0+p.t1)))
	
#### Minus Total Log-Likelihood (optim finds minimum, so negate LogLike)

outlist <- list(yhat=data.all$yhat, vhat=data.all$vhat, resid=data.all$resid, 
					mse=mse, rse=rse)

outlist	
}

###
##
#