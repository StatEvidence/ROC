################################################
## Function:	like.roc
## Version:		2.0		(old version = ROCLike.R)
## Purpose:		Compute Likelihood for ROC model
##	 			Parameters (a, b, theta.0, theta.1)
##
## Model:		Sens = G(a + b * F^(-1)(1-spec))
##				G,F ~ CDFs grp.1, grp.0 scores
##				H(scores.0; theta) ~ F  
##
## Current:		H() = theta.0 + theta.1*scores
##
## Author: 		JD Blume
## Date: 		April 2020 
################################################
#'
#' Compute Likelihood for ROC model
#'
#' @description This function will compute the likelihood for a given ROC model and data. Used for maximum likelihood computations.
#'
#' @param parms values for parameters \code{c(a, b, theta.0, theta.1)}
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param mm.list List of Model Matrix and numer of parameters for each formula. Obtained from \code{get.modelmatrix}.
#' @param model Model type. Options are \code{"binormal"}, \code{"bilogistic"}, \code{"other"}. The last option is defined by the next 4 parameters.
#' @param dist.0 Distribution F (e.g. negatives or group 0). Options are \code{"normal"}, \code{"logistic"}, \code{"gamma"}.
#' @param dpar.0 Parameter list for distribution F
#' @param dist.1 Distribution G (e.g. positives or group 1).
#' @param dpar.1 Parameter list for distribution G
#'
#' @details Computes log likelihood for ROC model.
#'
#' @return negative log-likeliihood
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

like.roc <- function(parms, roc.obj, mm.list,
						model="binormal", 
						dist.0=NA, dist.1=NA,
						dpar.0=NA, dpar.1=NA,...) {

#### Errors / Warnings

	if (is(roc.obj,"roc")!=TRUE) {stop("Object not ROC class.")}
	
	
	if (model=="binormal") {
		
		dist.0 <- dist.1 <- "normal" 
		dpar.0 <- dpar.1 <- c(0,1)
		
		}
		
		
	if (model=="bilogistic") {
	
		dist.0 <- dist.1 <- "logistic" 
		dpar.0 <- dpar.1 <- c(0,1)
	
		}


	if (!is.element(model,c("binormal", "bilogistic", "other")))	{
		
		stop("Model must be 'binormal', 'bilogistic' or 'other'")
		
		}


	if (model=="other")	{
	
		if (!is.element(dist.0,c("normal", "logistic", "gamma")))	{
	
			stop("dist.0 must be 'normal', 'logistic' or 'gamma'")
	
			}

		if (!is.element(dist.1,c("normal", "logistic", "gamma")))	{

			stop("dist.0 must be 'normal', 'logistic' or 'gamma'")

			}
		
		if (sum(is.na(dpar.0))>0) {
			
			stop("dpar.0 not properly specified; parameter vector needed")
			
			}
	
		if (sum(is.na(dpar.0))>0) {
		
			stop("dpar.1 not properly specified; parameter vector needed")
		
			}
	
		}
			
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

#### ROC Parameters (discrimination)
	
	a  <- mm.a  %*% parms[1:p.a]
	b  <- mm.b  %*% parms[(p.a+1):(p.a+p.b)]

## Transformation function H() Parameters (Calibration)
	
	theta.0 <- mm.t0 %*% parms[(p.a+p.b+1):(p.a+p.b+p.t0)]
	theta.1 <- mm.t1 %*% parms[(p.a+p.b+p.t0+1):(p.a+p.b+p.t0+p.t1)]

#### Model Kernal/Transformation

	H.all <- (theta.0 + theta.1*data.all$score)
	
	H.0  <- H.all[data.all$status==0]
	H.1  <- H.all[data.all$status==1]

	dH.0 <- theta.1[data.all$status==0]
	dH.1 <- theta.1[data.all$status==1]

#### Individual Log-Likelihood Contributions

	g.0 <- NULL
	if (dist.0=="normal")  {
		g.0 <- dnorm(H.0, mean=dpar.0[1], sd=dpar.0[2]) * dH.0
		}
		
	if (dist.0=="logistic") {
		g.0 <- dlogis(H.0, location=dpar.0[1], scale=dpar.0[2]) * dH.0
		}
		
	if (dist.0=="gamma") 	{
		g.0 <- dgamma(H.0, shape=dpar.0[1], scale=dpar.0[2]) * dH.0
		}

	g.1 <- NULL	
	if (dist.1=="normal") 	{
		g.1 <- dnorm(b*H.1 - a, mean=dpar.1[1], sd=dpar.1[2]) * b*dH.1
		}
		
	if (dist.1=="logistic") {
		g.1 <- dlogis(b*H.1 - a, location=dpar.1[1], scale=dpar.1[2]) * b*dH.1
		}
		
	if (dist.0=="gamma") 	{
		g.1 <- dgamma(b*H.1 - a, shape=dpar.1[1], scale=dpar.1[2]) * b*dH.1
		}
	
#### Minus Total Log-Likelihood (optim finds minimum, so negate LogLike)
	loss <- NULL
	
## Clean up probabilities  

	g.0[g.0 < 0] <- 0 
	g.0[g.0 > 1] <- 1
	
	g.1[g.1 < 0] <- 0  
	g.1[g.1 > 1] <- 1
	
## Compute total log-likelihood		

	loss <- sum(log(g.0)) + sum(log(g.1))
	loss <- -loss 	

return(loss)	
}

###
##
#