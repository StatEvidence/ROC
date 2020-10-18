################################################
## Function:	predict.roc
## Version:		1.0		
## Purpose:		Compute a.hat and b.hat for 
##				covariate pattern (with Var-Cov)
##
## Author: 		JD Blume
## Date: 		Sept 2020 
################################################
#'
#' Predict ROC intercept (a) and slope (b) for covariate pattern 
#'
#' @description This function computes the ROC intercept (a) and slop (b) for the given covaraite patterns. Fitted values for outcomes are already computed in rocfit object.
#'
#' @param fit.obj FITOBJ object from \code{fit.roc} function.
#' @param at dataframe of covariate profiles 
#' @param level Confidence Interval level
#'
#' @details Computes a contrast based on the a maxmum likelihood fit of the ROC model.
#'
#' @return A list with the following elements:
#' 
#' \describe{
#' \item{\code{a}}{Vector of a estimates}
#' \item{\code{b}}{Vector of b estimates}
#' \item{\code{vcov}}{Variance-covariance matrix for ROC coefficients}
#' \item{\code{auc}}{Vector of auc estimates}
#' \item{\code{auc.lo}}{Vector of lower auc CI estimates}
#' \item{\code{auc.hi}}{Vector of upper auc CI estimates}
#' \item{\code{t0}}{Vector of t0 estimates}
#' \item{\code{t1}}{Vector of t1 estimates}
#' }
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{www.statisticalevidence.com} \code{\link{set.roc}}
#' @keywords ROC curve, AUC, predict
#' 
#' @export
#'
#' @examples
#'

predict.fitobj <- function(fit.obj, at, level=0.95, ...) {


#### Errors / Warnings

	if (is(fit.obj, "fitobj")!=TRUE) {stop("Object not FITOBJ class.")}

	if (missing(at)) {stop("Predictors ('at') cannot be missing.")}

	
#### Set names and levels for covariate profile dataframe

	original <- fit.obj$data.all[, -c(1:2)]
		
	if ( NCOL(at)!=NCOL(original) ) {
		
		stop("Number of columns in 'at' does not match number of original predictors.") 
		
		}

	if ( any(names(at)!=names(original)) ) {
		
		warning("Predictors may be misaligned; 'at' column names do not match original data.")
	
		name.err <- rbind( at=c(names(at)), original=c(names(original)) )
		colnames(name.err) <- c(1:NCOL(at))
		name.err
		
		}

	if ( NCOL(original)>1 ) {

			for (i in 1:NCOL(original)) { 
		
				names(at)[i]  <- names(original)[i]		
				levels(at[,i]) <- levels(original[,i])
				
				}
				
		} else {
			
			names(at)  <- names(fit.obj$data.all)[3]
			levels(at[,1]) <- levels(original)
			
				}

#### Set design matrix for parameter (a, b, t0, t1) prediction 

	a.set <- model.matrix(fit.obj$a$formula[-2], data=at)

	b.set <- model.matrix(fit.obj$b$formula[-2], data=at)

	ab.end <- dim(a.set)[2] + dim(b.set)[2]

	t0.set <- model.matrix(fit.obj$t0$formula[-2], data=at)

	t1.set <- model.matrix(fit.obj$t1$formula[-2], data=at)


#### Compute predicted parameters (a, b, t0, t1) for each profile

	a.pred <- a.set%*%as.matrix(fit.obj$a$coef)

	b.pred <- b.set%*%as.matrix(fit.obj$b$coef)

	t0.pred <- t0.set%*%as.matrix(fit.obj$t0$coef)

	t1.pred <- t1.set%*%as.matrix(fit.obj$t1$coef)

	
#### Compute profile specific var-cov matrix for (a,b) predictions 

	M <- rbind( cbind( a.set, matrix(0, nrow=dim(a.set)[[1]], ncol=dim(b.set)[[2]] ) ),
			cbind( matrix(0, nrow=dim(b.set)[[1]], ncol=dim(a.set)[[2]] ), b.set ) )

	atlength <- dim(at)[[1]]

	rows.order <- c( rbind( c(1:atlength), c(1:atlength) + atlength ) )

	M <- M[rows.order,]


	## Compute global var-cov matrix for all (a,b)
	
	vcmat.all <- M %*% fit.obj$vcov[1:ab.end, 1:ab.end] %*% t(M)


	## Extract block diagonals vcov for each profile
	
	k <- 2 		# Block dimension (length of parameter vector)
	
	m.get <- diag(1:(atlength*2/k)) %x% matrix(1, k, k) 
	
	vcmat.each <- lapply(split(vcmat.all, m.get)[-1], matrix, k)

	if (atlength==1) {vcmat.each <- list(vcmat.all)}

#### Binormal AUC estimates for each (a,b) pair 

	auc.est <- auc.lo <- auc.hi <- NULL
	
	if (any(class(fit.obj)=="binormal")) {
	
		for (i in 1:atlength){

			hold <- NULL
			
			hold <- auc.binormal(	a=a.pred[i], b=b.pred[i],
									vcov=vcmat.each[[i]], level=level )

			auc.est[i] 	<- hold$auc
			auc.lo[i] 	<- hold$auc.ci[1]
			auc.hi[i] 	<- hold$auc.ci[2]
		
			}
	
		}
		
	if (any(class(fit.obj)=="bilogistic")) {

		for (i in 1:atlength){

			hold <- NULL
			
			hold <- list(model=fit.obj$model, a=data.frame(coef=a.pred[i]), b=data.frame(coef=b.pred[i]),
						vcov=vcmat.each[[i]])
						
			class(hold) <- c("fitobj", "bilogistic", "list")
			
			hold.crv 	<- curve(hold)		
			auc.est[i]  <- auc(hold.crv, level=level)$auc
						
			hold.crv$y  <- hold.crv$ci.lo
			auc.lo[i] 	<- auc(curve(hold), level=level)$auc
			
			hold.crv$y  <- hold.crv$ci.hi
			auc.hi[i] 	<- auc(curve(hold), level=level)$auc
	
			}

		}
		
	
#### Prepare output 

	outlist <- list(auc=data.frame(est=auc.est, lo=auc.lo, hi=auc.hi, at),
					parms=data.frame(a=a.pred, b=b.pred, t0=t0.pred, t1=t1.pred, at),
					vcov=vcmat.each, at=at,
					call=fit.obj$call, model=fit.obj$model, 
					records=NROW(at), features=NCOL(at), level=level)

	if (any(class(fit.obj)=="binormal")) {
	
	class(outlist) <- c("binormal", class(outlist)) 
		
		}
		
	if (any(class(fit.obj)=="bilogistic")) {

	class(outlist) <- c("bilogistic", class(outlist)) 
	
		}	 
	
	class(outlist) <- c("predict", class(outlist))

	invisible(outlist)
	
}







###
##
#