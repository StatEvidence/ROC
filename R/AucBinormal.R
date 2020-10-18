################################################
## Function:	auc.binormal
## Version:		2.0		
## Purpose:		Estimate AUC from binormal model
##
## Model:		Sens = Phi(a + b * F^(-1)(1-spec))
##				H(scores.0; theta) ~ N(0,1)  
##
## Author: 		JD Blume
## Date: 		May 2020 
################################################


##
## level		confidence Interval level
##
## pauc			partial AUC (from, to)
## rescale		flag for rescaling the PACU 
##				by maximum area in rectangle
##
## parms		starting values for
##				(a, b, theta.0, theta.1)
################################################
#'
#' Area Under the Binormal ROC Curve
#'
#' @description Comptues the analytical estimate of the AUC for a binormal model, along with a delta-method standard error.
#'
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param f.1 Regression formula for parameters. For example, \code{a~x+z} and \code{b~1}. Order does not matter. 
#' @param f.2 Same as \code{f.1}.
#' @param f.3 Same as \code{f.1}.
#' @param f.4 Same as \code{f.1}.
#' @param model Model specification. Options are \code{"binormal"}, \code{"Bilogistic"}, \code{"as.defined"}. The last option is defined by the next 4 parameters.
#' @param dist.0 Distribution F (e.g. negatives or group 0). Options are \code{"normal"}, \code{"logistic"}, \code{"gamma"}.
#' @param dpar.0 Parameter list for distribution F
#' @param dist.1 Distribution G (e.g. positives or group 1).
#' @param dpar.1 Parameter list for distribution G
#' @param parms values for parameters \code{c(a, b, theta.0, theta.1)}
#'
#' @details Performs a maxmum likelihood fit of the ROC model.
#'
#' @return A list with the following elements:
#' 
#' \describe{
#' \item{\code{coef}}{Maximum likelihood estimates of the parameter vector.}
#' \item{\code{fit.var}}{Variance-covariance matrix for coefficients}
#' }
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{www.statisticalevidence.com} \code{\link{set.roc}}
#' @keywords ROC curve, AUC
#' 
#' @export
#'
#' @examples
#'

auc.binormal <- function(fit.obj, a, b, vcov, level=0.95, ...) {

#### Errors / Warnings

	if (!missing(fit.obj)) {
		if (is(fit.obj, "fitobj")!=TRUE) {stop("Object not FITOBJ class.")}
		}

#### Set default parameters

	if (missing(a)) {a <- fit.obj$a$coef['(Intercept)']}
		
	if (missing(b)) {b <- fit.obj$b$coef['(Intercept)']}
	
	if (missing(vcov)) {
		vcov <- fit.obj$vcov[c(1,length(fit.obj$a$coef)+1), c(1,length(fit.obj$a$coef)+1)]
			}

#### Parametric AUC and delta-method CI limits for 'binormal' case

		auc.bin <- pnorm(a/sqrt(1+b^2))

		V		<- vcov

		delta.dir <- dnorm(a/sqrt(1+b^2))*c(
						(1/sqrt(1+b^2)),
						a*(-b)*((1+b^2)^(-3/2))
														)
		auc.se  <- sqrt(t(delta.dir)%*%V%*%(delta.dir))

		z.alpha =  -qnorm((1-level)/2)
		
		auc.lo  <- max(auc.bin - z.alpha*auc.se, 0)
		auc.hi  <- min(auc.bin + z.alpha*auc.se, 1)
	
### Results table
	
		binorm <- c("Binormal", format(round(c(auc.bin, auc.lo, auc.hi), 4), nsmall = 2))

		tab <- data.frame(rbind(binorm))
		colnames(tab) <- c("AUC       ", "Estimate ", "CI Lo  ", "CI Hi")
						
#### Results list  ($$ not correct...a and b are not always in that place)

	outlist <- list(auc=auc.bin, auc.ci=c(auc.lo, auc.hi), se=auc.se, tab=tab, level=level)
	class(outlist) <- c("auc", class(outlist))

return(outlist)
}

###
##
#