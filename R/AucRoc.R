################################################
## Function:	auc.roc
## Version:		2.0 
## Purpose:		Compute empirical area under ROC curve
##				Compute CIs for AUC (approx and robust)
## 
## Author: 		JD Blume
## Date: 		April 2020
################################################
#'
#' Empirical area under the ROC curve (AUC)
#'
#' @description This function computes the area under the empirical ROC curve. 
#'
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param level Confidence interval coverage level. Default is \code{95\%}.
#' @param digits Rounding accuracy for display. Default is \code{4}.
#'
#' @details Computes the empirical area under the ROC curve. This is also known as the Mann-Whitney-U-statistic and c-statistic. The reported AUC is actually \code{max(AUC, 1-AUC)}. This ensures the estimate is always greater than 0.5. Approximate and robust standard errors are used to obtain CIs.
#'
#' @return Returns the empirical estimate of the area under the ROC curve, along with an approximate and robust confidence intervals:
#' 
#' \describe{
#' \item{\code{auc}}{Area under the ROC curve (AUC).}
#' \item{\code{ci.mwu}}{Confidence interval based on Mann-Whitney U-statistic.}
#' \item{\code{ci.rob}}{Confidence interval based on distribution robust SEs.}
#' }
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Blume 2009 JSPI \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631183/} and \url{https://rdrr.io/cran/asht/man/wmwTest.html}
#' @seealso	\url{www.statisticalevidence.com} \code{\link{auc.curve}}
#' @keywords ROC curve, AUC
#' 
#' @export
#'
#' @examples
#'
#' ## Fake ROC data for examples 
#' y.0 <- round(rnorm(100, mean=0, sd=1)   , 1) 
#' y.1 <- round(rnorm(100, mean=1, sd=1.1) , 1)
#' myroc <- set.roc(score=y.0, grp=y.1, as.groups=TRUE)
#' 
#' ## Area under empirical ROC curve
#' auc(myroc)
#'
#' ## Empirical AUC from smoothing function
#' auc(smooth(myroc, adj=0.01))
#'

auc.roc <- function(roc.obj, level=0.95, digits=4, ...) {

#### Uses asht package from https://rdrr.io/cran/asht/man/wmwTest.html

#### Errors / Warnings

	if (is(roc.obj,"roc")!=TRUE) {stop("Object not ROC class.")}
		
#### Compute AUC
	y.0 <- roc.obj$data.0
	y.1 <- roc.obj$data.1

	w.min <- wilcox.test(y.0, y.1, exact=FALSE)$statistic/length(y.0)/length(y.1)
	area  <- max(w.min, 1-w.min)

#### Compute CIs 
##	'Robust' Varaince (Blume 2009): max varaince over all continuous distributions 

	var.upbd <- area*(1-area)/min(length(y.0),length(y.1))
	z.a		 <- -qnorm((1-level)/2)
	
	ci.rob <- area + c(-1,1) * z.a * sqrt(var.upbd)
	ci.rob <- pmin(ci.rob, 1)
	ci.rob <- pmax(ci.rob, 0)

##	Mann-Whitney-U based CIs for AUC (Asht package)

	wmw 	<- if (w.min > 0.5) {asht::wmwTest(y.1, y.0, conf.level=level)} else {asht::wmwTest(y.0, y.1, conf.level=level)}
	ci.mwu 	<- wmw$conf.int

	mwu <- c("Empirical", format(round(c(area, ci.mwu[1], ci.mwu[2]), 4), nsmall = 2))
	rob <- c("Robust", "  --", format(round(c(ci.rob[1], ci.rob[2]), 4), nsmall = 2))
	
	tab <- data.frame(rbind(mwu,rob))
	colnames(tab) <- c("AUC       ", "Estimate ", "CI Lo  ", "CI Hi")
	
	outlist <- list("auc"=area, "ci.mwu"=sort(ci.mwu), "ci.rob"=sort(ci.rob), tab=tab, level=level)
	class(outlist) <- c("auc", class(outlist))
	
	return(outlist) 
}

###
##
#