################################################
## Function:	get.formula
## Version:		1.0		
## Purpose:		Parse & Organize ROC forumlas 
##				for an ROC model
##	 			(a, b, theta.0, theta.1)
##
## Author: 		JD Blume
## Date: 		May 2020 
################################################
#'
#' Parse and organize the 4 ROC parameter formulas
#'
#' @description This function takes the 4 formula inputs and organizes and formats them for the ROC analysis.
#'
#' @param f.1 Regression formula for parameters where the dependent varaiable is one of a, b, t0, or t1. For example, \code{a~x+z} and \code{b~1}. 
#' @param f.2 Same as \code{f.1}.
#' @param f.3 Same as \code{f.1}.
#' @param f.4 Same as \code{f.1}.
#'
#' @details Formulas can be in any order. Dependent variable in each formula must be one of a, b, theta.0, or theta.1. Duplicate dependent variables not allowed.
#'
#' @return ROC model formulas for ROC parameters a, b, theta.0, and theta.1.
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{www.statisticalevidence.com} \code{\link{smooth.roc}}
#' @keywords ROC curve, AUC
#' 
#' @export
#'
#' @examples
#' f.1=a~1 ; f.2=t1~size ; f.3=NULL ; f.4=b~age
#' get.formula(f.1, f.2, f.3, f.4)
#' 
#' f.1=a~1 ; f.2=t1~size ; f.3=NULL ; f.4=b~age
#' get.formula(f.1, f.2, f.3, f.4)
#' 

get.formula <- function(f.1=NULL, f.2=NULL, f.3=NULL, f.4=NULL, ...) {

#### Check formula statements for validity
	
	lhs.parms <- c( f.1[[2]], f.2[[2]], f.3[[2]], f.4[[2]] )

	if (anyDuplicated(lhs.parms)!=0) {

			stop("Duplicate parameters. Only one formula per parameter (a, b, t0, t1) is allowed.")
		}

	check.parms <- match( lhs.parms, c("a", "b", "t0", "t1") )

	if (anyNA(check.parms)) {

			stop("Error in formula: Parameters must be a, b, t0 or t1.")
		}
		
#### Identify and parse formula statements (keep placeholder for NULL formulas)
	
	modeled.parms <- c( if (is.null(f.1)) {NA} else {all.vars(f.1[[2]])},
						if (is.null(f.2)) {NA} else {all.vars(f.2[[2]])},
						if (is.null(f.3)) {NA} else {all.vars(f.3[[2]])},
						if (is.null(f.4)) {NA} else {all.vars(f.4[[2]])} )
		
	f.a <- f.b <- f.t0 <- f.t1 <- NULL

	if (!is.element("a" , modeled.parms)) {f.a  <- a ~ 1} else { 
		f.a <- get(paste0("f.", which(modeled.parms=="a"), collapse=""))
		}	

	if (!is.element("b" , modeled.parms)) {f.b  <- b ~ 1} else {
		f.b <- get(paste0("f.", which(modeled.parms=="b"), collapse=""))
		}

	if (!is.element("t0", modeled.parms)) {f.t0 <- t0 ~ 1} else {
		f.t0 <- get(paste0("f.", which(modeled.parms=="t0"), collapse=""))
		}

	if (!is.element("t1", modeled.parms)) {f.t1 <- t1 ~ 1} else {
		f.t1 <- get(paste0("f.", which(modeled.parms=="t1"), collapse=""))
		}

#### Output Formula (Full and reduced)

	outlist <- list(f.a=f.a, f.b=f.b, f.t0=f.t0, f.t1=f.t1)

outlist
}

###
##
#