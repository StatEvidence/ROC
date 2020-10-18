################################################
## Function:	get.modelmatrix
## Version:		1.0		
## Purpose:		Obtain model matricies 
##				for an ROC model
##	 			(a, b, theta.0, theta.1)
##
## Author: 		JD Blume
## Date: 		May 2020 
################################################
#'
#' Obtain 4 ROC model matricies and ranks
#'
#' @description This function compute the model matrix for the 4 ROC formulas and computes rank etc.
#'
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param f.a Regression formula for a obtained from \code{get.formula} function.
#' @param f.b Regression formula for b obtained from \code{get.formula} function.
#' @param f.t0 Regression formula for t0 obtained from \code{get.formula} function.
#' @param f.t1 Regression formula for t1 obtained from \code{get.formula} function.
#'
#' @details Formulas are specific, can not be \code{NULL} nor \code{NA}.
#'
#' @return Model matrix for the 4 ROC formulas for a, b, theta.0, theta.1 and the parameter count for each formula statment. 
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

get.modelmatrix <- function(roc.obj, f.a, f.b, f.t0, f.t1,...) {

#### Errors / Warnings

	if (!is(roc.obj,"roc")) {
			stop("Object not ROC class.") 
		}
		
	if (any(is.null(f.a), is.null(f.b), is.null(f.t0), is.null(f.t1))) {
			stop("Model matrix: at least one formula is NULL.") 
		}

	is.formula <- function(x){inherits(x,"formula")}

	if (any(!is.formula(f.a), !is.formula(f.b), !is.formula(f.t0), !is.formula(f.t1))) {
			stop("Model matrix: one formula input is not of formula class.") 
		}

#### Set data for fitting

	data.all 	<- data.frame(roc.obj$data.all)	
	data.all.0 	<- data.frame(roc.obj$data.all[roc.obj$data.all$status==0,])
	data.all.1 	<- data.frame(roc.obj$data.all[roc.obj$data.all$status==1,])

#### remove response from formula

	f.a[[2]]  <- NULL
	f.b[[2]]  <- NULL
	f.t0[[2]] <- NULL
	f.t1[[2]] <- NULL
	
#### Set model matricies, get ranks and names

	mm.a  <- model.matrix(f.a, data.all.1)
	mm.b  <- model.matrix(f.b, data.all.1)
	mm.t0 <- model.matrix(f.t0, data.all) 
	mm.t1 <- model.matrix(f.t1, data.all) 

	p.a   <- dim(mm.a)[[2]]
	p.b   <- dim(mm.b)[[2]]
	p.t0  <- dim(mm.t0)[[2]]
	p.t1  <- dim(mm.t1)[[2]]

## Compute total log-likelihood		

	outlist <- list(p.a=p.a, p.b=p.b, p.t0=p.t0, p.t1=p.t1,
					mm.a=mm.a, mm.b=mm.b, mm.t0=mm.t0, mm.t1=mm.t1)
				
outlist
}

###
##
#