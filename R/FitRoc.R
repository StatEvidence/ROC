################################################
## Function:	fit.roc
## Version:		2.0		
## Purpose:		Maximum Likelihood Fit ROC curve
##	 			Parameters (a, b, theta.0, theta.1)
##				Return Estimates and SEs
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
#' Maximum Likelihood Fit for ROC model
#'
#' @description This function will maximize the likelihood function for ROC parameters. Used to estimate the ROC parameters.
#'
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param f.1 Regression formula for parameters. For example, \code{a~x+z} and \code{b~1}. Order does not matter. 
#' @param f.2 Same as \code{f.1}.
#' @param f.3 Same as \code{f.1}.
#' @param f.4 Same as \code{f.1}.
#' @param model Model specification. Options are \code{"binormal"}, \code{"bilogistic"}, \code{"other"}. The last option is defined by the next 4 parameters.
#' @param dist.0 Distribution F (e.g. negatives or group 0). Options are \code{"normal"}, \code{"logistic"}, \code{"gamma"}.
#' @param dpar.0 Parameter list for distribution F
#' @param dist.1 Distribution G (e.g. positives or group 1).
#' @param dpar.1 Parameter list for distribution G
#' @param parms Values for parameters \code{c(a, b, theta.0, theta.1)}
#' @param method Optimization algorithm passed to \code{optim}, options \code{"Nelder-Mead"} or \code{"BFGS"} are good options.
#'
#' @details Performs a maxmum likelihood fit of the ROC model.
#'
#' @return A list with the following elements:
#' 
#' \describe{
#' \item{\code{coef}}{Maximum likelihood estimates of the parameter vector.}
#' \item{\code{vcov}}{Variance-covariance matrix for coefficients}
#' \item{\code{more}}{more here..!}
#' \item{\code{adjr2}}{slightly modified adjusted R^2}
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

fit.roc <- function(roc.obj, f.1=NULL, f.2=NULL, f.3=NULL, f.4=NULL, 
						model="binormal",   
						dist.0=NA, dist.1=NA,
						dpar.0=NA, dpar.1=NA, 
						start.parms=NA, method="BFGS", maxit=500, ...) {
							
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

#### Capturing elements for display

	cl <- match.call()

	flist <- get.formula(f.1, f.2, f.3, f.4)
	
	mmats <- get.modelmatrix(roc.obj, f.a=flist$f.a, f.b=flist$f.b, 
						f.t0=flist$f.t0, f.t1=flist$f.t1)

	where.par 	<- with(mmats, 
						c(rep("a", mmats$p.a ), rep("b", mmats$p.b ), 
							rep("t0", mmats$p.t0 ), rep("t1", mmats$p.t1 ))
						)

#### Starting values (uses unbiased estimate of var and not MLE)

	if ( any(is.na(start.parms)) ) {
	
	a.start  <- with(roc.obj$data.all, (mean(score[status==1])-mean(score[status==0]))/sd(score[status==1]) )
	b.start  <- with(roc.obj$data.all, sd(score[status==0])/sd(score[status==1]) )
	t0.start <- with(roc.obj$data.all, -mean(score[status==0])/sd(score[status==0]))
	t1.start <- with(roc.obj$data.all, 1/sd(score[status==0]))

	start.parms <- c( a.start,  rep(0,    mmats$p.a-1 ),
					  b.start,	rep(0.01, mmats$p.b-1 ),
					  t0.start, rep(0,    mmats$p.t0-1),
					  t1.start, rep(0.01, mmats$p.t1-1) )
	}

#### Maximum likelihood Fit

	mlfit <- optim(par=start.parms, like.roc, roc.obj=roc.obj, 
					f.a=flist$f.a, f.b=flist$f.b, f.t0=flist$f.t0, f.t1=flist$f.t1,
					mm.list=mmats,
					model=model, 
					dist.0=dist.0, dpar.0=dpar.0,
					dist.1=dist.1, dpar.1=dpar.1, 
					hessian=TRUE, method=method, control=list(maxit=maxit), ...)
	
	fit.vcov 	<- solve(mlfit$hessian) 	# not negative because optim minimizing (-1)*log-like
	fit.se 		<- sqrt(diag(fit.vcov))
	
	fit.z		<- mlfit$par/fit.se 
	fit.pv		<- 2*pnorm(-abs(fit.z))
	
	n.0 <- with(roc.obj$data.all, length(score[status==0]))
	n.1 <- with(roc.obj$data.all, length(score[status==1]))
	
#### Add coefficient names

	names.a  <- dimnames(mmats$mm.a)[[2]]
	names.b  <- dimnames(mmats$mm.b)[[2]]
	names.t0 <- dimnames(mmats$mm.t0)[[2]]
	names.t1 <- dimnames(mmats$mm.t1)[[2]]
	names.all <- list(names.a=names.a, names.b=names.b, 
						names.t0=names.t0, names.t1=names.t1)

	names(mlfit$par) <- c(names.a, names.b, names.t0, names.t1)
	names(fit.se) <- c(names.a, names.b, names.t0, names.t1)
		
	names(fit.z) <- c(names.a, names.b, names.t0, names.t1)
	names(fit.pv) <- c(names.a, names.b, names.t0, names.t1)
	
	colnames(fit.vcov) <- c(names.a, names.b, names.t0, names.t1)
	rownames(fit.vcov) <- c(names.a, names.b, names.t0, names.t1)

#### Coefficient tables

	a.tab  <- data.frame(mlfit$par[where.par=="a"], fit.se[where.par=="a"],
							fit.z[where.par=="a"], fit.pv[where.par=="a"])
	colnames(a.tab) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
	
	b.tab  <- data.frame(mlfit$par[where.par=="b"], fit.se[where.par=="b"],
							fit.z[where.par=="b"], fit.pv[where.par=="b"])
	colnames(b.tab) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")

	t0.tab <- data.frame(mlfit$par[where.par=="t0"], fit.se[where.par=="t0"],
							fit.z[where.par=="t0"], fit.pv[where.par=="t0"])
	colnames(t0.tab) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
						
	t1.tab <- data.frame(mlfit$par[where.par=="t1"], fit.se[where.par=="t1"],
							fit.z[where.par=="t1"], fit.pv[where.par=="t1"])
	colnames(t1.tab) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
						
#### Get residuals

	resid.list <- fitted.roc(parms=mlfit$par, roc.obj=roc.obj, mm.list=mmats,...) 

	r.squared 		<- cor(roc.obj$data.all$score, resid.list$yhat)^2
	adj.r.squared 	<- 1-(1-r.squared)*(n.0+n.1-4)/(
						n.0+n.1-(mmats$p.a+mmats$p.b+mmats$p.t0+mmats$p.t1) )		#slightly modified
	
	resid.df		<- n.0+n.1 - (mmats$p.a+mmats$p.b+mmats$p.t0+mmats$p.t1)
	parms.df 		<- mmats$p.a+mmats$p.b+mmats$p.t0+mmats$p.t1
	
#### Organize results by formula

	a.list 	<- list(formula=flist$f.a,  coef=mlfit$par[where.par=="a"],  se=fit.se[where.par=="a"], tab=a.tab )
	b.list 	<- list(formula=flist$f.b,  coef=mlfit$par[where.par=="b"],  se=fit.se[where.par=="b"], tab=b.tab )
	t0.list <- list(formula=flist$f.t0, coef=mlfit$par[where.par=="t0"], se=fit.se[where.par=="t0"], tab=t0.tab )
	t1.list <- list(formula=flist$f.t1, coef=mlfit$par[where.par=="t1"], se=fit.se[where.par=="t1"], tab=t1.tab )

	model.list <- list(model=model, dist.0=dist.0, dist.1=dist.1, dpar.0=dpar.0, dpar.1=dpar.1)
#### Results list 

	outlist <- list(call=cl, loglike=-mlfit$value, data.all=roc.obj$data.all,
					a=a.list, b=b.list, t0=t0.list, t1=t1.list,
					coef=mlfit$par, se=fit.se, vcov=fit.vcov, 
					fitted.values=resid.list$yhat, residuals=resid.list$resid,
					fitvar=resid.list$vhat, mse=resid.list$mse, rse=resid.list$rse,
					n=c(n.0=n.0, n.1=n.1), r.squared=r.squared, adj.r.squared=adj.r.squared,
					resid.df=resid.df, parms.df=parms.df, model=model.list)

class(outlist) <- c("fitobj", model, class(outlist))
outlist
}

###
##
#