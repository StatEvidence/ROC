################################################
## Function:	set.roc
## Version:		2.0 	(old version = ROC.R)
## Purpose:		Compute empirical ROC points &
##				CI region based on exact CIs
##				for TPF and FPF
## 
## Author: 		JD Blume
## Date: 		April 2020
################################################
#'
#' Create an ROC object
#'
#' @description This function will create a ROC object that can be used for plotting empirical ROC point or estimating a smooth ROC curve or computing the area under the ROC curve (AUC).
#'
#' @param score	Outcome vector.
#' @param grp Grouping vector (must be 2 exactly levels).
#' @param x	Covariate matrix / dataframe (one column per covariate).
#' @param as.groups	Indicator for alternate data format (\code{TRUE} or \code{FALSE}). See Details.
#' @param level	Confidence interval level (default is \code{95\%}).
#' @param flip.grps Indicator if group with largest mean should be labeled as group 1. See details.
#' @param get.ci Flag to compute confidence intervals for the empirical ROC points.
#' @param center Flag to center all non-numeric and non-binary (0-1) variables. Default is \code{TRUE}.
#'
#' @details The group with the largest sample mean will be re-labeled as group '1'; group '0' has the smallest sample mean. This orientation often ensures ROC curve will lie above y=x line and the area under the ROC curve will be greater than 0.5. To reverse these default group assignments, set \code{flip.grps=FALSE}.
#'
#' To accomodate data structured in 2-group format, set \code{as.groups=TRUE}, and define \code{score} to be the outcome from one group and \code{grp} to be the outcome the other. This option is intended for quick and easy usage when there are no covariates data present.
#'
#' Centering of continuous variables is on by default. This facilitates the interpretation of the default AUC computations, which depend on the (a, b) intercepts.
#'
#' Set \code{get.ci=FALSE} to improve speed when confidence intervals are not needed.
#' 
#' @return An ROC object (list) is returned with the following elements:
#' 
#' \describe{
#' \item{\code{fp}}{Vector of the false positive points (1-specificity). Defined as P( score>=cutpoint | grp=0 ).}
#' \item{\code{tp}}{Vector of the true positive points (sensitivity). Defined as P( score>=cutpoint | grp=1 ).}
#' \item{\code{cutpoints}}{Vector of cutpoints corresponding to \code{(fp, tp)} pairs.}
#' \item{\code{fp.ci}}{Exact CI limits for fp estimates (Clopper and Pearson).}
#' \item{\code{tp.ci}}{Exact CI limits for tp estimates (Clopper and Pearson).}
#' \item{\code{data.0}}{Outcome vector for group 0.}
#' \item{\code{data.1}}{Outcome vector for group 1.}
#' \item{\code{data.all}}{Dataframe with \code{score} for outcome and \code{status} for group.}
#' }
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
#' ## Fake ROC data for examples
#' y.0 <- round( rnorm(100, mean=0, sd=1)   , 1)
#' y.1 <- round( rnorm(100, mean=1, sd=1.1) , 1)
#' age <- rnorm(200, mean=35, sd=8)
#' size <- sample(c("small", "medium", "big"), 200, replace=TRUE)
#'
#' ## Create an ROC object using 2-group structure
#' myroc <- set.roc(y.0, y.1, as.groups=TRUE)
#' 
#' ## Table of operating points
#' with(myroc, cbind(cutpoint, fp, fp.ci, tp, tp.ci))
#'
#' ## Fraction of group 1 scores at or above the overall median score
#' cc <- median(myroc$data.all$score)
#' with(myroc, cbind(cutpoint, tp, tp.ci))[myroc$cutpoint>cc,][1,]
#'
#' ## Create an ROC object using standard data structure
#' y <- c(y.0, y.1)
#' d <- c(rep(0,100),rep(1,100))
#' x <- data.frame(age, size)
#'
#' myroc <- set.roc(y, d, x)
#'

set.roc <- function(score, grp, x=NULL, as.groups=FALSE, level=0.95, flip.grps=TRUE, get.ci=TRUE, center=TRUE, ...){
	
#### Rearrange data from 2-group input

	if (as.groups==TRUE) {
		n 		<- 	c(length(score),length(grp))
		score 	<-	c(score,grp)
		grp		<-  c(rep("a",n[1]),rep("b",n[1]))
		}

#### Errors / Warnings

	if (length(levels(as.factor(grp)))!=2) {
		stop("Grouping variable is not dichotomous.")
		}

#### Set Data Frame

	status <- 1*(grp==levels(as.factor(grp))[2])

	all <- as.data.frame(cbind("score" = score, "status" = status))
	
	if (is.null(x)==FALSE) { 
						
		if (NROW(x)!=NROW(all)) {stop("Covariate and outcome record lengths do not match.")}
		
		all <- data.frame(cbind("score" = score, "status" = status), x) 
		}
	
	all <- all[complete.cases(all),]

	if (center==TRUE) {
		numeric <- sapply(all, is.numeric)
		binary 	<- sapply(all, function(x) { all(x %in% 0:1) })
		index <- numeric&!binary
		index[c('score','status')] <- FALSE
		all[index] <- lapply(all[index], scale, scale=FALSE, center=TRUE)
		}

#### Check label '1' is assigned to proper group (deafult grp w/ largest sample mean)

	grp.order <- (diff(aggregate(all$score, list(all$status), mean)$x)>=0)
	if (grp.order==FALSE & flip.grps==TRUE) { all$status <- 1-all$status }

#### Define group specific datasets
	
	data.0 <- with(all, score[status==0]) 
	data.1 <- with(all, score[status==1]) 

#### Get empirical empirical ROC points

	s.0 <- 1-ecdf(data.0)(sort(c(data.0,data.1)))
	s.1 <- 1-ecdf(data.1)(sort(c(data.0,data.1)))

	pts = unique(data.frame(rbind(c(1,1), cbind(fp=s.0,tp=s.1), c(0,0))))

#### Compute exact CIs

	x.ci <- y.ci <- matrix(NA, nrow=2, ncol=length(pts$fp))
	
	if (get.ci==TRUE) {
		
		## Get sample size for exact CIs

		FP.ct 	<- pts$fp*length(data.0)
		TP.ct 	<- pts$tp*length(data.1)

		##	Clopper-Pearson CIs for each ROC point (in x and y directions)

		x.ci <- sapply(FP.ct,
				function(x) binom.test(x, length(data.0), conf.level=level)$conf.int)
		dimnames(x.ci)[[1]] <- c("fp.lo", "fp.hi")
		
		y.ci <- sapply(TP.ct,
				function(x) binom.test(x, length(data.1), conf.level=level)$conf.int)
		dimnames(y.ci)[[1]] <- c("tp.lo", "tp.hi")
		}

#### Get unique cutpoint vector
	
	cutpt <- unique(c(sort(c(data.0, data.1))))

#### Return FPF, TPF, Cutpoint ; returns P(score > cutpoint | status)

	outlist <- list("fp" = pts$fp, "tp" = pts$tp, 
						"cutpoint" = c(cutpt, NA), 
						"fp.ci" = t(x.ci), "tp.ci" = t(y.ci),
						"data.0" = data.0, "data.1" = data.1,
						"data.all" = all)

class(outlist) <- "roc"
return(outlist)
}

###
##
#