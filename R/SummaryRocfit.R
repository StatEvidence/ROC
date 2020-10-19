#' Summary Fitobj
#' 
#' @export
#'

summary.fitobj <- function(object, digits=4) {
	cat("\nCall:\n")
	print(object$call)
	cat("\nObservations:\n")
	cat(paste(object$n[1], "class 0,", object$n[2] ,  "class 1" ))
	cat("\n")
    cat("\nMultiple R-squared: ", formatC(object$r.squared, digits=digits),
          ",\tAdjusted R-squared: ",formatC(object$adj.r.squared, digits=digits),
           "\n", sep="")
	cat("\nCoefficients:")
	cat("\n-->", deparse(object$a$formula),"\n")
	printCoefmat(object$a$tab, digits = digits, P.values = TRUE, signif.stars = FALSE, signif.legend = FALSE)
	cat("\n-->", deparse(object$b$formula),"\n")
	printCoefmat(object$b$tab, digits = digits, P.values = TRUE, signif.stars = FALSE, signif.legend = FALSE)
	cat("\n-->", deparse(object$t0$formula),"\n")
	printCoefmat(object$t0$tab, digits = digits, P.values = TRUE, signif.stars = FALSE, signif.legend = FALSE)
	cat("\n-->", deparse(object$t1$formula),"\n")
	printCoefmat(object$t1$tab, digits = digits, P.values = TRUE, signif.stars = FALSE, signif.legend = FALSE)
	cat("\nResidual standard error: ", formatC(object$rse, digits=digits), " on ",
       	 formatC(object$resid.df), " degrees of freedom\n", sep="")
 	cat("AIC: ", formatC(-2*object$loglike+2*object$parms.df, digits=digits), ", BIC: ",
        	 formatC(-2*object$loglike+object$parms.df*log(object$n[1]+object$n[2]), digits=digits), " \n", sep="")
	cat("\n")
#	cat("\nResiduals:\n")
#	print(data.frame(Min=min(ff$residuals), "1Q"=quantile(ff$residuals, 0.25), 
#						" Median"=median(ff$residuals) , '3Q'=quantile(ff$residuals, 0.75),
#						Max=min(ff$residuals), check.names=FALSE), digits=digits, row.names=FALSE)
#	cat("----------------------------------------------")
}


    



