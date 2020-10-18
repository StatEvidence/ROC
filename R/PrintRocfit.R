print.fitobj <- function(object) {
    cat("\nCall:\n")
	print(object$call)
	cat("\nCoefficients:\n")
	cat(deparse(object$a$formula),"\n")
	print(object$a$coef, digits=4)
	cat("\n")
	cat(deparse(object$b$formula),"\n")
	print(object$b$coef, digits=4)
	cat("\n")
	cat(deparse(object$t0$formula),"\n")
	print(object$t0$coef, digits=4)
	cat("\n")
	cat(deparse(object$t1$formula),"\n")
	print(object$t1$coef, digits=4)
	cat("\n")
    ## NextMethod(obj)	## To return complete object
}


