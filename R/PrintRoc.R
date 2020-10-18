print.roc <- function(object) {
	test=rbind(
	c(mean(object$data.0), mean(object$data.1)),
	c(sd(object$data.0), sd(object$data.1)),
	c(median(object$data.0), median(object$data.1)),
	c(IQR(object$data.0), IQR(object$data.1)),
	c(min(object$data.0), min(object$data.1)),
	c(max(object$data.0), max(object$data.1)) 
	)
	
	colnames(test) <- c("class 0", "class 1")
	rownames(test) <- c("Mean", "SD", "Median", "IQR", "Min", "Max")
	
	k=dim(object$data.all)[[2]]
	
	cat("\nObservations:\n")
	cat(paste(length(object$data.0), "class 0,", length(object$data.1) ,  "class 1\n" ))
	cat("\nCovariates:\n")
	cat(paste(shQuote(names(object$data.all)[3:k], type="cmd"), collapse=", "), "\n")
	cat("\nSummary:\n")
	print(round(test,4))
	cat("\nInvisible: fp, fp.ci, tp, tp.ci, \n") 
	cat("   cutpoint, data.0, data.1, data.all\n")
	cat("\n")
    ## NextMethod(obj)	## To return complete object
}

