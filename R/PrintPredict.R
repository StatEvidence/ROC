#' Print predict
#' 
#' @export
#'
print.predict <- function(object) {
    
	test=rbind(
		c(mean(object$auc$est), mean(object$auc$lo), mean(object$auc$hi)), 
		c(sd(object$auc$est), sd(object$auc$lo), sd(object$auc$hi)), 
		c(median(object$auc$est), median(object$auc$lo), median(object$auc$hi)), 
		c(IQR(object$auc$est), IQR(object$auc$lo), IQR(object$auc$hi)), 
		c(min(object$auc$est), min(object$auc$lo), min(object$auc$hi)), 
		c(max(object$auc$est), max(object$auc$lo), max(object$auc$hi))  
		)

	colnames(test) <- c("AUC", "Lower CI", "Upper CI")
	rownames(test) <- c("Mean", "SD", "Median", "IQR", "Min", "Max")
	
	cat("\nCall:\n")
	print(object$call)
	cat("\nModel:\n")
	cat(paste(object$model[1], "\n"))
	cat("\nPrediction summary: \n")
	cat(paste(object$records, "new records,", object$features ,  "covariates\n" ))
	cat("\n")
	print(round(test,4))
	cat(paste("\nConfidence Level:", object$level ,  "\n" ))
	cat("\nInvisible: all, a, b, vcov, t0, t1, \n") 
	cat("           auc, auc.lo, auc.hi\n")
	cat("\n")
	## NextMethod(obj)	## To return complete object
}

