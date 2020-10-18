print.auc <- function(object) {
	cat("\n")
	print(object$tab, row.names=FALSE, right=F)
	cat("----\n") 
	if (!is.null(object$level)) {
		cat("CI level: ", paste(round(object$level*100, 4)), "%\n",sep="")} else {
		cat("CI level: NA\n",sep="")}
	cat("\n")
    ## NextMethod(obj)	## To return complete object
}


