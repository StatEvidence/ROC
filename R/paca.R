#' Pancreatic Cancer Biomarker Data
#'
#' Cancer serum biomarker data from 141 individuals. 
#' 51 control patients and 90 cases with pancreatic cancer 
#' studied at the Mayo clinic. Carbohydrate antigen (CA19-9) and
#' cancer antigen (CA125) are in columns 1 and 2, respectively.
#' CA19-9 has higher accuracy.
#'
#' @docType data
#'
#' @usage data(paca)
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @keywords datasets
#'
#' @references Wieand et al. (1989) Biometrika 76(3):585-592
#' #(\href{https://doi.org/10.1093/biomet/76.3.585}{Biometrika})
#'
#' @source #\href{https://research.fhcrc.org/diagnostic-biomarkers-center/en/datasets.html}{Diagnostic and Biomarkers Statistical Center}
#'
#' @examples
#' data(paca)
#' auc(set.roc(log(paca$ca125), paca$case))
#' cutplot(set.roc(log(paca$ca199), paca$case))
"paca"


