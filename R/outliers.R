#' Outliers in ACE Evaluations
#' 
#' \code{ACEOutliers} returns outliers in MIREX ACE performance
#' according to Chauvenet's criterion.
#'
#' @param mirexace a \code{mirexace} object
#' @export
ACEOutliers <- function (mirexace) {
    chauvenet <- qnorm(1 / (4 * dim(mirexace$model$data)[1]), lower.tail=FALSE)
    resid.ace <- resid(mirexace$model)
    outliers <- mirexace$model$data[abs(resid.ace) > chauvenet * sd(resid.ace), ]
    outliers}


