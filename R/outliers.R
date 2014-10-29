ACEOutliers <- function (mirexace) {
    #' Identify outliers in MIREX ACE performance according to Chauvenet's 
    #' criterion.
    chauvenet <- qnorm(1 / (4 * dim(mirexace$model$data)[1]), lower.tail=FALSE)
    resid.ace <- resid(mirexace$model)
    outliers <- mirexace$model$data[abs(resid.ace) > chauvenet * sd(resid.ace), ]
    outliers}


