.PCorr <- function(m, n, adjust = NULL) {
    #' Replaces upper triangle of a correlation matrix with p-values, corrected for
    #' multiple comparisons.
    #'
    #' @param m the correlation matrix
    #' @param n the number of pairs for each entry in the matrix
    #' @param adjust the adjustment method for multiple comparisons with p.adjust
    #'
    upper <- row(m) < col(m)
    m[upper] <- m[upper] * sqrt((n - 2) / (1 - m[upper] ^ 2))
    m[upper] <- pt(abs(m[upper]), n - 2, lower.tail = FALSE) * 2
    m[upper] <- p.adjust(m[upper], adjust)
    m}

.WCorr <- function (a, b = a, w = rep(1, nrow(a))/nrow(a)) {
    #' Pearson's correlation matrix with weighted pairs. Adapted from
    #' http://stackoverflow.com/questions/9460664/weighted-pearsons-correlation.
    #'
    #' See also: corr() in the boot package
    #'
    ## Convert to matrices.
    if (is.data.frame(a)) a <- as.matrix(b)
    if (is.data.frame(b)) b <- as.matrix(b)
    ## Normalize weights.
    w <- w / sum(w)
    ## Centre matrices.
    a <- sweep(a, 2, colSums(a * w))
    b <- sweep(b, 2, colSums(b * w))
    ## Compute weighted correlation.
    t(w*a) %*% b / sqrt(colSums(w * a**2) %*% t(colSums(w * b**2)))}

ACECor <- function(mirexace) {
    #' Correlation matrix of algorithm performance for MIREX ACE.
    #'
    #' Correlations appear in the lower triangle and p-values appear in the 
    #' upper triangle, corrected for multiple comparisons.
    cor.ace <- with(mirexace$model$data, {
        if (mirexace$old.style) {
            .ranks. <- reshape(
                data      = data.frame(
                    song, 
                    algorithm, 
                    performance.rank, 
                    duration),
                idvar     = c("song"),
                timevar   = c("algorithm"),
                direction = "wide")
            rownames(.ranks.) <- .ranks.[, 1]
            .ranks. <- .ranks.[, c(2 * order(order(levels(algorithm))), 3)]
            colnames(.ranks.) <- c(levels(algorithm), "duration")
            .WCorr(
                .ranks.[, 1:nlevels(algorithm)], 
                w = .ranks.[, nlevels(algorithm) + 1])}
        else {
            .contrasts. <- contrasts(algorithm)
            .vcov. <- vcov(mirexace$model)
            cov2cor(.contrasts. %*% .vcov.[-1, -1] %*% t(.contrasts.))}})
    .PCorr(cor.ace, nlevels(mirexace$model$data$song), mirexace$adjust)}

PlotACEClust <- function(mirexace, ...) {
    #' Plots a hierarchical clustering of MIREX ACE algorithms based on their
    #' correlations in performance.
    #'
    plot(hclust(as.dist((1 - ACECor(mirexace))/2)),
         main="", sub="",
         ylab="Pearson's Distance", xlab="",
         ...)}
    
