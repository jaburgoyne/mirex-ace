ReadACEResults <- function (results.directory, analysis) {
    #' MIREX ACE results for a particular chord vocabulary or segmentation
    #' 
    #' @param results.directory directory containing the results in Johan
    #'                          Pauwels's output format
    #' @param analysis chord vocabulary or segmentation to load
    #'
    is.segmentation <- analysis == "Segmentation"
    path <- ""
    col.names <- c()
    col.classes <- c()
    if (is.segmentation) {
        path <- file.path(results.directory, "resultsSegmentation")
        col.names <- c("song", "", "under.segmentation", "over.segmentation")
        col.classes <- c("factor", "NULL", "numeric", "numeric")}
    else {
        path <- file.path(results.directory, paste0("resultsMirex", analysis))
        col.names <- c(
            "song", 
            "performance", 
            "duration",
            rep("", 6))
        col.classes <- c(
            "factor", 
            "numeric", 
            "numeric", 
            rep("NULL", 6))}
    results <- data.frame()
    for (algo.csv in list.files(path, ".*\\.csv", full.names = TRUE)) {
        algo <- substr(basename(algo.csv), 1, nchar(basename(algo.csv)) - 4)
        algo.results <- read.table(
            algo.csv,
            skip = 2, # Johan's header has two lines
            header = FALSE,
            sep = ",",
            dec = ".",
            col.names = col.names, 
            colClasses = col.classes)
        if (is.segmentation) {
            ## Replace over- and under-segmentation with the 
            ## harmonic mean of their arithmetic inverses. 
            algo.results$performance <- {
                1 / (0.5 
                     * (1/(1-algo.results$under.segmentation) 
                        + 1/(1-algo.results$over.segmentation)))}
            algo.results$over.segmentation <- c()}
        else {
            ## Rescale RCO results back to [0, 1].
            algo.results$performance <- algo.results$performance/100}
        results <- rbind(
            results, 
            cbind(algorithm = rep(algo, dim(algo.results)[1]), algo.results))}
    ## Johan's segmentation files have no column for duration, so add it.
    if (is.segmentation)
        results$duration <- ReadACEResults(results.directory, "Root")$duration
    ## Order algorithm factor by WCSR.
    results$algorithm <- with(results, { 
        factor(
            algorithm,
            levels = levels(algorithm)[
                order(
                    dummy.coef(
                        glm(performance ~ algorithm, 
                            weights = duration))
                    $algorithm,
                    decreasing = TRUE)])})
    ## Add rank information to support old-style analysis.
    results <- within(results, {
        performance.rank <- rep(NA, length(performance))
        for (group in levels(song))
            performance.rank[song==group] <- rank(performance[song==group])
        rm(group)})
    results}

EvaluateACE <- function (
    results.directory,
    analysis = c(
        "SeventhsBass",
        "Root", 
        "MajMin", "MajMinBass", 
        "Sevenths",
        "Segmentation"),
    old.style = FALSE,
    adjust = "fdr") {
    #' MIREX ACE evaluation for a particular chord vocabulary or segmentation
    #' 
    #' @param results.directory directory containing the results in Johan
    #'                          Pauwels's output format
    #' @param analysis chord vocabulary or segmentation to load
    #' @param old.style whether to use old-style (Friedman) evaluation
    #' @param adjust correction method for multiple comparisons (see 'p.adjust')
    #'
    analysis <- match.arg(analysis)
    results <- ReadACEResults(results.directory, analysis)
    if (old.style) {
        ace.eval <- summary(
            multcomp::glht(
                glm(formula = performance.rank ~ song + algorithm,
                    data    = results,
                    weights = duration),
                multcomp::mcp(algorithm = "Tukey")),
            test = multcomp::adjusted("fdr"))
        ## Use the standard post-hoc z-statistic for Friedman tests: see among
        ## other Salvador García, Alberto Fernández, Julián Luengo & Francisco
        ## Herrera (2010), “Advanced Nonparametric Tests for Multiple 
        ## Comparisons in the Design of Experiments in Computational 
        ## Intelligence and Data Mining: Experimental Analysis of Power”, 
        ## Information Sciences 180: 2044–64.
        ace.eval$test$sigma[1:length(ace.eval$test$sigma)] <- with(results, {
            sqrt(
                nlevels(results$algorithm) 
                * (nlevels(results$algorithm) - 1)
                / (6 * nlevels(results$song)))})
        ace.eval$test$tstat <- ace.eval$test$coefficients / ace.eval$test$sigma
        ace.eval$test$pvalues <- p.adjust(
            p      = 2 * pnorm(abs(ace.eval$test$tstat), lower.tail = FALSE),
            method = adjust)}
    else {
        ace.eval <- summary(
            multcomp::glht(
                geepack::geeglm(
                    formula   = performance ~ algorithm,
                    data      = results,
                    family    = binomial(link = "logit"),
                    corstr    = "exchangeable",
                    weights   = duration,
                    id        = song),
                multcomp::mcp(algorithm = "Tukey")),
            test = multcomp::adjusted(adjust))}
    ace.eval$old.style <- old.style
    ace.eval$adjust <- adjust
    class(ace.eval) <- c("mirexace", class(ace.eval))
    ace.eval}

summary.mirexace <- print # Summary would undo the choice of adjustment method.

plot.mirexace <- function(x, level = .005, ...) {
    #' Plot a compact letter display  with weighted response boxplots for
    #' MIREX ACE results.
    #'
    #' @param level significance level for compact letter display
    #'
    #' Adapted from plot.cld() in version 1.3 of the multcomp package (T.
    #' Hothorn, F. Bretz & P. Westfall, 2014).
    #'
    xcld <- cld(x, level = level)
    mcletters <- xcld$mcletters
    msletters <- mcletters$monospacedLetters
    vletters <- sapply(
        msletters,
        function(x) paste(strsplit(x, "")[[1]], "\n", collapse = ""))
    vletters <- vletters[gsub(" ", "", levels(xcld$x))]
    msletters <- msletters[gsub(" ", "", levels(xcld$x))]
    dat <- xcld[c("x", "y", "lp")]
    if (is.null(xcld$weights))
        dat$weights <- rep(1, NROW(x$y))
    else
        dat$weights <- xcld$weights
    dat <- as.data.frame(dat)
    
    if (!is.numeric(dat$y)) {
        if (is.integer(dat$y)) dat$y <- as.numeric(dat$y)
        else stop("Response must be numeric.")}
    if (!is.factor(dat$x)) stop("Covariate must be a factor.")
    
    bx <- list(
        stats = matrix(NA, 5, 0),
        n     = c(),
        conf  = matrix(NA, 2, 0),
        out   = c(),
        group = c(),
        names = c())
    for (i in 1:nlevels(dat$x)) {
        g <- levels(dat$x)[i]
        bx$names <- c(bx$names, g)
        ## Compute basic weighted quantiles.
        yg <- dat$y[!is.na(dat$y) & dat$x == g]
        wq <- Hmisc::wtd.quantile(
            x       = yg,
            weights = dat$weights[dat$x == g],
            normwt  = TRUE)
        iqr <- wq[4] - wq[2]
        ## Use the default outlier computation from boxplot.stats().
        wq[1] <- min(yg[yg >= wq[2] - 1.5 * iqr])
        wq[5] <- max(yg[yg <= wq[4] + 1.5 * iqr])
        bx$stats <- cbind(bx$stats, wq)
        out <- yg[yg < wq[1] | yg > wq[5]]
        bx$out <- c(bx$out, out)
        bx$group <- c(bx$group, rep(i, length(out)))
        n <- length(yg)
        bx$n <- c(bx$n, n)
        ## Use the default notch computation from boxplot.stats().
        s <- 1.58 * iqr / sqrt(n)
        bx$conf <- cbind(bx$conf, wq[3] + c(-s, s))}
    colnames(bx$stats) <- levels(dat$x)
    oldmar <- par("mar")
    par(mar = oldmar + c(0, 0, dim(mcletters$LetterMatrix)[2], 0)) # CLD space
    bxp(bx, 
        xlab = "Algorithm", 
        ylab = ifelse(
            x$old.style, 
            paste("Rank per Song (1 low; ", nlevels(dat$x), "high)"),
            "Performance"),
        las = 2,
        ...)
    axis(3, at = 1:nlevels(dat$x), labels = vletters)
    # Add points for means.
    points(
        1:nlevels(dat$x),
        {if (x$old.style) 
             dummy.coef(x$model)$"(Intercept)" + dummy.coef(x$model)$algorithm
         else
             plogis(coef(x$model)[1]
                    + contrasts(x$model$data$algorithm) %*% coef(x$model)[-1])},
        pch = 16)
    par(mar = oldmar)}
