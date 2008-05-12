lmer_ZStar <-
    function(formula, data, family = NULL, pre = NULL, method = c("REML", "ML"),
             control = list(), start = NULL, verbose = FALSE, subset,
             weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    mc <- match.call()
    if (is.null(pre)) {
        mc$pre <- NULL
        mv[[1]] <- as.name("lmer")
        return(eval(mc))
    }
    if (!is.null(family)) {             # call glmer_ZStar
        mc[[1]] <- as.name("glmer_ZStar")
        return(eval(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)

    fr <- lme4:::lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    FL <- lme4:::lmerFactorList(formula, fr$mf, 0L, 0L) # flist, Zt, cnames
    Y <- as.double(fr$Y)


    if (length(pre) != length(FL$fl))
        stop(paste("length of argument `pre' must be", length(FL)))
    for (i in seq_along(FL$fl))
        if (!is.null(pre[[i]]))
            FL$fl[[i]]$Zt <- FL$fl[[i]]$A <- pre[[i]] %*% FL$fl[[i]]$Zt
    lme4:::lmer_finalize(mc, fr, FL, start, match.arg(method), verbose)
}
