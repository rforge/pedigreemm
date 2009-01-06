lmer_ZStar <-
    function(formula, data, family = NULL, pre = NULL, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE, subset,
             weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    mc <- match.call()
    mc[[1]] <- as.name("lmer")
    if (is.null(pre)) {
        mc$pre <- NULL
        return(eval(mc))
    }
    mc$pre <- NULL
    mc$doFit <- FALSE
    lf <- eval.parent(mc)
    FL <- lf$FL

    if (length(pre) != length(FL$fl))
        stop(paste("length of argument `pre' must be", length(FL$fl)))
    ## FIXME: Should check for names.  Also check for dimensions matching if not NULL
    for (i in seq_along(FL$fl))
        if (!is.null(pre[[i]]))
            FL$trms[[i]]$Zt <- FL$trms[[i]]$A <- pre[[i]] %*% FL$trms[[i]]$Zt
    lf$FL <- FL
    ans <- do.call(if (!is.null(lf$glmFit)) lme4:::glmer_finalize else lme4:::lmer_finalize, lf)
    ans@call <- match.call()
    ans
}

