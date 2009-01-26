#### "pedigree" class methods

#' Constructor for pedigree objects
#'
#' A simple constructor for a pedigree object.  The main point for the
#' constructor is to use coercions to make the calls easier.
#' @param sire integer vector or factor representation of the sires
#' @param dam integer vector or factor representation of the dams
#' @param label character vector of labels
#! @return an pedigree object of class \linkS4class{pedigree}

#' @note \code{sire}, \code{dam} and \code{label} must all have the
#' same length and all labels in \code{sire} and \code{dam} must occur
#' in \code{label}
#' @export
pedigree <- function(sire, dam, label) {
    n <- length(sire)
    stopifnot(n == length(dam), n == length(label))
    sire <- as.integer(sire); dam <- as.integer(dam)
    sire[sire < 1 | sire > n] <- NA
    dam[dam < 1 | dam > n] <- NA
    new("pedigree", sire = sire, dam = dam,
	label = as.character(label))
}

#' Coerce a pedigree to a sparse triangular matrix
#'
#' Create a sparse, unit lower triangular matrix from a pedigree.  The
#' matrix the L factor in the LDL' Cholesky factorization of the
#' inverse of the relationship matrix.
#'
#' @export
setAs("pedigree", "sparseMatrix", # representation as T^{-1}
      function(from) {
	  sire <- from@sire
	  n <- length(sire)
	  animal <- seq_along(sire)
	  j <- c(sire, from@dam)
	  ind <- !is.na(j)
	  as(new("dtTMatrix", i = rep.int(animal, 2)[ind] - 1L,
		 j = j[ind] - 1L, x = rep.int(-0.5, sum(ind)),
		 Dim = c(n,n), Dimnames = list(from@label, NULL),
		 uplo = "L", diag = "U"), "dtCMatrix")
      })

## these data frames are now storage efficient but print less nicely
setAs("pedigree", "data.frame",
      function(from)
      data.frame(sire = from@sire, dam = from@dam,
		 row.names = from@label))

#' Convert a pedigree to a data frame
#'
#' Express a pedigree as a data frame with \code{sire} and
#' \code{dam} stored as factors.  If the pedigree is an object of
#' class \linkS4class{pedinbred} then the inbreeding coefficients are
#' appended as the variable \code{F}
#'
#' @param x a pedigree object of class \linkS4class{pedigree}
#' @return a data frame
#'
ped2DF <- function(x) {
    stopifnot(is(x, "pedigree"))
    lab <- x@label
    lev <- seq_along(lab)
    ans <- data.frame(sire = factor(x@sire, levels = lev, labels = lab),
                      dam  = factor(x@dam,  levels = lev, labels = lab),
                      row.names = lab)
    if (is(x, "pedinbred")) ans <- cbind(ans, F = x@F)
    ans
}

setMethod("show", signature(object = "pedigree"),
	  function(object) print(ped2DF(object)))

setMethod("head", "pedigree", function(x, ...)
	  do.call("head", list(x = ped2DF(x), ...)))

setMethod("tail", "pedigree", function(x, ...)
	  do.call("tail", list(x = ped2DF(x), ...)))


Linvt <- function(ped)
{
    stopifnot(is(ped, "pedigree"))
    t(as(Diagonal(x = sqrt(1/(1 +.Call("pedigree_inbreeding", ped)))) %*%
       as(ped, "sparseMatrix"), "triangularMatrix"))
}

#' Inbreeding coefficients from a pedigree
#'
#' Determine the diagonal factor in the inverse of the relationship
#' matrix from a pedigree.  Although in practice we use the
#' relationship matrix from a pedigree, it is easier to evaluate the
#' factors of the inverse of the relationship matrix.  This function
#' returns the diagonal matrix from the LDL' form of the Cholesky
#' factorization of the inverse of the relationship matrix.
#'
#' @param ped an object that inherits from class \linkS4class{pedigree}
#' @return a numeric vector
#' @export
inbreeding <- function(ped) {
    stopifnot(is(ped, "pedigree"))
    F <- .Call("pedigree_inbreeding", ped, PACKAGE = "pedigreemm")
    sire <- ped@sire
    dam <- ped@dam
    Fsire <- ifelse(is.na(sire), -1, F[sire])
    Fdam <-  ifelse(is.na(dam), -1, F[dam])
    ans <- 1 - 0.25 * (2 + Fsire + Fdam)
    names(ans) <- ped@label
    ans
}

#' Relationship factor from a pedigree
#'
#' Determine the right Cholesky factor of the relationship matrix for
#' the pedigree \code{ped}, possibly restricted to the specific labels
#' that occur in \code{labs}. 
#' 
#' @param ped a pedigree that includes the individuals who occur in svec
#' @param labs a character vector or a factor giving the labels to
#'    which to restrict the relationship matrix. If \code{labs} is a
#'    factor then the levels of the factor are used as the labels.
#'    Default is the complete set of labels in the pedigree.
#' @return a sparse, symmetric relationship matrix
#' @export

relfactor <- function(ped, labs = ped@label)
{
    stopifnot(is(ped, "pedigree"))
    labs <- factor(labs) # drop unused levels from a factor
    stopifnot(all(labs %in% ped@label))
    Diagonal(x = sqrt(inbreeding(ped))) %*%
        solve(t(as(ped, "sparseMatrix")),
              Matrix:::fac2sparse(factor(labs, levels = ped@label),
                                  drop = FALSE))
}

pedigreemm <-
    function(formula, data, family = NULL, REML = TRUE, pedigree = list(),
             control = list(), start = NULL, verbose = FALSE, 
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    mc <- match.call()
    lmerc <- mc
    lmerc[[1]] <- as.name("lmer")
    lmerc$pedigree <- NULL
                                        # check the pedigree argument
    if (!length(pedigree))              # call lmer instead
        return(eval.parent(lmerc))
    stopifnot(is.list(pedigree),
              length(names(pedigree)) == length(pedigree),
              all(sapply(pedigree, is, class2 = "pedigree")))
                                        # call lmer without pedigree and with doFit = FALSE
    lmerc$doFit <- FALSE
    lmf <- eval(lmerc, parent.frame())

    pnms <- names(pedigree)
    stopifnot(all(pnms %in% names(lmf$FL$fl)))
    asgn <- attr(lmf$FL$fl, "assign")
### FIXME: check here that the random effects term is of the form (1|factor)
    for (i in seq_along(pedigree)) {
        tn <- which(match(pnms[i], names(lmf$FL$fl)) == asgn)
        if (length(tn) > 1)
            stop("a pedigree factor must be associated with only one r.e. term")
        pedi <- pedigree[[i]]
        lit <- Linvt(pedi)
### FIXME: probably should use lmf$FL$trms[[tn]]$Zt instead.  Need to
### check that the number of rows is the number of unique levels of
### the factor.        
        Zt <- lmf$FL$trms[[tn]]$Zt
        rn <- rownames(Zt)
        stopifnot(all(rn %in% pedi@label))
        B <- solve(lit, Matrix:::fac2sparse(factor(rn,
                                                   levels = pedi@label),
                                            drop = FALSE))
        lmf$FL$trms[[tn]]$Zt <- lmf$FL$trms[[tn]]$A <- chol(crossprod(B)) %*% Zt
    }
    ans <- do.call(if (!is.null(lmf$glmFit)) lme4:::glmer_finalize else lme4:::lmer_finalize, lmf)
    ans@call <- match.call()
    ans
}

    
    
