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
setAs("pedigree", "dtCMatrix", # representation as T^{-1}
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


Linv <- function(ped)
{
    stopifnot(is(ped, "pedigree"))
    as(Diagonal(x = sqrt(1/(1 +.Call("pedigree_inbreeding", ped)))) %*%
       as(ped, "dtCMatrix"), "triangularMatrix")
}

pedmat <- function(Name, pedigree, type = c("id", "sire", "dam"))
### Should return the sparse forms of the Cholesky factor of the
### relationship matrix, reduced in the case of type = "sire" or "dam"    
{
    stopifnot(is(pedigree, "pedigree"))
    typ <- match.arg(type)
    nm <- as.name(Name)
    Tinv <- as(pedigree, "dtCMatrix")
    if (type == "sire") {
        ## Don't do this!
        ## generate the T matrix for the sires only
       ## ind <- as(diag(length(pedigree@label)),
         ##         "sparseMatrix")[, sort(unique(na.omit(pedigree@sire)))]
    }
}

setAs("pedigree", "pedinbred",
      function(from)
      new("pedinbred", sire = from@sire, dam = from@dam, label = from@label,
          F = .Call(pedigree_inbreeding, from)))

#' Diagonal factor of the relationship inverse
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
Ddiag <- function(ped) {
    stopifnot(is(ped, "pedigree"))
    F <- .Call("pedigree_inbreeding", ped, PACKAGE = "pedigreemm")
    sire <- ped@sire
    dam <- ped@dam
    Fsire <- ifelse(is.na(sire), -1, F[sire])
    Fdam <-  ifelse(is.na(dam), -1, F[dam])
    1 - 0.25 * (2 + Fsire + Fdam)
}

#' Relationship matrix from a pedigree
#'
#' Determine the relationship matrix for the pedigree \code{ped},
#' possibly restricted to the specific labels that occur in \code{labs}.
#' 
#' @param ped a pedigree that includes the individuals who occur in svec
#' @param labs a character vector or a factor giving the labels to
#'    which to restrict the relationship matrix. If \code{labs} is a
#'    factor then the levels of the factor are used as the labels.
#' @param drop.unused.levels logical scalar -- should unused levels in
#'    \code{labs} be dropped before matching to the pedigree.  Defaults
#'    to TRUE.  Only used when \code{labs} is a factor.
#' @return a sparse, symmetric relationship matrix
#' @export

relmat <- function(ped, labs = NULL, drop.unused.levels = TRUE)
{
    stopifnot(is(ped, "pedigree"))
    if (!is.null(labs)) {
        if (is(labs, "factor"))
            labs <- levels(labs[,drop = drop.unused.levels])
        stopifnot(all(levs %in% ped@label))
        Tt <-
            solve(t(as(ped, "dtCMatrix")),
                  Matrix:::fac2sparse(factor(levs, levels =
                                             ped@label), drop = FALSE))
    } else Tt <- solve(t(as(ped, "dtCMatrix")))
    crossprod(Diagonal(x = sqrt(Ddiag(ped))) %*% Tt)
}
