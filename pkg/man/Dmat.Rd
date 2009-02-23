\name{inbreeding}
\alias{inbreeding}
\title{Inbreeding coefficients from a pedigree...}
\usage{inbreeding(ped)}
\description{Inbreeding coefficients from a pedigree}
\details{
  Determine the diagonal factor in the inverse of the relationship
  matrix from a pedigree.  Although in practice we use the relationship
  matrix from a pedigree, it is easier to evaluate the factors of the
  inverse of the relationship matrix.  This function returns the diagonal
  matrix from the LDL' form of the Cholesky factorization of the inverse
  of the relationship matrix.
}
\value{a numeric vector}
\arguments{\item{ped}{an object that inherits from class \linkS4class{pedigree}}}
\examples{
ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
inbreeding(ped)
}
\keyword{misc}
