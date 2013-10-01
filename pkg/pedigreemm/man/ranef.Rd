\name{ranef}
\alias{ranef}
\title{Random Effects}
\description{
  A generic function to extract, and optionally accumulate, the random effects. 
}

\usage{
ranef( object, postVar = FALSE, drop = FALSE, whichel = names(ans), pedigree = TRUE, ...)
}

\arguments{
  \item{object}{a named object \code{\link[pedigreemm]{pedigreemm}}}
  \item{postVar}{an optional logical argument indicating if the conditional 
   variance covariance matrices, also called the posterior variances, of the 
    random effects should be included. Default is FALSE. \code{\link[lme4]{lmer}}}
  \item{drop}{a logical argument indicating whether random effects in the pedigree without data should be droped.}
  \item{pedigree}{a logical argument indicating whether the random effect should be linked to the pedigree}
  \item{\dots}{as in \code{\link[lme4]{lmer}}}
}
\details{
  All arguments to this function are the same as those to the function
  \code{\link[lme4]{lmer}} except pedigree which must be a logical argument and drop
}
\value{
  A list of data frames, one for each grouping factor for the random effects. 
  The number of rows in the data frame is the number of levels of the grouping factor. 
  The number of columns is the dimension of the random effect associated with each level of the factor.
  The random effect that is associated to a pedigree will be transformed by premultipling the
  random effect from lme4 by the cholesky factor of the relationship matrix.
}
\references{
   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel. 
   Technical Note: An R package for fitting generalized linear mixed models 
   in animal breeding. Journal of Animal Science, 88:497-504.
}
\seealso{\code{\linkS4class{pedigreemm}}, \code{\linkS4class{pedigree}},
  \code{\link{relfactor}}.}
\examples{

  p1 <- new("pedigree",
           sire = as.integer(c(NA,NA,1, 1,4,5)),
           dam  = as.integer(c(NA,NA,2,NA,3,2)),
           label = as.character(1:6))
  A<-getA(p1)
  cholA<-chol(A)  
  varU<-0.4; varE<-0.6; rep<-20
  n<-rep*6
  set.seed(108)
  bStar<- rnorm(6, sd=sqrt(varU))
  b<-crossprod(as.matrix(cholA),bStar)
  ID <- rep(1:6, each=rep)
  e0<-rnorm(n, sd=sqrt(varE))
  y<-b[ID]+e0
  fm1 <- pedigreemm(y ~ (1|ID) , pedigree = list(ID = p1))
  cor(ranef(fm1)$ID, as.vector(b))   
}
\keyword{Random effects}
