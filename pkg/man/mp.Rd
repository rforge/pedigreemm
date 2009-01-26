\name{mp}
\alias{mp}
\docType{data}
\title{Milk production}
\description{
  Records of the milk production of cattle.
}
\usage{data(mp)}
\format{
  A data frame with 3397 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{numeric identifier of cow}
    \item{\code{lact}}{number of lactation for which production is measured}
    \item{\code{herd}}{a factor indicating the herd}
    \item{\code{sire}}{a numeric vector}
    \item{\code{dim}}{a numeric vector}
    \item{\code{milk}}{a numeric vector}
    \item{\code{fat}}{a numeric vector}
    \item{\code{prot}}{a numeric vector}
    \item{\code{scs}}{a numeric vector}
  }
}
% \source{
% }
% \references{
% }
\examples{
str(mp)
}
\keyword{datasets}
