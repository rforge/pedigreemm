\name{dC}
\alias{dC}
\docType{data}
\title{Mastitis cases in dairy cattle}
\description{
  Records of the number of cases of clinical mastitis during the first
  lactation of a selection of dairly cattle, primarily Holsteins.
  The sire, herd and days in milk are also recorded for each
  cow.
}
%\usage{data(dC)}
\format{
  A data frame with 1675 observations on the following 8 variables.
  \describe{
    \item{\code{id}}{Identifier of the animal.}
    \item{\code{sire}}{Identifier of the animal's sire.}
    \item{\code{birth}}{year of birth of the animal (as a factor).}
    \item{\code{herd}}{herd id number (as a factor).}
    \item{\code{calvingYear}}{year of calving for this lactation.}
    \item{\code{DIM}}{total number of days in milk for the lactation.}
    \item{\code{mastitis}}{a factor indicating if the cow had any
      incidents of clinical mastitis during the lactation.}
    \item{\code{NCM}}{An ordered factor giving the number of clinical
      mastitis cases for the cow during this lactation.}
  }
}
\details{The pedigree of the sires is given in the companion
  \code{\link{pedSires}} data set.}
%\source{~~ reference to a publication or URL from which the data were obtained ~~}
%\references{~~ possibly secondary sources and usages ~~}
\seealso{\code{\link{pedSires}}, \code{\link{pedigree}}}
\examples{
str(dC)
summary(dC, maxsum = 10)
}
\keyword{datasets}
