\name{Holstein}
\alias{Holstein}
\docType{data}
\title{Mastitis cases in Holstein cows}
\description{
  Information on the pedigree and herd of Holstein cows.
}
\usage{data(Holstein)}
\format{
  A data frame with 42607 observations on the following 8 variables.
  \describe{
    \item{\code{maSum}}{An ordered factor giving the number of clinical
      mastitis cases for the cow during this lactation.}
    \item{\code{id}}{International code of the cow.}
    \item{\code{sire}}{International code of the cow's sire.}
    \item{\code{mastitis}}{a factor indicating if the cow had any
      incidents of clinical mastitis during the lactation.}
    \item{\code{birth}}{year of birth of the cow (as an ordered factor).}
    \item{\code{herd}}{herd id number.}
    \item{\code{fresh}}{year of calving for this lactation.}
    \item{\code{dim}}{total number of days in milk for the lactation.}
  }
}
%\details{~~ If necessary, more details than the __description__ above ~~}
%\source{~~ reference to a publication or URL from which the data were obtained ~~}
%\references{~~ possibly secondary sources and usages ~~}
\examples{
str(Holstein)
}
\keyword{datasets}
