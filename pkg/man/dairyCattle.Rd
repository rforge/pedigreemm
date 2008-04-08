\name{dairyCattle}
\alias{dairyCattle}
\docType{data}
\title{Mastitis cases in dairy cattle}
\description{
  Records of the number of cases of clinical mastitis during the first
  lactation of a selection of dairly cattle, primarily Holsteins.
  Pedigree (sire only), herd and days in milk are also recorded for each
  cow.
%  The pedigree of the sires is given in the companion
%  \code{\link{sirePed}} data set.
}
\usage{data(dairyCattle)}
\format{
  A data frame with 48055 observations on the following 8 variables.
  \describe{
    \item{\code{id}}{International code of the cow.}
    \item{\code{sire}}{International code of the cow sire.}
    \item{\code{mastitis}}{a factor indicating if the cow had any
      incidents of clinical mastitis during the lactation.}
    \item{\code{NCM}}{An ordered factor giving the number of clinical
      mastitis cases for the cow during this lactation.}
    \item{\code{birth}}{year of birth of the cow (as an ordered factor).}
    \item{\code{herd}}{herd id number (as a factor).}
    \item{\code{calvingYear}}{year of calving for this lactation.}
    \item{\code{DIM}}{total number of days in milk for the lactation.}
  }
}
%\details{~~ If necessary, more details than the __description__ above ~~}
%\source{~~ reference to a publication or URL from which the data were obtained ~~}
%\references{~~ possibly secondary sources and usages ~~}
\seealso{\code{\link{pedigree}}}
\examples{
str(dairyCattle)
}
\keyword{datasets}
