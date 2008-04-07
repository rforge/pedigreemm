\name{sirePed}
\alias{sirePed}
\docType{data}
\title{Pedigree of the sires in the Holstein data}
\description{
  A data frame with the pedigree of the sires mentioned in the
  \code{\link{Holstein}} data and auxillary information and the year of
  birth of the animal, where available.
}
\usage{data(sirePed)}
\format{
  A data frame with 38020 observations on the following 4 variables.
  \describe{
    \item{\code{label}}{International code of the animal.}
    \item{\code{sire}}{International code of the animal's sire (as a factor)}
    \item{\code{dam}}{International code of the animal's dam (as a factor).}
    \item{\code{yob}}{Year of birth (an ordered factor).}
  }
}
%\details{~~ If necessary, more details than the __description__ above ~~}
%\source{~~ reference to a publication or URL from which the data were obtained ~~}
%\references{~~ possibly secondary sources and usages ~~}
\seealso{\code{\link{dairyCattle}}}
\examples{
str(sirePed)
}
\keyword{datasets}
