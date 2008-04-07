\name{sirePedigree}
\alias{sirePedigree}
\docType{data}
\title{Pedigree of the sires from dairyCattle}
\description{
  A \code{\linkS4class{pedigree}} object giving (part of) the pedigree
  of the sires from the \code{\link{dairyCattle}} data frame.  The
  pedigree is traced back on sires only.
}
\usage{data(sirePedigree)}
\format{
  The format is:
Formal class 'pedigree' [package "pedigreemm"] with 3 slots
  ..@ sire : int [1:39591] NA NA NA NA NA NA NA NA NA NA ...
  ..@ dam  : int [1:39591] NA NA NA NA NA NA NA NA NA NA ...
  ..@ label: chr [1:39591] "AYCAN000000051910" "AYCAN000000052113" "AYCAN000000055131" "AYCAN000000055385" ...
}
%\details{}
%\source{}
%\references{}
\examples{
str(sirePedigree)
}
\keyword{datasets}
