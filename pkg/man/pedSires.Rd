\name{pedSires}
\alias{pedSires}
\docType{data}
\title{Pedigree of the sires from dC}
\description{
  A \code{\linkS4class{pedigree}} object giving (part of) the pedigree
  of the sires from the \code{\link{dC}} data frame.  The
  pedigree is traced back on sires only.
}
%\usage{data(pedSires)}
\format{
  The format is:
Formal class 'pedigree' [package "pedigreemm"] with 3 slots
  ..@ sire : int [1:352] NA NA NA NA NA NA NA NA NA NA ...
  ..@ dam  : int [1:352] NA NA NA NA NA NA NA NA NA NA ...
  ..@ label: chr [1:352] "1" "2" "3" "4" ...
}
\source{
  ~~ reference to a publication or URL from which the data were obtained ~~
}
\examples{
str(pedSires)
}
\keyword{datasets}
