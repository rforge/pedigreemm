\name{pedigreemm-class}
\docType{class}
\alias{pedigreemm-class}
\title{Pedigree-based mixed-effects model fits}
\description{A mixed-effects model fit by \code{\link{pedigreemm}}.
  This class extends class \code{"\link[lme4:merMod-class]{merMod}"} class and includes one
  additional slot, \code{relfac}, which is a list of (left) Cholesky
  factors of the relationship matrices derived from 
  \code{"\linkS4class{pedigree}"} objects.}
\section{Objects from the Class}{Objects are created by calls to the
  \code{\link{pedigreemm}} function.}
\section{Slots}{
  \describe{
    \item{\code{relfac}:}{A list of relationship matrix factors.  All
      other slots are inherited from class \code{"\link[lme4:merMod-class]{merMod}"}.}
  }
}
\section{Extends}{
Class \code{"\link[lme4:merMod-class]{merMod}"}, directly.
}
\section{Methods}{
No methods defined with class "pedigreemm" in the signature.
}
\seealso{\code{\link{pedigreemm}}
}
\examples{
showClass("pedigreemm")
}
\keyword{classes}
