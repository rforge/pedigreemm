\name{pedigreemm-class}
\docType{class}
\alias{pedigreemm-class}
\title{Pedigree-based mixed-effects model fits}
\description{A mixed-effects model fit by \code{\link{pedigreemm}}.
  This class extends class \code{"\linkS4class{mer}"} class and includes one
  additional slot, \code{relfac}, which is a list of (left) Cholesky
  factors of the relationship matrices derived from 
  \code{"\linkS4class{pedigree}"} objects.}
\section{Objects from the Class}{Objects are created by calls to the
  \code{\link{pedigreemm}} function.}
\section{Slots}{
	 \describe{
	   \item{\code{relfac}:}{A list of relationship matrix factors.  All
	     other slots are inherited from class \code{"\linkS4class{mer}"}.}
	   \item{\dots}{all other slots are inherited from the superclass
	     \code{"\linkS4class{mer}"}.
	   }
	 }
       }

\section{Extends}{
Class \code{"\linkS4class{mer}"}, directly.
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
