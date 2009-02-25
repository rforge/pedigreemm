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
    \item{\code{env}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{nlmodel}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{frame}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{call}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{flist}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{X}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{Xst}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{Zt}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{pWt}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{offset}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{y}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{Gp}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{dims}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{ST}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{V}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{A}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{Cm}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{Cx}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{L}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{deviance}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{fixef}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{ranef}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{u}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{eta}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{mu}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{muEta}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{var}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{resid}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{sqrtXWt}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{sqrtrWt}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{RZX}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{RX}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{ghx}:}{as in class \code{"\linkS4class{mer}"}.}
    \item{\code{ghw}:}{as in class \code{"\linkS4class{mer}"}.}
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
