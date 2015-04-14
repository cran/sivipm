\name{polyX-class}
\Rdversion{1.1}
\docType{class}
\alias{polyX-class}
%\alias{polyX}
%\alias{print,polyX-method}
%\alias{summary,polyX-method}
\alias{summary.polyX}
\alias{print.polyX}

\title{Class \code{"polyX"}}
\description{ An S4 class container for  a transformed dataset
  of inputs and its polynomial
  description. }
\section{Slots}{
  \describe{
    \item{\code{dataX.exp}}{  Transformed  X-inputs. A data.frame.
      The variables are
      the result of the calculation of the monomials of a
      polynomial on  raw X-inputs. }
    \item{\code{Pindic}}{ Object of class
      \code{\linkS4class{poly}}.}
 }}
 \section{Objects from the Class}{
   Objects from this class are created  by "new" or by the
   functions \code{\link{vect2polyX}}, \code{\link{vect2polyXT}},
   \code{\link{crpolyX}}, \code{\link{crpolyXT}}
 }
 \section{Methods}{
  \itemize{
    \item
    \code{summary, print,  show}: display functions
    \item
    \code{bind}: bind several  polyX objects.
    See \code{\link{bind.polyX}}
    }}
\keyword{classes}
