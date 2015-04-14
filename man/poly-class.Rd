\name{poly-class}
\Rdversion{1.1}
\docType{class}
\alias{poly-class}
%\alias{poly}
%\alias{print,poly-method}
%\alias{summary,poly-method}
\alias{print.poly}
\alias{summary.poly}


\title{Class \code{"poly"}}
\description{ An S4 class  container for a polynomial
  description. }
\section{Slots}{
  \describe{
    \item{\code{P}}{  List of length equal to the number of monomials.
     Each component is a vector of length equal to the
     polynomial degree. It contains the
     numbers of the variables in the monomial.}
   \item{\code{indic}}{  Matrix with as
        many rows as monomials and as many columns as X-inputs.
       The element \code{(i,j)} is 1 when the variable \code{j}
       is in the monomial, 0 otherwise.}
 }}
 \section{Objects from the Class}{
   Objects from this class are usually not intended to
   be accessed directly by the end-user. }
 \section{Methods}{
    \code{summary, print, or show}{: display functions.}
}
\keyword{classes}
 