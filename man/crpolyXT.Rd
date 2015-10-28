\name{crpolyXT}
\alias{crpolyXT}
\title{Create polyX Object from
  Standard  Polynomial, Transformed Inputs}
\description{Create an object of class \code{\linkS4class{polyX}}
   from transformed X-inputs and  standard  polynomial.}
\usage{
crpolyXT(varnames, dataXT, d, type = "full")
}
\arguments{
  \item{varnames}{ X-input names. A character vector of length equal to
    the number of X-inputs.}
   \item{dataXT}{ Transformed X-inputs. A
     data.frame with as many rows as observations, and which the number
     of columns is equal or greater than the number of monomials.
}
  \item{d}{ Polynomial degree.}
  \item{type}{ Polynomial type: see \code{\link{crpolyX}.}
}
}
\value{An object of class \code{\linkS4class{polyX}}.}
 \note{The first monomials are always  the X-input variables.}
\seealso{
 \code{\link{vect2polyXT}},   \code{\link{crpolyX}}
}
