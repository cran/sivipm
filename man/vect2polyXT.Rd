\name{vect2polyXT}
\alias{vect2polyXT}
\title{Create a polyX Object from  Polynomial Description}
\description{Create an object of class \code{\linkS4class{polyX}}
  from a dataset of transformed  inputs and a vector that
  describes the polynomial.}
\usage{vect2polyXT(varnames, dataXT, monomials)}
\value{An object of class \code{\linkS4class{polyX}}.}
\arguments{
  \item{varnames}{ X-input names. A character vector of length equal to
    the number of X-inputs.}
    \item{dataXT}{
 Transformed inputs. A
data.frame with as many rows as observations, and as many
columns as monomials.
}
\item{monomials}{Polynomial description. A character vector. Each element describes a monomial.
  The input variables are coded either by their numbers, or by their
  names. The character "*" denotes interaction between variables. See \code{\link{vect2polyX}} example.}
}
\seealso{
    \code{\link{crpolyXT}},
  \code{\link{vect2polyX}}
}
