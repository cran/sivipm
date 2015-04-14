\name{crpolyX}
\alias{crpolyX}
\title{Create a \code{polyX} Object from
   a Standard Polynomial }
\description{Create an object of class \code{\linkS4class{polyX}}
    from  raw dataset of  X-inputs and  standard polynomial.}
\usage{
crpolyX(dataX, d, type = "full")
}
\arguments{
  \item{dataX}{ Raw X-inputs. A
data.frame with as many rows as observations, and as many
columns as  X-input variables.}
  \item{d}{ Polynomial degree.}
  \item{type}{ Polynomial type. A character string among:
    \itemize{
      \item
      \code{full}: complete polynomial with all monomials
	of degree less or equal to \code{d},
	\item
	\code{power}: power terms 	of degree less or equal to \code{d},
	\item
	\code{interact}:  interactions of degree  less or equal to \code{d}.
}}}
\note{The first monomials are always  the X-input variables.}
\value{An object of class \code{\linkS4class{polyX}}.}
\examples{
X <- cornell0[,1:3]
# Creation of the polynomial of degree 2 including the power terms:
# P= X1 + X2 + X3 + X1*X1 + X2*X2 + X3*X3 
P <- crpolyX(X, 2, type="power")
}
 
\seealso{
 \code{\link{vect2polyX}},  \code{\link{crpolyXT}}
}
