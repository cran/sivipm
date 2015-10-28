\name{vect2polyX}
\alias{vect2polyX}
\title{Create  polyX Object from  Polynomial Description}
\description{Create an object of class \code{\linkS4class{polyX}}
  from  a dataset of X-inputs and a vector that
  describes the polynomial.}
\usage{vect2polyX(dataX, monomials)}
\value{An object of class \code{\linkS4class{polyX}}.}
\arguments{
    \item{dataX}{
 Raw X-inputs. A
data.frame with as many rows as observations, and as many
columns as  variables.
}
\item{monomials}{Polynomial description. A character vector. Each element describes a monomial.
  The input variables are coded either by their numbers, or by their
  names. The character "*" denotes interaction between variables. }
  }

\seealso{
    \code{\link{crpolyX}},
  \code{\link{vect2polyXT}}
}
\examples{
X <- cornell0[,1:3]
# Monomials expressed by variable numbers:
monomials <- c("1","2","3", "1*2*3")
polyXI <- vect2polyX (X, monomials)
# Monomials expressed by variable names:
monomials <- c("Distillation", "Reformat", "NaphthaT",
               "Distillation*Reformat*NaphthaT")
polyXI <- vect2polyX (X, monomials)
}
