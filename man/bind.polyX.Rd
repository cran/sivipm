\name{bind.polyX}
\docType{methods}
\alias{bind.polyX} % obligatoire
\alias{bind} % si on l' ote: warning
\alias{bind,polyX-method} % obligatoire pour S4method
\title{Bind Several  polyX Objects}
\description{
Bind several objects of class \code{\linkS4class{polyX}}.
}
\usage{
\S4method{bind}{polyX}(x, \ldots)
}
\arguments{
  \item{x}{
    Object of class \code{\linkS4class{polyX}}.
}
  \item{\dots}{Objects
  of class \code{\linkS4class{polyX}} to be binded to \code{x}. As many
  as required. It is  assumed that they have
      been  built on
the same dataset of inputs.}
}
\value{An object of class \code{\linkS4class{polyX}}, coding for
  a polynomial including the monomials of all
  the polynomials in the argument list. Its degree is the maximal
  degree of these polynomials.  Duplicated monomials are removed.}
\seealso{
  \code{\link{takeoff.polyX}}
  }

\examples{
# Create first polynomial
P <- vect2polyX(cornell0[,1:3],c("1","2","3", "3*3*3", "3*3"))
# Create second polynomial
P2 <- vect2polyX(cornell0[,1:3], c("1","2","3", "2*3"))
# Bind them: PP=X1 + X2 + X3 + X3*X3*X3 + X3*X3 + X2*X3
PP <- bind.polyX(P, P2)
}
\keyword{methods}
