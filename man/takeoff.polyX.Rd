\name{takeoff.polyX}
\docType{methods}
\alias{takeoff.polyX}
\alias{takeoff}
\alias{takeoff,polyX-method}
\title{Remove Monomials from a polyX Object}
\description{
  Remove monomials from an object  of class \code{\linkS4class{polyX}}.
}
\usage{
\S4method{takeoff}{polyX}(P, monomials)
}
\arguments{
  \item{P}{
Object of class \code{\linkS4class{polyX}}.
}
  \item{monomials}{
    Integer or character vector. The monomials to be removed.
    When integer,  the numbers of the
    monomials. When character,  their 
    expression using the variable names and the
character "*" to denote interaction.
}
}
\value{An object of class \code{\linkS4class{polyX}}.
}
\seealso{
\code{\link{bind.polyX}}
}
\examples{
X <- cornell0[,1:2] # X-inputs
# Polynomial creation
monomials <- c("Distillation","Reformat","Distillation*Reformat")
P <- vect2polyX (X, monomials)
# Remove the third  monomial:
P2 <- takeoff(P, c(3))
# Same expressed with using variable names:
P2 <- takeoff(P, c("Distillation*Reformat"))
}
