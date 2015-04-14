\name{factorsplit}
\alias{factorsplit}
\title{
Transform  Categorical Variables into Numeric Variables
}
\description{
  Transform each factor column of a dataframe into  numeric variables.
  \itemize{
    \item If the factor has more than 2 levels, it is split
   into as many
  0-1 columns as levels. The value 1 is set on
  the lines where the factor is equal to the level,
and 0 elsewhere.
\item If the factor has one or two levels, the first level is coded
 \code{-1} and the second one into \code{+1}.
}
The numeric columns are left unchanged.
}
\usage{
factorsplit(data)
}
\arguments{
  \item{data}{
 Data.frame.
}
}
\value{
  The data.frame \code{data}, where the factor columns are replaced
  by 0/1 or +1/-1 columns.
  The labels of the added columns are the concatenation
  of the name of the factor and the name of the level, with
  the underscore separator.
}
\note{
  When using the other functions of the package on a
  result of \code{factorsplit},
   be careful when describing the polynomial:
  the number and names of the variables in the returned data.frame
  may be different than the genuine ones.
}

\seealso{
\code{\link{crpolyX}},  \code{\link{vect2polyX}}
}
\examples{
a <- data.frame(V=1:3, color=c("red","green","black"),
                       temp=c("hot", "cold", "hot"))
b <- factorsplit(a)
# The columns in b are: V, color_red, color_green,  color_red, temp
}
