\name{sivipboot}
\alias{sivipboot}
\title{ Confidence Intervals for the Total Sensitivity Indices by  Bootstrap}
\description{Confidence intervals for the total sensitivity indices by a  bootstrap
  method.}
\usage{sivipboot(Y, XIndic, B, nc=2, graph=FALSE, alea=FALSE,
                 alpha=0.05)}
\value{A matrix with as many rows as input variables and two columns:
  the lower and upper bounds of
  the total sensitivity indices percentile bootstrap confidence intervals.}
\arguments{
  \item{Y}{ Outputs. A data.frame with as
    many rows as observations and
    as many columns  as response variables.}
 \item{XIndic}{  Object of class \code{\linkS4class{polyX}}
    which contains the polynomial
    description and the dataset of inputs.}
\item{B}{ Number of bootstrap replicates.}
\item{nc}{ Number of components.}
\item{graph}{If TRUE, boxplot display.}
\item{alea}{ If TRUE, an uniform random variable is included
	   in the analysis (see \code{\link{sivipm}}).}
	 \item{alpha}{ Level of the bootstrap confidence intervals.}}
\seealso{
  \code{\link{sivipm}}
}
