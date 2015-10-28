\name{sivipm-package}
\alias{sivipm-package}
\docType{package}
\title{ Sensitivity Indices with Dependent Inputs}
\description{Compute total and individual sensitivity indices, 
  significant components, and confidence
  intervals for the total sensitivity indices.}
\author{A. Bouvier, J.-P. Gauchi, E. Volatier

  Maintainer: Annie Bouvier <Annie.Bouvier@jouy.inra.fr>
  
}
\details{
The total and individual sensitivity indices are calculated
using a method based on the VIP of the PLS
regression, proposed by J.P. Gauchi.
The  significant components are calculated by the SIMCA software rule
  and by the Lazraq & Cléroux test.
The confidence
intervals for the total sensitivity indices are determined by
the bootstrap method.
A fast algorithm is implemented which allows big datasets analysis.
}
\references{
  \itemize{
    \item
    Gauchi, J.-P. and  Lehuta, S. and Mahévas, S. 2010.
    Optimal sensitivity analysis under constraints: Application to fisheries. In Procedia - Social and Behavioral Sciences, vol. 2. Elsevier,  pp. 7658-7659. Sixth International Conference on Sensitivity Analysis of Model Output (Milan, Italy), 19-22 July; co-organized by the ELEUSI research center of Bocconi University and by the Joint Research Center of the European Commission.
  \item
  Gauchi, J.-P. 2012.
    Global Sensitivity Analysis: The SIVIP method (SAS/IML language). 
  Rapport technique 2012-3. 
  INRA, UR1404, F-78350 Jouy-en-Josas, France.
  \item
  Gauchi, J.-P. 2015.
  A practical method of global sensitivity analysis
under constraints.
Rapport technique 2015-1.
INRA, UR1404, F-78350 Jouy-en-Josas, France.
\item
  Lazraq, A. and Cléroux, R. 2001.
  The PLS multivariate regression model: testing the significance of successive PLS components.
  Journal of Chemometrics.
  Vol. 15(6), pp 523-536.
  \item
  Lazraq, A. and Cléroux, R. and Gauchi, J.-P. 2003.
  Selecting both latent and explanatory variables in the PLS1 regression model. Chemometrics and Intelligent Laboratory Systems, Vol. 66(2),
  pp 117-126.
  Elsevier.         doi:10.1016/S0169-7439(03)00027-3.
  \item
  SIMCA Software.
\url{http://www.umetrics.com/products/simca}.
}
}

\keyword{ package }

