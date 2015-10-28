\name{sivipm}
\alias{sivipm}
\title{ Sensitivity Indices }
\usage{sivipm(Y,  XIndic,
        nc=2,  options= c("fo.isivip", "tsivip", "simca", "lazraq"),
        graph = FALSE, alea = FALSE, fast = FALSE,
        output=NULL)}
\description{Compute total and individual sensitivity indices, 
 significant components and auxiliary results.}
\arguments{
\item{Y}{ Outputs. A data.frame  with as
    many rows as observations and
    as many columns  as response variables.}
  \item{XIndic}{  Object of class \code{\linkS4class{polyX}}
    which contains the polynomial
  description and the dataset of inputs.}
\item{nc}{  Required number of components.}
\item{options}{Options to select what is calculated. A string vector.
  Valid values are:
  \itemize{
    \item{\code{fo.isivip}}{ first order individual sensitivity indices,}
    \item{\code{tsivip}}{ total sensitivity indices,}
    \item{\code{simca}}{ significant components calculated by
      the SIMCA software rule.
      See \code{Details}.}
    \item{\code{lazraq}}{ significant components calculated by
      the Lazraq and Cléroux test.   See \code{Details}.
}}}
	\item{graph}{If TRUE, a graph is drawn when \code{options} includes \code{tsivip}.}
	 \item{alea}{ If TRUE, an uniform random variable is included
	   in the analysis when
	   \code{options} includes \code{tsivip}.
	   Then, the non significant monomials are excluded from the
	   total
	   sensitivity indices calculation.}
	 \item{fast}{ If TRUE,  \code{Q2} and relative results
	   are calculated from  the Miller's formulae more adapted to
	   big datasets.} 
	 \item{output}{If non NULL, additional results are returned
	   in a  component named \code{output}. 
	   Character vector, which valid values are:
	   \itemize{
	     \item{\code{isivip}}{ to return \code{isivip}}
	     \item{\code{betaNat}}{ to return \code{betaNat} and
	       \code{betaNat0}}
	     \item{\code{VIP}}{ to return \code{VIP} and
	       \code{VIPind}}
%	     \item{\code{RSS}}{ to return \code{RSS}}
%	     \item{\code{PRESS}}{ to return \code{PRESS}}
	     \item{\code{Q2}}{ to return \code{Q2}, \code{Q2ckh} and
	       \code{Q2cum}}
	     \item{\code{PLS}}{ to return PLS results:
		 \code{mweights}, \code{weights}, \code{x.scores},
		 \code{x.loadings}, \code{y.scores},
		 \code{y.loadings}, \code{cor.tx}, \code{cor.ty},
		 \code{expvar}, \code{X.hat}, \code{Y.hat}
	       } } See "Value".\cr
	     It is advised to first determine the number of
	     significant components, by setting the options
	     \code{simca} or \code{lazraq}, before asking for
	     additional results.
	   }
	   }
	   \details{
	     \itemize{
	   \item    
  When the option \code{simca} or \code{lazraq}
  is set,  the significant components
  are computed   by the SIMCA software rule, or,
  by the Lazraq and Cléroux inferential test, at confidence level
  0.95, respectively.
 The option \code{simca} is ignored if there are
 missing values. The option \code{lazraq} is ignored if there are
 missing values and more than one response variables.
 \item
 When the option \code{alea} is set, 
the non significant monomials are those for which
the individual sensitivity indices is less or equal than
the one of the  random variable. These non significant monomials
are excluded from the  total sensitivity indices calculation.
\item
To analyze big datasets, the option \code{fast} is advised.
}
 }
\note{
 If the output is multivariate, \code{tsivip}
 are the generalized total sensitivity indices (GTSIVIP)
 and
 \code{isivip} are the generalized individual sensitivity
  indices (GISIVIP). }  
\value{An object of class \code{\linkS4class{sivip}},
  whose slots are:
 \item{\code{fo.isivip}}{When \code{options} includes \code{fo.isivip},
  the first order individual sensitivity indices.
 }
\item{\code{tsivip}}{When \code{options} includes \code{tsivip},
         total sensitivity indices for each input variable.}
       \item{\code{percent}}{When \code{options} includes \code{tsivip},
	 percentages of the total sensitivity indices.}

       \item{\code{monosignif}}{When \code{alea} is TRUE,
	 and \code{options} includes \code{tsivip},
	 logical vector which indicates the significant monomials.}
       \item{\code{correlalea}}{When \code{alea} is TRUE,
	 and \code{options} includes \code{tsivip},
	 the correlation matrix between the random variable
	 and the outputs.}
       
 \item{\code{simca.signifcomponents}}{When  \code{options} includes \code{simca},
   the significant components calculated by the S. Wold's rule
   (SIMCA software rule). Logical matrix with \code{nc} rows
   and as many columns as response variables plus one.
   Values are
 TRUE for the components where the test is positive at 95\% level, FALSE
   otherwise. The last column is the test applied to the Q2T (Q2
   cumulated on the response variables).
 }
  \item{\code{lazraq.signifcomponents}}{When  \code{options} includes \code{simca},
   the significant components calculated by the Lazraq and Cléroux inferential test. Logical
   vector of length \code{nc} with
 TRUE for the components where the test is positive at 95\% level, FALSE otherwise. }
 \item{\code{output}}{ When \code{output} is not NULL, a list with
   additional results, whose components depend on  \code{output}
   value. 
   \itemize{
     \item
     \code{isivip} Individual sensitivity indices for each
     monomial. Vector of length equal to the number
     of monomials.
 \item
 \code{betaNat} Natural \code{beta}.
 Matrix with as many
   rows as monomials and as many columns as response variables.
 \item
 \code{betaNat0} Natural \code{beta0} coefficient. Vector of length
     equal to the number of  response variables.
   
   \item  \code{VIP} Matrix of \code{nc} columns and as many rows as monomials.
   \item  \code{VIPind}  Matrix with as many rows as response variables and
 as many columns as monomials. 
%     \item  \code{RSS} Matrix  with as many columns as response variables
%     and \code{nc} +1 rows.
%     \item  \code{PRESS} Matrix  with as many columns as response variables
%      and \code{nc} rows.
     \item  \code{Q2} Matrix  with  as many columns as response
   variables +1 and \code{nc} rows. The last column is \code{Q2T}.
     \item  \code{Q2ckh} Matrix  with as many columns as response
     variables
     and \code{nc} rows (Q2 cumulated on the components).      
   \item  \code{Q2cum}, \code{Q2cum.percent} Matrix  with one column and  \code{nc} rows.

	    \item \code{PLS} PLS results. The dimension of the
	    components are indicated below in brackets. \code{nmono} denotes the number of monomials, \code{ny}, the number of response variables and \code{nobs}
the number of observations. 
\itemize{
\item   \code{betaCR} (\code{beta} centered and reduced. Vector \code{ny}),
\item     \code{mweights} (\code{nc} \code{X} \code{nmono}), 
\item      \code{x.scores} (\code{nc} \code{X} \code{nobs}),
\item     \code{x.loadings} (\code{nc} \code{X} \code{nmono}),
\item     \code{y.scores} (\code{nc} \code{X} \code{nobs}),
\item     \code{y.loadings} (\code{nc} \code{X} \code{ny}),
\item     \code{weights} (\code{nc} \code{X} \code{nmono}), 
\item     \code{cor.tx} (\code{nc} \code{X} \code{nmono}),
\item     \code{cor.ty} (\code{nc} \code{X} \code{ny}),
 \item    \code{expvar} (\code{4} \code{X} \code{nc}),
 \item      \code{x.hat}(\code{nobs} \code{X} \code{nmono}),
 \item    \code{y.hat} (\code{nobs} \code{X} \code{ny}).
 }}}
 }
 \examples{
X <- cornell0[,1:3] # X-inputs
Y <- as.data.frame( cornell0[,8]) # response variable
# Creation of the polynomial:
P <- vect2polyX(X, c("1", "2", "3", "3*3*3"))
# Compute total sensitivity indices:
A <- sivipm(Y, P, options=c("tsivip"))
# See the names of the returned components
getNames(A)
# The main results
summary(A)
# Calculation by using the fast algorithm:
B <- sivipm(Y, P, fast = TRUE, options=c("tsivip"))
}

