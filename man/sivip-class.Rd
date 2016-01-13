\name{sivip-class}
\docType{class}
\alias{sivip-class}
%\alias{sivip}
%\alias{print,sivip-method}
%\alias{summary,sivip-method}
\alias{print.sivip}
\alias{show.sivip}
\alias{summary.sivip}
\alias{getNames,sivip-method}
\alias{getNames}


\title{Class \code{"sivip"}}
\description{ An S4 class  container for the result of the function
\code{\link{sivipm}}. }
\section{Slots}{
  \describe{
     \item{\code{fo.isivip}}{ vector of length equal to the number
      of X-inputs. Individual sensitivity indices.}
    \item{\code{fo.isivip.percent}}{ sorted percentages of \code{fo.isivip}.}
    \item{\code{tsivip}}{ vector of length equal to the number
      of X-inputs. Total sensitivity indices. }
    \item{\code{tsivip.percent}}{  sorted percentages of \code{tsivip}.}
    \item{\code{monosignif}}{ a vector of length equal to the number
      of monomials.}
    \item{\code{correlalea}}{ a matrix. The number  of rows and columns
      is equal to the number
      of Y-inputs.}
    \item{\code{simca.signifcomponents}}{ a matrix. The number  of rows
      is equal to the number of components, and the number
      of columns is equal to  the number
      of Y-inputs plus one.}
    \item{\code{lazraq.signifcomponents}}{  a vector of length equal to the number
      of components.}
    \item{\code{output}}{ a list which contains the additional
      results.}}
See description in the function
\code{\link{sivipm}}.}
 \section{Objects from the Class}{
   Objects from this class are created by the function
\code{\link{sivipm}}. }
 \section{Methods}{
   \describe{
     \item{\code{ print}}{ display function.
     \code{print} has an optional argument, \code{all=FALSE}.
     When it is
     set to TRUE, all the components are displayed. Otherwise, the
     slots \code{fo.isivip}, \code{tsivip} and
     \code{output} are hidden.}
   \item{\code{ show}}{ equivalent to \code{ print(all=TRUE)}}
   \item{\code{ summary}}{ equivalent to \code{ print(all=FALSE)}}
   \item{\code{ getNames}}{ display the names of the non-null components}}
}

\keyword{classes}
 