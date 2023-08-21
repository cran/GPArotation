\name{print.GPArotation}
\alias{print.GPArotation}
\alias{summary.GPArotation}
\alias{print.summary.GPArotation}

\title{Print and Summary Methods for GPArotation}
\usage{
    \method{print}{GPArotation}(x, digits=3, sortLoadings=TRUE, rotateMat=FALSE, Table=FALSE, ...)
    \method{summary}{GPArotation}(object, digits=3, Structure=TRUE, ...)
    \method{print}{summary.GPArotation}(x, ...)
}
\arguments{
    \item{x}{a GPArotation object to summarize.}
    \item{object}{a summary.GPArotation to print.}
    \item{digits}{precision of printed numbers.}
    \item{sortLoadings}{display sorted loadings of a GPArotation object.}
    \item{rotateMat}{display the rotation matrix of a GPArotation object.}
    \item{Table}{display the convergence table of a GPArotation object.}
    \item{Structure}{display the structure matrix of a summary.GPArotation object of a oblique rotation.}
    \item{...}{further arguments passed to other methods.}
}

\value{The object printed or a summary object.}
\description{
   Print an object or summary of an object returned by \code{GPFRSorth},
   \code{GPFRSoblq}, \code{GPForth}, or \code{GPFoblq}.
   
   A GPArotation object by default to will print sorted loadings, Phi, and rotation matrix (if requested).
   
   Output includes contributions of factors \code{SS loadings} 
   (Sum of Squared loadings), (see e.g. Harman 1976, sections 2.4 and 12.4).
}

\references{
   Harman, H.H. (1976). \emph{Modern Factor Analysis.} The University of Chicago Press. 

}\seealso{
   \code{\link{GPForth}}, 
   \code{\link[base]{summary}}
}

\concept{rotation}
\keyword{internal}
