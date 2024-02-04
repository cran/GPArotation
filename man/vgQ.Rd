\name{vgQ}
\alias{vgQ}
\alias{vgQ.oblimin}
\alias{vgQ.quartimin}
\alias{vgQ.target}
\alias{vgQ.pst}
\alias{vgQ.oblimax}
\alias{vgQ.entropy}
\alias{vgQ.quartimax}
\alias{vgQ.varimax}
\alias{vgQ.simplimax}
\alias{vgQ.bentler}
\alias{vgQ.tandemI}
\alias{vgQ.tandemII}
\alias{vgQ.geomin}
\alias{vgQ.bigeomin}
\alias{vgQ.cf}
\alias{vgQ.infomax}
\alias{vgQ.mccammon}
\alias{vgQ.bifactor}
\alias{vgQ.varimin}

\title{Rotations}
\usage{
    vgQ.oblimin(L, gam=0)
    vgQ.quartimin(L)
    vgQ.target(L, Target=NULL)
    vgQ.pst(L, W=NULL, Target=NULL)
    vgQ.oblimax(L)
    vgQ.entropy(L)
    vgQ.quartimax(L)
    vgQ.varimax(L)
    vgQ.simplimax(L, k=nrow(L))
    vgQ.bentler(L)
    vgQ.tandemI(L)
    vgQ.tandemII(L)
    vgQ.geomin(L, delta=.01)
    vgQ.bigeomin(L, delta=.01)
    vgQ.cf(L, kappa=0)
    vgQ.infomax(L)
    vgQ.mccammon(L)
    vgQ.varimin(L)
    vgQ.bifactor(L)
}
\arguments{
    \item{L}{a factor loading matrix}
    \item{gam}{0=Quartimin, .5=Biquartimin, 1=Covarimin.}
    \item{Target}{rotation target for objective calculation.}
    \item{W}{weighting of each element in target.}
    \item{k}{number of close to zero loadings.}
    \item{delta}{constant added to Lambda^2 in objective calculation.}
    \item{kappa}{see details.}
}

\value{A list (which includes elements used by \code{GPForth} and \code{GPFoblq}) with:
  \item{f}{The value of the criterion at L.}
  \item{Gq}{The gradient at L.}
  \item{Method}{A string indicating the criterion.}
   }
\description{
   vgQ routines to compute value and gradient of the criterion (not exported from NAMESPACE)
}
\details{
 The \code{vgQ.*} versions of the code are called by the optimization routine and
 would typically not be used directly, so these methods are not exported
 from the package NAMESPACE. (They simply return the function value and gradient
 for a given rotation matrix.) You can print these functions, but the package 
 name needs to be specified since they are not exported. For example, use
 \code{GPArotation:::vgQ.oblimin} to view the function \code{vgQ.oblimin}.
 The T or Q ending on function names should be omitted for the \code{vgQ.*} 
 versions of the code so, for example, use \code{GPArotation:::vgQ.target} to view the target
 criterion calculation.
 \tabular{lll}{
\code{vgQ.oblimin}   \tab orthogonal or oblique \tab  oblimin family\cr
\code{vgQ.quartimin} \tab oblique	   \tab         				\cr
\code{vgQ.target}    \tab orthogonal or oblique \tab  target rotation			\cr
\code{vgQ.pst}	   \tab orthogonal or oblique \tab  partially specified target rotation	\cr
\code{vgQ.oblimax}   \tab oblique	   \tab         				\cr
\code{vgQ.entropy}   \tab orthogonal \tab  minimum entropy			\cr
\code{vgQ.quartimax} \tab orthogonal \tab 					\cr
\code{vgQ.varimax}   \tab orthogonal \tab 					\cr
\code{vgQ.simplimax} \tab oblique	   \tab         				\cr
\code{vgQ.bentler}   \tab orthogonal or oblique \tab  Bentler's invariant pattern simplicity criterion\cr
\code{vgQ.tandemI} \tab orthogonal	   \tab Tandem principle I criterion	\cr
\code{vgQ.tandemII} \tab orthogonal	   \tab Tandem principle II criterion	\cr
\code{vgQ.geomin}    \tab orthogonal or oblique  \tab 				   \cr
\code{vgQ.bigeomin}    \tab orthogonal or oblique  \tab 				   \cr
\code{vgQ.cf}		   \tab orthogonal or oblique \tab  Crawford-Ferguson family	   \cr
\code{vgQ.cubimax}   \tab orthogonal    \tab  	\cr
\code{vgQ.infomax}   \tab orthogonal or oblique   \tab 				   \cr
\code{vgQ.mccammon}   \tab orthogonal    \tab McCammon minimum entropy ratio  \cr
\code{vgQ.varimin} \tab orthogonal \tab varimin criterion \cr
\code{vgQ.bifactor} \tab orthogonal or oblique \tab bifactor/biquartimin rotation\cr
 }
 See \link{rotations} for use of arguments.
 
 New rotation methods can be programmed with a name "vgQ.newmethod". The 
 inputs are the matrix L, and optionally any additional arguments. The
 output should be a list with elements \code{f}, \code{Gq}, and \code{Method}.
 
 Gradient projection \emph{without} derivatives can be performed using the
 \code{GPArotateDF} package; type \code{vignette("GPArotateDF", package = "GPArotation")}
 at the command line.
 }
 \examples{
  GPArotation:::vgQ.oblimin
  getAnywhere(vgQ.oblimax)	
  
}
\seealso{
   \code{\link{rotations}}
}
\references{
  Bernaards, C.A. and Jennrich, R.I. (2005) Gradient Projection Algorithms 
     and Software for Arbitrary Rotation Criteria in Factor Analysis. 
     \emph{Educational and Psychological Measurement}, \bold{65}, 676--696.

  }

\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert.}
\concept{rotation}
\keyword{multivariate}


