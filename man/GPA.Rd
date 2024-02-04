\name{GPA}
\alias{GPFRSorth}
\alias{GPFRSoblq}
\alias{GPForth}
\alias{GPFoblq}

\title{Rotation Optimization}
\usage{
    GPFRSorth(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
       method="varimax", methodArgs=NULL, randomStarts=0)
    GPFRSoblq(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
       method="quartimin", methodArgs=NULL, randomStarts=0)

    GPForth(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
       method="varimax", methodArgs=NULL)
    GPFoblq(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
       method="quartimin", methodArgs=NULL)
    }
\arguments{
    \item{A}{initial factor loadings matrix for which the rotation 
             criterian is to be optimized.}
    \item{Tmat}{initial rotation matrix.}
    \item{normalize}{see details.}
    \item{eps}{convergence is assumed when the norm of the gradient is smaller
         than eps.}
    \item{maxit}{maximum number of iterations allowed in the main loop.}
    \item{method}{rotation objective criterian.}
    \item{methodArgs}{a list ofmethodArgs arguments passed to the rotation objective}
    \item{randomStarts}{number of random starts (GPFRSorth and GPFRSoblq)}
}

\description{
   Gradient projection rotation optimization routine used by various rotation
   objective.
   }
\value{A GPArotation object which is a list with elements 
  \item{loadings}{The rotated loadings, one column for each factor. If randomStarts
  		were requested then this is the rotated loadings matrix with the 
  		lowest criterion value.}
  \item{Th}{The rotation matrix, loadings \%*\% t(Th) = A.}
  \item{Table}{A matrix recording the iterations of the rotation optimization.}
  \item{method}{A string indicating the rotation objective function.}
  \item{orthogonal}{A logical indicating if the rotation is orthogonal.}
  \item{convergence}{A logical indicating if convergence was obtained.}
  \item{Phi}{t(Th) \%*\% Th.  The covariance matrix of the rotated factors.
  This will be the identity matrix for orthogonal
  rotations so is omitted (NULL) for the result from GPFRSorth and GPForth.}
  \item{Gq}{The gradient of the objective function at the rotated loadings.}
  \item{randStartChar}{A vector with characteristics of random starts 
  	(GPFRSorth and GPFRSoblq only; omitted if randomStarts =< 1).}
  }
\details{
 Gradient projection (GP) rotation optimization routines developed by
 Jennrich (2001, 2002) and Bernaards and Jennrich (2005).
 These functions can be used directly
 to rotate a loadings matrix, or indirectly through a rotation
 objective passed to a factor estimation routine such as \code{\link{factanal}}.
 A rotation of a matrix \code{A} is defined as \code{A \%*\% solve(t(Th))}. 
 In case of orthogonal rotation, the factors the rotation matrix
 \code{Tmat} is orthonormal, and the rotation simplifies 
 to \code{A \%*\% Th}. The rotation matrix \code{Th} is computed by GP rotation.
 
 The \code{GPFRsorth} and \code{GPFRSoblq} functions are the primary functions 
 for orthogonal and oblique rotations, respectively. These two functions 
 serve as wrapper functions for  \code{GPForth} and \code{GPFoblq}, with the
 added functionality of multiple random starts. 
 \code{GPForth} is the main GP algorithm for orthogonal rotation.
 \code{GPFoblq} is the main GP algorithm for oblique rotation.
 The \code{GPForth} and \code{GPFoblq} may be also be called directly. 
 
 Arguments in the wrapper functions \code{GPFRsorth} and \code{GPFRSoblq} 
 are passed to GP algorithms. Functions require an initial loadings matrix 
 \code{A} which fixes the equivalence class over which the optimization is done. 
 It must be the solution to the orthogonal factor analysis problem as obtained 
 from \code{factanal} or other factor estimation routines. 
 The initial rotation matrix is given by the \code{Tmat}. 
 By default the GP algorithm use the identity matrix as the initial rotation matrix. 
 
 For some rotation criteria local minima may exist. To start from random 
 initial rotation matrices, the \code{randomStarts} argument is available
 in \code{GPFRSorth} and \code{GPFRSoblq}. The returned object includes 
 the rotated loadings matrix with the lowest criterion value \code{f} 
 among attemnpted starts.Technically, this does not have to be the global 
 minimum. The \code{randomStarts} argument
 is not available \code{GPForth} and \code{GPFoblq}. However, for 
 \code{GPForth} and \code{GPFoblq} a single random
 initial rotation matrix may be set by \code{Tmat = \link{Random.Start}(ncol(A))}.
  
 The argument \code{method} can be used to specify a string indicating
 the rotation objective. Oblique rotation defaults to \code{"quartimin"} 
 and orthogonal rotation defaults to \code{"varimax"}. Available rotation objectives
 are \code{"oblimin"}, \code{"quartimin"}, \code{"target"}, \code{"pst"},
 \code{"oblimax"}, \code{"entropy"}, \code{"quartimax"}, \code{"varimax"}, 
 \code{"simplimax"}, \code{"bentler"}, \code{"tandemI"}, \code{"tandemII"}, 
 \code{ "geomin"}, \code{"cf"}, \code{"infomax"}, \code{"mccammon"}, \code{bifactor},
 and \code{"varimin"}. The string is prefixed with "vgQ." to give the actual function call. 
 See \code{\link{vgQ}} for details.
 
 Some rotation criteria (\code{"oblimin"}, \code{"target"}, \code{"pst"}, \code{"simplimax"},
 \code{"geomin"}, \code{"cf"}) require one or more additional arguments. See \code{link{rotations}}
 for details and default values, if applicable. 

 For examples of the indirect use see \code{\link{rotations}}. 
 
 The argument \code{normalize} gives an indication of if and how any
 normalization should be done before rotation, and then undone after rotation.
 If \code{normalize} is \code{FALSE} (the default) no normalization is done.
 If \code{normalize} is \code{TRUE} then Kaiser normalization is done. (So
 squared row entries of normalized \code{A} sum to 1.0. This is sometimes called
 Horst normalization.)
 If \code{normalize} is a vector of length equal to the number of indicators (=
 number of rows of \code{A}) then the colums are divided by \code{normalize} 
 before rotation and multiplied by \code{normalize} after rotation. 
 If \code{normalize} is a function then it should take \code{A} as an argument
 and return a vector which is used like the vector above. See Nguyen and Waller (2022)
 for detailed investigation of normalization on factor rotations, including 
 potential effect on qualitative interpretation of loadings.
}
\seealso{
   \code{\link{Random.Start}},
   \code{\link[stats]{factanal}},
   \code{\link{oblimin}},
   \code{\link{quartimin}},
   \code{\link{targetT}},
   \code{\link{targetQ}},
   \code{\link{pstT}},
   \code{\link{pstQ}},
   \code{\link{oblimax}},
   \code{\link{entropy}},
   \code{\link{quartimax}},
   \code{\link{Varimax}},
   \code{\link[stats]{varimax}},
   \code{\link{simplimax}},
   \code{\link{bentlerT}},
   \code{\link{bentlerQ}},
   \code{\link{tandemI}},
   \code{\link{tandemII}},
   \code{\link{geominT}},
   \code{\link{geominQ}},
   \code{\link{bigeominT}},
   \code{\link{bigeominQ}},
   \code{\link{cfT}},
   \code{\link{cfQ}},
   \code{\link{equamax}},
   \code{\link{parsimax}},
   \code{\link{infomaxT}},
   \code{\link{infomaxQ}},
   \code{\link{mccammon}},
   \code{\link{varimin}},
   \code{\link{bifactorT}},
   \code{\link{bifactorQ}},
   \code{\link[stats]{promax}}
}
\examples{
  # see rotations for more examples	
	
  data(Harman, package = "GPArotation")
  GPFRSorth(Harman8, method = "quartimax")
  quartimax(Harman8)
  GPFRSoblq(Harman8, method = "quartimin", normalize = TRUE)
  loadings( quartimin(Harman8, normalize = TRUE) )

  # using random starts
  data("WansbeekMeijer", package = "GPArotation")
  fa.unrotated  <- factanal(factors = 3, covmat=NetherlandsTV, normalize=TRUE, rotation="none")
  GPFRSoblq(loadings(fa.unrotated), normalize = TRUE, method = "oblimin", randomStarts = 100)
  oblimin(loadings(fa.unrotated), randomStarts=100)
  data(Thurstone, package = "GPArotation")
  geominQ(box26, normalize = TRUE, randomStarts=100)
  
  # displaying results of factor analysis rotation output
  origdigits <- options("digits")
  Abor.unrotated <- factanal(factors = 2, covmat = ability.cov, rotation = "none")
  Abor <- oblimin(loadings(Abor.unrotated), randomStarts = 20)
  Abor
  print(Abor)
  print(Abor, sortLoadings=FALSE) #this matches the output passed to factanal
  print(Abor, Table=TRUE)
  print(Abor, rotateMat=TRUE)
  print(Abor, digits=2)
  # by default provides the structure matrix for oblique rotation
  summary(Abor)
  summary(Abor, Structure=FALSE)
  options(digits = origdigits$digits)
  
  # GPArotation output does sort loadings, but use print to obtain if needed
  set.seed(334)
  xusl <- quartimin(Harman8, normalize = TRUE, randomStarts=100)
  # loadings without ordering (default)
  loadings(xusl)
  max(abs(print(xusl)$loadings - xusl$loadings)) == 0 # FALSE
  # output sorted loadings via print (not default)
  xsl <- print(xusl)
  max(abs(print(xsl)$loadings - xsl$loadings)) == 0 # TRUE
  
  # Kaiser normalization is used when normalize=TRUE
  factanal(factors = 2, covmat = ability.cov, rotation = "oblimin", 
  			control=list(rotate=list(normalize = TRUE)))
  # Cureton-Mulaik normalization can be done by passing values to the rotation
  # may result in convergence problems
  NormalizingWeightCM <- function (L) {
    Dk <- diag(sqrt(diag(L \%*\% t(L)))^-1) \%*\% L
    wghts <- rep(0, nrow(L))
    fpls <- Dk[, 1]
    acosi <- acos(ncol(L)^(-1/2))
    for (i in 1:nrow(L)) {
    	num <- (acosi - acos(abs(fpls[i])))
        dem <- (acosi - (function(a, m) ifelse(abs(a) < (m^(-1/2)), pi/2, 0))(fpls[i], ncol(L)))
        wghts[i] <- cos(num/dem * pi/2)^2 + 0.001
    	}
    Dv <- wghts * sqrt(diag(L \%*\% t(L)))^-1        
    Dv
  }
  quartimin(Harman8, normalize = NormalizingWeightCM(Harman8), randomStarts=100)
  quartimin(Harman8, normalize = TRUE, randomStarts=100)
	
  }
\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert}
\references{
   Bernaards, C.A. and Jennrich, R.I. (2005) Gradient Projection Algorithms 
     and Software for Arbitrary Rotation Criteria in Factor Analysis. 
     \emph{Educational and Psychological Measurement}, \bold{65}, 676--696.

   Jennrich, R.I. (2001). A simple general procedure for orthogonal rotation.
       \emph{Psychometrika}, \bold{66}, 289--306.

   Jennrich, R.I. (2002). A simple general method for oblique rotation. 
       \emph{Psychometrika}, \bold{67}, 7--19.
       
   Nguyen, H.V. and Waller, N.G. (2022). Local minima and factor rotations
   in exploratory factor analysis. \emph{ Psychological Methods}. 
   Advance online publication. https://doi.org/10.1037/met0000467

   }
\concept{rotation}
\keyword{multivariate}
