\name{Random.Start}
\alias{Random.Start}

\title{Generate a Random Orthogonal Rotation}
\usage{
    Random.Start(k)
    }
\arguments{
    \item{k}{An integer indicating the dimension of the square matrix.}
    }

\description{
   Random orthogonal rotation to use as Tmat matrix to start GPFRSorth, GPFRSoblq, GPForth, or GPFoblq.
   }
\value{An orthogonal matrix.}

\details{
  The random start function produces an orthogonal matrix with columns
  of length one based on the QR decompostion. This randomization procedures
  follows the logic of Stewart(1980) and Mezzari(2007), as of GPArotation version 2024.2-1.
  }
\seealso{
   \code{\link{GPFRSorth}}, 
   \code{\link{GPFRSoblq}},
   \code{\link{GPForth}}, 
   \code{\link{GPFoblq}},
   \code{\link{rotations}}
   }
\examples{
  # Generate a random ortogonal matrix of dimension 5 x 5
  Random.Start(5)
  
  # function for generating orthogonal or oblique random matrix
  Random.Start <- function(k = 2L,orthogonal=TRUE){
    mat <- matrix(rnorm(k*k),k)
    if (orthogonal){
      qr.out <- qr(matrix(rnorm(k * k), nrow = k, ncol = k))
      Q <- qr.Q(qr.out)
      R <- qr.R(qr.out)
      R.diag <- diag(R)
      R.diag2 <- R.diag/abs(R.diag)
      ans <- t(t(Q) * R.diag2)
      ans
      }
    else{
	  ans <- mat \%*\% diag(1/sqrt(diag(crossprod(mat))))
 	  }
    ans
    }
    	
  data("Thurstone", package="GPArotation")
  simplimax(box26,Tmat = Random.Start(3, orthogonal = TRUE))
  simplimax(box26,Tmat = Random.Start(3, orthogonal = FALSE))

  # covariance matrix is Phi = t(Th) \%*\% Th
  rms <- Random.Start(3, FALSE)
  t(rms) \%*\% rms # covariance matrix because oblique rms
  rms <- Random.Start(3, TRUE)
  t(rms) \%*\% rms # identity matrix because orthogonal rms
	
   }
   
\references{
Stewart, G. W. (1980). The Efficient Generation of Random Orthogonal Matrices 
with an Application to Condition Estimators. 
\emph{SIAM Journal on Numerical Analysis}, 
\bold{17}(3), 403--409. http://www.jstor.org/stable/2156882
 
Mezzadri, F. (2007). How to generate random matrices from the classical
compact groups. \emph{Notices of the American Mathematical Society}, \bold{54}(5), 592--604.
https://arxiv.org/abs/math-ph/0609050
}
\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert. Additional input from Yves Rosseel.
	}
\concept{rotation}
\keyword{multivariate}


