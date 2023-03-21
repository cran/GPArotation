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
  of length one based on the QR decompostion.
  }
\seealso{
   \code{\link{GPFRSorth}}, 
   \code{\link{GPFRSoblq}},
   \code{\link{GPForth}}, 
   \code{\link{GPFoblq}},
   \code{\link{oblimin}}
   }
\examples{
	Random.Start <- function(k,orthogonal=TRUE){
	   # routine for generating orthogonal or oblique random matrix
	   mat <- matrix(rnorm(k*k),k)
	   if (orthogonal){
 	      ans <- qr.Q(qr(mat))
 	   }
 	   else{
	      ans <- mat \%*\% diag(1/sqrt(diag(crossprod(mat))))
 	   }
	   ans
	}
	
	data("Thurstone", package="GPArotation")
	simplimax(box26,Tmat = Random.Start(3, TRUE))
	simplimax(box26,Tmat = Random.Start(3, FALSE))

	# covariance matrix is Phi = t(Th) \%*\% Th
	rms <- Random.Start(3, FALSE)
	t(rms) \%*\% rms # covariance matrix because oblique rms
	rms <- Random.Start(3, TRUE)
	t(rms) \%*\% rms # identity matrix because orthogonal rms
	
   }
\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert
	}
\concept{rotation}
\keyword{multivariate}


