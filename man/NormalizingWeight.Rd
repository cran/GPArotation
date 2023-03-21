\name{NormalizingWeight}
\alias{NormalizingWeight}

\title{Internal Utility for Normalizing Weights}
\usage{
    NormalizingWeight(A, normalize=FALSE)
    }
\arguments{
    \item{A}{A loading matrix.}
    \item{normalize}{An indication of if/how the matrix should be normalized.}
}

\description{
   See GPFRSoblq and GPFRSorth.
   }
\value{A matrix. This function is not exported in the NAMESPACE, and is only used by the GP rotation functions. See \link{GPFRSorth} for an example of its use.}

\keyword{internal}

