\name{00.GPArotation}
\alias{GPArotation}
\alias{GPArotation-package}
\alias{00.GPArotation-package}
\alias{GPArotation.Intro}
\docType{package}
\title{Gradient Projection Algorithms for Factor Rotation}
\description{GPA Rotation for Factor Analysis

The GPArotation package contains functions for the rotation of factor loadings
matrices. The functions implement Gradient Projection (GP) algorithms for orthogonal 
and oblique rotation. Additionally, a number of rotation criteria are provided.
The GP algorithms minimize the rotation criterion function, and provide the 
corresponding rotation matrix. For oblique rotation, the covariance / correlation matrix
of the factors is also provided. The rotation criteria implemented in this package 
are described in Bernaards and Jennrich (2005). Theory of the GP algorithm is 
described in Jennrich (2001, 2002) publications.

Additionally 2 rotation methods are provided that do not rely on GP (eiv and echelon)


\tabular{ll}{
Package: \tab GPArotation\cr
Depends: \tab R (>= 2.0.0)\cr
License: \tab GPL Version 2.\cr
URL: \tab https://optimizer.r-forge.r-project.org/GPArotation_www/\cr
}

Index of functions:

Wrapper functions that include random starts option \cr
\tabular{ll}{
\code{\link{GPFRSorth}} \tab Orthogonal rotation with random starts \cr
\code{\link{GPFRSorth}} \tab Oblique rotation with random starts \cr
}

Gradient Projection Rotation Algorithms (code unchanged since 2008)\cr
\tabular{ll}{
\code{\link{GPForth}} \tab Orthogonal rotation function \cr
\code{\link{GPForth}} \tab Oblique rotation function \cr
}

Utility functions\cr
\tabular{ll}{
\code{\link{Random.Start}} \tab Generate random a starting matrix \cr
\code{\link{NormalizingWeight}} \tab Kaiser normalization (not exported from NAMESPACE) \cr
}

Rotations\cr
\tabular{ll}{
\code{\link{oblimin}} \tab Oblimin rotation  \cr
\code{\link{quartimin}} \tab Quartimin rotation \cr
\code{\link{targetT}} \tab Orthogonal Target rotation \cr
\code{\link{targetQ}} \tab Oblique Target rotation \cr
\code{\link{pstT}} \tab Orthogonal Partially Specified Target rotation \cr
\code{\link{pstQ}} \tab Oblique Partially Specified Target rotation \cr
\code{\link{oblimax}} \tab Oblimax rotation \cr
\code{\link{entropy}} \tab Minimum Entropy rotation \cr
\code{\link{quartimax}} \tab Quartimax rotation \cr
\code{\link{Varimax}} \tab Varimax rotation \cr
\code{\link{simplimax}} \tab Simplimax rotation \cr
\code{\link{bentlerT}} \tab Orthogonal Bentler's Invariant Pattern Simplicity rotation \cr
\code{\link{bentlerQ}} \tab Oblique Bentler's Invariant Pattern Simplicity rotation \cr
\code{\link{tandemI}} \tab The Tandem Criteria Principle I rotation \cr
\code{\link{tandemII}} \tab The Tandem Criteria Principle II rotation \cr
\code{\link{geominT}} \tab Orthogonal Geomin rotation \cr
\code{\link{geominQ}} \tab Oblique Geomin rotation \cr
\code{\link{cfT}} \tab Orthogonal Crawford-Ferguson Family rotation \cr
\code{\link{cfQ}} \tab Oblique Crawford-Ferguson Family rotation \cr
\code{\link{equamax}} \tab Equamax rotation \cr
\code{\link{parsimax}} \tab Parsimax rotation \cr
\code{\link{infomaxT}} \tab Orthogonal Infomax rotation \cr
\code{\link{infomaxQ}} \tab Oblique Infomax rotation \cr
\code{\link{mccammon}} \tab McCammon Minimum Entropy Ratio rotation \cr
\code{\link{varimin}} \tab Varimin rotation \cr
\code{\link{bifactorT}} \tab Orthogonal Bifactor rotation \cr
\code{\link{bifactorQ}} \tab Oblique Bifactor rotation \cr
}


vgQ routines to compute value and gradient of the criterion (not exported from NAMESPACE)\cr
\tabular{ll}{
\code{\link{vgQ.oblimin}} \tab Oblimin vgQ  \cr
\code{\link{vgQ.quartimin}}\tab Quartimin vgQ \cr
\code{\link{vgQ.target}}\tab Target vgQ \cr
\code{\link{vgQ.pst}}\tab Partially Specified Target vgQ \cr
\code{\link{vgQ.oblimax}} \tab Oblimax vgQ \cr
\code{\link{vgQ.entropy}} \tab Minimum Entropy vgQ \cr
\code{\link{vgQ.quartimax}}\tab Quartimax vgQ \cr
\code{\link{vgQ.varimax}}\tab Varimax vgQ \cr
\code{\link{vgQ.simplimax}}\tab Simplimax vgQ \cr
\code{\link{vgQ.bentler}}\tab Bentler's Invariant Pattern Simplicity vgQ \cr
\code{\link{vgQ.tandemI}}\tab The Tandem Criteria Principle I vgQ \cr
\code{\link{vgQ.tandemII}}\tab The Tandem Criteria Principle II vgQ \cr
\code{\link{vgQ.geomin}}\tab Geomin vgQ \cr
\code{\link{vgQ.cf}} \tab Crawford-Ferguson Family vgQ \cr
\code{\link{vgQ.infomax}} \tab Infomax vgQ \cr
\code{\link{vgQ.mccammon}}\tab McCammon Minimum Entropy Ratio vgQ \cr
\code{\link{vgQ.varimin}}\tab Varimin vgQ \cr
\code{\link{vgQ.bifactor}}\tab Bifactor vgQ \cr
}

Data sets included in the GPArotation package \cr
\tabular{ll}{
\code{\link{Harman}}\tab Initial factor loading matrix for Harman's 8 physical variables \cr
\code{\link{Thurstone}} \tab box20 and box26 initial factor loadings matrices. \cr
\code{\link{WansbeekMeijer}} \tab Netherlands TV viewership. \cr
}



}



\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert.

Code is modified from original source \file{splusfunctions.net} available at 
  \url{https://optimizer.r-forge.r-project.org/GPArotation_www/}.

}
\references{
   The software reference is

   Bernaards, C.A. and Jennrich, R.I. (2005) Gradient Projection Algorithms 
     and Software for Arbitrary Rotation Criteria in Factor Analysis. 
     \emph{Educational and Psychological Measurement}, \bold{65}, 676--696.

   Theory of gradient projection algorithms may be found in:

   Jennrich, R.I. (2001). A simple general procedure for orthogonal rotation.
       \emph{Psychometrika}, \bold{66}, 289--306.

   Jennrich, R.I. (2002). A simple general method for oblique rotation. 
       \emph{Psychometrika}, \bold{67}, 7--19.
}
\keyword{ package }

\seealso{
   \code{\link{GPFRSorth}},
   \code{\link{GPFRSoblq}},
   \code{\link{rotations}},
   \code{\link{vgQ}}
}
