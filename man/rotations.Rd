\name{rotations}
\alias{rotations}
\alias{oblimin}
\alias{quartimin}
\alias{targetT}
\alias{targetQ}
\alias{pstT}
\alias{pstQ}
\alias{oblimax}
\alias{entropy}
\alias{quartimax}
\alias{Varimax}
\alias{simplimax}
\alias{bentlerT}
\alias{bentlerQ}
\alias{tandemI}
\alias{tandemII}
\alias{geominT}
\alias{geominQ}
\alias{bigeominT}
\alias{bigeominQ}
\alias{cfT}
\alias{cfQ}
\alias{equamax}
\alias{parsimax}
\alias{infomaxT}
\alias{infomaxQ}
\alias{mccammon}
\alias{varimin}
\alias{bifactorT}
\alias{bifactorQ}

\title{Rotations}
\usage{
    oblimin(A, Tmat=diag(ncol(A)), gam=0, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    quartimin(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    targetT(A, Tmat=diag(ncol(A)), Target=NULL, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0, L=NULL)
    targetQ(A, Tmat=diag(ncol(A)), Target=NULL, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0, L=NULL)
    pstT(A, Tmat=diag(ncol(A)), W=NULL, Target=NULL, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0, L=NULL)
    pstQ(A, Tmat=diag(ncol(A)), W=NULL, Target=NULL, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0, L=NULL)
    oblimax(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    entropy(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    quartimax(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5,maxit=1000,randomStarts=0)
    Varimax(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    simplimax(A, Tmat=diag(ncol(A)), k=nrow(A), normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    bentlerT(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    bentlerQ(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    tandemI(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    tandemII(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    geominT(A, Tmat=diag(ncol(A)), delta=.01, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    geominQ(A, Tmat=diag(ncol(A)), delta=.01, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    bigeominT(A, Tmat=diag(ncol(A)), delta=.01, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    bigeominQ(A, Tmat=diag(ncol(A)), delta=.01, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    cfT(A, Tmat=diag(ncol(A)), kappa=0, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    cfQ(A, Tmat=diag(ncol(A)), kappa=0, normalize=FALSE, eps=1e-5, 
    		maxit=1000, randomStarts=0)
    equamax(A, Tmat=diag(ncol(A)), kappa=ncol(A)/(2*nrow(A)), normalize=FALSE,
    		eps=1e-5, maxit=1000, randomStarts = 0)
    parsimax(A, Tmat=diag(ncol(A)), kappa=(ncol(A)-1)/(ncol(A)+nrow(A)-2), 
    		normalize=FALSE, eps=1e-5, maxit=1000, randomStarts = 0)
    infomaxT(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    infomaxQ(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    mccammon(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    varimin(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, randomStarts=0)
    bifactorT(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000,randomStarts=0)
    bifactorQ(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000,randomStarts=0)
    }
\arguments{
    \item{A}{an initial loadings matrix to be rotated.}
    \item{Tmat}{initial rotation matrix.}
    \item{gam}{0=Quartimin, .5=Biquartimin, 1=Covarimin.}
    \item{Target}{rotation target for objective calculation.}
    \item{W}{weighting of each element in target.}
    \item{k}{number of close to zero loadings.}
    \item{delta}{constant added to Lambda^2 in objective calculation.}
    \item{kappa}{see details.}
    \item{normalize}{parameter passed to optimization routine (GPForth or GPFoblq).}
    \item{eps}{parameter passed to optimization routine (GPForth or GPFoblq).}
    \item{maxit}{parameter passed to optimization routine (GPForth or GPFoblq).}
    \item{randomStarts}{parameter passed to optimization routine (GPFRSorth or GPFRSoblq).}
    \item{L}{provided for backward compatibility in target rotations only. Use A going forward.}
}

\value{A list (which includes elements used by \code{factanal}) with:
   \item{loadings}{Lh from \code{GPFRSorth} or \code{GPFRSoblq}.} 
   \item{Th}{Th from \code{GPFRSorth} or \code{GPFRSoblq}.}  
   \item{Table}{Table from \code{GPForth} or \code{GPFoblq}.}
   \item{method}{A string indicating the rotation objective function.}
   \item{orthogonal}{A logical indicating if the rotation is orthogonal.}
   \item{convergence}{Convergence indicator from \code{GPFRSorth} or \code{GPFRSoblq}.}
   \item{Phi}{t(Th) \%*\% Th.  The covariance matrix of the rotated factors.
   This will be the identity matrix for orthogonal
   rotations so is omitted (NULL) for the result from GPFRSorth and GPForth.}
   \item{randStartChar}{Vector indicating results from random starts 
   from \code{GPFRSorth} or \code{GPFRSoblq}}
   }
\description{
   Optimize factor loading rotation objective.
}
\details{
 These functions optimize a rotation objective. They can be used directly or the
 function name can be passed to factor analysis functions like \code{factanal}.
 Several of the function names end in T or Q, which indicates if they are
 orthogonal or oblique rotations (using \code{GPFRSorth} or \code{GPFRSoblq}
 respectively).  
    
 Rotations which are available are 
\tabular{lll}{
 \code{oblimin}   \tab oblique    \tab  oblimin family 			\cr
 \code{quartimin} \tab oblique	   \tab         				\cr
 \code{targetT}   \tab orthogonal \tab  target rotation			\cr
 \code{targetQ}   \tab oblique    \tab  target rotation			\cr
 \code{pstT}	   \tab orthogonal \tab  partially specified target rotation	\cr
 \code{pstQ}	   \tab oblique    \tab  partially specified target rotation	\cr
 \code{oblimax}   \tab oblique	   \tab         				\cr
 \code{entropy}   \tab orthogonal \tab  minimum entropy			\cr
 \code{quartimax} \tab orthogonal \tab 					\cr
 \code{varimax}   \tab orthogonal \tab 					\cr
 \code{simplimax} \tab oblique	   \tab         				\cr
 \code{bentlerT}  \tab orthogonal \tab  Bentler's invariant pattern simplicity criterion\cr
 \code{bentlerQ}  \tab oblique    \tab  Bentler's invariant pattern simplicity criterion\cr
 \code{tandemI}   \tab orthogonal \tab  Tandem principle I criterion		   \cr
 \code{tandemII}  \tab orthogonal \tab  Tandem principle II criterion		   \cr
 \code{geominT}   \tab orthogonal \tab 				   \cr
 \code{geominQ}   \tab oblique	   \tab         			   \cr
 \code{bigeominT} \tab orthogonal \tab 				   \cr
 \code{bigeominQ} \tab oblique	   \tab         			   \cr
 \code{cfT}	   \tab orthogonal \tab  Crawford-Ferguson family	   \cr
 \code{cfQ}	   \tab oblique	   \tab  Crawford-Ferguson family	   \cr
 \code{equamax}   \tab orthogonal \tab  Crawford-Ferguson family	   \cr
 \code{parsimax}  \tab orthogonal \tab  Crawford-Ferguson family	   \cr
 \code{infomaxT}  \tab orthogonal \tab 				   \cr
 \code{infomaxQ}  \tab oblique	   \tab         			   \cr
 \code{mccammon}  \tab orthogonal \tab  McCammon minimum entropy ratio    \cr
 \code{varimin}   \tab orthogonal \tab 				   \cr
 \code{bifactorT} \tab orthogonal \tab  Jennrich and Bentler bifactor rotation\cr
 \code{bifactorQ} \tab oblique    \tab  Jennrich and Bentler biquartimin rotation\cr
 }

 Note that \code{Varimax} defined here uses \code{vgQ.varimax} and 
 is not \code{varimax}
 defined in the \code{stats} package. \code{stats:::varimax} does Kaiser 
 normalization by default whereas \code{Varimax} defined here does not.
 
 The argument \code{kappa} parameterizes the family for the Crawford-Ferguson 
 method. If \code{m} is the number of factors and \code{p} is the number of 
 indicators then \code{kappa} values having special names are 0=Quartimax, 
 1/p=Varimax, m/(2*p)=Equamax, (m-1)/(p+m-2)=Parsimax, 1=Factor parsimony.

 }
\examples{
  # see GPFRSorth and GPFRSoblq for more examples
  
  # getting loadings matrices
  data("Harman", package="GPArotation")
  qHarman  <- GPFRSorth(Harman8, Tmat=diag(2), method="quartimax")
  qHarman <- quartimax(Harman8) 
  loadings(qHarman) - qHarman$loadings   #2 ways to get the loadings

  # factanal loadings used in GPArotation
  data("WansbeekMeijer", package="GPArotation")
  fa.unrotated  <- factanal(factors = 2, covmat=NetherlandsTV, normalize=TRUE, rotation="none")
  quartimax(loadings(fa.unrotated), normalize=TRUE)
  geominQ(loadings(fa.unrotated), normalize=TRUE, randomStarts=100)

  # passing arguments to factanal (See vignette for a caution)
  # vignette("GPAguide", package = "GPArotation")
  data(ability.cov)
  factanal(factors = 2, covmat = ability.cov, rotation="infomaxT")
  factanal(factors = 2, covmat = ability.cov, rotation="infomaxT", 
    control=list(rotate=list(normalize = TRUE, eps = 1e-6)))
  # when using factanal for oblique rotation it is best to use the rotation command directly
  # instead of including it in the factanal command (see Vignette).  
  fa.unrotated  <- factanal(factors = 3, covmat=NetherlandsTV, normalize=TRUE, rotation="none")
  quartimin(loadings(fa.unrotated), normalize=TRUE)

  # oblique target rotation of 2 varimax rotated matrices towards each other
  # See vignette for additional context and computation,
  trBritain <- matrix( c(.783,-.163,.811,.202,.724,.209,.850,.064,
    -.031,.592,-.028,.723,.388,.434,.141,.808,.215,.709), byrow=TRUE, ncol=2)
  trGermany <- matrix( c(.778,-.066, .875,.081, .751,.079, .739,.092,
    .195,.574, -.030,.807, -.135,.717, .125,.738, .060,.691), byrow=TRUE, ncol = 2)
  trx <- targetQ(trGermany, Target = trBritain)
  # Difference between rotated loadings matrix and target matrix 
  y <- trx$loadings - trBritain
  
  # partially specified target; See vignette for additional method
  A <- matrix(c(.664, .688, .492, .837, .705, .82, .661, .457, .765, .322, 
    .248, .304, -0.291, -0.314, -0.377, .397, .294, .428, -0.075,.192,.224,
    .037, .155,-.104,.077,-.488,.009), ncol=3)  
  SPA <- matrix(c(rep(NA, 6), .7,.0,.7, rep(0,3), rep(NA, 7), 0,0, NA, 0, rep(NA, 4)), ncol=3)
  targetT(A, Target=SPA)

  # using random starts
  data("WansbeekMeijer", package="GPArotation")
  fa.unrotated  <- factanal(factors = 3, covmat=NetherlandsTV, normalize=TRUE, rotation="none")
  # single rotation with a random start
  oblimin(loadings(fa.unrotated), Tmat=Random.Start(3))
  oblimin(loadings(fa.unrotated), randomStarts=1)
  # multiple random starts
  oblimin(loadings(fa.unrotated), randomStarts=100)

  # assessing local minima for box26 data
  data(Thurstone, package = "GPArotation")
  infomaxQ(box26, normalize = TRUE, randomStarts = 150)
  geominQ(box26, normalize = TRUE, randomStarts = 150)
  # for detailed investigation of local minima, consult package 'fungible' 
  # library(fungible)
  # faMain(urLoadings=box26, rotate="geominQ", rotateControl=list(numberStarts=150))
  # library(psych) # package 'psych' with random starts:
  # faRotations(box26, rotate = "geominQ", hyper = 0.15, n.rotations = 150)

  }
\seealso{
   \code{\link{GPFRSorth}}, 
   \code{\link{GPFRSoblq}}, 
   \code{\link{vgQ}},
   \code{\link{eiv}}, 
   \code{\link{echelon}}, 
   \code{\link{WansbeekMeijer}}, 
   \code{\link[stats]{factanal}}, 
   \code{\link[stats]{varimax}}
}
\references{
  Bernaards, C.A. and Jennrich, R.I. (2005) Gradient Projection Algorithms 
     and Software for Arbitrary Rotation Criteria in Factor Analysis. 
     \emph{Educational and Psychological Measurement}, \bold{65}, 676--696.

   Bifactor rotation, bifactorT and bifactorQ are called bifactor and
      biquartimin in Jennrich, R.I. and Bentler, P.M. 
      (2011) Exploratory bi-factor analysis. \emph{Psychometrika}, \bold{76}. 

  }

\author{Coen A. Bernaards and Robert I. Jennrich
        with some R modifications by Paul Gilbert.}
\concept{rotation}
\keyword{multivariate}


