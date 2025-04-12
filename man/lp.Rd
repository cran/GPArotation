\encoding{UTF-8}
\name{lp}
\alias{Lp rotation}
\alias{GPForth.lp}
\alias{GPFoblq.lp}

\title{\eqn{L^p}{Lp} Rotation}
\usage{
     GPForth.lp(A, Tmat=diag(rep(1, ncol(A))), p=1, normalize=FALSE, eps=1e-05, 
            maxit=10000, gpaiter=5) 
     GPFoblq.lp(A, Tmat=diag(rep(1, ncol(A))), p=1, normalize=FALSE, eps=1e-05, 
            maxit=10000, gpaiter=5)  
}

\arguments{
    \item{A}{Initial factor loadings matrix to be rotated.}
    \item{Tmat}{Initial rotation matrix.}
    \item{p}{Component-wise \eqn{L^p}{Lp} where 0 < p \eqn{=<}{=<} 1.}
    \item{normalize}{Not recommended for \eqn{L^p}{Lp} rotation.}
    \item{eps}{Convergence is assumed when the norm of the gradient is smaller than eps.}
    \item{maxit}{Maximum number of iterations allowed in the main loop.}
    \item{gpaiter}{Maximum iterations for GPA rotation. 
    The goal is to decrease the objective value, not optimize the inner loop. 
    Warnings may appear, but they can be ignored if the main loop converges.}
}

\value{
 A \code{GPArotation} object, which is a list containing:
  \item{loadings}{Rotated loadings matrix, with one column per factor. If \code{randomStarts} 
  were used, this contains the loadings with the lowest criterion value.}
  \item{Th}{Rotation matrix, satisfying \code{loadings \%*\% t(Th) = A}.}
  \item{Table}{Matrix recording iteration details during optimization.}
  \item{method}{String indicating the rotation objective function.}
  \item{orthogonal}{Logical indicating whether the rotation is orthogonal.}
  \item{convergence}{Logical indicating whether convergence was achieved. 
  Convergence is controlled element-wise by tolerance.}
  \item{Phi}{Covariance matrix of rotated factors, \code{t(Th) \%*\% Th}.}
}

\description{
   Performs \eqn{L^p}{Lp} rotation to obtain sparse loadings. 
}

\details{
These functions optimize an \eqn{L^p}{L^p} rotation objective, where \code{0 < p =< 1}. A smaller \code{p} promotes sparsity in the loading matrix but increases computational difficulty. For guidance on choosing \code{p}, see the Concluding Remarks in the references.

Since the \eqn{L^p}{Lp} function is nonsmooth, a different optimization method is required compared to smooth rotation criteria. Two new functions, \code{GPForth.lp} and \code{GPFoblq.lp}, replace \code{GPForth} and \code{GPFoblq} for orthogonal and oblique \eqn{L^p}{Lp} rotations, respectively.

The optimization method follows an iterative reweighted least squares (IRLS) approach. It approximates the nonsmooth objective function with a smooth weighted least squares function in the main loop and optimizes it using GPA in the inner loop.

Normalization is not recommended for \eqn{L^p}{Lp} rotation. Its use may have unexpected results.

}


\examples{
  data("WansbeekMeijer", package = "GPArotation")
  fa.unrotated <- factanal(factors = 2, covmat = NetherlandsTV, rotation = "none")
  
  options(warn = -1)
  
  # Orthogonal rotation
  # Single start from random position
  fa.lpT1 <- GPForth.lp(loadings(fa.unrotated), p = 1)
  # 10 random starts
  fa.lpT <- lpT(loadings(fa.unrotated), Tmat=Random.Start(2), p = 1, randomStarts = 10)
  print(fa.lpT, digits = 5, sortLoadings = FALSE, Table = TRUE, rotateMat = TRUE)
  
  p <- 1
  # Oblique rotation
  # Single start
  fa.lpQ1 <- GPFoblq.lp(loadings(fa.unrotated), p = p)
  # 10 random starts
  fa.lpQ <- lpQ(loadings(fa.unrotated), p = p, randomStarts = 10)
  summary(fa.lpQ, Structure = TRUE)


  # this functions ensures consistent ordering of factors of a
  # GPArotation object for cleaner comparison
  # Inspired by fungible::orderFactors and fungible::faSort functions
  sortFac <- function(x){
    # Only works for object of class GPArotation 
    if (!inherits(x, "GPArotation")) {stop("argument not of class GPArotation")}
    cln <- colnames(x$loadings) 
    # ordering for oblique slightly different from orthogonal
    ifelse(x$orthogonal, vx <- colSums(x$loadings^2), 
      vx <- diag(x$Phi \%*\% t(x$loadings) \%*\% x$loadings) )
    # sort by squared loadings from high to low
    vxo <- order(vx, decreasing = TRUE)
    vx <- vx[vxo]
    # maintain the right sign
    Dsgn <- diag(sign(colSums(x$loadings^3))) [ , vxo]
    x$Th <- x$Th \%*\% Dsgn 
    x$loadings <- x$loadings \%*\% Dsgn
    if (match("Phi", names(x))) { 
            # If factor is negative, reverse corresponding factor correlations
            x$Phi <- t(Dsgn) \%*\% x$Phi \%*\% Dsgn
      }
    colnames(x$loadings) <- cln 
    x 
  } 
  
  # seed set to see the results of sorting
  set.seed(1020)
  fa.lpQ1 <- lpQ(loadings(fa.unrotated),p=1,randomStarts=10)
  fa.lpQ0.5 <- lpQ(loadings(fa.unrotated),p=0.5,randomStarts=10)
  fa.geo <- geominQ(loadings(fa.unrotated), randomStarts=10)
  
  # with ordered factor loadings
  res <- round(cbind(sortFac(fa.lpQ1)$loadings, sortFac(fa.lpQ0.5)$loadings, 
    sortFac(fa.geo)$loadings),3)
  print(c("oblique-  Lp 1           Lp 0.5         geomin")); print(res)

  # without ordered factor loadings
  res <- round(cbind(fa.lpQ1$loadings, fa.lpQ0.5$loadings, fa.geo$loadings),3)
  print(c("oblique-  Lp 1           Lp 0.5         geomin")); print(res)

  options(warn = 0)
  
}

\seealso{
  \code{\link{lpT}},
  \code{\link{lpQ}},
  \code{\link{vgQ.lp.wls}}
}

\references{
  Liu, X., Wallin, G., Chen, Y., & Moustaki, I. (2023). Rotation to sparse 
  loadings using \eqn{L^p}{Lp} losses and related inference problems. 
  \emph{Psychometrika}, \bold{88}(2), 527--553.
}

\author{Xinyi Liu, with minor modifications for GPArotation by C. Bernaards}

\concept{rotation}
\keyword{multivariate}