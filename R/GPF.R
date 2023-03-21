# legacy functions that contain the actual GP algorithms
# these function shall not be changed
# functions have not changes since 2008

GPForth <- function(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
		    method="varimax", methodArgs=NULL){
 if((!is.logical(normalize)) || normalize) {
     W <- NormalizingWeight(A, normalize=normalize)
     normalize <- TRUE
     A <- A/W
     }
 if(1 >= ncol(A)) stop("rotation does not make sense for single factor models.")
 al <- 1
 L <- A %*% Tmat
 #Method <- get(paste("vgQ",method,sep="."))
 #VgQ <- Method(L, ...)
 Method <- paste("vgQ",method,sep=".")
 VgQ <- do.call(Method, append(list(L), methodArgs))
 G <- crossprod(A,VgQ$Gq)
 f <- VgQ$f
 Table <- NULL
 #set initial value for the unusual case of an exact initial solution 
 VgQt <- do.call(Method, append(list(L), methodArgs))   
 for (iter in 0:maxit){
   M <- crossprod(Tmat,G)
   S <- (M + t(M))/2
   Gp <- G - Tmat %*% S
   s <- sqrt(sum(diag(crossprod(Gp))))
   Table <- rbind(Table, c(iter, f, log10(s), al))
   if (s < eps)  break
   al <- 2*al
   for (i in 0:10){
     X <- Tmat - al * Gp
     UDV <- svd(X)
     Tmatt <- UDV$u %*% t(UDV$v)
     L <- A %*% Tmatt
     #VgQt <- Method(L, ...)
     VgQt <- do.call(Method, append(list(L), methodArgs))
     if (VgQt$f < (f - 0.5*s^2*al)) break
     al <- al/2
     }
   Tmat <- Tmatt
   f <- VgQt$f
   G <- crossprod(A,VgQt$Gq)
   }
 convergence <- (s < eps)
 if ((iter == maxit) & !convergence)
     warning("convergence not obtained in GPForth. ", maxit, " iterations used.")
 if(normalize) L <- L * W
 dimnames(L) <- dimnames(A)
 r <- list(loadings=L, Th=Tmat, Table=Table, 
        method=VgQ$Method, orthogonal=TRUE, convergence=convergence, Gq=VgQt$Gq)
 class(r) <- "GPArotation"
 r
}

GPFoblq <- function(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
		    method="quartimin", methodArgs=NULL){
 if(1 >= ncol(A)) stop("rotation does not make sense for single factor models.")
 if((!is.logical(normalize)) || normalize) {
     W <- NormalizingWeight(A, normalize=normalize)
     normalize <- TRUE
     A <- A/W
     }
 al <- 1
 L <- A %*% t(solve(Tmat))
 #Method <- get(paste("vgQ",method,sep="."))
 #VgQ <- Method(L, ...)
 Method <- paste("vgQ",method,sep=".")
 VgQ <- do.call(Method, append(list(L), methodArgs))
 G <- -t(t(L) %*% VgQ$Gq %*% solve(Tmat))
 f <- VgQ$f
 Table <- NULL
 #Table <- c(-1,f,log10(sqrt(sum(diag(crossprod(G))))),al)
 #set initial value for the unusual case of an exact initial solution 
 VgQt <- do.call(Method, append(list(L), methodArgs))   
 for (iter in 0:maxit){
   Gp <- G - Tmat %*% diag(c(rep(1,nrow(G)) %*% (Tmat*G)))
   s <- sqrt(sum(diag(crossprod(Gp))))
   Table <- rbind(Table,c(iter,f,log10(s),al))
   if (s < eps) break
   al <- 2*al
   for (i in 0:10){
     X <- Tmat - al*Gp
     v <- 1/sqrt(c(rep(1,nrow(X)) %*% X^2))
     Tmatt <- X %*% diag(v)
     L <- A %*% t(solve(Tmatt))
     #VgQt <- Method(L, ...)
     VgQt <- do.call(Method, append(list(L), methodArgs))
     improvement <- f - VgQt$f 
     if (improvement >  0.5*s^2*al) break
     al <- al/2
     }
   Tmat <- Tmatt
   f <- VgQt$f
   G <- -t(t(L) %*% VgQt$Gq %*% solve(Tmatt))
   }
 convergence <- (s < eps)
 if ((iter == maxit) & !convergence)
     warning("convergence not obtained in GPFoblq. ", maxit, " iterations used.")
 if(normalize) L <- L * W
 dimnames(L) <- dimnames(A)
 # N.B. renaming Lh to loadings in specificific rotations 
 #   uses fact that  Lh is first.
 r <- list(loadings=L, Phi=t(Tmat) %*% Tmat, Th=Tmat, Table=Table,
      method=VgQ$Method, orthogonal=FALSE, convergence=convergence, Gq=VgQt$Gq)
 class(r) <- "GPArotation"
 r
}
