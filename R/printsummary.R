# summaries for GPArotation
# S3 class functions
# pint.GPArotation
# summary.GPArotation
# print.summary.GPArotation

print.GPArotation <- function (x, digits=3L, sortLoadings=TRUE, rotateMat=FALSE, Table=FALSE, ...){
   cln <- colnames(x$loadings) # Based on Fungible faSort and orderFactors 
   ifelse(x$orthogonal, vx <- colSums(x$loadings^2), vx <- diag(x$Phi %*% t(x$loadings) %*% x$loadings))
   if (sortLoadings){
     vxo <- order(vx, decreasing = TRUE)
     vx <- vx[vxo]
     Dsgn <- diag(sign(colSums(x$loadings^3))) [ , vxo]
     x$Th <- x$Th %*% Dsgn
     x$loadings <- x$loadings %*% Dsgn
     if ("Phi" %in% names(x)) {
       x$Phi <- t(Dsgn) %*% x$Phi %*% Dsgn
     }
     colnames(x$loadings) <- cln 
   }
   cat(if(x$orthogonal)"Orthogonal" else "Oblique")
   cat(" rotation method", x$method)
   cat((if(!x$convergence)" NOT" ), "converged") 
   cat(if ("randStartChar" %in% names(x)) " at lowest minimum.\n" else ".\n")
   
   if ("randStartChar" %in% names(x)){
     cat("Of ",x$randStartChar[1]," random starts ",round(100 * x$randStartChar[2] / x$randStartChar[1]),"% converged, ", sep="")
     cat(round(100 * x$randStartChar[3] / x$randStartChar[1]),"% at the same lowest minimum.\n", sep="")
     if (x$randStartChar[4] > 1){
         cat("Random starts converged to ",x$randStartChar[4], " different local minima\n", sep="")
       }
	}
	
   cat(if ("randStartChar" %in% names(x)) "Loadings at lowest minimum:\n" else "Loadings:\n")
   print(x$loadings, digits=digits)

    varex <- rbind(`SS loadings` = vx)
    colnames(varex) <- cln
    if (is.null(attr(x, "covariance"))) {
        varex <- rbind(varex, `Proportion Var` = vx/nrow(x$loadings))
        varex <- rbind(varex, `Cumulative Var` = cumsum(vx/nrow(x$loadings)))
    }
    cat("\n")
    print(round(varex, digits))

   if(!x$orthogonal){
   	 dimnames(x$Phi) <- list(cln, cln)
     cat("\nPhi:\n")
     print(x$Phi, digits=digits)
     }

   if(rotateMat){
     cat("\nRotating matrix:\n")
     print(t(solve(x$Th)), digits=digits)
     }
     
   if(Table){
     cat("\nIteration table:\n")
     print(x$Table, digits=digits)
     }
   invisible(x)
   }	

summary.GPArotation <- function (object, digits=3L, Structure=TRUE, ...){ 
   r <- list(loadings=object$loadings, Phi=object$Phi, method=object$method, 
   		orthogonal=object$orthogonal, convergence=object$convergence, 
   		iters= rev(object$Table[, 1])[1], Structure=Structure, digits=digits) 
   class(r) <- "summary.GPArotation"
   r
   }

print.summary.GPArotation <- function (x, ...){
   cat(ifelse(x$orthogonal, "Orthogonal", "Oblique"))
   cat(" rotation method", x$method)
   if(!x$convergence) cat(" NOT" )
   cat(" converged in ", x$iters, " iterations.\n", sep="")

   prstr <- !x$orthogonal && x$Structure
   cat(ifelse(prstr, "Pattern (loadings):\n", "Loadings:\n"))
   print(x$loadings, digits=x$digits)
   if (prstr){
     cat("\nStructure:\n")
     print(x$loadings %*% x$Phi, digits=x$digits)
     }
  }
