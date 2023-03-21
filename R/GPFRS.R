# Main functions that serve as a wrapper function to GPForth and GPFoblq
# these functions have the added functionality of random starts
# Functions added in version 2023-1.1

GPFRSorth <- function(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
       method="varimax", methodArgs=NULL, randomStarts = 0){
	if (0 < randomStarts){ Tmat <- Random.Start(ncol(A)) }
	r <- GPForth(A, Tmat = Tmat, normalize = normalize, eps = eps, maxit = maxit, 
				method = method, methodArgs = methodArgs)
	if (randomStarts > 1){
		Qvalues <- rev(r$Table[,2])[1]  
		Qmin <- Qvalues
		Qconverged <- r$convergence
		for (inum in 2:randomStarts){
			gpout <- GPForth(A, Tmat = Random.Start(ncol(A)), normalize = normalize, eps = eps, maxit = maxit, 
				method = method, methodArgs = methodArgs)
	  		Qvalues <- c(Qvalues, rev(gpout$Table[,2])[1]) 
	  		Qconverged <- c(Qconverged, gpout$convergence)
	  		if (rev(gpout$Table[,2])[1] < Qmin){ 
	  			r <- gpout 
	  			Qmin <- rev(gpout$Table[,2])[1] 
	  			}
    	    }
		Qmin <- eps * round(Qmin * 1/eps,0)
		Qvalues <- eps * round(Qvalues * 1/eps,0)
		Qvaluessame <- Qvalues == Qmin
		randStartChar <- c(randomStarts, sum(Qconverged), sum(Qvaluessame), length(unique(Qvalues)))
		names(randStartChar) <- c("randomStarts","Converged","atMinimum", "localMins")
	    r <- list(loadings = r$loadings, Th = r$Th, Table = r$Table, method = r$method,  
    	    orthogonal = TRUE, convergence = r$convergence, Gq=r$Gq, randStartChar = randStartChar)
		class(r) <- "GPArotation"
	}
	colnames(r$Table) <- c("iter", "f", "log10(s)", "alpha")
	r
}


GPFRSoblq <- function(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
       method="quartimin", methodArgs=NULL, randomStarts=0){
	if (0 < randomStarts){ Tmat <- Random.Start(ncol(A)) }
	r <- GPFoblq(A, Tmat = Tmat, normalize = normalize, eps = eps, maxit = maxit, 
				method = method, methodArgs = methodArgs)
	if (randomStarts > 1){
		Qvalues <- rev(r$Table[,2])[1]  
		Qmin <- Qvalues
		Qconverged <- r$convergence
		for (inum in 2:randomStarts){
			gpout <- GPFoblq(A, Tmat = Random.Start(ncol(A)), normalize = normalize, eps = eps, maxit = maxit, 
				method = method, methodArgs = methodArgs)
	  		Qvalues <- c(Qvalues, rev(gpout$Table[,2])[1]) 	  		
	  		Qconverged <- c(Qconverged, gpout$convergence)
	  		if (rev(gpout$Table[,2])[1] < Qmin){ 
	  			r <- gpout 
	  			Qmin <- rev(gpout$Table[,2])[1] 
	  			}
    	    }	
		Qmin <- eps * round(Qmin * 1/eps,0)
		Qvalues <- eps * round(Qvalues * 1/eps,0)
		Qvaluessame <- Qvalues == Qmin
		randStartChar <- c(randomStarts, sum(Qconverged), sum(Qvaluessame), length(unique(Qvalues)))
		names(randStartChar) <- c("randomStarts","Converged","atMinimum", "localMins")
	    r <- list(loadings = r$loadings, Phi = r$Phi, Th = r$Th, Table = r$Table, method = r$method, 
    	    orthogonal = FALSE, convergence = r$convergence,  Gq=r$Gq, randStartChar = randStartChar)
		class(r) <- "GPArotation"
	}
	dimnames(r$Phi) <- list(colnames(r$loadings),colnames(r$loadings))
	colnames(r$Table) <- c("iter", "f", "log10(s)", "alpha")
	r
}
