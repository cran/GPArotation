# The functions for iterative re-weighted least squares for Lp rotation were provided by 
# Xinyi Liu. The functions include an inner and an outer loop for each iteration.
# Adding randomstarts then adds another loop for random starts inclusion.

# lpQ and lpT are wrapper function for random starts.  This function has to exist separately from the
# other RS functions b/c of the call to the lpT.GPFoblq and lpT.GPForth functions.  
# Generalizing this function to work with the existing GPFRSoblq would involve a
# change to the legacy code for GPFoblq to include an if-then calling for lp.
# Since I did not want to change the legacy code, the RS function has to exist
# separately.   Coen Bernaards, February 2025.

# 4 functions in the lp.R file:
# lpQ: the wrapper function for oblique Lp rotation
# GPFoblq.lp: the main Lp GPFoblq function
# lpT: the wrapper for orthogonal rotation
# GPForth.lp: the main Lp GPForth function


# lpQ is oblique rotation
lpQ<-function (A, Tmat = diag(ncol(A)),p=1, normalize = FALSE, eps = 1e-05, 
    maxit = 1000, randomStarts = 0, gpaiter = 5) 
{
    # message("Warnings may appear, but as long as the convergence in the output value is TRUE, the warnings can be ignored.")
    if (0 < randomStarts) {
        Tmat <- Random.Start(ncol(A))
    }
    if (normalize){
       message("Normalization is not recommended for Lp rotation. Use may have unexpected results")
    }
    r <- GPFoblq.lp(A,Tmat = Tmat, p,normalize = normalize, eps = eps, maxit = maxit, gpaiter = gpaiter)
    if (randomStarts > 1) {
        	#3 Qvalues <- sum(abs(r$loadings)^p)
        Qvalues <- sum(abs(r$loadings)^p)/nrow(A)
        Qmin <- Qvalues
        Qconverged <- r$convergence
        for (inum in 2:randomStarts) {
            gpout <-GPFoblq.lp(A,Tmat = Random.Start(ncol(A)), p,normalize = normalize, eps = eps, maxit = maxit, gpaiter = gpaiter)
            Qvalues <- c(Qvalues, sum(abs(gpout$loadings)^p)/nrow(A))
            Qconverged <- c(Qconverged, gpout$convergence)
            if (sum(abs(gpout$loadings)^p)/nrow(A) < Qmin) {
                r <- gpout
                #3 Qmin <- sum(abs(gpout$loadings)^p)
                Qmin <- sum(abs(gpout$loadings)^p)/nrow(A)
            }
        }
        Qmin <- eps * round(Qmin * 1/eps, 0)
        Qvalues <- eps * round(Qvalues * 1/eps, 0)
        Qvaluessame <- Qvalues == Qmin
        randStartChar <- c(randomStarts, sum(Qconverged), sum(Qvaluessame), 
            length(unique(Qvalues)))
        names(randStartChar) <- c("randomStarts", "Converged", 
            "atMinimum", "localMins")
        r <- list(loadings = r$loadings, Phi = r$Phi, Th = r$Th, 
            Table = r$Table, method = r$method, orthogonal = FALSE, 
            convergence = r$convergence, Gq = r$Gq, randStartChar = randStartChar, 
            Qvalues = Qvalues)
        class(r) <- "GPArotation"
    }
    dimnames(r$Phi) <- list(colnames(r$loadings), colnames(r$loadings))
    return(r)
}
GPFoblq.lp <-function (A, Tmat = diag(rep(1, 
    ncol(A))),p=1,  normalize = FALSE, eps = 1e-05,maxit = 10000, gpaiter = 5) 
{
    j_sum <- 0
    start_time <- proc.time()
    convergence <- FALSE
    #2 W <- (A^2 + eps)^(p/2 - 1)
    W <- ((A%*%t(solve(Tmat)))^2 + eps)^(p / 2 - 1)# corrected initial weights 
    for (it in 1:maxit) {
        # r <- GPFoblq(A, Tmat = Tmat, normalize = normalize, eps = eps, 
        #     maxit = gpaiter, method = "lp.wls", methodArgs = list(W = W))
 		r <- suppressWarnings(GPFoblq(A, Tmat = Tmat, normalize = normalize, eps = eps, 
            maxit = gpaiter, method = "lp.wls", methodArgs = list(W = W))) 
        T_new <- r$Th
        L <- r$loadings
        k <- nrow(r$Table)
        ft <- r$Table[k, 2]
        j <- r$Table[k, 1]
        j_sum <- j_sum + j
        W <- (L^2 + eps)^(p/2 - 1)
        if (max(abs(L - A %*% solve(t(Tmat)))) < eps) {
            Tmat <- T_new
            convergence <- TRUE
            break
        }
        Tmat <- T_new
    }
    Phi <- t(Tmat) %*% Tmat
    dimnames(Phi) <- list(colnames(A), colnames(A))
    Table <- data.frame(iter = it, f = sum(abs(L)^p)/nrow(L), time = proc.time()[3] - 
        start_time[3])
    r <- list(loadings = L, Phi = Phi, Th = Tmat, Table = Table, 
        method = paste0("Lp rotation, p=", p), orthogonal = FALSE, 
        convergence = convergence)  
    class(r) <- "GPArotation"
    return(r)
}



# lpT for Orthogonal Rotation
lpT<-function (A, Tmat = diag(ncol(A)),p=1, normalize = FALSE, eps = 1e-05, 
    maxit = 1000, randomStarts = 0,gpaiter = 5) 
{
    # message("Warnings may appear, but as long as the convergence in the output value is TRUE, the warnings can be ignored.")
    if (0 < randomStarts) {
        Tmat <- Random.Start(ncol(A))
    }
    if (normalize){
       message("Normalization is not recommended for Lp rotation. Use may have unexpected results")
    }
    r <- GPForth.lp(A,Tmat = Tmat, p,normalize = normalize, eps = eps, maxit = maxit, gpaiter = gpaiter)
    if (randomStarts > 1) {
        Qvalues <- sum(abs(r$loadings)^p)/nrow(A)
        Qmin <- Qvalues
        Qconverged <- r$convergence
        for (inum in 2:randomStarts) {
            gpout <-GPForth.lp(A,Tmat = Random.Start(ncol(A)), p,normalize = normalize, eps = eps, maxit = maxit, gpaiter = gpaiter)
            Qvalues <- c(Qvalues, sum(abs(gpout$loadings)^p)/nrow(A))
            Qconverged <- c(Qconverged, gpout$convergence)
            #3 if (sum(abs(gpout$loadings)^p) < Qmin) {
            if (sum(abs(gpout$loadings)^p)/nrow(A)< Qmin) {
                r <- gpout
                #3 Qmin <- sum(abs(gpout$loadings)^p)
                Qmin <- sum(abs(gpout$loadings)^p)/nrow(A)
            }
        }
        Qmin <- eps * round(Qmin * 1/eps, 0)
        Qvalues <- eps * round(Qvalues * 1/eps, 0)
        Qvaluessame <- Qvalues == Qmin
        randStartChar <- c(randomStarts, sum(Qconverged), sum(Qvaluessame), 
            length(unique(Qvalues)))
        names(randStartChar) <- c("randomStarts", "Converged", 
            "atMinimum", "localMins")
        r <- list(loadings = r$loadings, Phi = r$Phi, Th = r$Th, 
            Table = r$Table, method = r$method, orthogonal = TRUE, 
            convergence = r$convergence, Gq = r$Gq, randStartChar = randStartChar, 
            Qvalues = Qvalues)
        class(r) <- "GPArotation"
    }
    dimnames(r$Phi) <- list(colnames(r$loadings), colnames(r$loadings))
    return(r)
}
GPForth.lp <-function (A, Tmat = diag(rep(1, 
    ncol(A))),p=1,  normalize = FALSE, eps = 1e-05,maxit = 10000, gpaiter = 5) 
{
    j_sum <- 0
    start_time <- proc.time()
    convergence <- FALSE
    #2 W <- (A^2 + eps)^(p/2 - 1)
    W <- ((A%*%t(solve(Tmat)))^2 + eps)^(p / 2 - 1)# corrected initial weights 
    for (it in 1:maxit) {
        #r <- GPForth(A, Tmat = Tmat, normalize = normalize, eps = eps, 
        #    maxit = gpaiter, method = "lp.wls", methodArgs = list(W = W))
        r <- suppressWarnings(GPForth(A, Tmat = Tmat, normalize = normalize, eps = eps, 
            maxit = gpaiter, method = "lp.wls", methodArgs = list(W = W)))
        T_new <- r$Th
        L <- r$loadings
        k <- nrow(r$Table)
        ft <- r$Table[k, 2]
        j <- r$Table[k, 1]
        j_sum <- j_sum + j
        W <- (L^2 + eps)^(p/2 - 1)
        if (max(abs(L - A %*% solve(t(Tmat)))) < eps) {
            Tmat <- T_new
            convergence <- TRUE
            break
        }
        Tmat <- T_new
    }
    Phi <- t(Tmat) %*% Tmat
    dimnames(Phi) <- list(colnames(A), colnames(A))
    Table <- data.frame(iter = it, f = sum(abs(L)^p)/nrow(L), time = proc.time()[3] - 
        start_time[3])
    r <- list(loadings = L, Phi = Phi, Th = Tmat, Table = Table, 
        method = paste0("Lp rotation, p=", p), orthogonal = TRUE, 
        convergence = convergence) 
    class(r) <- "GPArotation"
    return(r)
}


