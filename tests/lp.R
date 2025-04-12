
 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()



#sortFac <- function(x){  # Based on Fungible faSort 
#  vx <- order(colSums(x$loadings^2), decreasing = TRUE)
#  Dsgn <- diag(sign(colSums(x$loadings^3))) [ , vx]
#  x$Th <- x$Th %*% Dsgn 
#  x$loadings <- x$loadings %*% Dsgn
#  if ("Phi" %in% names(x)) { 
#    x$Phi <- diag(1/diag(Dsgn)) %*% x$Phi %*% Dsgn
#    }
#  x 
#} 


sortFac <- function(x){  # Based on Fungible faSort 
  vx <- order(colSums(x$loadings^2), decreasing = TRUE)
  Dsgn <- diag(sign(colSums(x$loadings^3))) [ , vx]
  x$Th <- x$Th %*% Dsgn 
  x$loadings <- x$loadings %*% Dsgn
  if ("Phi" %in% names(x)) { 
    x$Phi <- diag(1/diag(Dsgn)) %*% x$Phi %*% Dsgn
    }
  x 
} 
fuzz <- 1e-5 # using  eps=1e-5 these tests do not do better than this
all.ok <- TRUE  

#######################################
#######################################
#test 1

L <- rbind(diag(3),diag(3),diag(3),diag(3),diag(3),diag(3),diag(3),diag(3),diag(3),diag(3))
True_rot<-matrix(c(1,0.02079577,0.5024378,
                   0,0.99978374,0.2635086,
                   0,0,0.8234801),3,3,byrow=TRUE)

L1 <- L%*%t(True_rot)

r1 <- lpQ(L1,diag(3),maxit=1000)

tst <- t(matrix(c(
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881,
 9.990889e-01, 6.082765e-05, 0.0018180081,
 5.911473e-05, 9.997556e-01, 0.0008914537,
 1.819385e-03, 8.932991e-04, 0.9988458881
 ), 3, 30))

 if( fuzz < max(abs(r1$loadings - tst ))) {
    cat("irls: Calculated value is not the same as test value in test test 1. Value:\n")
    print(r1$loadings, digits=18)
    cat("difference:\n")
    print(r1$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 

# test 2
set.seed(1001)
r2 <- lpQ(L1,p=1,randomStarts=10)

tst <- t(matrix(c(
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399,
 6.068202e-05, 9.990904e-01, 0.0018150388,
 9.997560e-01, 5.852987e-05, 0.0008900377,
 8.933085e-04, 1.819065e-03, 0.9988460399
 ), 3, 30))
 

 if( fuzz < max(abs(sortFac(r2)$loadings - tst ))) {
    cat("irls Calculated value is not the same as test value in test 2. Value:\n")
    print(r2$loadings, digits=18)
    cat("difference:\n")
    print(r2$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 



#test 3
data("WansbeekMeijer", package="GPArotation")
fa.unrotated  <- factanal(factors = 2, covmat=NetherlandsTV, rotation="none")

set.seed(100102)
r3 <- lpQ(fa.unrotated$loadings, p=0.75, randomStarts = 50)

tst <- t(matrix(c(
-0.002173496,  0.792259797,
 0.100607850,  0.781465993,
 0.002154186,  0.772058937,
 0.617379876,  0.138657906,
 0.707274693,  0.090566478,
 0.822845923, -0.009242103,
 0.725228135,  0.002691264
), 2, 7))

 if( fuzz < max(abs(sortFac(r3)$loadings - tst ))) {
    cat("irls: Calculated value is not the same as test value in test 3. Value:\n")
    print(r3$loadings, digits=18)
    cat("difference:\n")
    print(r3$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 



#test 4
data("WansbeekMeijer", package="GPArotation")
fa.unrotated  <- factanal(factors = 3, covmat=NetherlandsTV, rotation="none")

set.seed(100102)
r4 <- lpT(fa.unrotated$loadings, p=0.75, randomStarts = 50)

tst <- t(matrix(c(
0.4231105,  0.669080641,  0.0078944327,
0.5205654,  0.659289508, -0.0008080892,
0.4202373,  0.648462232, -0.0127761636,
0.7171426,  0.100648821, -0.0888999851,
0.7436880,  0.086465344,  0.0308002902,
0.8111649, -0.001129993, -0.0001462372,
0.7373292, -0.000421333,  0.6718224963
), 3, 7))

 if( fuzz < max(abs(sortFac(r4)$loadings - tst ))) {
    cat("irls: Calculated value is not the same as test value in test 4. Value:\n")
    print(r4$loadings, digits=18)
    cat("difference:\n")
    print(r4$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 


if (! all.ok) stop("some tests FAILED")
