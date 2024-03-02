# Random start routine that outputs 1 orthogonal random matrix

#Random.Start <- function(k){
#  qr.Q(qr(matrix(rnorm(k*k),k)))
#  }

# Stewart, G. W. (1980). The Efficient Generation of Random Orthogonal Matrices 
# with an Application to Condition Estimators. SIAM Journal on Numerical Analysis, 
# 17(3), 403-409. http://www.jstor.org/stable/2156882
# 
# Mezzadri, F. (2007). How to generate random matrices from the classical
# compact groups. Notices of the American Mathematical Society, 54(5), 592-604.
# see https://arxiv.org/abs/math-ph/0609050
#
# This updated version changed as of GPArotation 2024-2.1
# Thanks to Yves Rosseel for pointing this out
#
Random.Start <- function(k = 2L) {
    qr.out <- qr(matrix(rnorm(k * k), nrow = k, ncol = k))
    Q <- qr.Q(qr.out)
    R <- qr.R(qr.out)
    R.diag <- diag(R)
    R.diag2 <- R.diag/abs(R.diag)
    out <- t( t(Q) * R.diag2 )
    out
}  
