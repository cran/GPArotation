# Random start routine that outputs 1 orthogonal random matrix

Random.Start <- function(k){
  qr.Q(qr(matrix(rnorm(k*k),k)))
  }
