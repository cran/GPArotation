# testing that the print.GPArotation output is identical 
# for 2 runs of quartimin rotation, that have 2
# different looking loadings matrices wrt sign and order
# the print.GPArotation should look identical


 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()

require("stats")  
require("GPArotation")  

athl <- matrix(c(
  .73, -.07, .50, .82, -.01, .27,  .77, -.46, -.22,  .78, .17, .03,
  .77, .41, .13,  .81, -.01, .27,  .71, -.45, -.30,  .82, .12, -.11,
  .66, -.15, -.45,  .39, .76, -.40), byrow=T,  ncol =3)
## z1 gives the results that have the right ordering and sign of the factors
## z2 is a random other order and sign
set.seed(238)
z1 <- quartimin(athl, Tmat = Random.Start(3))
head(z1$loadings)
set.seed(46)
z2 <- quartimin(athl, Tmat = Random.Start(3))
head(z2$loadings)


#> z1
#Oblique rotation method Quartimin converged.
#Loadings:
#         [,1]    [,2]     [,3]
# [1,]  0.9451 -0.0535 -0.18033
# [2,]  0.7725  0.1431  0.01187
# [3,]  0.1323  0.8617 -0.12847
# [4,]  0.5377  0.2050  0.28753
# [5,]  0.6888 -0.0555  0.44072
# [6,]  0.7665  0.1386  0.00987
# [7,]  0.0150  0.8967 -0.08931
# [8,]  0.4047  0.3792  0.32647
# [9,] -0.1056  0.7915  0.24071
#[10,] -0.0155 -0.0165  0.94994
#
#                [,1]  [,2]  [,3]
#SS loadings    3.034 2.405 1.401
#Proportion Var 0.303 0.240 0.140
#Cumulative Var 0.303 0.544 0.684
#
#Phi:
#      [,1]  [,2]  [,3]
#[1,] 1.000 0.554 0.259
#[2,] 0.554 1.000 0.186
#[3,] 0.259 0.186 1.000
#> z2
#Oblique rotation method Quartimin converged.
#Loadings:
#         [,1]    [,2]     [,3]
# [1,]  0.9451 -0.0535 -0.18033
# [2,]  0.7725  0.1431  0.01187
# [3,]  0.1323  0.8617 -0.12847
# [4,]  0.5377  0.2050  0.28753
# [5,]  0.6888 -0.0555  0.44072
# [6,]  0.7665  0.1386  0.00987
# [7,]  0.0150  0.8967 -0.08930
# [8,]  0.4047  0.3792  0.32647
# [9,] -0.1056  0.7915  0.24071
#[10,] -0.0155 -0.0165  0.94994
#
#                [,1]  [,2]  [,3]
#SS loadings    3.034 2.405 1.401
#Proportion Var 0.303 0.240 0.140
#Cumulative Var 0.303 0.544 0.684
#
#Phi:
#      [,1]  [,2]  [,3]
#[1,] 1.000 0.554 0.259
#[2,] 0.554 1.000 0.186
#[3,] 0.259 0.186 1.000
