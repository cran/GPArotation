### R code from vignette source 'Guide.Stex'

###################################################
### code chunk number 1: Guide.Stex:7-8
###################################################
 options(continue="  ")


###################################################
### code chunk number 2: Guide.Stex:13-14
###################################################
library("GPArotation")  


###################################################
### code chunk number 3: Guide.Stex:30-33
###################################################
  data(ability.cov)
  z <- factanal(factors = 2, covmat = ability.cov, rotation="oblimin")
  loadings(z)


