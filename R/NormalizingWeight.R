# Kaiser normalization

#NormalizingWeight <- function(A, normalize=FALSE){
# if ("function" == mode(normalize)) normalize <- normalize(A)
# if (is.logical(normalize)){
#    if (normalize) normalize <- sqrt(rowSums(A^2))
#    else return(array(1, dim(A)))
#    }
# if (is.vector(normalize)) 
#    {if(nrow(A) != length(normalize))
#        stop("normalize length wrong in NormalizingWeight")
#     return(array(normalize, dim(A)))
#    }
# stop("normalize argument not recognized in NormalizingWeight")
#}
# 

#
# Version below submitted by Kim-Laura Speck, Uni Kassel, 25 October 2023
#
# avoid NaNs in matrix A by adding machine precision values to zeros

NormalizingWeight <- function(A, normalize=FALSE){
  # custom function to normalize; input from user
  if ("function" == mode(normalize)) normalize <- normalize(A)
  # Kaiser normalization
  if (is.logical(normalize)){
    if (normalize) {
      normalize <- sqrt(rowSums(A^2)) # this is only a vector
      # avoid division by zero exceptions by checking that double != exactly zero
      idxZero <- which(normalize == 0) 
      # add machine precision to values that are exactly zero
      normalize[idxZero] <- normalize[idxZero] + .Machine$double.eps
    } else return(array(1, dim(A)))
  }
  if (is.vector(normalize)) 
  {if(nrow(A) != length(normalize))
    stop("normalize length wrong in NormalizingWeight")
    return(array(normalize, dim(A)))
  }
  stop("normalize argument not recognized in NormalizingWeight")
}
