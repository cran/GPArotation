# Kaiser normalization

NormalizingWeight <- function(A, normalize=FALSE){
 if ("function" == mode(normalize)) normalize <- normalize(A)
 if (is.logical(normalize)){
    if (normalize) normalize <- sqrt(rowSums(A^2))
    else return(array(1, dim(A)))
    }
 if (is.vector(normalize)) 
    {if(nrow(A) != length(normalize))
        stop("normalize length wrong in NormalizingWeight")
     return(array(normalize, dim(A)))
    }
 stop("normalize argument not recognized in NormalizingWeight")
}
