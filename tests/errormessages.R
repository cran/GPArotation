# test for to understand if there is breaking in the code
# when an error is produced.

 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()


all.ok <- TRUE

# 1-factor model loadings vector

xv <- runif(5)

# Testing if single factor models will break when error is called

#test 1
y <- try(GPArotation::quartimin(xv), TRUE)
  if (!inherits(y, "try-error")) {
	print("error messages: test 1 failed")
	all.ok <- FALSE  
	}


#test 2
y <- try(GPForth(xv, method = "quartimax"), TRUE)
  if (!inherits(y, "try-error")) {
	print("error messages: test 2 failed")
	all.ok <- FALSE  
	}

#test 3
y <- try(GPFoblq(xv, method = "quartimin"), TRUE)
  if (!inherits(y, "try-error")) {
	print("error messages: test 3 failed")
	all.ok <- FALSE  
	}

# same but with matrix instead of vector	
xw <- matrix(xv)

#test 4
y <- try(GPForth(xw, method = "quartimax"), TRUE)
if (! grep("rotation does not make sense for single factor models", attr(y, "condition")$message) ) 
{
	print("error messages: test 4 failed")
	all.ok <- FALSE  
	}

#test 5
y <- try(GPFoblq(xw, method = "quartimin"), TRUE)
if (! grep("rotation does not make sense for single factor models", attr(y, "condition")$message) ) 
{
	print("error messages: test 5 failed")
	all.ok <- FALSE  
	}

if (! all.ok) stop("some tests FAILED")

