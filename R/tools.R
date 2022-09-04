# 
# 
# 
# image() does not work with matrices of factors, so we define a method here for the 
# camodel_initmat class. 
# 
#'@export
image.camodel_initmat <- function(x, ...) { 
  x <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
  graphics::image(x, ...)
}
