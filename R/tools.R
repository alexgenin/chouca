# 
# 
# 
#'@export
image.camodel_initmat <- function(x, ...) { 
  x <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
  image(x)
}
