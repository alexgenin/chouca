
# Define some classif neighborhood kernels
nb_kernel_von_neumann <- matrix(c(0, 1, 0, 
                                  1, 0, 1, 
                                  0, 1, 0), 
                                byrow = TRUE, ncol = 3, nrow = 3)  
nb_kernel_moore <- matrix(c(1, 1, 1, 
                            1, 0, 1, 
                            1, 1, 1), 
                          byrow = TRUE, ncol = 3, nrow = 3)  

# Build a neighborhood kernel
# 
#'@title Define the neighborhood used in a stochastic cellular automaton 
#' 
#' @description This function builds standard neighborhood kernels for use in a 
#'   cellular automaton
#' 
#' @param type The type of neighborhood kernel, either "von_neumann", "moore", or
#'   "circular" 
#' 
#' @param distance The distance (sometimes called "order") that the neighborhood should
#'   span. A value of 1 for example will only consider neighbors directly adjacent to a
#'   focal cell, and a value of 3 will consider neighbors that are up to three cells away 
#'   from a focal cell. 
#' 
#' @details 
#' 
#' A key concept in a stochastic cellular automaton is the neighborhood, which defines
#' which cells are used to compute \eqn{q_i}, the proportion of cells around a focal 
#' cell in state \eqn{i}. 
#' 
#' Typical neighborhoods include the Moore and Von Neumann neighborhoods of order 1, 
#' which means that \eqn{q} will be computed using the 8 cells all around a focal cell, 
#' or only its 4 nearest neighbors, respectively. These types of neighborhoods can be 
#' extended to larger distances, i.e. considering neighbors that are not only adjacent 
#' to a cell, but at a given integer distance \code{distance}. See Examples section to 
#' visualize different types of neighborhoods. 
#' 
#'@examples
#'
#'# Classic moore and von_neumann neighborhoods
#'print( nb_kernel("moore", 1) )
#'
#'print( nb_kernel("von_neumann", 1) )
#'
#'# Moore kernel of order 2
#'ker <- nb_kernel("moore", 2)
#'image(ker)
#' 
#'# Different kernel types and distances 
#'opar <- par(no.readonly = TRUE)
#'par(mfcol = c(3, 5)) 
#'lapply(c(1, 3, 6, 9, 12), function(distance) { 
#'  image( nb_kernel("von_neumann", distance) )
#'  title(sprintf("Von Neumann kernel, distance %s", distance))
#'  image( nb_kernel("moore", distance) )
#'  title(sprintf("Moore kernel, distance %s", distance))
#'  image( nb_kernel("circular", distance) )
#'  title(sprintf("Circular kernel, distance %s", distance))
#'})
#'
#'par(opar)
#' 
#'# Circular kernel with a given distance
#'opar <- par(no.readonly = TRUE)
#'par(mfrow = c(3, 3)) 
#'lapply(seq(1, 9*2, by = 2), function(distance) { 
#'  image( nb_kernel("circular", distance) )
#'  title(sprintf("Circular kernel, distance %s", distance))
#'})
#'par(opar)
#'
#'
#'@export
nb_kernel <- function(type, distance) { 
  if ( length(type) != 1 ) { 
    stop("The kernel type must be a length-1 character vector")
  }
  
  cr <- distance + 1 # center
  
  if ( type == "moore" ) { 
    k <- matrix(TRUE, nrow = 2*distance + 1, ncol = 2*distance + 1)
    k[1+distance, 1+distance] <- FALSE
    
  } else if ( type %in% c("diamond", "von_neumann") ) { 
    
    k <- matrix(FALSE, nrow = 2*distance + 1, ncol = 2*distance + 1)
    for ( i in seq(-distance, distance) ) { 
      for ( j in seq(-distance, distance) ) { 
        is_in <- ( abs(i) + abs(j) ) <= distance 
        k[cr+i, cr+j] <- is_in
      }
    }
  
  } else if ( type == "circular" ) { 
    
    k <- matrix(FALSE, nrow = 2*distance + 1, ncol = 2*distance + 1)
    for ( i in seq(-distance, distance) ) { 
      for ( j in seq(-distance, distance) ) { 
        is_in <- ( i^2 + j^2 ) <= distance^2
        k[cr+i, cr+j] <- is_in
      }
    }
    
  }
  
  # Center of kernel is always false
  k[cr, cr] <- FALSE
  
  return(k)
}

# Build the neighborhood kernel depending on neighbor specification 
build_neighbor_kernel <- function(neighbors) { 
  
  if ( length(neighbors) == 1 && neighbors == 4 ) { 
    kernel <- nb_kernel("von_neumann", 1)
  }
  
  if ( length(neighbors) == 1 && neighbors == 8 ) { 
    kernel <- nb_kernel("moore", 1)
  }
  
  if ( is.matrix(neighbors) ) { 
    kernel <- neighbors
  }
  
  if ( ncol(kernel) != nrow(kernel) ) { 
    stop("The neighbors kernel must be a square matrix")
  }
  
  if ( ncol(kernel) %% 2 == 0 || nrow(kernel) %% 2 == 0  ) { 
    stop("The neighbors matrix must have an odd number of rows and columns")
  }
  
  # Is this needed ?
  # Make sure center is FALSE
  # c <- 1 + (ncol(kernel) -1)/ 2
  # kernel[c, c] <- 0
  
  # Enforce logical type
  kernel <- matrix(kernel > 0.5, nrow = nrow(kernel), ncol = ncol(kernel))
  
  kernel
}
