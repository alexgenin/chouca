# 
# Functions used to create/handle CA grids
# 

#' @title Generate an initial matrix for a \code{chouca} model
#'
#' @description Helper function to create a spatially-random initial landscape (matrix)
#'   with specified covers for a cellular automaton
#'
#' @param mod A stochastic cellular automaton model created by \code{\link{camodel}}
#'
#' @param pvec A numeric vector of covers for each state in the initial configuration,
#'   possibly with named elements.
#'
#' @param nrow The number of rows of the output matrix
#'
#' @param ncol The number of columns of the output matrix
#'
#' @details
#'
#'   This function is a helper to build a starting configuration (matrix) for a
#'     stochastic cellular automaton based on the definition of the model and the
#'     specified starting covers (in \code{pvec}). It will produce a landscape with
#'     expected global cover of each state equal to the covers in \code{pvec}, and a
#'     completely random spatial structure.
#'
#'   The length of the \code{pvec} vector must match the number of possible cell states
#'     in the model. If present, the names of \code{pvec} must match the states
#'     defined in the model. In this case, they will be used to determine which state
#'     gets which starting cover instead of the order of the values.
#'
#'   The \code{pvec} will be normalized to sum to one, and the function will produce a
#'     warning if this produces a meaningful change in covers.
#'
#'   If you already have a matrix you want to use as a starting configuration, we
#'     recommend you to use \code{\link{as.camodel_initmat}} to convert it to an
#'     object that \code{\link{run_camodel}} can use.
#' 
#' @returns 
#' 
#'  This function returns a matrix containing values as factors, with levels
#'    corresponding to the model states (defined in the \code{mod} argument) and
#'    dimensions set by \code{nrow} and \code{ncol}. This matrix has the
#'    \code{\link{class}} \code{camodel_initmat} so that it can be displayed with the
#'    \code{image} generic function. 
#' 
#' @seealso as.camodel_initmat
#' 
#' @examples
#'
#' # Run the Game of Life starting from a random grid
#' game_of_life <- ca_library("gameoflife")
#' grid <- generate_initmat(game_of_life, c(LIVE = .1, DEAD = .9), nrow = 64)
#' out <- run_camodel(game_of_life, grid, times = seq(0, 128))
#' image(out) # final configuration
#' 
#' # Logistic growth of plants
#' mod <- camodel(
#'   transition(from = "empty", to = "plant", ~ r * p["plant"]),
#'   transition(from = "plant", to = "empty", ~ m),
#'   parms = list(r = 1, m = .03),
#'   wrap = TRUE,
#'   neighbors = 8
#' )
#' grid <- generate_initmat(mod, c(empty = .99, plant = .01), nrow = 128)
#' image(grid) # initial state
#' out <- run_camodel(mod, grid, times = seq(0, 30))
#' image(out) # final state
#' plot(out) #
#' 
#'@export
generate_initmat <- function(mod, pvec, nrow, ncol = nrow) {

  if ( any(is.na(pvec)) ) {
    stop("NAs in pvec are not supported")
  }

  ns <- mod[["nstates"]]
  if ( length(pvec) != ns ) {
    stop("the length of state covers does not match the number of states in the model definition")
  }

  if ( ! is.null(names(pvec)) ) {
    diff <- setdiff(names(pvec), mod[["states"]])
    if ( length(diff) > 0 ) {
      stop("State names in 'pvec' do not match the states defined in the model")
    }
    # Make sure order is right. /!\ indexing with factors == indexing with integers !!!
    pvec <- pvec[ as.character(mod[["states"]]) ]
  }

  spvec <- sum(pvec)
  if ( sum(pvec) > 1.000001 | sum(pvec) < 0.99999 ) {
    warning("The initial covers do not sum to one, they will be rescaled")
  }
  pvec <- pvec / spvec

  # Generate the matrix
  m <- matrix(sample(mod[["states"]],
                     replace = TRUE, prob = pvec, size = nrow * ncol),
              nrow = nrow, ncol = ncol)

  # Adjust and set class
  m <- as.camodel_initmat(m, levels = mod[["states"]])

  return(m)
}

#'@title Convert a matrix to a CA model landscape
#'
#'@description Convert a matrix to a CA model landscape for later use with
#'  \link{run_camodel} or \link{run_meanfield}.
#'
#'@param m The input matrix (numeric, character or factor)
#'
#'@param levels The levels to use in the resulting landscape. If \code{NULL}, the unique
#'  values of the input matrix are used as levels. Set this manually if you want the
#'  resulting landscape to have extra levels that are not present in the original matrix.
#'
#' @returns This function returns a matrix containing values as factors, with levels
#'   corresponding to the \code{levels} argument. This matrix has the
#'   \code{\link{class}} \code{camodel_initmat} so that it can be displayed with the
#'   \code{image} generic function and works well with CA-related functions 
#'   (such as \code{\link{run_camodel}}). 
#'
#'@seealso generate_initmat, run_camodel, run_meanfield
#'
#'@examples
#' 
#' # Simple conversion of a matrix with regular patterns
#' x <- seq(0, 2 * pi, l = 256)
#' z <- outer(x, x, function(x, y) as.numeric(sin(10*x) + cos(10*y) > 0.8))
#' mat <- as.camodel_initmat(z)
#' summary(mat)
#' image(mat)
#' 
#' # This is a character matrix. We need to convert it to use it as input to
#' # run_camodel()
#' size <- 64
#' m <- matrix(ifelse(runif(size^2) < .5, "TREE", "EMPTY"), nrow = size, ncol = size)
#' m <- as.camodel_initmat(m)
#' summary(m) # this is a landscape object
#' image(m)
#'
#' # Start a simulation using this matrix
#' mod <- ca_library("forestgap")
#' out <- run_camodel(mod, m, seq(0, 256))
#' plot(out)
#'
#' \donttest{
#' # Run a glider in the game of life
#' mod <- ca_library("gameoflife")
#' init <- matrix(c(0, 0, 1, 0, 0, 0, 0,
#'                  0, 0, 0, 1, 0, 0, 0,
#'                  0, 1, 1, 1, 0, 0, 0,
#'                  0, 0, 0, 0, 0, 0, 0,
#'                  0, 0, 0, 0, 0, 0, 0,
#'                  0, 0, 0, 0, 0, 0, 0),
#'                 nrow = 6, ncol = 7, byrow = TRUE)
#' init[] <- ifelse(init == 1, "LIVE", "DEAD")
#' # image() does not work on init here without conversion by as.camodel_initmat
#' init <- as.camodel_initmat(init)
#' image(init)
#' 
#' # Run the model and display simulation output as it is running
#' ctrl <- list(custom_output_fun = landscape_plotter(mod, fps_cap = 5),
#'              custom_output_every = 1)
#' out <- run_camodel(mod, init, times = seq(0, 32), control = ctrl)
#' }
#'
#'@export
as.camodel_initmat <- function(m, levels = NULL) {

  if ( inherits(m, "camodel_initmat") ) {
    return(m)
  }

  if ( ! is.matrix(m) ) {
    stop("This object cannot be converted into a camodel_initmat object")
  }
  ml <- m

  if ( is.null(levels) ) {
    levels <- unique(as.vector(m))
  }

  # If not a factor, convert to it -> this will create levels based on m
  if ( ! is.factor(m) ) {
    ml <- factor(m, levels = levels)
    dim(ml) <- dim(m)
  }

  class(ml) <- c("camodel_initmat", "factor", "matrix")
  ml
}
