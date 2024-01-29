#
# This file contains the function that compiles a declarative model definition into
# something that is acceptable for the CA engines.
#

#' @title Definition of a stochastic cellular automaton
#'
#' @description  High-level definition of a stochastic cellular automaton
#'
#' @param ... A number of transition descriptions, as built by the
#'   \code{\link{transition}} function (see Details and Examples)
#'
#' @param neighbors The number of neighbors to use in the cellular automaton (4 for 4-way
#'   or von-Neumann neigborhood, or 8 for an 8-way or Moore neighborhood)
#'
#' @param wrap If \code{TRUE}, then the 2D grid on which the model is run wraps around
#'   at the edges (the top/leftmost cells will be considered neighbors of the
#'   bottom/rightmost cells)
#'
#' @param parms A named list of parameters, which should be all numeric, single values
#'
#' @param all_states The complete set of states of the model (a character vector). If
#'   unspecified, it will be guessed from the transition rules, but it is a good idea
#'   to pass it here to make sure the model definition is correct.
#'
#' @param verbose Whether information should be printed when parsing the model
#'   definition.
#'
#' @param check_model A check of the model definition is done to make sure there
#'   are no issues with it (e.g. probabilities outside the [0,1] interval, or an
#'   unsupported model definition). A quick check that should catch most problems is
#'   performed if check_model is "quick", an extensive check that tests all possible
#'   neighborhood configurations is done with "full", and no check is performed with
#'   "none".
#'
#' @param epsilon A small value under which the internal model coefficients values are
#'   considered to be equal to zero. The default value should work well here, except
#'   if you run models that have extremely small transition probabilities (<1e-8).
#'
#' @param fixed_neighborhood When not using wrapping around the edges 
#'   (\code{wrap = FALSE}), the number of neighbors per cell is variable, which can 
#'   slow down the simulation. Set this option to \code{TRUE} to consider that the number
#'   of neighbors is always four or eight, regardless of the position of the cell in the
#'   landscape, at the cost of approximate dynamics at the edges of the landscape.
#'
#' @details
#' 
#' This help page describes in detail technical points related to the definition of 
#' models in \code{chouca}. If this is your first time working with \code{chouca}, 
#' you may like the longer introduction in the vignette, accessible using 
#' \code{vignette("chouca-package")}.
#' 
#' \code{camodel} allows defining a stochastic cellular automaton model by its set of
#' transition rules. These are defined by a set of calls to the \code{transition()}
#' function. Each of these calls defines the two states of the transition, and the
#' probability as a one-sided formula involving constants and the special vectors 
#' p and q.
#' 
#' \code{transition()} calls takes three arguments: the state from which the transition
#' is defined, the state to which the transition goes, and a transition probability,
#' defined as a one-sided formula. This formula can include numerical constants,
#' parameters defined in the named list \code{parms}, and any combination of
#' \code{p['a']} and \code{q['b']}, which respectively represent the proportion 
#' of cells in a landscape in state 'a', and the proportion of neighbors of a given cell
#' in state 'b' ('a', and 'b' being, of course, any of the possible states defined in the
#' model). Such formula could typically look like \code{~ 0.2 + 0.3 * p["a"] + q["b"]}. 
#' See below for examples of model definitions. 
#' 
#' It is important to remember when using this function that \code{chouca} only
#' supports models where the probabilities depend on constant parameters, the global
#' proportion of each state in the landscape, and the local proportion of cells around
#' a given cell. In other words, all transition probabilities should have the following 
#' functional form: 
#' 
#' \deqn{a_0 + \sum_{k=1}^S g_k(q_k) + s(q, q) + s(p, q) + s(q, q)}{
#'       a_0 + \sum_{k=1}^S g_k(q_k) + s(q, q) + s(p, q) + s(q, q)}
#'
#' where \eqn{a_0} is a constant, \eqn{g_k} are univariate functions of \eqn{q_k}, 
#' the proportions of neighbors of a cell in state k, and \eqn{q} is the vector 
#' containing all the \eqn{q_k} for k between \eqn{1} and \eqn{S}, the total number of
#' states in the model. Similarly, \eqn{p} is the length-\eqn{S} vector containing the
#' proportion of cells in each state in the whole grid. \eqn{s} above is the sum, defined 
#' for two vectors \eqn{x = (x_1, ..., x_S)} and \eqn{y = (y_1, ..., y_S)} as 
#' 
#' \deqn{ 
#'   a_1 x_1^{\alpha_1} y_1^{\beta_1} + 
#'   a_2 x_1^{\alpha_2} y_2^{\beta_2} + 
#'   a_3 x_1^{\alpha_3} y_3^{\beta_3} + 
#'   a_4 x_2^{\alpha_3} y_1^{\beta_3} + 
#'   a_4 x_2^{\alpha_3} y_2^{\beta_3} + 
#'   \dots + 
#'   a_K x_S^{\alpha_K} y_S^{\beta_K}
#' }
#' 
#' where the \eqn{a_k}, \eqn{\alpha_k} and \eqn{\beta_k} are constants for all \eqn{k}, 
#' and \eqn{K} is the total number of terms (equal to \eqn{S^2}). Note that 
#' \eqn{\alpha_K} and \eqn{\beta_K} are capped to 5. This can be overriden using
#' \code{options(chouca.degmax = n)}, but we do not recommend changing it as higher 
#' values typically make the package slow and/or leads to numerical instabilities.
#' The functions \eqn{g_k} above can be any univariate functions of \eqn{q_k}, so
#' \code{chouca} effectively supports any type of transition rule involving the 
#' neighborhood of a cell, including some 'threshold' rules that involve a single state 
#' (and only one). For example, a rule such as "more than
#' 5 neighbors in a given state make a cell switch from state A to B" is OK, but
#' combining states may not be supported, such as "more than 5
#' neighbors in state A *and* 2 in state B means a cell switches from A to B". When in
#' doubt, just write your model, and \code{chouca} will tell you if it cannot run it
#' accurately by running model checks.
#' 
#' Model checks are controlled by the argument \code{check_model}. When set to "quick" 
#' or "full", a check is performed to make sure the functional form above is able to
#' accurately represent probabilities of transitions in the model, with "full" enabling
#' more extensive testing, and "none" removing it entirely. Coefficients in the formula 
#' above are rounded down to zero when below \code{epsilon}. This may be an issue if 
#' your transition probabilities are close to zero: consider reducing \code{epsilon} to 
#' a smaller value in this case, or adjusting your model parameters. 
#' 
#' When space does not wrap around (\code{wrap = FALSE}), cells in the corners
#' or in the edges will have a lower number of neighbors. The proportions
#' of cells in a given state \eqn{k}, \eqn{q_k}, will thus be computed with a reduced
#' number of cells. For example, a cell in a corner will have only 2 neighbors
#' when using a 4x4 neighborhood, so \eqn{q_k} is computed using only two cells, and 
#' can be only equal to 0, 0.5 or 1. 
#' 
#' To run a model once it is defined, the function \code{\link{run_camodel}} can be
#' used, or \code{\link{run_meanfield}} for a mean-field approximation. An initial
#' landscape for a simulation can be created using \code{\link{generate_initmat}}.
#' 
#' You can update a model definition with new parameters (all of them or a subset)
#' using the \code{\link[=update.ca_model]{update}} method. The model graph with the 
#' different states and transitions can be displayed using the \code{plot} method 
#' (this requires the package igraph).
#' 
#' @returns 
#'
#NOTE: Also update the same section in ca_library() when changing this
#'
#' This function returns a \code{\link{list}} object with class \code{ca_model}, with 
#' the following named components. Please note that most are for internal use and may 
#' change with package updates. 
#' 
#' \describe{
#'   \item{\code{transitions}}{the list of transitions of the model, as returned 
#'       by \code{\link{transition}} }
#' 
#'   \item{\code{nstates}}{the number of states of the model}
#' 
#'   \item{\code{parms}}{the parameter values used for the model}
#' 
#'   \item{\code{beta_0},\code{beta_q}, \code{beta_pp}, \code{beta_pq}, \code{beta_qq}}{ 
#'     internal tables used to represent probabilities of transitions when running
#'     simulations, these tables are for internal use and probably not interesting for 
#'     end users, but more information is provided in the package source code}
#'   
#'   \item{\code{wrap}}{Whether the model uses a toric space that wraps around the edge}
#' 
#'   \item{\code{neighbors}}{The type of neighborhood (4 or 8)}
#'   
#'   \item{\code{epsilon}}{The \code{epsilon} values used in the model definition, below 
#'     which transition probabilities are assumed to be zero}
#'   
#'   \item{\code{xpoints}}{(for internal use only) The number of values used to 
#'      represent the proportion of neighbors of a cell in each state}
#'   
#'   \item{\code{max_error}, \code{max_rel_error}}{vector of numeric values containing 
#'     the maximum error and maximum relative error on each transition probability}
#' 
#'   \item{\code{fixed_neighborhood}}{flag equal to \code{TRUE} when cells have
#'     a fixed number of neighbors}
#' 
#' }
#' 
#' 
#' @seealso run_camodel, run_meanfield, generate_initmat, run_meanfield,
#'   update.ca_model, ca_library
#'
#' @examples
#'
#' # Redefine Kubo's 1996 forest gap model
#' kubo <- camodel(
#'   transition(from = "TREE",
#'              to   = "EMPTY",
#'              prob = ~ d + delta * q["EMPTY"] ),
#'   transition(from = "EMPTY",
#'              to   = "TREE",
#'              prob = ~ alpha),
#'   parms = list(d = 0.125,
#'                delta = 0.5,
#'                alpha = 0.2),
#'   all_states = c("EMPTY", "TREE"),
#'   neighbors = 4,
#'   wrap = TRUE
#' )
#'
#' # Display it as a graph
#' plot(kubo)
#'
#' # A fun plant model
#' mod <- camodel(
#'   transition("plant", "empty", ~ death * ( 1 - (2*q["plant"]-1)^2) ),
#'   transition("empty", "plant", ~ q["plant"]^2 ),
#'   all_states = c("empty", "plant"),
#'   wrap = TRUE,
#'   neighbors = 4,
#'   parms = list(death = 0.2496)
#' )
#'
#'
#' # Conway's Game of Life
#' mod <- camodel(
#'   transition("LIVE", "DEAD", ~ q["LIVE"] < (2/8) | q["LIVE"] > (3/8)),
#'   transition("DEAD", "LIVE", ~ q["LIVE"] == (3/8)),
#'   wrap = TRUE,
#'   neighbors = 8,
#'   all_states = c("DEAD", "LIVE")
#' )
#'
#' # A spiral-generating rock-paper-scissor model
#' mod <- camodel(
#'   transition(from = "r", to = "p", ~ q["p"] > 0.25 ),
#'   transition(from = "p", to = "c", ~ q["c"] > 0.25 ),
#'   transition(from = "c", to = "r", ~ q["r"] > 0.25 ),
#'   parms = list(prob = 1),
#'   wrap = TRUE,
#'   neighbors = 8
#' )
#'
#' # Display the model as a graph
#' plot(mod)
#'
#' # Running the above model (see also the help files for the relevant functions)
#' init <- generate_initmat(mod, c(r = 1/3, p = 1/3, c = 1/3), nrow = 128)
#' out <- run_camodel(mod, init, times = seq(0, 128))
#' plot(out)
#'
#' # Update a model definition using update()
#' mod <- camodel(
#'   transition("plant", "empty", ~ m),
#'   transition("empty", "plant", ~ r * q["plant"] * ( 1 - q["plant"] ) ),
#'   all_states = c("empty", "plant"),
#'   wrap = TRUE,
#'   neighbors = 4,
#'   parms = list(m = 0.35, r = 0.4)
#' )
#'
#' mod_updated <- update(mod, parms = list(m = 0.05, r = 1))
#' init <- generate_initmat(mod_updated, c(plant = 0.8, empty = 0.2), nrow = 128)
#' out <- run_camodel(mod_updated, init, times = seq(0, 128))
#' plot(out)
#' image(out)
#'
#' # You can also specify only part of the parameters, the others will be
#' # kept to their original values
#' mod_updated <- update(mod, parms = list(m = 0.035))
#'
#'@export
camodel <- function(...,
                    neighbors,
                    wrap,
                    parms = list(),
                    all_states = NULL,
                    check_model = "quick",
                    verbose = FALSE,
                    epsilon = sqrt(.Machine[["double.eps"]]),
                    fixed_neighborhood = FALSE) {

  if ( ( ! identical(parms, list()) ) &&
        ( ! is.list(parms) || is.null(names(parms)) || any(names(parms) == "") ) ) {
    stop("parms must be a list all elements named")
  }

  if ( any( c("p", "q") %in% names(parms) ) ) {
    stop("no parameters must be named 'q' or 'p', these are reserved to refer to densities")
  }

  if ( is.logical(check_model) && isFALSE(check_model) ) {
    check_model <- "none"
  }

  if ( ! check_model %in% c("none", "full", "quick") ) {
    stop("'check_model' must be one of 'none', 'quick' or 'full'")
  }

  # Read all possible states
  msg <- function(txt) if ( verbose ) message(txt, appendLF = FALSE)
  msg("Parsing CA model definition\n")

  # Read transition objects
  transitions <- list(...)

  # Check that they are all transitions
  if ( ! all(sapply(transitions, inherits, "camodel_transition")) ) {
    stop("One of the model transition was not processed by the transition() function. Please check your function arguments and the ?camodel help on how to define model state transitions.")
  }

  states <- plyr::ldply(transitions, function(o) as.data.frame(o[c("from", "to")]) )
  uniqstates <- unique(c(states[ ,"from"], states[ ,"to"]))
  if ( is.null(all_states) ) {
    nstates <- length(uniqstates)
    states <- uniqstates
  } else {
    nstates <- length(all_states)
    states <- all_states
    if ( ! all(uniqstates %in% states) ) {
      stop("Some cell states defined in transitions are not defined in 'all_states'")
    }
  }

  # Make sure transitions are unique
  trcodes <- plyr::laply(transitions, function(tr) {
    paste(tr[["from"]], tr[["to"]], sep = " -> ")
  })
  if ( any(duplicated(trcodes)) ) {
    msg <- paste("Duplicated transition(s) in model definition:\n",
                 paste("  ", trcodes[duplicated(trcodes)], collapse = "\n"))
    stop(msg)
  }

  msg(sprintf("Using %s cell states\n", nstates))
  msg(sprintf("Using %s transitions\n", length(transitions)))

  # Parse transitions
  # In the worst case scenario, we need 25 points, so that the number of neighbors is
  # always divisible by 2, 3, 4 or 8. (25 because zero included)
  xpoints <- 1 + ifelse(wrap, neighbors, ifelse(neighbors == 8, 120, 24))

  # If we want a fixed neighborhood, regardless of whether if wrap or not, then
  # overwrite xpoints with the fixed value.
  xpoints <- ifelse(fixed_neighborhood, 1+neighbors, xpoints)

  transitions_parsed <- lapply(transitions, parse_transition, states, parms, xpoints,
                               epsilon, neighbors, check_model)

  # Pack transitions into matrices
  beta_q  <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_q")
  beta_pp <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_pp")
  beta_pq <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_pq")
  beta_qq <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_qq")
  beta_0  <- plyr::ldply(transitions_parsed, pack_table_fromto, "beta_0")

  # Extract mean/max errors
  max_error <- plyr::laply(transitions_parsed, function(o) o[["max_error"]])
  max_rel_error <- plyr::laply(transitions_parsed, function(o) o[["max_rel_error"]])


  caobj <- list(transitions = transitions,
                nstates = nstates,
                states = factor(states, levels = states), # convert explicitely to factor
                parms = parms,
                beta_0 = beta_0,
                beta_q = beta_q,
                beta_pp = beta_pp,
                beta_pq = beta_pq,
                beta_qq = beta_qq,
                wrap = wrap,
                neighbors = neighbors,
                epsilon = epsilon,
                xpoints = xpoints,
                max_error = max_error,
                max_rel_error = max_rel_error,
                fixed_neighborhood = fixed_neighborhood)

  class(caobj) <- c("ca_model", "list")
  return(caobj)
}

#'@describeIn camodel
#'
#' @param from The state from which the transition is defined
#'
#' @param to The state to which the transition is defined
#'
#' @param prob a one-sided formula describing the probability of transition between the two
#'  states (see Details section for more information).
#'
#'@export
transition <- function(from, to, prob) {
  if ( length(from) != 1 ) {
    stop("from is not a single-length vector")
  }
  if ( length(to) != 1 ) {
    stop("from is not a single-length vector")
  }
  if ( ! inherits(prob, "formula") ) {
    stop("prob must be a one-sided formula ('~ <expr>')")
  }

  o <- list(from = from,
            to   = to,
            prob = prob)
  class(o) <- c("camodel_transition", "list")

  return(o)
}

# Helper function to pack transitions into model-level tables
pack_table_fromto <- function(tr, table) {

  if ( nrow(tr[[table]]) == 0 ) {
    data.frame(from = integer(0), to = integer(0), tr[[table]])
  } else {
    data.frame(from = tr[["from"]], to = tr[["to"]], tr[[table]])
  }
}

# Update a ca_model with new arguments
# first argument needs to be 'object' to respect S3 method naming
#' @title Update a cellular automaton
#'
#' @description Update the definition of a stochastic cellular automaton
#'   (SCA), using new parameters, type of wrapping, or any other parameters
#'   entering in the definition of the model.
#'
#' @param object The SCA object (returned by \code{\link{camodel}})
#'
#' @param parms a named list of parameters, which should be all numeric,
#'   single values. If this list contains only a subset of model parameters, the
#'   old parameter values will be reused for those not provided.
#'
#' @param neighbors The number of neighbors to use in the cellular automaton
#'   ('4' for 4-way or von-Neumann neghborhood, or '8' for an 8-way or Moore
#'   neighborhood)
#'
#' @param wrap If \code{TRUE}, then the 2D grid on which the model is run wraps
#'   around at the edges (the top/leftmost cells will be considered neighbors
#'   of the bottom/rightmost cells)
#'
#' @param fixed_neighborhood When not using wrapping around the edges (\code{wrap = TRUE},
#'   the number of neighbors per cell is variable, which can slow down the simulation.
#'   Set this option to \code{TRUE} to consider that the number of neighbors is always
#'   four or eight, regardless of the position of the cell in the landscape, at the cost
#'   of approximate dynamics on the edge of the landscape.
#'
#' @param check_model A check of the model definition is done to make sure there
#'   are no issues with it (e.g. probabilities outside the [1,0] interval, or an
#'   unsupported model definition). A quick check that should catch most problems is
#'   performed if check_model is "quick", an extensive check that tests all
#'   neighborhood configurations is done with "full", and no check is performed with
#'   "none".
#'
#' @param verbose whether information should be printed when parsing the model
#'   definition.
#'
#' @param ... extra arguments are ignored
#'
#' @details
#'  This function updates some aspects of a pre-defined stochastic celullar
#'    automaton, such as parameter values, the type of neighborhood, whether
#'    to wrap around the edge of space, etc. It is handy when running multiple
#'    simulations, and only a few aspects of the model needs to be changed, such as
#'    parameter values. Note that this function cannot add or remove states to a model.
#'
#'  Note that the \code{parms} list may only specify a subset of the model parameters
#'    to update. In this case, old parameter values not specified in the call to
#'    \code{update} will be re-used.
#'
#' @seealso camodel, run_camodel...
#'
#' @returns 
#'  
#'  This function returns a list with class \code{ca_model} with the changes applied to 
#'  the original model (see \code{\link{camodel}} for details about this type of
#'  object). 
#'
#'@examples
#'
#'
#' # Update the parameters of a model
#' mussels <- ca_library("musselbed")
#' mussels[["parms"]] # old parameters
#' mussels_new <- update(mussels, parms = list(d = 0.2, delta = 0.1, r = 0.8))
#' mussels_new[["parms"]] # updated parameters
#'
#' # Update the type of neighborhood, wrapping around the edges, and
#' # the parameters
#' mussels_new <- update(mussels,
#'                       parms = list(d = 0.2, delta = 0.1, r = 0.8),
#'                       wrap = TRUE,
#'                       neighbors = 8)
#' mussels_new
#'
#' \donttest{
#' 
#' # Run the model for different values of d, the death rate of mussels
#' ds <- seq(0, 0.25, length.out = 12)
#' initmat <- generate_initmat(mussels, c(0.5, 0.5, 0), nrow = 64, ncol = 64)
#' results <- lapply(ds, function(this_dvalue) {
#'   musselmod <- update(mussels, parms = list(d = this_dvalue))
#'   run <- run_camodel(musselmod, initmat, times = seq(0, 128))
#'   data.frame(d = this_dvalue,
#'              as.data.frame(tail(run[["output"]][["covers"]], 1)))
#' })
#' 
#' results <- do.call(rbind, results)
#' plot(results[ ,"d"], results[ ,"MUSSEL"], type = "b",
#'      xlab = "d", ylab = "Mussel cover")
#'
#' }
#'
#'@export
update.ca_model <- function(object,
                            parms = NULL,
                            neighbors = NULL,
                            wrap = NULL,
                            fixed_neighborhood = NULL,
                            check_model = "quick",
                            verbose = FALSE,
                            ...) {

  if ( is.null(wrap) ) {
    wrap <- object[["wrap"]]
  }
  if ( is.null(neighbors) ) {
    neighbors <- object[["neighbors"]]
  }
  if ( is.null(parms) ) {
    parms <- object[["parms"]]
  }
  if ( is.null(fixed_neighborhood) ) {
    fixed_neighborhood <- object[["fixed_neighborhood"]]
  }

  if ( ! is.null(parms) ) {
    parms_orig <- object[["parms"]]
    if ( is.null(names(parms)) || any( names(parms) == "" ) ) {
      stop("Not all elements are named in the parameter list")
    }
    if ( ! all( names(parms) %in% names(parms_orig) ) ) {
      stop("New parameter names do not match original parameters")
    }
    parms_orig[names(parms)] <- parms
    parms <- parms_orig
  }

  # Extract model parameters, and do the call
  newcall <- c(object[["transitions"]], # always a list, so result of c() is a list
               neighbors = neighbors,
               wrap = wrap,
               fixed_neighborhood = fixed_neighborhood,
               parms = list(parms),
               all_states = list(object[["states"]]),
               check_model = check_model,
               verbose = verbose,
               epsilon = object[["epsilon"]])

  do.call(camodel, newcall)
}
