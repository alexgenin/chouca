# 
# Advanced methods for defining ca models
# 
# Simple interface for models of the following form using explicit matrices
# 
# P(i -> j) = beta0 + sum_i { beta_p[i] * p[i] } + sum_i { beta_q[i] * q[i] }
#

#'@title Definition of a stochastic cellular automaton based on numerical arrays
#'
#'@description Definition of a stochastic cellular automaton using arrays and matrices 
#'  of numbers instead of symbolic expressions (as in \code{\link{camodel}})
#'
#'@inheritParams camodel
#'
#'@param beta_0 the constant term of the probabilities of transition (see Details 
#'  section), given as a square matrix of real numbers
#'
#'@param beta_p the coefficient for p[i] for i between 1 and the number of states 
#'  (see Details section), as a cubic array of numbers
#'
#'@param beta_q the coefficient for q[i] for i between 1 and the number of states 
#'  (see Details section), as a cubic array of numbers
#'
#'@param build_transitions Transition definitions 
#'  (as if produced by \code{\link{transition}}) will be included in the resulting 
#'  object if \code{TRUE}. Set to \code{FALSE} to disable this, which may improve 
#'  speed at the cost of loss of some of the package functionality
#'
#'@details 
#'
#' \code{\link{camodel}} will perform badly when defining models with a large number of 
#'   transitions or states. This function provides a much faster alternative, but it 
#'   assumes that all transitions from state \eqn{i} to state \eqn{j} follow the 
#'   following form: 
#'   
#' \deqn{ 
#'   P(i \rightarrow j) = \beta^0_{i, j} + 
#'   \beta^p_{i, j, 1} p_1 + \beta^p_{i, j, 2} p_2 + ... + \beta^p_{i, j, K} p_K + 
#'   \beta^q_{i, j, 1} q_1 + \beta^p_{i, j, 2} q_2 + ... + \beta^p_{i, j, K} q_K
#' }
#' 
#' where \eqn{\beta_0}, \eqn{\beta^p} and \eqn{\beta_q} are constant model paramaters, 
#' and \eqn{p}, \eqn{q} are the proportion of cells in the landscape and neighborhood 
#' in each state, respectively. \eqn{K} is the number of states in the model. 
#' 
#' \eqn{\beta^0} is a K x K square matrix, and \eqn{\beta^p} and 
#' \eqn{\beta^q} are a K x K x K cubic arrays (with again, K the number of states).
#' The function will check that this is valid and throw an error if this is not true. 
#' Make sure the values in those arrays are in the right order (it is 
#' \code{array[from, to, coef_for_state_k]}). To make sure the model definition is 
#' correct, you can print the model on the R command line to see if transition 
#' definitions were correctly specified. 
#' 
#' The resulting object will be similar to what would be produced by 
#' \code{\link{camodel}}, but model-building will be much quicker for complicated models 
#' (lots of states and transitions). You can further improve 
#' performance by setting \code{build_transitions} to \code{FALSE}. In this case, the 
#' generation of symbolic model transitions will be skipped entirely. They are not needed
#' to run the model, but can be handy to make sure you specified your model correctly. 
#' Disable the generation of symbolic model transitions only when you are sure your model 
#' is correctly-defined. 
#' 
#' The diagonal entries in the matrix \code{beta_0} are zero by definition, and the 
#' function will warn if that is not the case. Note that this function will not check 
#' whether the model can yield probabilities below zero or above one. 
#' 
#'@references 
#' 
#' Durrett, Richard, and Simon A. Levin. 1994. “Stochastic Spatial Models: A User’s 
#' Guide to Ecological Applications.” Philosophical Transactions of the Royal Society
#' of London. Series B: Biological Sciences 343 (1305): 329–50. 
#' https://doi.org/10.1098/rstb.1994.0028.
#' 
#'@seealso camodel 
#'
#'@examples
#'
# # Implement Kubo's model of forest dynamics using camodel() and camodel_mat()
#'mod_classic <- camodel(
#'  transition(from = "TREE",
#'              to   = "EMPTY",
#'              prob = ~ 0.125 + 0.05 * q["EMPTY"] ),
#'  transition(from = "EMPTY",
#'              to   = "TREE",
#'              prob = ~ 0.2 * p["TREE"]),
#'  neighbors = 4,
#'  wrap = TRUE,
#'  all_states = c("EMPTY", "TREE"),
#'  check_model = "quick"
#')
#'
#'states <- c("EMPTY", "TREE")
#'beta_0_mat <- diag(2) * 0
#'colnames(beta_0_mat) <- rownames(beta_0_mat) <- states
#'beta_0_mat["TREE", "EMPTY"] <- 0.125 
#'
#'beta_p_array <- array(0, dim = rep(2, 3), 
#'                      dimnames = list(states, states, states))
#'beta_p_array["EMPTY", "TREE", "TREE"] <- 0.2 
#'
#'beta_q_array <- array(0, dim = rep(2, 3), 
#'                      dimnames = list(states, states, states))
#'beta_q_array["TREE", "EMPTY", "EMPTY"] <- 0.05
#'
#'mod_mat <- camodel_mat(beta_0 = beta_0_mat, 
#'                       beta_p = beta_p_array, 
#'                       beta_q = beta_q_array, 
#'                       all_states = states, 
#'                       neighbors = 4, 
#'                       wrap = TRUE)
#' 
#' 
#'# Voter model (defined as in Durrett & Levin, 1994, p342). See also ?ca_library
#'gamma <- 0.1
#'kappa <- 30
#'coefs <- array(0, dim = rep(kappa, 3))
#'for ( state in seq.int(kappa) ) { 
#'  coefs[ , state, state] <- gamma
#'  coefs[state, state, ] <- 0
#'}
#' 
#'mod <- camodel_mat(
#'  beta_q = coefs, 
#'  wrap = TRUE, 
#'  all_states = as.character(seq.int(kappa)), 
#'  neighbors = 8, 
#'  build_transitions = TRUE
#')
#'
#'# This model takes a long time to run for high values of kappa
#'\donttest{ 
#' nrows <- ncols <- 256
#' initmm <- generate_initmat(mod, rep(1/kappa, kappa), 
#'                             nrow = nrows, ncol = ncols)
#' iters <- seq(0, 128)
#' run_camodel(mod, initmm, iters, control = list(engine = "compiled"))
#'}
#'
#'@export
camodel_mat <- function(beta_0 = NULL, 
                        beta_p = NULL, 
                        beta_q = NULL, 
                        all_states, 
                        neighbors,
                        wrap, 
                        fixed_neighborhood = FALSE,
                        epsilon = sqrt(.Machine[["double.eps"]]), 
                        build_transitions = TRUE) { 
  
  ns <- unique(c(dim(beta_0), dim(beta_p), dim(beta_q), length(all_states)))
  
  if ( ! length(ns) == 1 ) { 
    msg <- paste0("Mismatch in array dimensions:\n", 
                  "  all_states:", length(all_states), 
                  "  beta_0:", paste(dim(beta_0), collapse = " "), 
                  "  beta_q:", paste(dim(beta_q), collapse = " "), 
                  "  beta_p:", paste(dim(beta_p), collapse = " "), "\n", 
                  "All dimensions should be equal")
    stop(msg)
  }
  
  if ( is.null(beta_0) ) { 
    beta_0 <- zero_matrix(ns)
  }
  
  if ( is.null(beta_p) ) { 
    beta_p <- zero_array(ns)
  }
  
  if ( is.null(beta_q) ) { 
    beta_q <- zero_array(ns) 
  }
  
  if ( ! is_square(beta_0) ) { 
    stop("beta_0 must be a square matrix")
  }
  
  if ( ! is_cube(beta_p) ) { 
    stop("beta_p must be a cubic array")
  }
  
  if ( ! is_cube(beta_q) ) { 
    stop("beta_q must be a cubic array")
  }
  
  if ( any(diag(beta_0) != 0) ) { 
    warning("Diagonal entries in beta_0 are non-zero, they will be ignored")
  }
  
  warn_p <- warn_q <- FALSE
  for ( from in seq.int(ns) ) { 
    if ( ! warn_p && any(beta_p[from, from, ] != 0) ) { 
      warning("Diagonal entries in beta_p are non-zero, they will be ignored")
      warn_p <- TRUE
    }
    if ( ! warn_q && any(beta_q[from, from, ] != 0) ) { 
      warning("Diagonal entries in beta_q are non-zero, they will be ignored")
      warn_q <- TRUE
    }
  }
  
  if ( is.null(all_states) ) { 
    all_states <- as.character(seq.int(ns))
  }
  
  s_seq <- seq.int(ns) - 1
  
  # Handle beta_0 component
  beta_0_tab <- data.frame(expand.grid(from = all_states, 
                                       to = all_states, 
                                       stringsAsFactors = FALSE), 
                           coef = as.vector(beta_0))
  beta_0_tab <- beta_0_tab[beta_0_tab[ ,"coef"] >= epsilon, ]
  
  # Handle beta_p component 
  beta_pp_list <- list()
  for ( from in seq.int(ns) ) { 
    for ( to in seq.int(ns) ) { 
      for ( s in seq.int(ns) ) { 
        if ( beta_p[from, to, s] != 0 ) { 
          beta_pp_list <- c(beta_pp_list, 
                            list(data.frame(from = all_states[from], 
                                            to = all_states[to], 
                                            state_1 = all_states[s], 
                                            state_2 = all_states[1], 
                                            expo_1 = 1, 
                                            expo_2 = 0, 
                                            coef = beta_p[from, to, s])))
        }
      }
    }
  }
  
  if ( length(beta_pp_list) > 0 ) { 
    beta_pp_tab <- do.call(rbind, beta_pp_list)
#     beta_pp_tab[ ,"from"] <- factor(beta_pp_tab[ ,"from"], levels = all_states)
#     beta_pp_tab[ ,"to"] <- factor(beta_pp_tab[ ,"to"], levels = all_states)
  } else { 
    beta_pp_tab <- dummy_beta_pp # defined in transition.R
  }
  
  # Handle beta_qq component 
  beta_qq_list <- list()
  for ( from in seq.int(ns) ) { 
    for ( to in seq.int(ns) ) { 
      for ( s in seq.int(ns) ) { 
        if ( beta_q[from, to, s] != 0 ) { 
          beta_qq_list <- c(beta_qq_list, list(
            data.frame(from = all_states[from], 
                       to = all_states[to], 
                       state_1 = all_states[s], 
                       state_2 = all_states[1], 
                       expo_1 = 1, 
                       expo_2 = 0, 
                       coef = beta_q[from, to, s])
          ))
        }
      }
    }
  }
  
  if ( length(beta_qq_list) > 0 ) { 
    beta_qq_tab <- do.call(rbind, beta_qq_list)
#     beta_q_tab[ ,"from"] <- factor(beta_q_tab[ ,"from"], levels = all_states)
#     beta_q_tab[ ,"to"] <- factor(beta_q_tab[ ,"to"], levels = all_states)
  } else { 
    beta_qq_tab <- dummy_beta_qq # defined in transition.R
  }
  
  # Build transitions 
  if ( build_transitions ) { 
    transitions <- list()
    for ( from in seq.int(ns) ) { 
      for ( to in seq.int(ns) ) { 
        str <- c()
        if ( abs(beta_0[from, to]) > 0 ) { 
          str <- c(as.character(beta_0[from, to]))
        }
        if ( any(abs(beta_p[from, to, ]) > 0) ) { 
          str <- c(str, format_coefs(beta_p[from, to, ], all_states, "p"))
        }
        if ( any(abs(beta_q[from, to, ]) > 0) ) { 
          str <- c(str, format_coefs(beta_q[from, to, ], all_states, "q"))
        }
        if ( length(str) > 0 ) { 
          str <- paste0("~ ", paste0(str, collapse = " + "))
          transitions <- append(transitions, 
                                list(transition(from = all_states[from], 
                                                to   = all_states[to], 
                                                prob = stats::as.formula(str))))
        }
      }
    }
  } else { 
    transitions <- list(transition(from = NA, to = NA, 
                                   prob = ~ NA))
  }
  
  camod_object <- list(
    transitions = transitions, 
    nstates = ns, 
    states = factor(all_states, levels = all_states), 
    parms = list(), 
    beta_0 = ord_by_state(beta_0_tab, all_states), 
    beta_q = dummy_beta_q, 
    beta_pp = ord_by_state(beta_pp_tab, all_states), 
    beta_pq = dummy_beta_pq, 
    beta_qq = beta_qq_tab, 
    wrap = wrap, 
    neighbors = neighbors, 
    kernel = build_neighbor_kernel(neighbors), 
    epsilon = epsilon, 
    xpoints = get_q_npoints(wrap, neighbors, fixed_neighborhood), # unused, but for completion
    max_error = rep(NA, 4), 
    max_rel_error = rep(NA, 4), 
    fixed_neighborhood = fixed_neighborhood
  )
  class(camod_object) <- c("ca_model", "list")
  
  return(camod_object)
}

is_square <- function(mat) { 
  is.matrix(mat) && diff(dim(mat)) == 0
}

is_cube <- function(cub) { 
  is.array(cub) && all(diff(dim(cub)) == 0)
}

zero_matrix <- function(s) { 
  diag(s) * 0
}

zero_array <- function(s) { 
  array(0, dim = c(s, s, s))
}

# Get a series of coefficients and states, and make a nice formula out of them
format_coefs <- function(coefs, all_states, prefix) {
  kept_coefs <- abs(coefs) > 0
  coefs_str <- paste0(coefs[kept_coefs], "*")
  coefs_str <- ifelse(coefs[kept_coefs] == 1.0, "", coefs_str) # do not show a leading one
  paste0(coefs_str, paste0(prefix, "['", all_states[kept_coefs], "']"))
}
