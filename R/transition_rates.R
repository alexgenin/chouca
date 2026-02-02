#
# This file contains a utility function that returns transition rates to next
# iteration given a model.
#

next_state_probs <- function(model, mat, log = FALSE) {

  # Convert mat to internal representation
  d <- dim(mat)
  mat <- as.integer(mat) - 1
  dim(mat) <- d

  # Prepare data for internal representation. Here we always use substeps = 1 because
  # we want the raw transition probabilities
  betas <- pack_betas(model[["transitions_parsed"]])
  betas <- adjust_states_and_probs(betas, model[["states"]], substeps = 1)
  model[names(betas)] <- betas

  # Split tabs for cpp code that cannot handle mixed-type values in a single table
  model <- split_tabs(model)

  # Make the transition matrix (ns*ns matrix with TRUE when there is a transition).
  # This improves performance most of the time because transition matrices are often
  # quite sparse.
  # TODO: factorize code with run_camodel()
  transition_mat <- matrix(FALSE, model[["nstates"]], model[["nstates"]])
  diag(transition_mat) <- TRUE
  colnames(transition_mat) <- rownames(transition_mat) <- model[["states"]]
  for ( tr in model[["transitions"]] ) {
    transition_mat[tr[["from"]], tr[["to"]]] <- TRUE
  }
  model[["transition_matrix"]] <- transition_mat

  # Handle the normalization function argument before passing to cpp
  if ( model[["normfun"]] == "identity" ) {
    model[["normfun"]] <- 0
  } else if ( model[["normfun"]] == "softmax" ) {
    model[["normfun"]] <- 1
  }

  get_transition_probas_cpp(mat, model, return_log = log)
}
