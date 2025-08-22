#
# This file contains a utility function that returns transition rates to next
# iteration given a model.
#

next_state_probs <- function(model, mat) {

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

  get_transition_probas_cpp(mat, model)
}
