#
# This file contains a utility function that returns transition rates to next
# iteration given a model.
#

next_state_probs <- function(mod, mat) {

  # Convert mat to internal representation
  d <- dim(mat)
  mat <- as.integer(mat) - 1
  dim(mat) <- d

  # Convert state factors to internal representation
  fix <- function(x) {
    as.integer(factor(as.character(x), levels = mod[["states"]])) - 1
  }

  # Prepare data for internal representation, and take substeps into account
  betas <- mod[c("beta_0", "beta_q", "beta_pp", "beta_qq", "beta_pq")]
  betas <- lapply(betas, function(tab) {
    # Convert references to state to internal representation
    tab[ ,"from"] <- fix(tab[ ,"from"])
    tab[ ,"to"] <- fix(tab[ ,"to"])
    if ( "state_1" %in% colnames(tab) ) {
      tab[ ,"state_1"] <- fix(tab[ ,"state_1"])
    }
    if ( "state_2" %in% colnames(tab) ) {
      tab[ ,"state_2"] <- fix(tab[ ,"state_2"])
    }

    as.matrix(tab)
  })
  mod[c("beta_0", "beta_q", "beta_pp", "beta_qq", "beta_pq")] <- betas

  mod <- split_tabs(mod)

  get_transition_probas_cpp(mat, mod)
}
