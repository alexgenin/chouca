#
# This files holds functions that are common to all SCA engines and help preparing data
# before it is passed to compiled code.
#

# In C/C++ code, we handle the coefficients as arrays of integers or doubles, this mean
# we have to split the tables of coefficients into arrays containing different data
# types. This function does just that, spliting beta_xx into beta_xx_index for ints, and
# beta_xx_vals for floats.
split_tabs <- function(ctrl) {

  for ( tab in c("beta_0", "beta_q", "beta_pp", "beta_qq", "beta_pq") ) {
    tabix <- ctrl[[tab]]
    tabix <- tabix[ ,intersect(colnames(tabix), c("from", "to", "state_1", "state_2",
                                                  "qs", "expo_1", "expo_2")),
                   drop = FALSE]
    ctrl[[paste0(tab, "_index")]] <- intmat(tabix)

    tabfl <- ctrl[[tab]]
    tabfl <- tabfl[ ,intersect(colnames(tabfl), "coef"), drop = FALSE]
    ctrl[[paste0(tab, "_vals")]] <- tabfl
  }

  return(ctrl)
}

# Construct the matrix of from/to states. This matrix holds TRUE when a transition
# from this state to another exists, FALSE otherwise.
make_transition_matrix <- function(states, betas_list) {

  transition_mat <- matrix(FALSE, nrow = length(states), ncol = length(states))
  colnames(transition_mat) <- rownames(transition_mat) <- states
  all_transitions <- lapply(betas_list, function(o) unique(o[ ,c("from", "to")]))
  all_transitions <- unique(do.call(rbind, all_transitions))

  # Mind the +1 because states are coded from 0 to (ns-1) at this stage
  if ( nrow(all_transitions) > 0 ) {
    for ( i in seq.int(nrow(all_transitions)) ) {
      transition_mat[ all_transitions[i, "from"] + 1,
                      all_transitions[i, "to"]   + 1 ] <- TRUE
    }
  }

  return(transition_mat)
}

# Convert betas tables so that they store states as integers, and take the number
# of substeps into account. This function takes the list of betas tables, and returns
# the same list, transformed.
adjust_states_and_probs <- function(betas_list, states, substeps) {

  fix <- function(x) {
    as.integer(factor(as.character(x), levels = states)) - 1
  }

  lapply(betas_list, function(tab) {
    # Convert references to state to internal representation
    tab[ ,"from"] <- fix(tab[ ,"from"])
    tab[ ,"to"] <- fix(tab[ ,"to"])
    if ( "state_1" %in% colnames(tab) ) {
      tab[ ,"state_1"] <- fix(tab[ ,"state_1"])
    }
    if ( "state_2" %in% colnames(tab) ) {
      tab[ ,"state_2"] <- fix(tab[ ,"state_2"])
    }

    # We divide all probabilities by the number of substeps
    tab[ ,"coef"] <- tab[ ,"coef"] / substeps

    as.matrix(tab)
  })
}

# Make sure the kernel holds a zero or FALSE value at its center
adjust_kernel <- function(kern) {

  if ( is.logical(kern) ) {
    kern[ceiling(nrow(kern)/2), ceiling(ncol(kern)/2)] <- FALSE
  } else {
    kern[ceiling(nrow(kern)/2), ceiling(ncol(kern)/2)] <- 0
  }

  return(kern)
}

# Convert mat to integers
intmat <- function(m) {
  mn <- matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
  colnames(mn) <- colnames(m)
  mn
}

