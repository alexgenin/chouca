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

# Convert betas tables so that they store states as integers, and take the number
# of substeps into account. This function takes the list of betas tables, and returns
# the same list, transformed.
adjust_states_and_probs <- function(betas_list, states, substeps) {

  # Make lookup table
  lookup <- seq_along(states) - 1
  names(lookup) <- states

  purrr::map(betas_list, function(tab) {

    # Convert references to state to internal representation
    tab[ ,"from"] <- lookup[ tab[ ,"from"] ]
    tab[ ,"to"] <- lookup[ tab[ ,"to"] ]

    if ( "state_1" %in% colnames(tab) ) {
      tab[ ,"state_1"] <- lookup[ tab[ ,"state_1"] ]
    }
    if ( "state_2" %in% colnames(tab) ) {
      tab[ ,"state_2"] <- lookup[ tab[ ,"state_2"] ]
    }

    # Order factors so that it is more cache friendly, i.e. state 0 comes first in the
    # table, etc.
    if ( nrow(tab) > 0 ) {
      ord <- order(tab[ ,"from"], tab[ ,"to"])
      tab <- tab[ord, ]
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

ord_by_state <- function(tbl, states) {
  ord <- with(tbl, order(factor(from, levels = states),
                         factor(to, levels = states)))
  tbl <- tbl[ord, ]
  row.names(tbl) <- NULL
  tbl
}

# Helper function to pack transitions coefficients into tables that contain coefficients
# for all transitions
pack_betas <- function(tr_parsed) {
  beta_names <- c("beta_0", "beta_q", "beta_pp", "beta_pq", "beta_qq")

  # Sometimes the transitions are already packed in the component transitions_parsed.
  # This happens when creating a model with camodel_mat for example
  if ( all(beta_names %in% names(tr_parsed)) ) {
    betas <- tr_parsed[beta_names]
  } else {
    betas <- purrr::map(beta_names, function(tbl) {
      l <- purrr::map(tr_parsed, function(tr) {
        tbl <- tr[[tbl]]
        if ( nrow(tbl) == 0 ) {
          data.frame(from = character(0), to = character(0), tbl)
        } else {
          data.frame(from = tr[["from"]], to = tr[["to"]], tbl)
        }
      })
      purrr::list_rbind(l)
    })
  }

  names(betas) <- beta_names
  betas
}
