#
# Fit a polynomial function
#
#

# Maximum degree for bivariate polynomial fitting
DEGMAX <- 5
DEGINIT <- 1
ERROR_REL_MAX <- 0.001 # 0.1 % error is the limit to produce a warning

# 
# This function parses the definition of a transition, and computers the internal 
# coefficients used to describe the probability as a function of q, p, and the products 
# qp, pp, qq (see full description in the vignette vignette("chouca-package")). This is 
# very much an internal function. 
# 
parse_transition <- function(tr, state_names, parms, xpoints, epsilon, neighbors,
                             check_model) {

  if ( ! inherits(tr, "camodel_transition") ) {
    m <- paste("The transition definition has not been defined using transition().",
               "Please do so using transition(from = ..., to = ..., prob = ~ ...)")
    stop(m)
  }
  ns <- length(state_names)

  envir_expr <- environment(tr[["prob"]])
  pexpr <- as.expression( as.list(tr[["prob"]])[[2]] )
  zero <- rep(0, ns)
  names(zero) <- state_names

  # We accumulate warnings, and print them later if any arose.
  warn <- character()
  
  # Compute the probability of transition from its symbolic expression using the 
  # values of p and q
  prob_with <- function(p, q) {
    prob <- as.numeric( eval(pexpr,
                             envir = c(parms, list(p = p, q = q)),
                             enclos = envir_expr) )
    if ( is.na(prob) ) {
      m <- "NA/NaNs in model probabilities, please check your model definition."
      stop(m)
    }

    if ( prob < 0 ) {
      m <- paste0("Transition probabilities may go below zero. Problematic transition:\n",
                 "  ", tr[["from"]], " -> ", tr[["to"]])
      warn <<- c(warn, m)
    }

    prob
  }

  # Get intercept for this transition. Note that beta_0 always has one line, but
  # we use a df here because we are going to pack multiple values into a df later
  # (using plyr::ldply).
  beta0dbl <- prob_with(p = zero, q = zero)
  beta_0 <- data.frame(coef = beta0dbl)
  
  # Get coefficient/exponent table for q. Note that we allow any function of q, i.e. 
  # the probability can depend on any univariate function f(q_i) for all possible states 
  # i. This is ok because there is a finite number of values for f(q_i) as q_i can only 
  # take values 0, 1/nb, 2/nb, ... nb/nb (=1) where nb is the number of neighbors. 
  # 
  # We evaluate q at 'xpoints' points for vectors of q with zeros everywhere except one
  # state.
  beta_q <- plyr::ldply(state_names, function(s) {
    q <- zero
    xs <- seq(0, 1, length.out = xpoints)
    nqs <- seq_along(xs) - 1
    ys <- sapply(xs, function(x) {
      q[s] <- x
      prob_with(p = zero, q = q) - beta0dbl
    })

    data.frame(state_1 = s, qs = nqs, coef = ys)
  })

  # We need to fit a sparse polynomial to get alpha * q_i^k * q_j^l. Level of sparsity 
  # is obtained by cross-validation in fitprod
  beta_qq <- fitprod(function(q) {
    prob_with(p = zero, q = q) - beta0dbl -
      sum(sapply(seq_along(q), function(i) {
        q[-i] <- 0
        prob_with(p = zero, q = q) - beta0dbl
      }))
  }, state_names, epsilon, tr, parms, no_zero_exponent = TRUE)

  # We need to fit a sparse polynomial to get alpha * p_i^k * p_j^l. Level of sparsity 
  # is obtained by cross-validation in fitprod
  beta_pp <- fitprod(function(p) {
    prob_with(p = p, q = zero) - beta0dbl
  }, state_names, epsilon, tr, parms, no_zero_exponent = FALSE)
  
  
  # Fit a sparse polynomial to get alpha * p_i^k * q_j^l. Level of sparsity 
  # is obtained by cross-validation in fitprod2. The function is different because 
  # fitprod() can make a few simplifying assumptions as it uses only the p vector. 
  beta_pq <- fitprod2(function(p, q) {
    prob_with(p = p, q = q) -
      prob_with(p = zero, q = q) -
      prob_with(p = p, q = zero) +
      beta0dbl
  }, state_names, epsilon, tr, parms)

  # Print all the warnings accumulated so far
  if ( length(warn) > 0 ) {
    # Some of them may be repeated
    warn <- unique(warn)
    lapply(warn, warning, call. = FALSE)
  }

  # Subset matrices where coefficients are zero
  notzero <- function(X) abs(X) > epsilon
  beta_0   <- beta_0[notzero(beta_0[ ,"coef"]),   , drop = FALSE]
  beta_pp  <- beta_pp[notzero(beta_pp[ ,"coef"]), , drop = FALSE]
  beta_pq  <- beta_pq[notzero(beta_pq[ ,"coef"]), , drop = FALSE]
  beta_qq  <- beta_qq[notzero(beta_qq[ ,"coef"]), , drop = FALSE]

  # For beta_q, we do not subset things that are zero, because this allows doing
  # pointer/line arithmetic to directly fetch the corresponding value in the table, without
  # having to 'scan' for values.
  beta_q <- plyr::ddply(beta_q, ~ state_1, function(df) {
    if ( any(notzero(df[ ,"coef"])) ) {
      df
    } else {
      df[FALSE, ] # return zero-line df
    }
  })

  # When we use the same state on both side, we put all the exponents in
  # expo_1, and set expo_2 to zero. This is clearer and probably faster.
  beta_pp <- normalize_coefs(beta_pp)
  beta_qq <- normalize_coefs(beta_qq)

  # Check if we reconstruct the probabilities correctly
  max_error <- mean_error <- max_rel_error <- NA_real_
  if ( check_model %in% c("full", "quick") ) {

    # filter = 2 -> only keep values summing to neighbors
    # cap is the maximum number of lines in all_qss, as a full test can test some
    # time since it scales horribly with the number of states
    cap <- ifelse(check_model == "quick", 256, 0)
    all_qss <- generate_all_qs(neighbors, ns, filter = 2,
                               line_cap = cap)[ ,1:ns, drop = FALSE]

    colnames(all_qss) <- state_names
    ps <- matrix(stats::runif(ns * nrow(all_qss)), 
                 nrow = nrow(all_qss), 
                 ncol = ns)
    ps <- t(apply(ps, 1, function(X) X / sum(X)))
    colnames(ps) <- state_names

    p_refs <- plyr::laply(seq.int(nrow(all_qss)), function(i) {
      prob_with(q = all_qss[i, ] / neighbors, p = ps[i, ])
    })

    p_preds <- plyr::laply(seq.int(nrow(all_qss)), function(i) {

      # alpha component (always length one, as there is only one intercept for a
      # transition).
      beta <- ifelse(nrow(beta_0) == 0, 0, beta_0[ ,"coef"])

      # TODO: the join is very slow, try to do something about it
      qssum <- 0
      if ( nrow(beta_q) > 0 ) {
        # q component.
        qslong <- data.frame(state_1 = state_names,
                             # this gets the qs point corresponding to this proportion
                             # of neighbors.
                             qs = all_qss[i, ] / neighbors * (xpoints - 1))

        qslong <- plyr::join(qslong, beta_q, type = "left", by = names(qslong),
                             match = "first")
        # If nomatch, then coef is NA, which will count as zero here
        qssum <- sum(qslong[ ,"coef"], na.rm = TRUE)
      }

      # pp component
      beta_pp[ ,"ps1"] <- ps[i, beta_pp[ ,"state_1"]]
      beta_pp[ ,"ps2"] <- ps[i, beta_pp[ ,"state_2"]]
      ppssum <- sum(
        beta_pp[ ,"coef"] *
          beta_pp[ ,"ps1"] ^ beta_pp[ ,"expo_1"] *
          beta_pp[ ,"ps2"] ^ beta_pp[ ,"expo_2"]
      )

      # pq component
      beta_pq[ ,"ps"] <- ps[i, beta_pq[ ,"state_1"]]
      beta_pq[ ,"qs"] <- all_qss[i, beta_pq[ ,"state_2"]]
      pqssum <- sum(
        beta_pq[ ,"coef"] *
          beta_pq[ ,"ps"] ^ beta_pq[ ,"expo_1"] *
          (beta_pq[ ,"qs"] / neighbors) ^ beta_pq[ ,"expo_2"]
      )

      # qq component
      beta_qq[ ,"qs1"] <- all_qss[i, beta_qq[ ,"state_1"]]
      beta_qq[ ,"qs2"] <- all_qss[i, beta_qq[ ,"state_2"]]
      qqssum <- sum(
        beta_qq[ ,"coef"] *
          (beta_qq[ ,"qs1"] / neighbors) ^ beta_qq[ ,"expo_1"] *
          (beta_qq[ ,"qs2"] / neighbors) ^ beta_qq[ ,"expo_2"]
      )
      
      # Total probability expression is the sum of all components
      beta + qssum + ppssum + pqssum + qqssum
    })
    
    # Compute the errors in probability reconstruction
    mean_error <- mean(abs(p_refs - p_preds))
    max_error  <- max(abs(p_refs - p_preds))
    max_rel_error <- max(ifelse((p_refs - p_preds) != 0,
                                abs((p_refs - p_preds)/p_refs), 0))
    
    
    if ( max_rel_error > ERROR_REL_MAX || max_error > ERROR_REL_MAX ) {
      msg <- paste(
        "Residual error in computed probabilities",
        paste0("  max error: ", format(max_error, digits = 3)),
        paste0("  max rel error: ", format(max_rel_error, digits = 3)),
        "Problematic probability expression: ",
        paste(utils::capture.output(print(tr)), collapse = "\n"),
        sep = "\n"
      )
      warning(msg, call. = FALSE)
    }
  }

  list(from = tr[["from"]],
       to   = tr[["to"]],
       beta_0 = beta_0,
       beta_q = beta_q,
       beta_pp = beta_pp,
       beta_pq = beta_pq,
       beta_qq = beta_qq,
       prob = pexpr,
       mean_error = mean_error,
       max_error = max_error,
       max_rel_error = max_rel_error)
}






fitprod2 <- function(yfun, state_names, epsilon, tr, parms) {

  ns <- length(state_names)

  new_max_abs_error <- 1e10
  old_max_abs_error <- Inf
  deg <- DEGINIT
  final_mat <- NULL

  # If new_max_abs_error is very very small, we can bail right away, we
  # most probaby found the polynomial
  while ( new_max_abs_error > epsilon &&
          new_max_abs_error < old_max_abs_error &&
          deg <= getOption("chouca.degmax", default = DEGMAX) ) {
    old_max_abs_error <- new_max_abs_error

    vals <- expand.grid(state_1 = seq.int(ns),
                        state_2 = seq.int(ns),
                        expo_1  = seq(0, deg),
                        expo_2  = seq(0, deg))
    vals <- vals[vals[ ,"expo_1"] > 0 & vals[ ,"expo_2"] > 0, ]

    # Clamp to low degree
    # vals <- subset(vals, expo_1 + expo_2 <= deg)

    # Check the uniqueness of products
    hash <- with(vals, primes[state_1]^expo_1 * primes[ns+state_2]^expo_2)
    vals <- vals[which(! duplicated(hash)), ]
    vals <- as.matrix(vals)

    # Compute function derivative for each line of vals, check if that changes the
    # probability. If it does not. we set those coefs to zero
    p0 <- rep(1/ns, ns)
    names(p0) <- state_names
    q0 <- rep(1/ns, ns)
    names(q0) <- state_names
    y0 <- yfun(p0, q0)

    # Differentiate over p, over q, if both are zero, then the coef is zero.
    delta_yp <- lapply(seq.int(nrow(vals)), function(i) {
      p1 <- q1 <- p0
      p1[ vals[i, "state_1"] ] <- 1 - (0.01 + 1/ns )
      yp <- yfun(p1, q1)
      dp <- yp - y0

      q1[ vals[i, "state_2"] ] <- 1 - (0.01 + 1/ns )
      yq <- yfun(p0, q1)
      dq <- yq - y0

      c(dp, dq)
    })
    delta_yp <- do.call(rbind, delta_yp)

    keep_lines <- apply(delta_yp, 1, function(X) all(abs(X) > epsilon))
    vals <- vals[keep_lines, , drop = FALSE]

    if ( nrow(vals) == 0 ) {
      best_final_mat <- cbind(vals, coef = numeric(0))
      break
    }

    ps <- lapply(seq.int(ns), function(state) {
      if ( state %in% vals[ ,"state_1"] || state %in% vals[ ,"state_2"] ) {
        x <- seq(0, 1, l = 3/2 * (2 + log(nrow(vals)) )) ^ ( 1 / deg )
      } else {
        0
      }
    })

    ps <- as.matrix(do.call(expand.grid, c(ps, ps)))
    qs <- ps[ ,seq(1, ns), drop = FALSE]
    ps <- ps[ ,seq(ns+1, ncol(ps)), drop = FALSE]

  #   ps <- subset(ps, apply(ps, 1, sum) == 1)
    colnames(ps) <- colnames(qs) <- state_names
    
    test_set <- which(seq(1, nrow(ps)) %% 3 == 0)
    ps_test <- ps[test_set, , drop = FALSE]
    qs_test <- qs[test_set, , drop = FALSE]
    ps <- ps[-test_set, , drop = FALSE]
    qs <- qs[-test_set, , drop = FALSE]

    ys <- sapply(seq.int(nrow(ps)), function(i) { yfun(ps[i, ], qs[i, ]) })

    coefs <- matrix(rep(0, nrow(vals)), ncol = 1)
    colnames(coefs) <- "coef"

    # If all ys are zero, then all coefficients are zero, then no need to optimize
    # anything, we just return the coefficients
    if ( ! all( abs(ys) < .Machine$double.eps ) ) {

      mframe <- do.call(cbind, plyr::llply(seq.int(nrow(vals)), function(i) {
        ps[ , vals[i, "state_1"]] ^ vals[i, "expo_1"] *
          qs[ , vals[i, "state_2"]] ^ vals[i, "expo_2"]
      }))

      coefs[] <- stats::lm.fit(mframe, ys, tol = 1e-16)$coefficients
    }

    final_mat <- cbind(vals, coefs)

    ypred <- quick_pred_cpp(coefs, ps_test, qs_test, vals)
    ys_test <- sapply(seq.int(nrow(ps_test)), function(i) {
      yfun(ps_test[i, ], qs_test[i, ])
    })

#     plot(sort(ps[ ,2]*qs[ ,2]),
#          ys[order(ps[ ,2]*qs[ ,2])], type = "l")
#     points(ps[ ,2]*qs[ ,2], mframe %*% coefs, col = "red", pch = 20)
#     lines(sort(ps[ ,2]*qs[ ,2]),
#           (mframe %*% coefs)[order(ps[ ,2]*qs[ ,2])], col = "red", pch = 21)
#     points(ps_test[ ,2]*qs_test[ ,2], ypred, col = "darkgreen", pch = 21)
#     lines(sort(ps_test[ ,2]*qs_test[ ,2]),
#           (ys_test)[order(ps_test[ ,2]*qs_test[ ,2])], col = "darkgreen", pch = 21)


    new_max_abs_error <- max(abs(ys_test - ypred))
#
    # par(mfrow = c(1, 2))
    # plot(ps[ ,2]*qs[ ,2], ys, xlim = c(0, 1), ylim = range(c(ys, ys_test)))
    # points(ps[ ,2]*qs[ ,2], quick_pred_cpp(coefs, ps, qs, vals), col = "red")
    # plot(ps_test[ ,2]*qs_test[ ,2], ys_test,
    #      xlim = c(0, 1), ylim = range(c(ys, ys_test)))
    # points(ps_test[ ,2]*qs_test[ ,2], ypred, col = "red")
#   #   plot(ps_test[ ,2]*qs_test[ ,2], ys_test - ypred, col = "red")

#     print(deg)
#     print(new_max_abs_error)
#     if ( is.na(new_max_abs_error) ) browser()

    if ( new_max_abs_error < old_max_abs_error ) {
      best_final_mat <- final_mat
    }

    deg <- deg + 1
  }

  final_mat <- as.data.frame(best_final_mat)
  for ( col in c("state_1", "state_2") ) {
    final_mat[ ,col] <- state_names[final_mat[ ,col]]
  }

  return(final_mat)
}


fitprod <- function(yfun, state_names, epsilon, tr, parms,
                    no_zero_exponent = FALSE) {

  ns <- length(state_names)

  new_max_abs_error <- 1e10
  old_max_abs_error <- Inf
  deg <- DEGINIT
  final_mat <- NULL

  while ( new_max_abs_error > epsilon &&
          new_max_abs_error < old_max_abs_error &&
          deg <= getOption("chouca.degmax", default = DEGMAX) ) {
    old_max_abs_error <- new_max_abs_error

    vals <- expand.grid(state_1 = seq.int(ns),
                        state_2 = seq.int(ns),
                        expo_1  = seq(0, deg),
                        expo_2  = seq(0, deg))
    
    # sometimes we do not want to include zero exponents, i.e. terms of the 
    # sum that depend only on one q_i and not the product q_i q_j. This is because these
    # single terms are already captured by the f(q_i) terms. 
    if ( no_zero_exponent ) {
      vals <- vals[vals[ ,"expo_1"] > 0 & vals[ ,"expo_2"] > 0, ]
    } else {
      vals <- vals[vals[ ,"expo_1"] > 0 | vals[ ,"expo_2"] > 0, ]
    }

    # Check the uniqueness of products.
    primes <- primes[1:ns]
    hash <- with(vals, primes[state_1]^expo_1 * primes[state_2]^expo_2)
    vals <- vals[which(!duplicated(hash)), , drop = FALSE]
    vals <- as.matrix(vals)

    # Compute function derivative for each line of vals, check if that changes the
    # probability. If it does not. we set those coefs to zero
    p0 <- rep(1/ns, ns)
    names(p0) <- state_names
    p1 <- p0
    y0 <- yfun(p0)

    # Differentiate over p, over q, if both are zero, then the coef is zero.
    delta_yp <- lapply(seq.int(nrow(vals)), function(i) {
      p1 <- p0
      p1[ vals[i, "state_1"] ] <- 1 - (0.01 + 1/ns )
      dpa <- yfun(p1) - y0

      p1 <- p0
      p1[ vals[i, "state_2"] ] <- 1 - (0.01 + 1/ns )
      dpb <- yfun(p1) - y0

      c(dpa, dpb)
    })
    delta_yp <- do.call(rbind, delta_yp)

    # We keep lines which, if the two values are non-constant, have both derivatives
    # that are non zero, or if one of the values has exponent zero, keep the first
    # derivative above zero.
    keep_lines <- ifelse(vals[ ,"expo_2"] == 0,
                         abs(delta_yp[ ,1]) > epsilon,
                         abs(delta_yp[ ,1]) > epsilon & abs(delta_yp[ ,2]) > epsilon)

    vals <- vals[keep_lines, , drop = FALSE]

    if ( nrow(vals) == 0 ) {
      best_final_mat <- cbind(vals, coef = numeric(0))
      break
    }

    ps <- lapply(seq.int(ns), function(state) {
      if ( state %in% vals[ ,"state_1"] || state %in% vals[ ,"state_2"] ) {
        x <- seq(0, 1, l = 3/2 * (4 + nrow(vals) )) ^ ( 1 / deg )
      } else {
        0
      }
    })

    ps <- as.matrix(do.call(expand.grid, ps))
    colnames(ps) <- state_names
    
    # We decide on a test set to evalute the polynomial fit. We will stop adding terms 
    # when the test error increases. 
    test_set <- which(seq(1, nrow(ps)) %% 3 == 0)
    ps_test <- ps[test_set, , drop = FALSE]
    ps <- ps[-test_set, , drop = FALSE]

    # We use this instead of apply because we need to preserve column names for yfun
    ys <- sapply(seq.int(nrow(ps)), function(i) yfun(ps[i, ]))

    coefs <- matrix(rep(0, nrow(vals)), ncol = 1)
    colnames(coefs) <- "coef"

    # If all ys are zero, then all coefficients are zero, then no need to optimize
    # anything, we just return the coefficients

    if ( ! all( abs(ys) < .Machine$double.eps ) ) {
      # Create model frame
      mframe <- do.call(cbind, plyr::llply(seq.int(nrow(vals)), function(i) {
        ps[ , vals[i,"state_1"]] ^ vals[i, "expo_1"] *
          ps[ , vals[i,"state_2"]] ^ vals[i, "expo_2"]
      }))
      coefs[] <- stats::lm.fit(mframe, ys, tol = 1e-16)$coefficients

#     # if ( length(coefs) > 20 ) browser()
      # optimans <- stats::optim(coefs, ssq,
      #                          method = "BFGS",
      #                          control = list(maxit = 1000L,
      #                                         trace = TRUE,
      #                                         parscale = rep(1e-5, length(coefs)),
      #                                         fnscale  = rep(1e-16, length(coefs)),
      #                                         abstol = epsilon/10))
      #
      # coefs[] <- optimans[["par"]]
#       plot(coefs, coefs2)
    }

    final_mat <- cbind(vals, coefs)

    # Try to turn off coefficients and see if the max_abs_error changes
    ypred <- quick_pred_cpp(coefs, ps_test, ps_test, vals)
    ys_test <- sapply(seq.int(nrow(ps_test)), function(i) yfun(ps_test[i, ]))

    new_max_abs_error <- max(abs(ys_test - ypred))

    # par(mfrow = c(1, 2))
    # plot(ps[ ,2] * ps[ ,2], ys,
    #      xlim = c(0, 1), ylim = range(c(ys, ys_test)))
    # points(ps[ ,2]*ps[ ,2], quick_pred_cpp(coefs, ps, ps, vals), col = "red")
    # plot(ps_test[ ,2]*ps_test[ ,2], ys_test,
    #      xlim = c(0, 1), ylim = range(c(ys, ys_test)))
    # points(ps_test[ ,2]*ps_test[ ,2], ypred, col = "red")
    #   plot(ps_test[ ,2]*ps_test[ ,2], ys_test - ypred, col = "red")

#     print(deg)
#     print(new_max_abs_error)
#     if ( length(coefs) > 30 ) browser()

    if ( new_max_abs_error < old_max_abs_error ) {
      best_final_mat <- final_mat
    }

    deg <- deg + 1
  }

  final_mat <- as.data.frame(best_final_mat)
  for ( col in c("state_1", "state_2") ) {
    final_mat[ ,col] <- state_names[final_mat[ ,col]]
  }

  return(final_mat)
}


# This function reorganizes a matrix of coefficients (for beta_pp and beta_qq), i.e.
#  * it fixes the exponents: when we have p[state]^a * p[state]^b, we simplify this 
#      into p[state]^(a+b)
#  * it reorders the rows of states so it is more cache-friendly (hopefully)
normalize_coefs <- function(vals) {

  # Normalize coefficients
  same_state <- vals[ ,"state_1"] == vals[ ,"state_2"]
  vals[same_state, "expo_1"] <- vals[same_state, "expo_1"] + vals[same_state, "expo_2"]
  vals[same_state, "expo_2"] <- 0

  vals <- vals[order(vals[ ,"state_1"], vals[ ,"state_2"],
                     vals[ ,"expo_2"], vals[ ,"expo_1"]), ]
  vals
}


# A list of the first 1000 prime numbers. This is more than enough to handle up to 256
# states, which is the max we support in compiled code anyway.
primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379,
383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487,
491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727,
733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853,
857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977,
983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187,
1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427,
1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613,
1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867,
1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087,
2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333,
2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557,
2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789,
2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037,
3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307,
3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539,
3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769,
3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019,
4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261,
4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523,
4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793,
4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039,
5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323,
5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569,
5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827,
5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091,
6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337,
6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619,
6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871,
6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151,
7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459,
7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687,
7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919)
