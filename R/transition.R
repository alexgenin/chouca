# 
# Fit a polynomial function
# 
# 

# Maximum degree for bivariate polynomial fitting
DEGMAX <- 5
DEGINIT <- 1
ERROR_REL_MAX <- 0.001 # 0.1 % error is the limit to produce a warning

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
  # we use a df anyway because we are going to pack multiple values into a df later 
  # (using plyr::ldply).
  beta_0 <- data.frame(coef = prob_with(p = zero, q = zero))
  
  # Get coefficient/exponent table for q. We evaluate q at 'xpoints' points for vectors 
  # of q with zeros everywhere except one state. 
  beta_q <- plyr::ldply(state_names, function(s) { 
    q <- zero
    xs <- seq(0, 1, length.out = xpoints) 
    nqs <- seq_along(xs) - 1
    ys <- sapply(xs, function(x) { 
      q[s] <- x
      prob_with(p = zero, q = q) - prob_with(p = zero, q = zero)
    })
    
    data.frame(state_1 = s, qs = nqs, coef = ys)
  })
  
  beta_qq <- fitprod(function(q) { 
    prob_with(p = zero, q = q) - prob_with(p = zero, q = zero) - 
      sum(sapply(seq_along(q), function(i) { 
        q[-i] <- 0 
        prob_with(p = zero, q = q) - prob_with(p = zero, q = zero)
      }))
  }, state_names, epsilon, tr, parms, nosingle = TRUE)
  
  # We need to fit a sparse polynomial to get alpha * p_a^k * p_b^l
  beta_pp <- fitprod(function(p) { 
    prob_with(p = p, q = zero) - prob_with(p = zero, q = zero)
  }, state_names, epsilon, tr, parms, nosingle = FALSE)
  
  beta_pq <- fitprod2(function(p, q) { 
    prob_with(p = p, q = q) - 
      prob_with(p = zero, q = q) - 
      prob_with(p = p, q = zero) + 
      prob_with(p = zero, q = zero)
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
  beta_q   <- beta_q[notzero(beta_q[ ,"coef"]),   , drop = FALSE]
  beta_pp  <- beta_pp[notzero(beta_pp[ ,"coef"]), , drop = FALSE]
  beta_pq  <- beta_pq[notzero(beta_pq[ ,"coef"]), , drop = FALSE]
  beta_qq  <- beta_qq[notzero(beta_qq[ ,"coef"]), , drop = FALSE]
  
  # When we use the same state on both side, we put all the exponents in 
  # expo_1, and set expo_2 to zero. This is clearer and probably faster. 
  beta_pp <- normalize_coefs(beta_pp)
  beta_qq <- normalize_coefs(beta_qq)
  
  # Check if we reconstruct the probabilities correctly 
  max_error <- mean_error <- max_rel_error <- NA_real_
  if ( check_model ) { 
    all_qss <- generate_all_qs(neighbors, ns, filter = TRUE)[ ,1:ns]
    all_qss <- all_qss[apply(all_qss, 1, sum) == neighbors, ]
    colnames(all_qss) <- state_names
    ps <- matrix(stats::runif(ns*nrow(all_qss)), nrow = nrow(all_qss), ncol = ns)
    ps <- t(apply(ps, 1, function(X) X / sum(X)))
    colnames(ps) <- state_names
    
    p_refs <- plyr::laply(seq.int(nrow(all_qss)), function(i) { 
      prob_with(q = all_qss[i, ]/neighbors, p = ps[i, ])
    })
    
    p_preds <- plyr::laply(seq.int(nrow(all_qss)), function(i) { 
      
      # alpha component (always length one, as there is only one intercept for a 
      # transition).
      beta <- ifelse(nrow(beta_0) == 0, 0, beta_0[ ,"coef"])
      
      # q component. 
      qslong <- data.frame(state_1 = state_names, 
                           # this gets the qs point corresponding to this proportion 
                           # of neighbors. 
                           qs = all_qss[i, ] / neighbors * (xpoints - 1))
      
      # TODO: the join is very slow, try to do something 
      qslong <- plyr::join(qslong, beta_q, type = "left", by = names(qslong), 
                           match = "first")
      # If nomatch, then coef is NA, which should count as zero here
      qssum <- sum(qslong[ ,"coef"], na.rm = TRUE)
      
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
      
      beta + qssum + ppssum + pqssum + qqssum
    })
    
    mean_error <- mean(abs(p_refs - p_preds))
    max_error  <- max(abs(p_refs - p_preds))
    max_rel_error <- max(ifelse((p_refs - p_preds) != 0, 
                                abs((p_refs - p_preds)/p_refs), 0))
    
    if ( max_rel_error > ERROR_REL_MAX ) { 
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
    vals <- subset(vals, expo_1 > 0 & expo_2 > 0)
    
    # Clamp to low degree
    # vals <- subset(vals, expo_1 + expo_2 <= deg)
    
    # Check the uniqueness of products
    prods <- sapply(seq.int(nrow(vals)), function(i) { 
      paste(sort(c(rep(paste0("p", vals[i, "state_1"]), vals[i, "expo_1"]), 
                   rep(paste0("q", vals[i, "state_2"]), vals[i, "expo_2"]))), 
            collapse = "")
    }) 
    vals <- vals[which(! duplicated(prods)), ]
    vals <- as.matrix(vals)
    
    # Compute function derivative for each line of vals, check if that changes the 
    # probability. If it does not. we set those coefs to zero
    p0 <- rep(1/ns, ns)
    names(p0) <- state_names
    q0 <- rep(1/ns, ns)
    names(q0) <- state_names
    y0 <- yfun(p0, q0)
    
    # Differentiate over p, over q, if both are zero, then the coef is zero.
    delta_yp <- plyr::laply(seq.int(nrow(vals)), function(i) { 
      p1 <- q1 <- p0
      p1[ vals[i, "state_1"] ] <- 1 - (0.01 + 1/ns )
      yp <- yfun(p1, q1)
      dp <- yp - y0
      
      q1[ vals[i, "state_2"] ] <- 1 - (0.01 + 1/ns )
      yq <- yfun(p0, q1)
      dq <- yq - y0
      
      c(dp, dq)
    })
    
    keep_lines <- apply(delta_yp, 1, function(X) all(abs(X) > epsilon))
#     cat("Removed ", sum(!keep_lines), " coefs\n")
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
    qs <- ps[ ,seq(1, ns)]
    ps <- ps[ ,seq(ns+1, ncol(ps))]
    
  #   ps <- subset(ps, apply(ps, 1, sum) == 1)
    colnames(ps) <- colnames(qs) <- state_names
    
    test_set <- which(seq(1, nrow(ps)) %% 3 == 0)
    ps_test <- ps[test_set, ]
    qs_test <- qs[test_set, ]
    ps <- ps[-test_set, ]
    qs <- qs[-test_set, ]
    
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
                    nosingle = FALSE) { 
  
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
    
    if ( nosingle ) { 
      vals <- subset(vals, expo_1 > 0 & expo_2 > 0)
    } else { 
      vals <- subset(vals, expo_1 > 0 | expo_2 > 0)
    }
    
    # Clamp to low degree
    # vals <- subset(vals, expo_1 + expo_2 <= deg)
    
    # Check the uniqueness of products
    prods <- sapply(seq.int(nrow(vals)), function(i) { 
      paste(sort(c(rep(paste0("p", vals[i, "state_1"]), vals[i,"expo_1"]), 
                   rep(paste0("p", vals[i, "state_2"]), vals[i,"expo_2"]))), 
            collapse = "")
    }) 
    vals <- vals[which(!duplicated(prods)), , drop = FALSE]
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
    
#     cat("Removed ", sum(!keep_lines), " coefs\n")
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
    
    test_set <- which(seq(1, nrow(ps)) %% 3 == 0)
    ps_test <- ps[test_set, ]
    ps <- ps[-test_set, ]
    
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
      
#     #   if ( length(coefs) > 20 ) browser()
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
#   
    # par(mfrow = c(1, 2))
    # plot(ps[ ,2] * ps[ ,2], ys, 
    #      xlim = c(0, 1), ylim = range(c(ys, ys_test)))
    # points(ps[ ,2]*ps[ ,2], quick_pred_cpp(coefs, ps, ps, vals), col = "red")
    # plot(ps_test[ ,2]*ps_test[ ,2], ys_test, 
    #      xlim = c(0, 1), ylim = range(c(ys, ys_test)))
    # points(ps_test[ ,2]*ps_test[ ,2], ypred, col = "red")
#   #   plot(ps_test[ ,2]*ps_test[ ,2], ys_test - ypred, col = "red")
    
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
#  * it fixes the exponents: when we have p[state]^a * p[state]^b, we simplify this into 
# p[state]^(a+b) 
#  * it reorders the rows of states so it is more cache-friendly (hopefully)
normalize_coefs <- function(vals) { 
  
  # Normalize coefficients
  same_state <- vals[ ,"state_1"] == vals[ ,"state_2"] 
  vals[same_state, "expo_1"] <- vals[same_state,"expo_1"] + vals[same_state,"expo_2"]
  vals[same_state, "expo_2"] <- 0 
  
  vals <- vals[order(vals[ ,"state_1"], vals[ ,"state_2"], 
                      vals[ ,"expo_2"], vals[ ,"expo_1"]), ]
  vals
}
