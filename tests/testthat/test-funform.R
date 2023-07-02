# 
# These tests are there to make sure chouca can approximate functions up to degree 4 
# 


trans_exact <- transition(from = "0", to = "+", ~ 1 + exp(- p["+"] * q["+"] ))

exp_approx <- function(x, d) { 
  if ( length(x) > 1 ) { 
    return(sapply(exp_approx, x, d = d))
  }
  
  sum(sapply(seq(0, d), function(k) { 
    x^k / factorial(k)
  }))
}

mkmod <- function(tr) { 
  camodel(tr, 
          wrap = TRUE, 
          neighbors = 4, 
          check_model = "full")
}


err <- function(deg) { 
  
  m <- mkmod(transition(from = "0", to = "+", 
                        ~ 1 + exp_approx( - p["+"] * q["+"], deg)))
  
  data.frame(err = m[["max_error"]], 
             relerr = m[["max_rel_error"]])
}

allerrs <- plyr::ldply(seq.int(0, 9), function(d) { 
  plyr::ldply(seq.int(5), function(n) { 
    data.frame(nrep = n, deg = d, err(d)) 
  })
}, .progress = "none")

# library(ggplot2)
# 
# ggplot(allerrs, aes(x = deg, y = relerr)) + 
#   geom_point() + 
#   scale_y_continuous(trans = "log10")
# 





m <- mkmod(trans_exact)
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 1)))
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 2)))
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 3)))
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 4)))

options(chouca.degmax = 10)
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 5)))
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 6)))
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 7)))
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 8)))
m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 9)))
options(chouca.degmax = NULL) # unset


# Expect warning if not enough degrees to work with 
options(chouca.degmax = 1)
expect_warning({ 
  m <- mkmod(transition(from = "0", to = "+", ~ 1 + exp_approx(-p["+"]*q["+"], 5)))
})
options(chouca.degmax = NULL)



# Check the reconstruction of coefs
m <- mkmod(
  transition(from = "0", to = "+", 
             ~ 1 - 0.1*p["+"]*p["0"] + 0.2*(p["+"]*p["0"])^2 + 0.03*(p["+"]*p["0"])^3 +
                0.1*(p["+"]*p["0"])^4)
)

expect_true({ 
  abs(m[["beta_0"]][ ,"coef"] - 1) < 1e-8
})

expect_true({ 
  all(abs(m[["beta_pp"]][ ,"coef"] - c(-0.1, 0.2, 0.03, 0.10)) < 1e-8)
})




m <- mkmod(
  transition(from = "0", to = "+", 
             ~ 1 - 0.1*p["+"]*q["0"] + 0.2*(p["+"]*q["0"])^2 + 0.03*(p["+"]*q["0"])^3 +
                0.1*(p["+"]*q["0"])^4)
)

expect_true({ 
  abs(m[["beta_0"]][ ,"coef"] - 1) < 1e-8
})

expect_true({ 
  all(abs(m[["beta_pq"]][ ,"coef"] - c(-0.1, 0.2, 0.03, 0.10)) < 1e-8)
})








sin_approx <- function(x, d) { 
  
  if ( length(x) > 1 ) { 
    return( sapply(x, sin_approx, d = d) )
  }
  
  sum(sapply(seq(0, d), function(k) { 
    (-1)^k / factorial(2 * k + 1) * x^(2 * k + 1)
  }))
}

err <- function(deg) { 
  
  # We suppress warnings because we can get something saying it goes below zero 
  suppressWarnings({ 
    m <- mkmod(transition(from = "0", to = "+", 
                          ~ 1.01 + sin_approx( - p["+"] * q["+"] * pi / 1.2, deg)))
  })
  
#   m <- mkmod(transition(from = "0", to = "+", 
#                         ~ 1.01 + sin_approx( - p["+"] * q["+"] * pi / 1.2, deg)))
  
  xy <- expand.grid(x = seq(0, 1, l = 32), 
                    y = seq(0, 1, l = 32))
  xy[ ,"z"] <- with(xy, sin_approx(- x * y * pi / 1.2, 2))
  xy[ ,"pr"] <- sapply(seq.int(nrow(xy)), function(i) { 
    with(as.data.frame(m$beta_pq), { 
      sum( coef * xy[i,"x"]^expo_1 * xy[i,"y"]^expo_2 )
    })
  })
#   
#   ggplot(xy, aes(x = x, y = y, fill = z)) + 
#     geom_raster() + 
#     geom_contour(aes(z = z)) 
#   
#   ggplot(xy, aes(x = x, y = z, group = y, color = y)) + 
#     geom_line() + 
#     geom_line(aes(y = pr), col = "red")
#   
  data.frame(err = m[["max_error"]], 
             relerr = m[["max_rel_error"]])
}

x <- seq(0, 1, l = 12)
# plot(x, sin(x*pi))
# lines(x, sin_approx(x*pi, 0))


allerrs <- plyr::ldply(seq.int(0, 2), function(d) { 
  plyr::ldply(seq.int(5), function(n) { 
    data.frame(nrep = n, deg = d, err(d)) 
  })
}, .progress = "none")

# We expect to be able to reconstruct sinus mclaurin series up to degree 2 (which means 
#  polynomial of degree 5 total)
expect_true({
  all( allerrs[ ,"relerr"] < 1e-3 ) 
})

# 
# ggplot(allerrs, aes(x = deg, y = relerr)) + 
#   geom_point() + 
#   scale_y_continuous(trans = "log10")
# 


# m <- mkmod(transition(from = "0",to = "+",
#                       ~ 1.1 + sin(-p["+"]*q["+"] * pi / 1.2)))

# m <- mkmod(transition(from = "0",to = "+",
#                       ~ 1.1 + sin(-p["+"]*p["0"] * pi / 1.2)))

# m <- mkmod(transition(from = "0",to = "+",
#                       ~ 1.1 + sin(-q["+"]*q["0"] * pi / 1.2)))


# 1.1 - pq[+] * pi / 1.2 + (pi/1.2)^3/(3*2) - (pi/1.2)^5/(5*4*3*2)
#  = 1.1 + p[+]q[+] (- 2.6 + 2.99 - 1.02 )
expect_warning({ 
  m <- mkmod(transition(from = "0",to = "+", 
                        ~ 1.1 + sin_approx(-p["0"]*q["+"] * pi / 1.2, 0)))
}, regexp = "Transition probabilities may go below zero")


# We expect to recover the first two terms of the Mac Laurin series for sin
m <- mkmod(transition(from = "0",to = "+", 
                      ~ 1.1 + sin_approx(-p["+"]*q["0"] * pi / 1.2, 1)))
expect_true({ 
  all( abs(m[["beta_pq"]][["coef"]] - c(-2.617994, 2.990575)) < 1e-5 )
})

# Here things get more complicated because some non-zero coefficients become very 
# small but creep in nevertheless
# m <- mkmod(transition(from = "0",to = "+",
#                       ~ 1.1 + sin_approx(-p["+"]*q["+"] * pi / 1.2, 2)))






