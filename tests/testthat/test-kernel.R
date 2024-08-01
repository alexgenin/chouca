# 
# Test the kernel-building function
# 


moore_1x1 <- structure(c(1, 1, 1, 1, 0, 1, 1, 1, 1), dim = c(3L, 3L))
von_neumann_1x1 <- structure(c(0, 1, 0, 1, 0, 1, 0, 1, 0), dim = c(3L, 3L))
von_neumann_2x2 <- structure(c(0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 
                               1, 0, 0, 0, 1, 0, 0), dim = c(5L, 5L))
circular_3x3 <- structure(c(0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 
                            1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0), dim = c(7L, 7L))


# Simple kernels
test_that("Kernel-generation works correctly", { 
  expect_true({ 
    all( nb_kernel("moore", 1) == moore_1x1 ) 
  })
  expect_true({ 
    all( nb_kernel("von_neumann", 1) == von_neumann_1x1 ) 
  })

  expect_true({ 
    all( nb_kernel("von_neumann", 2) == von_neumann_2x2  )
  })
  expect_true({ 
    all( nb_kernel("von_neumann", 2) == nb_kernel("diamond", 2)  )
  })

  expect_true({ 
    all( nb_kernel("circular", 3) == circular_3x3 )
  })
})

