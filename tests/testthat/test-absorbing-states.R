# 
# This file tests various things related to absorbing states (states we never get out 
# of)
# 
test_that("Absorbing states are correctly handled", { 

  mod <- camodel(
    transition(from = "a", to = "b", ~ 1), 
    transition(from = "a", to = "c", ~ 1), 
    transition(from = "b", to = "a", ~ 1), 
    transition(from = "b", to = "d", ~ 1), 
    wrap = TRUE, 
    neighbors = 8
  )
  
  expect_true({
    all( mod[["absorbing_states"]] == c("c", "d") ) 
  })
  
  
  
  
})
