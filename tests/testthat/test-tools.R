# 
# Test the image() function, just by making sure the code runs
# 

mod <- ca_library("gameoflife")
initmm <- generate_initmat(mod, c(0.5, 0.5), 10, 10)
image(initmm)
image(initmm)
