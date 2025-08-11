# 
# This function does some data wrangling before we pass things to the c++ engine. For 
# example, it will split matrices of integers and doubles so that types are matched in 
# the cpp engine. 
# 

camodel_cpp_engine_wrap <- function(ctrl) { 
  ctrl <- split_tabs(ctrl)
  camodel_cpp_engine(ctrl)
}
