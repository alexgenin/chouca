# # 
# # This tests chouca's results against code used to run Danet's model of plants 
# # 
# # Danet, Alain, Florian Dirk Schneider, Fabien Anthelme, and Sonia Kefi. 2021. "Indirect
# # Facilitation Drives Species Composition and Stability in Drylands." Theoretical 
# # Ecology 14 (2): 189–203. https://doi.org/10.1007/s12080-020-00489-0.
# # 
# 
# danet_covers <- readRDS("~/work/2023/chouca/tests/testthat/danet_output/covers.rds")
# 
# # Alain's model
expect_warning({ 
  mod <- ca_library("aridvege-danet")
})

# 
# # Excerpt of code
# # 
# #     parms = c(z = 4, del = .1, b = .8, c = .2, g = .2, m = .2,
# #       gamma1 = .1, r = .01, f = .9, d = .1, protection_type = list("first_protect"),
# #       u = 5),
# #     init = matrix(sample.int(4, size = 100*100, replace = TRUE, prob = c(.4, .4, .1, .1)), nrow = 100, ncol = 100),
# 
# 
# initmm <- generate_initmat(mod, c(N = .4, P = .4, `0` = .1, `-` = .1), 
#                            nr = 100, nc = 100)
# 
# ctrl <- list(engine = "compiled", 
#              precompute_probas = TRUE)
# 
# chouca_simus <- plyr::ldply(seq.int(nrow(danet_covers)), function(i) { 
#   ps <- mod[["parms"]]
#   ps[["g"]] <- danet_covers[i, "g"]
#   ps[["b"]] <- danet_covers[i, "b"]
#   print(danet_covers[i,c("g", "b")])
#   mod <- update(mod, parms = ps)
#   print(danet_covers[i,c("g", "b")])
#   
#   output <- run_camodel(mod, initmm, seq.int(1000), control = ctrl)
# }, .progress = "time")
# 
# 
# 
# 
# 
