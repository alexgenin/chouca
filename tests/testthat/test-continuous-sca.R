# # 
# # This tests that continuous SCA are correct
# # 
# 
# # A continuous CA is equivalent to a discrete CA under some circumnstances 
# #  - one state out of each transition 
# #  - very small probabibilities of transition 
# #  - delta_t = 1
# #  => in this case, P ~= rate for continuous CA (because proba = 1 - exp(-dtx)), 
# #  which is the def of discrete CA 
# 
# base_rate <- 0.1
# mod <- camodel(transition("e", "a", ~ base_rate * p["a"]), 
#                transition("a", "e", ~ base_rate * 9/10), 
#                wrap = TRUE, 
#                continuous = FALSE, 
#                neighbors = 8)
# 
# init <- generate_initmat(mod, c(0.5, 0.5), nr = 1024)
# 
# control_args <- list(engine = "compiled", 
#                      console_output_every = 128, 
#                      delta_t = 1)
# 
# run <- run_camodel(mod, init, seq(1200, 1300, by = 1), control = control_args)
# 
# sims_cont <- plyr::ldply(seq.int(3), function(n) { 
#   mod <- update(mod, continuous = TRUE)
#   run <- run_camodel(mod, init, seq.int(10*1024), control = control_args)
#   data.frame(nrep = n, as.data.frame(run[["output"]][["covers"]]))
# }, .progress = "time")
# 
# sims_disc <- plyr::ldply(seq.int(3), function(n) { 
#   mod <- update(mod, continuous = FALSE)
#   run <- run_camodel(mod, init, 10*1024, control = control_args)
#   data.frame(nrep = n, as.data.frame(run[["output"]][["covers"]]))
# }, .progress = "time")
# 
# sims <- rbind(data.frame(type = "discrete",   sims_disc), 
#               data.frame(type = "continuous", sims_cont))
# 
# library(ggplot2)
# ggplot(sims, aes(x = t, y = a, color = type, group = paste(type, nrep))) + 
#   geom_line()
# 
# means <- plyr::ddply(sims, ~ t + type, plyr::summarise, a = mean(a))
# ac <- subset(means, type == "continuous")[ ,"a"]
# ad <- subset(means, type == "discrete")[ ,"a"]
# 
# expect_true({ 
#   all(abs(ac - ad) < 0.01)
# })
# 
# 
# 
# sims_cont <- plyr::ldply(10^seq(log10(2), log10(10*10*1024), 
#                                 l = 16), 
#                          function(npts) { 
#   npts <- round(npts)
#   mod <- update(mod, continuous = TRUE)
#   run <- run_camodel(mod, init, 
#                      seq(0, 10*1024, l = npts), 
#                      control = control_args)
#   data.frame(npts = npts, dt = (10*1024)/npts, 
#              as.data.frame(run[["output"]][["covers"]]))
# }, .progress = "time")
# 
# 
# library(ggplot2)
# library(plyr)
# sims_cont_summ <- plyr::ddply(sims_cont, ~ npts, function(df) df[nrow(df), ])
# sims_cont_summ <- mutate(sims_cont_summ, 
#                          pmax = 1 - exp(-0.01 * dt))
# ggplot(sims_cont_summ, aes(x = pmax, y = a)) + 
#   geom_point() + 
#   scale_x_continuous(trans = "log10")
# 
# 
# ggplot(sims_cont, aes(x = t, y = a, 
#                       color = 1 - exp(-base_rate * dt), 
#                       group = paste(npts))) + 
#   geom_line() + 
#   geom_text(aes(x = t, y = a, label = 1 - exp(-0.01 * dt)), 
#             color = "black", 
#             data = sims_cont_summ) + 
#   scale_color_distiller(palette = "Viridis")
# 
# 
# 
