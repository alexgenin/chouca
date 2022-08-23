# 
# Run a benchmark of chouca 
# 

library(plyr)
library(devtools)
library(ggplot2)
library(lubridate)

GIT_ORIG <- "git@github.com:alexgenin/chouca.git"
TEST_COMMITS <- c("fea6ff41c3cda84e138012bc6a719327a8aba56f", 
                  "06a7e291ada8d810de2be580282c5331dc983da2", 
                  "c9cebaa15784cfc989f5dd5f1549d8efa09fd133", 
                  "2994a802fda6d7b65052bbb9b911b61c223f9f11")
COMMIT_LAST <- tail(TEST_COMMITS, 1)

# Download latest chouca package in directory, compile and load it 
PKGDIR <- file.path(tempdir(), "choucabench")
dir.create(PKGDIR)
system(paste0("git clone -b master ", GIT_ORIG, " ", PKGDIR))

BENCH_SIZES <- 2^seq(4, 10)
NREPS       <- 3
CXXF <- "-O2 -Wall"
ENGINES <- c("cpp", "compiled") 
# ENGINES <- c("compiled") 
ALL_MODELS <- c("forestgap", "musselbed", "gameoflife", "rockpaperscissor")
# ALL_MODELS <- c("musselbed")

time_mod <- function(mod, init, control, niter) { 
  timings <- system.time({ 
    a <- try( run_camodel(mod, init, niter = niter, control = control) )
  })
  list(error = inherits(a, "try-error"), timings = timings)
}

mkbench <- function(sizes, nrep, cxxf, commit, tmax = "auto") { 
  
  # Checkout correct commit in PKGDIR
  system(sprintf("cd '%s' && git stash save && git checkout %s", PKGDIR, commit))
  
  # cd into pkgdir 
  owd <- getwd()
  on.exit({ setwd(owd) })
  setwd(PKGDIR)
  
  # Load devtools, compile and load package 
  clean_dll(PKGDIR)
  withr::with_makevars(c(CXX11FLAGS = CXXF), load_all(PKGDIR))
  
  control <- list(console_output_every = 0)
  ldply(ALL_MODELS, function(model) { 
    ldply(ENGINES, function(engine) { 
      ldply(rev(sizes), function(size) { 
        cat(sprintf("Benching engine %s with size %s\n", engine, size))
        
        # Model used for benchmark
        # NOTE: we use ::: because at some point ca_library() was not exported. And we 
        # use try() because the model list has changed. 
        mod <- try({ chouca:::ca_library(model) })
        if ( inherits(mod, "try-error") ) { 
          return( NULL ) 
        }
        
        # Initialize matrix and parameters
        inits <- sapply(mod$states, function(o) { 
          mean(sample(mod$states, size = 1024, replace = TRUE) == o)
        })
        init <- generate_initmat(mod, inits / sum(inits), size, size)
        if ( tmax == "auto" ) { 
          tmax <- round(1000 * 1 / log(size))
        }
        control[["engine"]] <- engine
        ldply(c(TRUE, FALSE), function(precompute_probas) { 
          control[["precompute_probas"]] <- precompute_probas
        
          # We first run a small simulation to warm up the engine (= compile the code)
          warmup <- time_mod(mod, init, control, 1)
          
          ldply(seq.int(NREPS), function(nrep) { 
            
            inits <- sapply(mod$states, function(o) { 
              mean(sample(mod$states, size = 1024, replace = TRUE) == o)
            })
            # Generate new matrix for each iteration
            init <- generate_initmat(mod, inits / sum(inits), size, size)
            
            timings <- time_mod(mod, init, control, tmax)
            
            mcells_per_s <- (tmax * size^2) / timings[["timings"]][["elapsed"]] / 1e6
            data.frame(nrep = nrep, size = size, tmax = tmax, 
                       precompute_probas = precompute_probas, 
                       model = model, 
                       mcells_per_s = mcells_per_s, 
                       engine = engine, 
                       finished = ! (timings[["error"]] | warmup[["error"]] ), 
                       warmup = warmup[["timings"]]["elapsed"], 
                       elapsed = timings[["timings"]]["elapsed"])
          })
        })
      })
    })
  }, .progress = "time")
}



# Benchmark last commit  
if ( FALSE ) { 
  bench_engines <- mkbench(BENCH_SIZES, NREPS, CXXF, COMMIT_LAST)
  
  # Take averages 
  bench_engines_summ <- subset(bench_engines, finished)
  bench_engines_summ <- ddply(bench_engines_summ, 
                              ~ engine + size + precompute_probas, 
                              function(df) { 
    if ( df[1, "engine"] == "cpp" & df[1, "precompute_probas"] == TRUE ) { 
      return(NULL) 
    }
    with(df, data.frame(mcells_per_s = mean(mcells_per_s), 
                        kiter_per_s  = mean(tmax / elapsed / 1e3), 
                        simu_type = paste0(engine[1], 
                                           ifelse(precompute_probas[1], "+memoise", ""))))
  })
  
  iter_per_s_plot <- ggplot(bench_engines_summ, 
        aes(x = size, y = mcells_per_s, color = simu_type)) + 
    geom_point() + geom_line(aes(group = simu_type)) + 
    scale_x_continuous(trans = "log", 
                      breaks = BENCH_SIZES) + 
    scale_y_continuous(trans = "log10", 
                       labels = function(X) sprintf("%sM", format(X, digits = 2))) + 
    scale_color_brewer(palette = "Set2", 
                       name = "Engine type") + 
    annotation_logticks(base = 10, 
                        sides = "l", 
                        long  = unit(0.2, "cm"), 
                        mid   = unit(0.1, "cm"), , 
                        short = unit(0.05, "cm"), 
                        alpha = 0.5) + 
    labs(x = "Matrix size", 
         y = "Cells evaluations.s⁻¹") 
  
  speed_plot <- ggplot(bench_engines_summ, 
                       aes(x = size, y = kiter_per_s, color = simu_type)) + 
    geom_point() + geom_line(aes(group = simu_type)) + 
    scale_x_continuous(trans = "log10", 
                       breaks = BENCH_SIZES) + 
    scale_y_continuous(trans = "log10", 
                       labels = function(X) sprintf("%sk", format(X, digits = 2))) + 
    scale_color_brewer(palette = "Set2", 
                       guide = "none", 
                       name = "Engine type") + 
    annotation_logticks(base = 10, 
                        sides = "l", 
                        long  = unit(0.2, "cm"), 
                        mid   = unit(0.1, "cm"), , 
                        short = unit(0.05, "cm"), 
                        alpha = 0.5) + 
    labs(x = "Matrix size", 
         y = "Iterations.s⁻¹")
  
  linuxcmd <- with(as.list(Sys.info()), paste(sysname, release))
  proctype <- system("cat /proc/cpuinfo|grep 'model name'|head -n1|cut -d':' -f2", 
                     intern = TRUE)
  
  library(patchwork)
  ggsave({ 
    speed_plot + iter_per_s_plot + 
      plot_layout(nrow = 1, widths = c(.5, .5)) + 
      plot_annotation(caption = paste(linuxcmd, proctype, sep = "\n"))
  }, width = 10, height = 3, dpi = 200, file = "./benchmarks_last_commit.png")
}


# 
# Many small simulations setting
# 

bench_engines <- mkbench(c(64, 128), 12, CXXF, COMMIT_LAST, tmax = 128)

ggplot(subset(bench_engines, finished), 
       aes(x = size, y = elapsed, color = engine)) + 
  geom_point() + 
  geom_line(aes(group = paste(nrep, engine, model, precompute_probas), 
                linetype = precompute_probas)) + 
  facet_grid( ~ model ) + 
  scale_x_continuous(trans = "log", 
                    breaks = BENCH_SIZES) + 
  scale_y_continuous(trans = "log10") + 
  labs(x = "Matrix size", 
       y = "Elapsed time")




bench_commits <- ldply(TEST_COMMITS, function(commit) { 
  
  this_results <- mkbench(BENCH_SIZES, NREPS, CXXF, commit)
  
  # Get some info about the commit 
  commit_msg <- system(paste0("cd ", PKGDIR, " && git log -1 --pretty=%B | cat | head -n1"), 
                       intern = TRUE)
  commit_datetime <- as_datetime(as.numeric({ 
    system(paste0("cd ", PKGDIR, " && git log -1 --pretty=%at | cat | head -n1"), 
           intern = TRUE)
  }))
  
  data.frame(commit = commit, 
             commit_msg = commit_msg, 
             commit_datetime = commit_datetime, 
             this_results, row.names = NULL)
})

ggplot(subset(bench_commits, finished), 
       aes(x = size, y = mcells_per_s, 
           color = paste0(year(commit_datetime), "-", month(commit_datetime), 
                          "-", day(commit_datetime), " ", 
                          substr(commit, 1, 6), " ", commit_msg))) + 
  geom_point() + 
  geom_line(aes(group = paste(nrep, commit, model))) + 
  facet_grid( engine ~ model, scales = "free_y") + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  labs(x = "Matrix size", 
       y = "Million cells evaluated per second")


ggplot(subset(bench_commits, finished), 
       aes(x = size, y = tmax / elapsed / 1e3, color = commit)) + 
  geom_point() + 
  geom_line(aes(group = paste(nrep, commit, engine, model))) + 
  facet_grid( engine ~ model, scales = "free_y") + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_y_continuous(trans = "log10") + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  
  labs(x = "Matrix size", 
       y = "kIter/s")




# Check the effect of compiling native with O3, or native with Ofast 
cxxfs <- c(CXXF) # , 
#            "-O3 -mtune=native -march=native", 
#            "-Ofast -mtune=native -march=native")
COMMIT_LAST <- tail(TEST_COMMITS, 1)
bench_optim <- ldply(cxxfs, function(cxxf) { 
  data.frame(cxxf = cxxf, mkbench(BENCH_SIZES, NREPS, CXXF, COMMIT_LAST))
})

ggplot(bench_optim, aes(x = size, y = mcells_per_s, color = engine, shape = engine)) + 
  geom_point() + geom_line(aes(group = paste(engine, nrep))) + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  labs(x = "Matrix size", 
       y = "Million cells evaluated per second")

ggplot(bench_optim, aes(x = size, y = mcells_per_s, color = engine, shape = engine)) + 
  geom_point() + geom_line(aes(group = paste(engine, nrep))) + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  labs(x = "Matrix size", 
       y = "Million cells evaluated per second")


ggplot(bench_optim, aes(x = size, y = (warmup + elapsed) / tmax, 
                        color = engine, shape = engine)) + 
  geom_point() + geom_line(aes(group = paste(engine, nrep))) + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_y_continuous(trans = "log10") + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  labs(x = "Matrix size")



