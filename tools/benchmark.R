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
                  "c9cebaa15784cfc989f5dd5f1549d8efa09fd133")

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

mkbench <- function(sizes, nrep, cxxf, commit) { 
  
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
        tmax <- round(1000 * 1 / log(size))
        control[["engine"]] <- engine
        control[["precompute_probas"]] <- TRUE
        
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
                     model = model, 
                     mcells_per_s = mcells_per_s, 
                     engine = engine, 
                     finished = ! (timings[["error"]] | warmup[["error"]] ), 
                     warmup = warmup[["timings"]]["elapsed"], 
                     elapsed = timings[["timings"]]["elapsed"])
        })
      })
    })
  }, .progress = "time")
}



# Benchmark last commit  
if ( FALSE ) { 
  COMMIT_LAST <- tail(TEST_COMMITS, 1)
  bench_engines <- mkbench(BENCH_SIZES, NREPS, CXXF, COMMIT_LAST)
  
  ggplot(subset(bench_engines, finished), 
        aes(x = size, y = mcells_per_s, color = engine, shape = engine)) + 
    geom_point() + geom_line(aes(group = paste(engine, nrep))) + 
    scale_x_continuous(trans = "log", 
                      breaks = BENCH_SIZES) + 
    scale_color_brewer(palette = "Set2") + 
    facet_wrap( ~ model ) + 
    labs(x = "Matrix size", 
         y = "Million cells evaluated per second")
  
  ggplot(subset(bench_engines, finished), 
        aes(x = size, y = tmax / elapsed / 1e3, color = engine)) + 
    geom_point() + 
    geom_line(aes(group = paste(nrep, engine, model))) + 
  #   facet_grid( ~ engine ) + 
    scale_x_continuous(trans = "log", 
                      breaks = BENCH_SIZES) + 
    scale_y_continuous(trans = "log10") + 
    labs(x = "Matrix size", 
        y = "kIter/s")
}








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



# TODO: make a short simulation scenario for benchmarking


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



