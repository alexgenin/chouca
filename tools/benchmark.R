# 
# Run a benchmark of chouca 
# 

library(plyr)
library(devtools)
library(ggplot2)
library(lubridate)

GIT_ORIG <- "git@github.com:alexgenin/chouca.git"
TEST_COMMITS <- c("58959efce8b69831d78103beee6d4f9eb8477718", # 2022-04-13
                  "ec1a941a73a87b8525660ed3bf855c6e0bffb798", 
                  "1b889584e49d8af12a7ce11059445b5d987a0d00", # 2022-04-30
                  "4e35aa5e408fc1ae29ffc7ba3c2803b0d9ef1510")
CTRL_LIST_VERS <- c(1, 1, 1, 2)

# Download latest chouca package in directory, compile and load it 
PKGDIR <- file.path(tempdir(), "choucabench")
dir.create(PKGDIR)
system(paste0("git clone -b master ", GIT_ORIG, " ", PKGDIR))

BENCH_SIZES <- 2^seq(4, 8)
NREPS       <- 3
CXXF <- "-g -O2 -Wall"
ENGINES <- c("cpp", "compiled") 

time_mod <- function(mod, init, control, niter, ctrl_ver) { 
  if ( ctrl_ver == 1 ) { 
    timings <- system.time({ 
      # Translate control list into old syntax
      control[["ca_engine"]] <- control[["engine"]]
      control[["engine"]] <- NULL
      control[["console_output_every"]] <- 1
      control[["console_output"]] <- FALSE
      a <- try(run_camodel(mod, init, niter = niter, control = control), 
                silent = FALSE)
    })
  } else { 
    timings <- system.time({ 
      a <- try( run_camodel(mod, init, niter = niter, control = control) )
    })
  }
  list(error = inherits(a, "try-error"), timings = timings)
}

mkbench <- function(sizes, nrep, cxxf, commit) { 
  
  # Checkout correct commit in PKGDIR
  system(sprintf("cd '%s' && git checkout %s", PKGDIR, commit))
  
  # cd into pkgdir 
  owd <- getwd()
  on.exit({ setwd(owd) })
  setwd(PKGDIR)
  
  # Load devtools, compile and load package 
  clean_dll(PKGDIR)
  withr::with_makevars(c(CXX11FLAGS = CXXF), load_all(PKGDIR))
  
  control <- list(console_output_every = 0)
  
  ldply(ENGINES, function(engine) { 
    ldply(rev(sizes), function(size) { 
      cat(sprintf("Benching engine %s with size %s\n", engine, size))
      
      # Model used for benchmark
      mod <- forestgap()
      
      # We use rectangles because we always want to test the package on rectangles
      init <- generate_initmat(mod, c(0.5, 0.5), size, size)
      tmax <- round(1000 * 1 / log(size))
      control[["engine"]] <- engine
      
      # Get the control list version 
      ctrl_ver <- CTRL_LIST_VERS[which(TEST_COMMITS == commit)]
      
      # We first run a small simulation to warm up the engine (= compile the code)
      warmup <- time_mod(mod, init, control, 10, ctrl_ver)
      
      ldply(seq.int(NREPS), function(nrep) { 
        
        # Generate new matrix for each iteration
        init <- generate_initmat(mod, c(0.5, 0.5), size, size)
        
        timings <- time_mod(mod, init, control, tmax, ctrl_ver)
        
        mcells_per_s <- (tmax * size^2) / timings[["timings"]][["elapsed"]] / 1e6
        data.frame(nrep = nrep, size = size, tmax = tmax, 
                   mcells_per_s = mcells_per_s, 
                  engine = engine, 
                  finished = ! (timings[["error"]] | warmup[["error"]] ), 
                  warmup = warmup[["timings"]]["elapsed"], 
                  elapsed = timings[["timings"]]["elapsed"])
      })
    })
   }, .progress = "time")
  
}



# Benchmark chouca 
if ( FALSE ) { 
  COMMIT_LAST <- tail(TEST_COMMITS, 1)
  bench_engines <- mkbench(BENCH_SIZES, NREPS, CXXF, COMMIT_LAST)
  
  ggplot(subset(bench_engines, finished), 
        aes(x = size, y = mcells_per_s, color = engine, shape = engine)) + 
    geom_point() + geom_line(aes(group = paste(engine, nrep))) + 
    scale_x_continuous(trans = "log", 
                      breaks = BENCH_SIZES) + 
    scale_color_brewer(palette = "Set2") + 
    labs(x = "Matrix size", 
         y = "Million cells evaluated per second")

  ggplot(subset(bench_engines, finished), 
        aes(x = size, y = tmax / elapsed / 1e3, color = engine)) + 
    geom_point() + 
    geom_line(aes(group = paste(nrep, engine))) + 
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
  geom_line(aes(group = paste(nrep, commit))) + 
  facet_grid( ~ engine ) + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  labs(x = "Matrix size", 
       y = "Million cells evaluated per second")


ggplot(subset(bench_commits, finished), 
       aes(x = size, y = tmax / elapsed / 1e3, color = engine)) + 
  geom_point() + 
  geom_line(aes(group = paste(nrep, commit, engine))) + 
#   facet_grid( ~ engine ) + 
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



