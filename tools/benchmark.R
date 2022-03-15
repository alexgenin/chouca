# 
# Run a benchmark of chouca 
# 

library(plyr)
library(devtools)
library(ggplot2)
library(lubridate)


GIT_ORIG <- "git@github.com:alexgenin/chouca.git"
TEST_COMMITS <- c("c80479a4a12016226c259a65402ee22085b5a985", # 2022-03-13
                  "c4106a38479b24220950f273ebb2ad96d89cee5b", # 2022-03-13 (after 2..5)
                  "3558b9de5119e32c3bd7c8fa5f3a2cdd7cb5f64b") # 2022-03-14 

# Download latest chouca package in directory, compile and load it 
PKGDIR <- file.path(tempdir(), "choucabench")
dir.create(PKGDIR)
system(paste0("git clone -b master ", GIT_ORIG, " ", PKGDIR))

BENCH_SIZES <- 2^seq(1, 10)
NREPS       <- 3
CXXF <- "-g -O2 -Wall"
ENGINES <- c("cpp", "compiled") 

mkbench <- function(sizes, nrep, cxxf, commit) { 
  
  # Checkout correct commit in PKGDIR
  system(sprintf("cd '%s' && git checkout %s", PKGDIR, commit))
  
  # Load devtools, compile and load package 
  clean_dll(PKGDIR)
  withr::with_makevars(c(CXX11FLAGS = CXXF), load_all(PKGDIR))
  
  # Model used for benchmark
  mod <- forestgap()
  
  control <- list(console_output = FALSE)
  
  ldply(rev(sizes), function(size) { 
    cat(sprintf("Benching size %s\n", size))
    ldply(seq.int(NREPS), function(nrep) { 
      ldply(ENGINES, function(engine) { 
        # We use rectangles because we always want to test the package on rectangles
        init <- generate_initmat(mod, c(0.5, 0.5), size, size)
        tmax <- round(1000 * 1 / log(size))
        control[["ca_engine"]] <- engine
        timings <- system.time({ 
          a <- try(run_camodel(mod, init, niter = tmax, control = control), 
                   silent = FALSE)
        })
        error <- FALSE
        if ( inherits(a, "try-error") ) { 
          error <- TRUE
        }
        timings <- as.list(timings)
        mcells_per_s <- (tmax * size^2) / timings[["elapsed"]] / 1e6
        data.frame(nrep = nrep, size = size, tmax = tmax, mcells_per_s = mcells_per_s, 
                   engine = engine, finished = ! error, timings)
      })
    })
  }, .progress = "time")
}





bench_results <- ldply(TEST_COMMITS, function(commit) { 
  
  this_results <- mkbench(BENCH_SIZES, NREPS, CXXF, commit)
  
  # Get some info about the commit 
  commit_msg <- system(paste0("cd ", PKGDIR, " && git log -1 --pretty=%B | cat | head -n1"), 
                       intern = TRUE)
  commit_datetime <- as_datetime(as.numeric({ 
    system(paste0("cd ", PKGDIR, " && git log -1 --pretty=%at | cat | head -n1"), 
           intern = TRUE)
  }))
  
  data.frame(commit = commit, 
             commit_msg = commit_msg, commit_datetime = commit_datetime, 
             this_results, row.names = NULL)
})



# Print Mcells/s
bench <- mutate(bench_results, 
                mcells_per_s = (tmax * size^2) / elapsed / 1e6)

ggplot(bench_results, aes(x = size, y = mcells_per_s, 
                          color = paste(substr(commit, 1, 6), commit_msg, sep = " "))) + 
  geom_point() + 
  facet_grid( ~ engine ) + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  labs(x = "Matrix size", 
       y = "Million cells evaluated per second")



# Check the effect of compiling native with O3, or native with Ofast 
cxxfs <- c(CXXF, 
           "-O3 -mtune=native -march=native", 
           "-Ofast -mtune=native -march=native")
COMMIT_LAST <- tail(TEST_COMMITS, 1)
bench_optim <- ldply(cxxfs, function(cxxf) { 
  data.frame(cxxf = cxxf, mkbench(BENCH_SIZES, NREPS, CXXF, COMMIT_LAST))
})

ggplot(bench_optim, aes(x = size, y = mcells_per_s, color = cxxf)) + 
  geom_smooth() + geom_point() + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  labs(x = "Matrix size", 
       y = "Million cells evaluated per second")



