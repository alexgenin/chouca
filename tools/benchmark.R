# 
# Run a benchmark of chouca 
# 

library(plyr)
library(devtools)
library(ggplot2)
library(lubridate)


GIT_ORIG <- "git@github.com:alexgenin/chouca.git"
TEST_COMMITS <- c("c80479a4a12016226c259a65402ee22085b5a985", # 2022-03-14
                  "9c9fef1b9562cd16860354cf7da6c6ba0ea89258", # 2022-03-14 (2)
                  "3b1e4c04d979808233dbda077cf1fb51f5d868c3", # 2022-03-14 (3)
                  "688135410e79f4313cd736fd5139e1f2b5ef08ad") # 2022-03-14 (4)

# Download latest chouca package in directory, compile and load it 
PKGDIR <- tempdir()
system(paste0("git clone -b master ", GIT_ORIG, " ", PKGDIR))

BENCH_SIZES <- 2^seq(1, 10)
NREPS       <- 3
CXXF <- "-g -O2 -Wall"

bench_results <- ldply(TEST_COMMITS, function(commit) { 
  
  # Checkout correct commit in PKGDIR
  system(sprintf("cd '%s' && git checkout %s", PKGDIR, commit))
  
  # Load devtools, compile and load package 
  clean_dll(PKGDIR)
  withr::with_makevars(c(CXX11FLAGS = CXXF), load_all(PKGDIR))
  
  control <- list(console_output = FALSE)
  
  # Model used for benchmark
  mod <- forestgap()
  
  this_results <- ldply(rev(BENCH_SIZES), function(size) { 
    cat(sprintf("Benching size %s\n", size))
    ldply(seq.int(NREPS), function(nrep) { 
      # We use rectangles because we always want to test the package on rectangles
      init <- generate_initmat(mod, c(0.5, 0.5), size, size)
      tmax <- round(1000 * 1 / log(size))
      timings <- system.time({ 
        run_camodel(mod, init, niter = tmax, control = control)
      })
      data.frame(nrep = nrep, size = size, tmax = tmax, as.list(timings))
    })
  }, .progress = "time")
  
  # Get some info about the commit 
  commit_msg <- system(paste0("cd ", PKGDIR, " && git log -1 --pretty=%B | cat | head -n1"), 
                       intern = TRUE)
  commit_datetime <- as_datetime(as.numeric({ 
    system(paste0("cd ", PKGDIR, " && git log -1 --pretty=%at | cat | head -n1"), 
           intern = TRUE)
  }))
  
  data.frame(commit = commit, commit_msg = commit_msg, commit_datetime = commit_datetime, 
             this_results, row.names = NULL)
})

# Print Mcells/s
bench <- mutate(bench_results, 
                mcells_per_s = (tmax * size^2) / elapsed / 1e6)

ggplot(bench, aes(x = size, y = mcells_per_s, 
                  color = paste(substr(commit, 1, 6), commit_msg, sep = " "))) + 
  geom_smooth() + geom_point() + 
  scale_x_continuous(trans = "log", 
                     breaks = BENCH_SIZES) + 
  scale_color_brewer(palette = "Set2", name = "commit") + 
  labs(x = "Matrix size", 
       y = "Million cells evaluated per seconds")



