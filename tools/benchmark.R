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
                  "a07407ea3b9d0bfef5658641c4db5c4dad83a92d", # 2022-03-17 
                  "989dde0d1318319aa940f23b391cabf6ec753791", # 2022-03-19
                  "245c05834cc67394cf731f47f328aabec033ed5c", # 2022-04-05
                  "2480a01b2b253c4ae1793f63914411c12d444ed7", # 2022-04-12
                  "9fcaf03c07da9251eb39b9cc127483177ee25ae0")

# Download latest chouca package in directory, compile and load it 
PKGDIR <- file.path(tempdir(), "choucabench")
dir.create(PKGDIR)
system(paste0("git clone -b master ", GIT_ORIG, " ", PKGDIR))

BENCH_SIZES <- 2^seq(4, 10)
NREPS       <- 3
CXXF <- "-g -O2 -Wall"
ENGINES <- c("compiled", "cpp") 

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
  
  # Model used for benchmark
  mod <- musselbed()
  
  control <- list(console_output = FALSE)
  
  ldply(ENGINES, function(engine) { 
    ldply(rev(sizes), function(size) { 
      cat(sprintf("Benching engine %s with size %s\n", engine, size))
      
      # We use rectangles because we always want to test the package on rectangles
      init <- generate_initmat(mod, c(0.5, 0.4, 0.1), size, size)
      tmax <- round(1000 * 1 / log(size))
      control[["ca_engine"]] <- engine
      
      # We first run a small simulation to warm up the engine (= compile the code)
      warmup <- system.time({ 
        try( run_camodel(mod, init, niter = 10, control = control) )
      })
        
      ldply(seq.int(NREPS), function(nrep) { 
        
        # Generate new matrix for each iteration
        init <- generate_initmat(mod, c(0.5, 0.4, 0.1), size, size)
        
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
                  engine = engine, finished = ! error, warmup = warmup["elapsed"], 
                  timings)
      })
    })
   }, .progress = "time")
  
}



# Benchmark chouca 

# Check the effect of compiling native with O3, or native with Ofast 
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
           color = paste0(year(commit_datetime), "-", month(commit_datetime), " ", 
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



