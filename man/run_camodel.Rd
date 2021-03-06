% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{run_camodel}
\alias{run_camodel}
\title{Running a cellular automata}
\usage{
run_camodel(mod, initmat, niter, control = list())
}
\arguments{
\item{mod}{A stochastic cellular automaton model defined using \code{\link{camodel}}}

\item{initmat}{An initial matrix to use for the simulation, possibly created using 
\code{\link{generate_initmat}}}

\item{niter}{The number of iterations for which to run the model}

\item{control}{a named list with settings to alter how the simulation is run (see 
full list of settings in secion 'Details')}
}
\value{
A \code{ca_model_result} objects, which is a list with the following 
  components: 

\enumerate{ 
  
  \item \code{model} The original model used for the model run 
  
  \item \code{initmat} The initial landscape (matrix) used for the model run 
  
  \item \code{niter} The number of iterations used 
  
  \item \code{control} The control list used for the model run, containing the options
    used for the run 
  
  \item \code{output} A named list containing the model outputs. The 'covers' 
    component contains a matrix with the first column containing the time step, and the 
    other columns the proportions of cells in a given state. The 'snapshots' componentµ
    contains the landscapes recorded as matrices, with a 't' attribute indicating 
    the corresponding time step of the model run. 
}
}
\description{
Run a pre-defined stochastic cellular automaton
}
\details{
\code{run_camodel()} is the function to run a pre-defined cellular automaton. It loads 
  the model definition, and runs the simulation for a pre-defined number of 
  iterations. It will run the simulation for \code{niter} iterations, starting from 
  the initial landscape matrix \code{initmat}.

The \code{control} list must have named elements, and allows altering the 
  way the simulation is run. The possible options are the following: 
 
 \enumerate{ 
   \item \code{substeps} The number of substeps within each iteration. Probabilities 
     defined in the model definition may go above one for some models, which produces 
     results that are approximate. To avoid this problem, each time step of the 
     model can be divided into several 'substeps'. Transitions between cell 
     probabilities may occur between each substep, but with a probability divided 
     by the number of substeps. This makes sure that every probability evaluated 
     during the model run is below one. Default is 1.
   
   \item \code{save_covers_every} The period of time between which the global covers 
     of each state in the landscape is saved. Set to zero to turn off saving 
     them. Default is 1 (save covers at each iteration).
   
   \item \code{save_snapshots_every} Period of time between which a snapshot of the
     2D grid is saved. Set to zero to turn off the saving of snapshots (default
     option).
    
   \item \code{console_callback_every} Sets the number of iterations between which 
     progress report is printed on the console. Set to zero to turn off progress 
     report. Default is to print progress every ten iterations.
   
   \item \code{neighbors} The number of neighbors to use. Use the integer value 4 to 
     use a four-way (von-Neumann) neighborhood, or 8 for an 8-way (Moore) 
     neighborhood. Any other values will produce an error. Default is to use 4
     neighbors.
   
   \item \code{wrap} Set to \code{TRUE} to use a toric space in which edges wrap 
     around, i.e. that cells on a side of the landscape are considered neighbors of
     cells on the other side. 
   
   \item \code{engine} The engine to use to run the simulations. Accepted values 
     are 'r', to use the pure-R engine, 'cpp' to use the C++ engine, or 'compiled', to
     compile the model code on the fly. Default is to use the C++ engine. Note that 
     the 'compiled' engine uses its own random number generator, and for this reason
     may produce results that are different from the two other engines.
   
   \item \code{olevel} (Compiled engine only) The optimization level to use when
     compiling the model code (default, O2, O3 or Ofast). This requires compiling with 
     gcc. By default, \code{\link[Rcpp]{sourceCpp}} options are used 
     (option 'default').
   
   \item \code{unroll_loops} (Compiled engine only) Set to \code{TRUE} to unroll loops 
     or not when compiling the model. Default is \code{FALSE}. Requires compiling with
     gcc.  
   
   \item \code{precompute_probas} (Compiled engine only) Set to \code{TRUE} to 
     precompute probabilities of transitions for all possible combinations of 
     neighborhood. When working with a model with a low number of states, this 
     can increase simulation speed dramatically. Default is to do it when the 
     possible total number of neighborhood combinations is below the number of cells 
     in the landscape (so that it is effectively less work to pre-compute transitions 
     than to compute them at each iteration). 
   
   \item \code{verbose_compilation} (Compiled engine only) Set to \code{TRUE} to print 
     Rcpp messages when compiling the model. Default is \code{FALSE}. 
   
}
}
\examples{

# Run a model with default parameters
mod <- ca_library("forestgap")
im  <- generate_initmat(mod, c(0.4, 0.6), nr = 100, nc = 50)
run_camodel(mod, im, niter = 100) 

# Set some options and use the compiled engine
ctrl <- list(engine = "compiled", save_covers_every = 1, save_snapshots_every = 100)
run <- run_camodel(mod, im, niter = 200)

covers <- run[["output"]][["covers"]]
matplot(covers[ ,1], covers[ ,-1], type = "l")

}
\seealso{
camodel
}
