## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib chouca, .registration = TRUE
## usethis namespace: end

#'@name chouca
#'
#'@title chouca: a package for stochastic cellular automata
#'
#'@description \code{chouca} is a package that can be used to implement and run stochastic
#'  cellular automata (SCA). SCA are model that describe dynamics over a grid of cells.
#'  Each cell can be in one of a set of states, and at each time step, can switch to
#'  another one with a given transition probabilities. These transitions probabilities
#'  typically depend on the neighborhood of a focal cell.
#'
#'  The \code{chouca} package is able to simulate SCA models efficiently, including
#'  by emitting and compiling the required C++ code at runtime. It does not support all
#'  cellular automata, but typically only those where transition probabilities only depend
#'  on (1) constant parameters, (2) the proportion of cells in a given state in the
#'  landscape, and (3) the proportion of neighboring cells in a given state. More
#'  information on the types of supported model can be found in our publication, but in
#'  any case \code{chouca} is able to identify and warn you if your model is unsupported.
#'
#'  The package workflow has typically four steps (1) define the
#'  model, (2) create an initial landscape (grid of cells), (3) run the model, and
#'  (4) display/extract the results. We describe each step below, and give examples at the
#'  end of this document. A more complete description of the workflow and the package 
#'  is available in the vignette \code{vignette("chouca-package", package = "chouca")}.
#'
#'  (1) Model definition
#'
#'    Models can be defined using the \code{\link{camodel}} function. A typical call
#'    would be something looking like this, for a simple model of plants growing over
#'    space:
#'
#'    \preformatted{
#'    mod <- camodel(
#'      transition(from = "bare", to = "plant", ~ r1 * p["plant"] + r2 * q["plant"]),
#'      transition(from = "plant", to = "bare", ~ m),
#'      parms = list(r1 = 0.1, r2 = 0.1, m = 0.05),
#'      wrap = TRUE,
#'      neighbors = 8
#'    )
#'    }
#'
#'    This model defines two transitions (between the "bare" and the
#'    "plant" state and vice versa). These transitions are defined using the
#'    transition() function, which arguments define where the transition
#'    goes from, to, and an expression on how to compute its probability. In this model,
#'    the probability that the first transition occurs depends on the proportion of
#'    "plant" cells in the landscape, \code{p["plant"]}, and the proportion of neighbors
#'    in the "plant" state, \code{q["plant"]}. The model has three parameters,
#'    \code{r1}, \code{r2} and \code{m} whose value is passed through the named list 
#'    \code{parms}. The \code{neighbors} argument defines the type of neighborhood 
#'    (here 8-way, or Moore neighborhood), and we specify that the model should run over 
#'    a toric space that wraps around the edges of the grid (\code{wrap = TRUE}).
#'
#'    More information about the creation of models is available at \code{\link{camodel}}.
#'
#'  (2) Creation of the initial landscape
#'
#'    An initial grid of cell (or landscape) can be created using
#'      \code{\link{generate_initmat}}, which will fill a grid with the states provided
#'      by a model object created using \code{\link{camodel}} above, and respecting the
#'      specified proportions:
#'
#'    \code{
#'    init_grid <- generate_initmat(mod, c(bare = 0.4, plant = 0.6),
#'                                  nrow = 128, ncol = 90)
#'    }
#'
#'    Here, we create a 128x90 rectangular grid that contains 40% bare area and 60% plant,
#'      distributed randomly through space.
#'
#'    If you already have a specific grid of cells (as an R \link{matrix})
#'    you would like to use, we recommend to process it first through
#'    \code{\link{as.camodel_initmat}} so it will play nicely with the rest of the package
#'    functions.
#'
#'  (3) Running the model
#'
#'    You can feed the model and the initial landscape to \code{\link{run_camodel}},
#'    which will run the model and output the results at the time step specified by
#'    the \code{times} argument:
#'
#'    \code{
#'      out <- run_camodel(mod, init_grid, times = seq(0, 1024))
#'    }
#'
#'    \code{\link{run_camodel}} has many options to control how the simulation is run
#'    and the model outputs are saved, and are documented in the function help page.
#'
#'  (4) Extracting results
#'
#'    The results of the simulation run can be extracted from the resulting object
#'    using the `[[` operators, typically using \code{out[["output"]][["covers"]]}
#'    or \code{out[["output"]][["snapshots"]]} to extract the global proportions of
#'    cells in each state and the landscapes, respectively. These data can then be used
#'    for further analyses. Standard methods can be used to display the results (e.g.
#'    \code{plot()} or \code{image()}).
#'
#' Each step of this workflow can be adjusted as necessary, either to display the results
#' as the simulation is run (see e.g. \code{\link{landscape_plotter}}), or use different
#' simulation backends (see options in \code{\link{run_camodel}}). You can also run the
#' equivalent mean-field model of the SCA using \code{\link{run_meanfield}}, which
#' assumes that local and global proportion of states are equal.
#'
#' \code{chouca} comes with a set of pre-implemented models, which you can access using
#' \code{\link{ca_library}}.
#'
#' If you use and like \code{chouca}, we would appreciate you to cite our corresponding
#' publication:
#' 
#'  Genin A, Dupont G, Valencia D, Zucconi M, Avila-Thieme M, Navarrete
#'  S, Wieters E (2023). "Easy, fast and reproducible Stochastic Cellular
#'  Automata with 'chouca'." \doi{10.1101/2023.11.08.566206}
#'
#'@examples
#'
#' # The above example in full
#' mod <- camodel(
#'   transition(from = "bare", to = "plant", ~ r1 * p["plant"] + r2 * q["plant"]),
#'   transition(from = "plant", to = "bare", ~ m),
#'   parms = list(r1 = 0.1, r2 = 0.1, m = 0.05),
#'   wrap = TRUE,
#'   neighbors = 8
#' )
#'
#' # Display the structure of the model
#' plot(mod)
#' 
#' init_grid <- generate_initmat(mod, c(bare = 0.4, plant = 0.6),
#'                               nrow = 128, ncol = 90)
#' out <- run_camodel(mod, init_grid, times = seq(0, 128))
#' 
#' # Display results
#' plot(out)
#' image(out)
#' 
#' # Run the meanfield model (uses deSolve internally)
#' if ( requireNamespace("deSolve", quietly = TRUE) ) {
#'   out <- run_meanfield(mod, init_grid, times = seq(0, 1024))
#'   plot(out)
#' }
#'
#'@seealso camodel, generate_initmat, run_camodel, run_meanfield, ca_library
NULL
