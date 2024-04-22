#
# This file contains functions that will initialize and run a model
#

#' @title Run a cellular automata
#'
#' @description Run a pre-defined stochastic cellular automaton
#'
#' @param mod A stochastic cellular automaton model defined using \code{\link{camodel}}
#'
#' @param initmat An initial matrix to use for the simulation, possibly created using
#'   \code{\link{generate_initmat}}
#'
#' @param times A numeric vector describing the time sequence for which output is
#'   wanted. Time will always start at zero but output will only be saved at the time
#'   steps specified in this vector.
#'
#' @param control a named list with settings to alter the way the simulation is run (see
#'   full list of settings in 'Details' section)
#'
#' @seealso camodel, generate_initmat, trace_plotter, landscape_plotter, run_meanfield
#'
#' @details
#'
#' \code{run_camodel()} is the workhorse function to run cellular automata. It runs the
#'  simulation and outputs the results at the time steps specified by the \code{times}
#'  argument, starting from the initial landscape \code{initmat} (a matrix typically
#'  created by \code{\link{generate_initmat}}).
#'
#' Note that the simulation is run for all time steps, but output is only
#' provided for the time steps specified in \code{times}.
#'
#' The \code{control} list must have named elements, and allows altering the
#'   way the simulation is run, including the live display of covers or
#'   landscapes (see \code{\link{trace_plotter}} or
#'   \code{\link{landscape_plotter}}).
#'
#' Possible options are the following:
#'
#'  \enumerate{
#'
#'    \item \code{save_covers_every} By default, global covers are saved for each time step
#'      specified in \code{times}. Setting this argument to values higher than one will
#'      skip some time steps (thinning). For example, setting it to 2 will make
#'      \code{run_camodel} save covers only every two values specified in \code{times}.
#'      Set to 0 to skip saving covers. This value must be an integer.
#'
#'    \item \code{save_snapshots_every} In the same way as covers, landscape snapshots
#'      can be saved every set number of values in \code{times}. By default, only the
#'      initial and final landscape are saved. Set to one to save the landscape for each
#'      value specified in \code{times}. Higher values will skip elements in \code{times}
#'      by the set number. Set to zero to turn off the saving of snapshots. This
#'      value must be an integer.
#'
#'    \item \code{console_output_every} Set the number of iterations between which
#'      progress report is printed on the console. Set to zero to turn off progress
#'      report. The default option is to print progress five times during the simulation.
#'
#'    \item \code{custom_output_fun} A custom function can be passed using this
#'      argument to compute something on the landscape as the simulation is being run.
#'      This function can return anything, but needs to take two arguments, the first
#'      one being the current time in the simulation (single numeric value), and the
#'      other one the current landscape (a matrix). This can be used to plot the
#'      simulation results as it is being run, see \code{\link{landscape_plotter}} and
#'      \code{\link{trace_plotter}} for such use case.
#'
#'    \item \code{custom_output_every} If \code{custom_output_fun} is specified, then
#'      it will be called for every time step specified in the \code{times} vector.
#'      Increase this value to skip some time points, in a similar way to covers
#'      and snapshots above.
#'
#'    \item \code{substeps} Stochastic CA can run into issues where the probabilities
#'      of transitions are above one. A possible solution to this is to run the model
#'      in 'substeps', i.e. an iteration is divided in several substeps, and
#'      the substeps are run subsequently with probabilities divided by this amount. For
#'      example, a model run with 4 substeps means that each iteration will be divided
#'      in 4 'sub-iterations', and probabilities of transitions are divided by 4 for
#'      each of those sub-iterations.
#'
#'    \item \code{engine} The engine used to run the simulations. Accepted values
#'      are 'cpp' to use the C++ engine, or 'compiled', to emit and compile the model
#'      code on the fly. Default is to use the C++ engine. Note that the 'compiled'
#'      engine uses its own random number generator, and for this reason may produce
#'      simulations that are different from the C++ engine (it does respect the R seed
#'      however). You will need a compiler to use the 'compiled' engine, which
#'      may require you to install
#'      \href{https://cran.r-project.org/bin/windows/Rtools/}{Rtools} on Windows
#'      systems.
#'
#'    \item \code{precompute_probas} (Compiled engine only) Set to \code{TRUE} to
#'      precompute probabilities of transitions for all possible combinations of
#'      neighborhood. When working with a model with a low number of states
#'      (typically 3 or 4), this can increase simulation speed dramatically.
#'      By default, a heuristic is used to decide whether to enable
#'      precomputation or not.
#'
#'    \item \code{verbose_compilation} (Compiled engine only) Set to \code{TRUE} to print
#'      Rcpp messages when compiling the model. Default is \code{FALSE}.
#'
#'    \item \code{force_compilation} (Compiled engine only) \code{chouca} has a
#'      cache system to avoid recompiling similar models. Set this
#'      argument to \code{TRUE} to force compilation every time the model is run.
#'
#'    \item \code{write_source} (Compiled engine only) A file name to which
#'      the C++ code used to run the model will be written (mostly for
#'      debugging purposes).
#'
#'    \item \code{cores} (Compiled engine only) The number of threads to use
#'      to run the model. This provides a moderate speedup in most cases, and
#'      is sometimes counter-productive on small landscapes. If you plan on
#'      running multiple simulations, you are probably better off parallelizing
#'      at a higher level. See also the 'Performance' section in the vignette, 
#'      accessible using the command \code{vignette("chouca-package")}.
#'      
#' }
#'
#' @returns A \code{ca_model_result} objects, which is a list with the following
#'   components:
#'
#' \enumerate{
#'
#'   \item \code{model} The original model used for the model run (see 
#'     return value of \code{\link{camodel}} for more details about these objects).
#'
#'   \item \code{initmat} The initial landscape (matrix) used for the model run, 
#'     such as what is returned by \code{\link{generate_initmat}}. 
#'
#'   \item \code{times} The time points vector at which output is saved
#'
#'   \item \code{control} The control list used for the model run, containing the options
#'     used for the run
#'
#'   \item \code{output} A named list containing the simulation outputs. The 'covers'
#'     component contains a matrix with the first column containing the time step, and the
#'     other columns the proportions of cells in a given state. The 'snapshots' component
#'     contains the landscapes recorded as matrices(\code{camodel_initmat} objects), with 
#'     a 't' attribute indicating the corresponding time step of the model run. The 
#'     'custom' component contains the results from calling a custom function provided 
#'     as 'custom_output_fun' in the control list (see examples below).
#' }
#'
#' @examples
#'
#' # Run a model with default parameters
#' mod <- ca_library("musselbed")
#' im  <- generate_initmat(mod, c(0.4, 0.6, 0), nrow = 100, ncol = 50)
#' out <- run_camodel(mod, im, times = seq(0, 100))
#' plot(out)
#' 
#' # Disable console output 
#' opts <- list(console_output_every = 0) 
#' out <- run_camodel(mod, im, times = seq(0, 100), control = opts)
#' 
#' # 
#' 
#' \donttest{
#'
#' # Run the same model with the 'compiled' engine, and save snapshots. This
#' # requires a compiler on your computer (typically available by installing
#' # 'Rtools' on Windows)
#' ctrl <- list(engine = "compiled", save_covers_every = 1, save_snapshots_every = 100)
#' run <- run_camodel(mod, im, times = seq(0, 100), control = ctrl)
#' plot(run)
#' 
#' # 
#' oldpar <- par(mfrow = c(1, 2))
#' image(run, snapshot_time = 0)
#' image(run, snapshot_time = 100)
#' par(oldpar)
#' 
#' # Disable console output
#' ctrl <- list(console_output_every = 0)
#' run <- run_camodel(mod, im, times = seq(0, 100), control = ctrl)
#' plot(run)
#'
#' # Very verbose console output (display compilation information, etc.)
#' ctrl <- list(console_output_every = 1,
#'              verbose_compilation = TRUE,
#'              engine = "compiled",
#'              force_compilation = TRUE)
#' run <- run_camodel(mod, im, times = seq(0, 100), control = ctrl)
#'
#' # Turn on or off the memoisation of transition probabilities (mind the speed
#' # difference)
#' ctrl <- list(engine = "compiled", precompute_probas = FALSE)
#' run <- run_camodel(mod, im, times = seq(0, 256), control = ctrl)
#' ctrl2 <- list(engine = "compiled", precompute_probas = TRUE)
#' run2 <- run_camodel(mod, im, times = seq(0, 256), control = ctrl2)
#'
#' # Use a custom function to compute statistics while the simulation is running
#' fun <- function(t, mat) {
#'   # Disturbed cell to mussel cell ratio
#'   ratio <- mean(mat == "DISTURB") / mean(mat == "MUSSEL")
#'   data.frame(t = t, ratio = ratio)
#' }
#' ctrl <- list(custom_output_fun = fun, custom_output_every = 1)
#'
#' run <- run_camodel(mod, im, times = seq(0, 256), control = ctrl)
#' stats <- do.call(rbind, run[["output"]][["custom"]])
#' plot(stats[ ,1], stats[ ,2], ylab = "DISTURB/MUSSEL ratio", xlab = "time", type = "l")
#'
#' }
#'@export
run_camodel <- function(mod, initmat, times,
                        control = list()) {

  ns <- mod[["nstates"]]
  states <- mod[["states"]]

  if ( ! all(levels(initmat) %in% states) ) {
    stop("States in the initial matrix do not match the model states")
  }

  if ( ! all(round(times) == times) && all(times >= 0) ) {
    stop("Time steps must be non-negative integers")
  }
  times <- as.integer(round(times))

  # Make sure the levels of the init landscape are all there, and in the right order,
  # as from now on we are going to use their integer representation
  if ( length(levels(initmat)) != length(states) ||
       ! all( levels(initmat) == states ) ) {
    finitmat <- levels(initmat)[as.integer(initmat)]
    dim(finitmat) <- dim(initmat)
    initmat <- as.camodel_initmat(finitmat, levels = states)
  }

  # Read parameters
  control <- load_control_list(control, max(times))

  # NOTE: callbacks defined below will modify things in the current environment, so that
  # this function returns the output of the simulation.

  # Handle cover-storage callback
  cover_callback <- function(t, ps, n) { }
  cover_callback_active <- is_positive(control[["save_covers_every"]])
  if ( cover_callback_active ) {
    n_exports <- length(unique(floor(times)))

    if ( control[["save_covers_every"]] > 1 ) {
      n_exports <- 1 + n_exports %/% control[["save_covers_every"]]
    }

    global_covers <- matrix(NA_real_, ncol = 1+ns, nrow = n_exports)
    colnames(global_covers) <- c("t", as.character(states))
    cur_line <- 1

    cover_callback <- function(t, ps, n) {
      global_covers[cur_line, ] <<- c(t, ps / n)
      cur_line <<- cur_line + 1
    }
  }

  # Handle snapshot-storage callback
  snapshot_callback <- function(t, m) { }
  snapshot_callback_active <- is_positive(control[["save_snapshots_every"]])
  if ( snapshot_callback_active ) {
    snapshots <- list()

    snapshot_callback <- function(t, m) {
      d <- dim(m)
      m <- factor(states[m+1], levels = states)
      dim(m) <- d
      attr(m, "t") <- t
      m <- as.camodel_initmat(m)
      snapshots <<- c(snapshots, list(m))
    }
  }

  # Handle console output callback
  console_callback <- function(iter, ps, n) { }
  console_callback_active <- is_positive(control[["console_output_every"]])

  if ( console_callback_active ) {
    last_time <- proc.time()["elapsed"]
    first_time <- last_time
    last_iter <- 0
    tmax <- max(times)
    niter <- max(times)

    console_callback <- function(iter, ps, n) {
      new_time <- proc.time()["elapsed"]

      iter_per_s <- (iter - last_iter) / (new_time - last_time)

      cover_string <- paste(states[seq.int(length(ps))],
                            format(ps/n, digits = 3, width = 4),
                            sep = ":", collapse = " ")

      perc <- paste0(format(100 * (iter / niter), digits = 1, width = 3), " %")

      speed <- ifelse(is.nan(iter_per_s) | iter_per_s < 0, "",
                      paste("[", sprintf("%0.2f", iter_per_s), " iter/s]", sep = ""))
      speed <- ifelse(iter == 0, "", speed)

      tstring <- paste0("iter = ", format(iter, width = ceiling(log10(niter))))
      cat(paste0(tstring, " (", perc, ") ", cover_string, " ", speed, "\n"))
      last_time <<- new_time
      last_iter <<- iter
    }
  }

  # Custom callback
  custom_callback <- function(t, mat) { }
  custom_callback_active <- is_positive(control[["custom_output_every"]])
  if ( custom_callback_active ) {
    custom_output <- list()
    custom_callback <- function(t, mat) {
      # Transform back to factor from internal representation
      fmat <- states[mat+1]
      dim(fmat) <- dim(mat)
      fmat <- as.camodel_initmat(fmat) # make sure the right class is attached
      this_output <- list(control[["custom_output_fun"]](t, fmat))
      if ( ! is.null(this_output) ) {
        custom_output <<- append(custom_output, this_output)
      }
    }
  }

  # Convert initmat to internal representation
  d <- dim(initmat)
  initmat <- as.integer(initmat) - 1
  dim(initmat) <- d

  # Convert state factors to internal representation
  fix <- function(x) {
    as.integer(factor(as.character(x), levels = states)) - 1
  }

  # Prepare data for internal representation, and take substeps into account
  betas <- mod[c("beta_0", "beta_q", "beta_pp", "beta_qq", "beta_pq")]
  betas <- lapply(betas, function(tab) {
    # Convert references to state to internal representation
    tab[ ,"from"] <- fix(tab[ ,"from"])
    tab[ ,"to"] <- fix(tab[ ,"to"])
    if ( "state_1" %in% colnames(tab) ) {
      tab[ ,"state_1"] <- fix(tab[ ,"state_1"])
    }
    if ( "state_2" %in% colnames(tab) ) {
      tab[ ,"state_2"] <- fix(tab[ ,"state_2"])
    }

    # We divide all probabilities by the number of substeps
    tab[ ,"coef"] <- tab[ ,"coef"] / control[["substeps"]]

    as.matrix(tab)
  })

  # Construct the matrix of from/to states
  transition_mat <- matrix(FALSE, nrow = ns, ncol = ns)
  colnames(transition_mat) <- rownames(transition_mat) <- states
  for (i  in seq_along(mod[["transitions"]]) ) {
    transition_mat[ mod[["transitions"]][[i]][["from"]],
                    mod[["transitions"]][[i]][["to"]] ] <- TRUE
  }

  # Adjust the control list to add some components
  control_list <- c(control,
                    mod[c("wrap",
                          "neighbors", "xpoints", "fixed_neighborhood")],
                    betas,
                    list(init     = initmat,
                         times    = times,
                         nstates  = ns,
                         transition_mat = transition_mat,
                         console_callback       = console_callback,
                         cover_callback         = cover_callback,
                         snapshot_callback      = snapshot_callback,
                         custom_callback        = custom_callback))

  engine <- control[["engine"]][1]
  if ( tolower(engine) %in% c("cpp", "c++") ) {
    camodel_cpp_engine_wrap(control_list)
  } else if ( tolower(engine) %in% c("compiled") ) {
    camodel_compiled_engine_wrap(control_list)
  } else {
    stop(sprintf("%s is an unknown CA engine", engine))
  }

  # Store artefacts and return result
  results <- list(model = mod,
                  initmat = initmat,
                  times = times,
                  control = control,
                  output = list())

  if ( cover_callback_active ) {
    results[["output"]][["covers"]] <- global_covers
  }

  if ( snapshot_callback_active ) {
    results[["output"]][["snapshots"]] <- snapshots
  }

  if ( custom_callback_active ) {
    results[["output"]][["custom"]] <- custom_output
  }

  class(results) <- list("ca_model_result", "list")
  return(results)
}


load_control_list <- function(l, tmax) {

  if ( length(l) > 0 && ( is.null(names(l)) || any( names(l) == "" ) ) ) {
    stop("Some elements of the control list are not named")
  }

  control_list <- list(
    save_covers_every = 1,
    save_snapshots_every = NULL,
    console_output_every = NULL,
    custom_output_every = NULL,
    custom_output_fun = NULL,
    substeps = 1,
    engine = "cpp",
    # Compiled engine options
    precompute_probas = "auto",
    verbose_compilation = FALSE,
    force_compilation = FALSE,
    write_source = NULL,
    cores = 1
  )

  for ( nm in names(l) ) {
    if ( nm %in% names(control_list) ) {
      control_list[[nm]] <- l[[nm]]
    } else {
      stop(sprintf("%s is not a CA model option", nm))
    }
  }

  if ( is.null(control_list[["save_snapshots_every"]]) ) {
    control_list[["save_snapshots_every"]] <- tmax
  }

  if ( is.null(control_list[["console_output_every"]]) ) {
    if ( max(tmax) > 4 ) {
      control_list[["console_output_every"]] <- ceiling(tmax/4)
    } else {
      control_list[["console_output_every"]] <- 0
    }
  }

  # If we passed a custom function, but did not specify the number of times we want to
  # call it, assume we do it for all values in the times vector.
  if ( is.null(control_list[["custom_output_every"]]) &&
       ! is.null(control_list[["custom_output_fun"]]) ) {
    control_list[["custom_output_every"]] <- 1
  }

  # When custom_output_fun is unspecified, custom_output_every defaults to 0
  if ( is.null(control_list[["custom_output_every"]]) ) {
    control_list[["custom_output_every"]] <- 0
  }

  if ( any(duplicated(names(l))) ) {
    stop("Duplicated elements in control list")
  }

  check_length1_integer(control_list[["save_covers_every"]], "save_covers_every", 0)
  check_length1_integer(control_list[["save_snapshots_every"]], "save_snapshots_every", 0)
  check_length1_integer(control_list[["console_output_every"]], "console_output_every", 0)
  check_length1_integer(control_list[["custom_output_every"]], "custom_output_every", 0)
  check_length1_integer(control_list[["substeps"]], "substeps", 1)
  check_length1_integer(control_list[["cores"]], "cores", 1)

  if ( ! control_list[["engine"]] %in% c("cpp", "compiled", "r") ) {
    stop(sprintf("Engine must be either 'cpp' or 'compiled'"))
  }

  if ( ! is.logical(control_list[["verbose_compilation"]]) ) {
    stop("'verbose_compilation' option must be TRUE or FALSE")
  }

  if ( ! is.logical(control_list[["force_compilation"]]) ) {
    stop("'force_compilation' option must be TRUE or FALSE")
  }

  if ( ! ( is.logical(control_list[["precompute_probas"]]) ||
          control_list[["precompute_probas"]] == "auto") ) {
    stop("precompute probas must be TRUE, FALSE or 'auto'")
  }

  if ( control_list[["custom_output_every"]] > 0 &&
       ! is.function(control_list[["custom_output_fun"]]) ) {
    stop("Custom output was turned on but no custom function was given")
  }

  if ( control_list[["custom_output_every"]] == 0 &&
       is.function(control_list[["custom_output_fun"]]) ) {
    warning("A custom output function was provided, but it will not be used")
  }
  
  # Some platforms do not support OpenMP, so we disable the support of multiple 
  # cores
  control_list[["cores"]] <- openmp_platform_check(control_list[["cores"]])
  
  control_list
}

check_length1_integer <- function(x, str, minx = 0) {
  err <- FALSE
  if ( is.null(x) || is.na(x) || length(x) != 1 || ( ! is.numeric(x) ) ) {
    err <- TRUE
  } else if ( x < minx ) {
    err <- TRUE
  }
  if ( err ) {
    msg <- sprintf("Option '%s' must be an integer >= %s", str, minx)
    stop(msg)
  }
  invisible(TRUE)
}

is_positive <- function(x) {
  if ( is.null(x) || is.na(x) || length(x) == 0 ) {
    return(FALSE)
  }
  if ( length(x) == 1 && x >= 1 ) {
    return(TRUE)
  }
  # zero
  return(FALSE)
}


openmp_platform_check <- function(ncores) { 
  # not implemented as things appear to run smoothly so far
  return(ncores) 
}
