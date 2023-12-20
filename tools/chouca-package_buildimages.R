
library(chouca)

if ( ! dir.exists("vignettes") ) { 
  stop("Wrong working directory, are you in the chouca directory?")
} else { 
  setwd("vignettes")
}
on.exit({ 
  setwd("..")
})

set.seed(123)

mod <- ca_library("aridvege")
init <- generate_initmat(mod, rep(1/3, 3), nr = 64, nc = 64)
export_dir <- tempdir()

export_landscapes <- function(t, mat) {
  if ( t <= 128 ) {
    return(TRUE)
  }
  
  png(file.path(export_dir, sprintf("%04d.png", t)))
  par(mar = rep(0, 4))
  image(as.camodel_initmat(mat),
        xaxt = "n", yaxt = "n",
        col = c("white", "#F4B300", "#8A0332"))
  dev.off()
  return(TRUE)
}

run_camodel(mod, init, times = seq(0, 256),
            control = list(custom_output_every = 1,
                           custom_output_fun = export_landscapes,
                           console_output_every = 0))

# Copy one frame 
files <- tail(dir(export_dir, full = TRUE), 1)
file.copy(files, "chouca_sca_example.png", overwrite = TRUE)

# Assemble into video 
cmd <- paste0("ffmpeg -y -framerate 16 -pattern_type glob -i '", export_dir, "/*.png' ", 
              "-c:v libx264 -pix_fmt yuv420p chouca_sca_gifcreate.mp4")
system(cmd)

cmd <- paste0("ffmpeg -y -i chouca_sca_gifcreate.mp4 -s 256:256 chouca_sca_gifcreate.gif")
system(cmd)

# Upload
system({ 
  "scp chouca_sca_gifcreate.mp4 chouca_sca_gifcreate.gif alex@lecairn.org:~/www/public/files/"
})

# Cleanup 
file.remove(c("chouca_sca_gifcreate.mp4", "chouca_sca_gifcreate.gif"))







# Landscape plotting example

# Cleanup
file.remove(dir(export_dir, full.names = TRUE))

mod <- ca_library("rock-paper-scissor")

ctrl <- list(custom_output_every = 1,
             custom_output_fun = export_landscapes)

init <- generate_initmat(mod, rep(1, 3)/3, 100, 178)

out <- run_camodel(mod, init, seq(0, 512), control = ctrl)

# Copy one frame 
files <- tail(dir(export_dir, full = TRUE), 1)
file.copy(files, "chouca_landscape_plotting_example.png", overwrite = TRUE)

# Assemble into video 
cmd <- paste0("ffmpeg -y -framerate 16 -pattern_type glob -i '", export_dir, "/*.png' ", 
              "-c:v libx264 -pix_fmt yuv420p chouca_landscape_plotting_example.mp4")
system(cmd)

# Upload
system({ 
  "scp chouca_landscape_plotting_example.mp4 alex@lecairn.org:~/www/public/files/"
})

# Cleanup 
file.remove(c("chouca_landscape_plotting_example.mp4"))






# Trace plotting example

# Cleanup
file.remove(dir(export_dir, full.names = TRUE))

mod <- ca_library("rock-paper-scissor")
init <- generate_initmat(mod, rep(1, 3)/3, 100, 178)

ctrl <- list(custom_output_every = 1,
             custom_output_fun = trace_plotter(mod, init, 
                                               new_window = FALSE))

init <- generate_initmat(mod, rep(1, 3)/3, 100, 178)
pdf(file.path(export_dir, "plots.pdf"))
out <- run_camodel(mod, init, seq(0, 128), control = ctrl)
dev.off()

system(sprintf("cd %s && pdftk %s burst", export_dir, file.path(export_dir, "plots.pdf")))
pngs <- dir(export_dir, full.names = TRUE, pattern = "*.pdf")
plyr::l_ply(seq_along(pngs), function(i) { 
  system(sprintf("convert -density 96 -flatten %s %s/%04d.png", pngs[i], export_dir, i))
}, .progress = "time")

# Copy one frame 
files <- tail(dir(export_dir, full = TRUE, pattern = "*.png"), 1)
file.remove(files) # the last one is buggy for some reason
files <- tail(dir(export_dir, full = TRUE, pattern = "*.png"), 1)
file.copy(files, "chouca_trace_plotting_example.png", overwrite = TRUE)

# Assemble into video 
cmd <- paste0("ffmpeg -y -framerate 16 -pattern_type glob -i '", export_dir, "/*.png' ", 
              "-c:v libx264 -pix_fmt yuv420p chouca_trace_plotting_example.mp4")
system(cmd)


# Upload
system({ 
  "scp chouca_trace_plotting_example.mp4 alex@lecairn.org:~/www/public/files/"
})

# Cleanup 
file.remove(c("chouca_trace_plotting_example.mp4"))










