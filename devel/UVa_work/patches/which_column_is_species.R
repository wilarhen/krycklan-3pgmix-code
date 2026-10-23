#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Which column holds dominant tree species: DTS or SPM?
#
# Peichl's workbook labels DTS "Dominant tree species" and SPM "Soil parent
# material". Their marginal counts match Peichl et al. (2023) CROSSWISE
# (DTS 16/34 vs the paper's soil 15/35; SPM 27/23 vs the paper's species
# 28/22), and SPM shows no elevation separation, which soil must show here.
# This settles it against the NMD land cover map, which classifies species
# and knows nothing about soil.
#
# Reads only a window around the plots, not the national raster.
#
# Usage:
#   Rscript which_column_is_species.R                  # default T7 path
#   Rscript which_column_is_species.R /some/other/nmd.tif
# ---------------------------------------------------------------------------

suppressMessages({library(readxl); library(terra)})
t0 <- Sys.time()

XL  <- path.expand("~/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx")
BP  <- "/run/media/wilar/T7 Shield/slu_data_bases/3PG"
DEF <- file.path(BP, "nationalndcoverdata", "nmd2018bas_ogeneraliserad_v1_0.tif")
OUT <- path.expand("~/Documents/SLU/3PG/devel/validation_data/nmd_krycklan_window.tif")

args <- commandArgs(trailingOnly = TRUE)
nmd_path <- if (length(args) >= 1 && nzchar(args[1])) path.expand(args[1]) else DEF

# reuse a previously saved window if it exists
if (file.exists(OUT)) {
  cat("using the saved window:", OUT, "\n\n"); nmd_path <- OUT; saved <- TRUE
} else {
  saved <- FALSE
  if (!file.exists(nmd_path)) {
    cat("Raster not found:\n  ", nmd_path, "\n")
    if (!dir.exists(BP)) cat("The T7 Shield does not appear to be mounted.\n")
    quit(status = 1)
  }
  cat("NMD raster:", nmd_path, "\n\n")
}

# ---- plot attributes ------------------------------------------------------
x <- suppressMessages(read_excel(XL, sheet = "50 main plots", skip = 2,
                                 .name_repair = "unique"))
x <- x[!is.na(x$Plot_ID), ]
u <- x[!duplicated(x$Plot_ID), c("Plot_ID", "DTS", "SPM", "LAImax", "Elevation")]
u$Plot_ID <- trimws(as.character(u$Plot_ID))

cc <- suppressMessages(read_excel(XL, sheet = "plot coordinates"))
names(cc) <- trimws(names(cc))
id_col <- grep("plot",      names(cc), ignore.case = TRUE, value = TRUE)[1]
e_col  <- grep("^E_|east",  names(cc), ignore.case = TRUE, value = TRUE)[1]
n_col  <- grep("^N_|north", names(cc), ignore.case = TRUE, value = TRUE)[1]
if (any(is.na(c(id_col, e_col, n_col)))) {
  cat("Could not identify columns. Sheet has:\n"); print(names(cc)); quit(status = 1)
}
co <- data.frame(Plot_ID = trimws(as.character(cc[[id_col]])),
                 E = suppressWarnings(as.numeric(cc[[e_col]])),
                 N = suppressWarnings(as.numeric(cc[[n_col]])))
co <- co[!is.na(co$E) & !is.na(co$N), ]
d <- merge(u, co, by = "Plot_ID")
cat("plots with attributes and coordinates:", nrow(d), "of", nrow(u), "\n")
if (nrow(d) < 10) { cat("Merge failed.\n"); print(head(u$Plot_ID)); print(head(co$Plot_ID)); quit(status = 1) }

# ---- window the raster to the plots --------------------------------------
pts <- vect(d[, c("E", "N")], geom = c("E", "N"), crs = "EPSG:3006")
e <- ext(pts); pad <- 1000
win <- ext(xmin(e) - pad, xmax(e) + pad, ymin(e) - pad, ymax(e) + pad)
cat("window:", round(xmin(win)), round(xmax(win)), round(ymin(win)), round(ymax(win)), "\n")

r <- try(rast(nmd_path, win = win), silent = TRUE)
if (inherits(r, "try-error")) {
  cat("win= not supported, cropping instead (slower)\n")
  r <- crop(rast(nmd_path), win)
}
cat("window read in", round(difftime(Sys.time(), t0, units = "secs")), "s;",
    ncol(r), "x", nrow(r), "cells\n")

if (!saved) {
  ok <- try(writeRaster(r, OUT, overwrite = TRUE, datatype = "INT2U"), silent = TRUE)
  if (!inherits(ok, "try-error")) cat("window saved to", OUT, "- reruns are instant\n")
}
cat("\n")

if (!identical(crs(r, proj = TRUE), crs(pts, proj = TRUE))) pts <- project(pts, crs(r))
d$nmd <- as.integer(terra::extract(r, pts)[, 2])
d$nmd_f <- ifelse(!is.na(d$nmd) & d$nmd >= 121 & d$nmd <= 127, d$nmd - 10, d$nmd)

cat("NMD classes at the plots:\n"); print(table(d$nmd_f, useNA = "ifany"))
cat("  111 pine | 112 spruce | 113,114 mixed | 115-117 broadleaf | 118 non-forest\n\n")

d$nmd_sp <- ifelse(d$nmd_f == 111, "Pine", ifelse(d$nmd_f == 112, "Spruce", NA))
cmp <- d[!is.na(d$nmd_sp), ]
cat("plots the map calls pure pine or pure spruce:", nrow(cmp), "\n")
if (nrow(cmp) < 5) {
  cat("\nToo few for a verdict. Cross-tabs against all classes:\n")
  cat("\nDTS:\n"); print(table(NMD = d$nmd_f, DTS = d$DTS))
  cat("\nSPM:\n"); print(table(NMD = d$nmd_f, SPM = d$SPM)); quit(status = 0)
}

score <- function(col, nm) {
  lab <- ifelse(cmp[[col]] == 1, "Pine", "Spruce")
  tb <- table(map = cmp$nmd_sp, column = lab)
  if (nrow(tb) < 2 || ncol(tb) < 2) { cat("\n--", nm, "-- degenerate\n"); print(tb); return(NA) }
  a <- sum(diag(tb)) / sum(tb); s <- sum(tb[cbind(c(1, 2), c(2, 1))]) / sum(tb)
  cat("\n--", nm, "as labelled (1=Pine, 2=Spruce) --\n"); print(tb)
  cat(sprintf("  agreement: %.1f%%   with 1/2 reversed: %.1f%%\n", 100 * a, 100 * s))
  max(a, s)
}

a_dts <- score("DTS", "DTS"); a_spm <- score("SPM", "SPM")

cat("\n=== VERDICT ===\n")
cat(sprintf("  DTS best: %.1f%%   SPM best: %.1f%%\n", 100 * a_dts, 100 * a_spm))
if (is.na(a_dts) || is.na(a_spm) || abs(a_dts - a_spm) < 0.10) {
  cat("  Too close to call.\n")
} else if (a_spm > a_dts) {
  cat("  SPM tracks the species map -> LABELS ARE SWAPPED.\n")
  cat("  The analysis has been grouping by soil, not species.\n")
} else {
  cat("  DTS tracks the species map -> labels are correct as given.\n")
  cat("  The paper's marginal counts then need a separate explanation.\n")
}
cat("\nNMD assigns a species class only above 70% crown cover, so neither\n")
cat("will score near 100%. What matters is which scores higher.\n")
cat("\ntotal:", round(difftime(Sys.time(), t0, units = "secs")), "s\n")
