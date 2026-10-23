## ---------------------------------------------------------------------------
## 42_fig3_annotations.R
##
## The values printed on figure 3 (30_fig3_initial.R, run on the server), so
## that its caption can cite a results file. Same computation as script 30:
## validation_pairs.csv (script 03, the initial cube; 123 plot-years, 41 plots), plot means.
##   panel (a)  r(simulated, measured tree GPP)
##   panel (c)  r(residual = simulated - measured, measured LAI)
## The plot-year values quoted in Section 3.1 are given alongside.
##
##   Rscript 42_fig3_annotations.R      (from UVa_work/)
## Output: results/fig3_annotations.csv
## ---------------------------------------------------------------------------

d <- read.csv('validation_out/validation_pairs.csv', stringsAsFactors = FALSE)
d <- d[is.finite(d$sim) & is.finite(d$obs_tree), ]   # all 123 rows are the main set
a <- aggregate(cbind(sim, obs_tree, LAImax) ~ Plot_ID, d, mean, na.rm = TRUE)
a$resid <- a$sim - a$obs_tree
d$resid_py <- d$sim - d$obs_tree
R <- data.frame(
  quantity = c('r(sim, obs)', 'r(residual, measured LAI)'),
  plot_means = c(cor(a$sim, a$obs_tree), cor(a$LAImax, a$resid, use = 'complete.obs')),
  n_plots = nrow(a),
  plot_years = c(cor(d$sim, d$obs_tree), cor(d$LAImax, d$resid_py, use = 'complete.obs')),
  n_plot_years = nrow(d))
print(R, row.names = FALSE, digits = 3)
write.csv(R, 'results/fig3_annotations.csv', row.names = FALSE)
