interpolate_flux <- function(df, r, pred_grid = rasterToPoints(r, spatial = TRUE), n_posterior = 1000) {
  grid_coords <- rasterToPoints(r)

  gdf <- as.geodata(df, coords.col = 1:2, data.col = 3)
  maxdist <- variog(gdf)$max.dist
  seqphi_sparse <- seq(0, 2 * maxdist)
  # Now we do the kriging analysis:
  bsp <- krige.bayes(
    gdf,
    loc = grid_coords,
    prior = prior.control(
            beta.prior = "flat",
            sigmasq.prior = "reciprocal",
            phi.prior = "uniform", phi.discrete = seqphi_sparse,
            tausq.rel.prior = "fixed", tausq.rel = 0
    ),
    output = output.control(n.posterior = n_posterior, messages = FALSE)
  )
  return(bsp)
}

plot_mean_prediction <- function(bsp, r, df) {
  # write mean prediction to raster
  r_flux <- setValues(r, bsp$predictive$mean)
  names(r_flux)[] <- "flux"
  dfp <- as.data.frame(r_flux, xy = TRUE)
  p <- ggplot(dfp, aes(x, y, z = flux))
  p <- p + scale_fill_viridis()
  p <- p + geom_tile(aes(fill = flux))
  p <- p + geom_contour(colour = "pink")
  p <- p + geom_point(data = df, aes(x = easting, y = northing, size = flux), colour = "red")
  p <- p + geom_text(data = df, aes(x = easting, y = northing,
    label = round(flux, 0)), size = 4)
  p
  return(p)
}

plot_sigma_prediction <- function(bsp, r, df) {
  # plot the standard deviations (square root of kriging variance)
  r_flux <- setValues(r, sqrt(bsp$predictive$variance))
  names(r_flux)[] <- "flux"
  dfp <- as.data.frame(r_flux, xy = TRUE)
  p <- ggplot(dfp, aes(x, y, z = flux))
  p <- p + scale_fill_viridis()
  p <- p + geom_tile(aes(fill = flux))
  p <- p + geom_contour(colour = "pink")
  p <- p + geom_point(data = df, aes(x = easting, y = northing, size = flux), colour = "red")
  p
  return(p)
}

plot_realisations <- function(bsp, r, df, n_realisations = 4) {
  s_flux <- stack()
  for (i in seq_len(n_realisations)) {
    # write one sample prediction to raster
    r_flux <- setValues(r, bsp$predictive$simulations[, i])
    names(r_flux)[] <- "sample"
    s_flux <- stack(s_flux, r_flux)
  }
  dfp <- as.data.frame(s_flux, xy = TRUE, long = TRUE)
  dfp <- setNames(dfp, c("x", "y", "realisation", "flux"))
  p <- ggplot(dfp, aes(x, y, z = flux))
  p <- p + scale_fill_viridis()
  p <- p + geom_tile(aes(fill = flux))
  p <- p + geom_contour(colour = "pink")
  p <- p + geom_point(data = df, aes(x = easting, y = northing, size = flux), colour = "red")
  p <- p + facet_wrap(~ realisation)
  p
  return(p)
}
