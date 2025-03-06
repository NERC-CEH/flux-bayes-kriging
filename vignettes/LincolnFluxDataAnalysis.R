## ----rendering, eval=FALSE----------------------------------------------------
## library(rmarkdown)
## library(here)
## knitr::purl("vignettes/LincolnFluxDataAnalysis.Rmd",
##   output = "vignettes/LincolnFluxDataAnalysis.R")
system.time(rmarkdown::render("vignettes/LincolnFluxDataAnalysis.Rmd"))
## tools::texi2pdf("LincolnFluxDataAnalysis.tex")
## source("LincolnFluxDataAnalysis.R")
## knitr::knit("LincolnFluxDataAnalysis.Rnw")


## ----clean up and load required packages--------------------------------------
rm(list=ls(all=TRUE))
here::i_am("vignettes/LincolnFluxDataAnalysis.Rmd")
# install.packages(c("geoR"))
library(ggplot2)
library(raster)
library(geoR)


## ----set_constants------------------------------------------------------------
######### Set up prediction grid. #############
# Define size and location of grid:
# extentField <- extent(494100, 494350, 380750, 381050)
extentField <- extent(494247.6, 494291.6, 380905.1, 380945.1)
# Smaller resolution = higher computation time
res <- 3 # resolution for prediction grid, m
# Define the coordinate rotation system (crs) for the grid
crs_OSGB <- raster::crs("EPSG:27700")  # a terra crs object

# define a raster object with the above extent, resolution and crs
r <- raster(extentField, resolution = res, crs = crs_OSGB)


## ----read_chamber_data, eval=TRUE---------------------------------------------
df <- read.csv(file = here::here("data/df_flux.csv"))
names(df)
summary(df)


## ----plot_histogram-----------------------------------------------------------
ggplot(df, aes(x = flux_linear)) + geom_histogram()


## ----plot_spatially-----------------------------------------------------------
r_obs <- raster::rasterize(x = df[, 3:4], r, field = df[, 5])
plot(r_obs)
points(df[, 3:4])
# coordinates(df) <- ~ lon + lat
# r_obs <- raster::rasterize(df, r, field = "chamberID")
# plot(r_obs)
# plot(df)
title(main = "N2O flux from 16 chamber locations rasterised on 3-m grid (nmol/m2/s)")


## ----interpolateChamberFlux---------------------------------------------------

interpolate_flux <- function(df, r, pred_grid = rasterToPoints(r, spatial = TRUE), n_realisations = 10) {
        grid_coords <- rasterToPoints(r)

        # geo_Fn2o <- as.geodata(df, coords.col = 3:4, data.col = 2, covar.col = 7)
        geo_Fn2o <- as.geodata(df, data.col = 2, covar.col = 7)
        maxdist <- variog(geo_Fn2o)$max.dist
        seqphi_sparse <- seq(0, 2 * maxdist)
        # Now we do the kriging analysis:
        geo_Fn2o_bayes <- krige.bayes(
                geo_Fn2o,
                loc = grid_coords,
                model = model.control(
                        trend.d = "cte", trend.l = "cte", cov.model = "matern",
                        # model=model.control( trend.d=trend.spatial(~ndvi, geo_Fn2o),
                        #                     trend.l=trend.spatial(~ndvi, pred_grid),
                        #                     cov.model="matern",
                        kappa = 0.5, aniso.pars = NULL, lambda = 1
                ),
                prior = prior.control(
                        beta.prior = "flat",
                        sigmasq.prior = "reciprocal",
                        phi.prior = "uniform", phi.discrete = seqphi_sparse,
                        tausq.rel.prior = "fixed", tausq.rel = 0
                ),
                output = output.control(n.posterior = n_realisations, messages = FALSE)
        )

        # write mean prediction to raster
        r_flux <- setValues(r, geo_Fn2o_bayes$predictive$mean)
        # r_flux <- flip(r_flux, direction = "y")

        # We create a map of the mean interpolated fluxes
        image(geo_Fn2o_bayes, useRaster=TRUE,
                val = geo_Fn2o_bayes$predictive$mean,
                xlim = c(r_flux@extent@xmin, r_flux@extent@xmax), ylim = c(r_flux@extent@ymin - 10, r_flux@extent@ymax),
                col = terrain.colors(21), xlab = "lon", ylab = "lat"
        )
        points(geo_Fn2o, cex.max = 1, col = "black", add = T, pt.divide = "equal")
        contour(geo_Fn2o_bayes, add = T, nlev = 11)
        legend.krige(
                val = geo_Fn2o_bayes$predictive$mean,
                col = terrain.colors(21),
                x.leg = c(r_flux@extent@xmin, r_flux@extent@xmax), y.leg = c(r_flux@extent@ymin - 10, r_flux@extent@ymin - 5)
        )


        n_realisations <- 1 # number of realisations to plot
        # par(mfrow = c(2, 2)) # sqrt(n_realisations) to fit in square plot
        for (i in (1:n_realisations)) {
                image(geo_Fn2o_bayes, useRaster=TRUE,
                        val = geo_Fn2o_bayes$predictive$simulations[, i],
                        xlim = c(r_flux@extent@xmin, r_flux@extent@xmax), ylim = c(r_flux@extent@ymin - 10, r_flux@extent@ymax),
                        col = terrain.colors(21), xlab = "lon", ylab = "lat"
                )
                points(geo_Fn2o, cex.max = 0.5, col = "black", add = T, pt.divide = "equal")
                legend.krige(
                        val = geo_Fn2o_bayes$predictive$mean,
                        col = terrain.colors(21),
                        x.leg = c(r_flux@extent@xmin, r_flux@extent@xmax), y.leg = c(r_flux@extent@ymin - 10, r_flux@extent@ymin - 5)
                )
        }
        names(r_flux)[] <- "flux_interpolated"
        return(list(geo_Fn2o_bayes = geo_Fn2o_bayes, r = r_flux))
}

plot_mean_prediction <- function(df, r, pred_grid = rasterToPoints(r, spatial = TRUE), n_realisations = 10) {

}

plot_realisations <- function(df, r, pred_grid = rasterToPoints(r, spatial = TRUE), n_realisations = 10) {

}

bsp <- interpolate_flux(df, r)
plot(bsp$r)
colMeans(bsp$posterior$sample)

