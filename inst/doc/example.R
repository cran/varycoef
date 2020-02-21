## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, 
                      fig.width=7, fig.height=4)

## ----install------------------------------------------------------------------
# install from CRAN
# install.packages("varycoef")

# attach package
library(varycoef)

# general package help file (check out examples with manual vignette!!)
help("varycoef")

## ----SVC example--------------------------------------------------------------
# setting seed
set.seed(123)
# number of SVC
p <- 3
# sqrt of total number of observations
m <- 20
# covariance parameters
(pars <- data.frame(var = c(0.1, 0.2, 0.3), 
                    scale = c(0.3, 0.1, 0.2)))
nugget.var <- 0.05

# function to sample SVCs
sp.SVC <- fullSVC_reggrid(m = m, p = p, 
                          cov_pars = pars, 
                          nugget = nugget.var)

library(sp)
# visualization of sampled SVC
spplot(sp.SVC, colorkey = TRUE)

## ----meuse intro--------------------------------------------------------------
# attach sp and load data
library(sp)
data("meuse")

# documentation
help("meuse")

# overview
summary(meuse)
dim(meuse)

## ----hist plot----------------------------------------------------------------
par(mfrow = 1:2)
# histogram of (log) cadmium 
hist(meuse$cadmium); hist(log(meuse$cadmium))
par(mfrow = c(1, 1))

meuse$l_cad <- log(meuse$cadmium)

## ----spatial plot-------------------------------------------------------------
# construct spatial object
sp.meuse <- meuse
coordinates(sp.meuse) <- ~x+y
proj4string(sp.meuse) <- CRS("+init=epsg:28992")

# using package tmap
library(tmap)
# producing an interactive map
tmap_leaflet(tm_shape(sp.meuse) + tm_dots("l_cad", style = "cont"))

## ----linear regression--------------------------------------------------------
# linear model (LM)
lm.fit <- lm(l_cad ~ 1+dist+lime+elev, data = meuse)
summary(lm.fit)

## ----LM residuals-------------------------------------------------------------
oldpar <- par(mfrow = c(1, 2))
plot(lm.fit, which = 1:2)
par(oldpar)

## ----LM spatial residuals-----------------------------------------------------
# add LM residuals to data frame
sp.meuse$LM_res <- resid(lm.fit)
head(sp.meuse)
# plot residuals at corresponding locations
tmap_leaflet(tm_shape(sp.meuse) + tm_dots("LM_res", style = "cont"))

## ----geostatistical-----------------------------------------------------------
library(gstat)
# empirical variogram
eV <- variogram(LM_res ~ 1, sp.meuse)
# define variogram model with initial values
mV <- vgm(0.2, "Exp", 300, 0.4)
# fit model 
(fV <- fit.variogram( eV, mV))
# plot empirical and fitted
plot(eV, model=fV)

## ----ordinary kriging---------------------------------------------------------
# study area
data("meuse.grid")
coordinates(meuse.grid) <- ~x+y
proj4string(meuse.grid) <- proj4string(sp.meuse)
# kriging
GS.fit <- krige(LM_res ~ 1, sp.meuse, 
                newdata = meuse.grid, model = fV)
# output
tmap_leaflet(tm_shape(GS.fit) + tm_dots("var1.pred", style = "cont"))

## ----varycoef preparation-----------------------------------------------------
# response variable
y <- meuse$l_cad
# covariates for fixed effects
X <- model.matrix(~1+dist+lime+elev, data = meuse)
# locations
locs <- as.matrix(meuse[, 1:2])/1000

## ----varycoef random effect---------------------------------------------------
# covariates for SVC (random effects)
W <- model.matrix(~1+dist+lime, data = meuse)

## ----varycoef control---------------------------------------------------------
# construct initial value (recall transformation from meters to kilometers)
(init <- c(
  # 3 times for 3 SVC
  rep(c(
    # range
    fV$range[2]/1000,
    # variance
    fV$psill[2]),     
  3), 
  # nugget
  fV$psill[1]
))
# control settings vor MLE
control <- SVC_mle_control(
  # profile likelihood optimization
  profileLik = TRUE,
  # initial values
  init = init
)

## ----varycoef fit-------------------------------------------------------------
# MLE
VC.fit <- SVC_mle(y = y, X = X, W = W, locs = locs,
                  control = control)
# outcome
summary(VC.fit)

# residuals
oldpar <- par(mfrow = c(1, 2))
plot(VC.fit, which = 1:2)
par(mfrow = c(1, 1))
plot(VC.fit, which = 3)
par(oldpar)

## ----varycoef predict locations-----------------------------------------------
newlocs <- coordinates(meuse.grid)/1000
# prediciton
VC.pred <- predict(VC.fit, newlocs = newlocs)
# outcome
head(VC.pred)
# transformation to a spatial points data frame
colnames(VC.pred)[1:3] <- c("Intercept", "dist", "lime")
sp.VC.pred <- VC.pred
coordinates(sp.VC.pred) <- coordinates(meuse.grid)
proj4string(sp.VC.pred) <- proj4string(sp.meuse)

# points of interest (POI), we come back to them later
POI1.id <- 235
POI2.id <- 2016
POI1 <- meuse.grid[POI1.id, ]
POI2 <- meuse.grid[POI2.id, ]

tm_POIs <- 
  tm_shape(POI1) +
  # square (pch = 15)
  tm_symbols(shape = 15, col = "black") +
  tm_shape(POI2) +
  # triangle (pch = 17)
  tm_symbols(shape = 17, col = "black")

## ----fig.height=6-------------------------------------------------------------
# tm.GS: GS fit
tm.GS <- tm_shape(GS.fit) + tm_dots("var1.pred", style = "cont")

# tm1: Intercept
tm1 <- tm_shape(sp.VC.pred) + 
  tm_dots("Intercept", style = "cont") +
  tm_POIs
# tm2: dist
tm2 <- tm_shape(sp.VC.pred) + 
  tm_dots("dist", style = "cont") +
  tm_POIs
# tm1: Intercept
tm3 <- tm_shape(sp.VC.pred) + 
  tm_dots("lime", style = "cont") +
  tm_POIs

tmap_arrange(list(tm.GS, tm1, tm2, tm3), ncol = 2, nrow = 2)

## ----coef means---------------------------------------------------------------
coef(VC.fit)

## ----POI equations------------------------------------------------------------
# POI1:
# SVC deviations at POI1
VC.pred[POI1.id, 1:3]
# combined
(coefs.poi1 <- coef(VC.fit) + c(as.numeric(VC.pred[POI1.id, 1:3]), 0))
# POI2
# SVC deviations at POI2
VC.pred[POI2.id, 1:3]
# combined
(coefs.poi2 <- coef(VC.fit) + c(as.numeric(VC.pred[POI2.id, 1:3]), 0))

options(digits = 3)

