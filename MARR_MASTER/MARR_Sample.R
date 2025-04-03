# ==============================================================================

# Keras installation ------------------------------------------------------
# If keras is not installed run the following lines to install
# once installed this can be removed or commented out
# 
# install.packages("keras")
# library(keras)
# install_keras()
# 
# -------------------------------------------------------------------------
###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 1:                             ###
###                    DATA INPUT AND INITIALIZATION                    ###
###                                                                     ###
###########################################################################
###########################################################################
# Set directory -----------------------------------------------------------
dir = "path\\to\\directory"
setwd(dir)

# Load required packages and Functions and basically everything
# if packages not installed (except for keras) they will be installed
# through the source file.
library(keras)
library(tensorflow)
suppressWarnings(source("MARR_Source_Functions.R")) # Do not change
# -------------------------------------------------------------------------
# Load the master model saved in your working directory -------------------
model <- load_model_hdf5("MARR_Master_Model.h5") # Do not change
# -------------------------------------------------------------------------
# =========================================================================
# Universal name (recommended: siteyr) and other universal variables ------
fn.Uni <- "siteYR" # Filename
sea.Uni <- "seaDir"      # Direction of sea ("N", "E", "S", or "W")
mhw.Uni <-  0.8    # Mean high water for shoreline
d2s <- 100          # Distance from shoreline to back of dune mask parameter
# -------------------------------------------------------------------------

# Create output folders: everything is saved in these folders -------------
# DO NOT CHANGE THIS
nm.Dirs <- c(paste(fn.Uni, "_Preprocess", sep=""),
             paste(fn.Uni, "_ModelResults_Raw", sep=""),
             paste(fn.Uni, "_ModelResults_Intermediate", sep=""),
             paste(fn.Uni, "_ModelResults_Final", sep=""))
# Check if folder already exists, if not create it
for(d in 1:length(nm.Dirs)){
  ifelse(!dir.exists(nm.Dirs[d]), dir.create(nm.Dirs[d]),
         "Folder exists already")
}
# -------------------------------------------------------------------------
# Load in data ------------------------------------------------------------
# DEM
M <- raster("rasterNameAnd.ext")  # Raster load
MT <- rast(M)                     # SpatRaster conversion

# OPTIONAL: if the dataset contains areas you are not interested in
# processing, select the region you don't want to process. (ie, a bluff
# section, removing bathymetry, etc.), not always applicable.
# {
#   x11()
#   plot(MT, col=turbo(200))
#   subMT <- sel(MT)
#   dev.off()
#   MT <- subMT
#   M <- raster(MT)
#   D <- as.matrix(M)
# }
# -------------------------------------------------------------------------

# Make hillshade for visualization adjust angle and direction
HS <- shade(terrain(MT * 2, "slope", unit="radians"),
            terrain(MT * 2, "aspect", unit="radians"),
            angle = 30, direction = 65)

# Create a copy of the raster as a matrix and set elevations <= -1 to NA
D <- as.matrix(M, wide = T)
D[D <= -1] <- NA

# Optional: Check the data makes sense visually
plot(rast(D), col=turbo(200))
plot(HS, col=grey(1:200/200))

subMT <- MT
# OPTIONAL: Subset part of the DEM for visualizing results as you go
#           select subset by clicking upper left and lower right corners
#           of rectangle on map.
# {
#   x11()
#   plot(MT, col=turbo(200))
#   subMT <- sel(MT)
#   dev.off()
# }

# Global variables for SpatRaster creation
crs.Uni <- crs(MT)
ext.Uni <- ext(MT)

# Metadata of the raster used for SurfNorm --------------------------------
NAvalue(M) <- -9999
nil <- NAvalue(M)

res <- res(M)
l = res[1]
# -------------------------------------------------------------------------

# OPTIONAL: smooth or resample data if needed -----------------------------
# Smooth the data
# M2 <- focal(M, w=matrix(1, 5, 5), mean)
# D <- as.matrix(M2)
# M <- M2
# MT <- rast(M)
# rm(M2)
#
# Resample the data
# M2 <- disaggregate(M, fact=3, method='bilinear',
#                    filename="Data_Resample.tif")
# M2 <- raster("Data_Resample.tif")
# M2 <- focal(M2, w=matrix(1, 5, 5), mean)
# D <- as.matrix(M2)
# M <- M2
# MT <- rast(M)
# rm(M2)
# -------------------------------------------------------------------------

# OPTIONAL: Manually loading different model and environment --------------
# load(file.choose()) # environment
# model <- load_model_hdf5(file.choose()) # model
# -------------------------------------------------------------------------

# 1A: Scales --------------------------------------------------------------
# Find the scale best suited for the upper end of the avg RR calc ---------
# perc = percentage of nr/nc to sample 10% is the default
#
# Global environment variables:
# Id.Sc = max scale
# mn.Sc = min scale
# shore = shoreline matrix
# -------------------------------------------------------------------------
suppressWarnings(scale.Find(perc=0.1))

# OPTIONAL: set your scales manually
# Id.Sc <- 27 # Max, or 'Ideal' scale
# mn.Sc <- 7  # Min scale
# -------------------------------------------------------------------------


# 1B: Variables -----------------------------------------------------------
# Preprocess the data to get the input variables SN vectors, etc. ---------
# Saves into the preprocessing folder
#
# Global environment variables:
# rr         = relative relief matrix
# RX, RY     = ridge feature matrices
# SX, SY     = swale feature matrices
# UX, UY, UZ = surface normal feature matrices
# XE, YN     = easting and northing matrices
# tpi, tri   = tpi and tri feature matrices
# shoreMask  = shoreline mask matrix

dirPre <- paste(dir,"\\",fn.Uni, "_Preprocess", sep="")
setwd(dirPre)
system.time(suppressWarnings(preProcess.Morph(makeRR=T, makeSN=T, inMn=mn.Sc, snWn=5, shore.Ext=d2s)))
plot.Inputs()
setwd(dir)

# Check DEM alongshore ridges
plot(subMT, col=turbo(200), main='Ridges')
plot(rast(as.matrix(R*shoreMask), ext=ext.Uni, crs=crs.Uni),
     col='black', add=T, legend = F)

# Check RR alongshore swales
plot(subMT, col=turbo(200), main='Swale')
plot(rast(as.matrix(S*shoreMask), ext=ext.Uni, crs=crs.Uni),
     col='black', add=T, legend = F)
# ------------------------------------------------------------------------------
# Create i and j mats for cross reference --------------------------------------
# Global environment variables:
# iMat = matrix filled with row values
# Mat = matrix filled with col values

createIJMat()
# ------------------------------------------------------------------------------

# -------------------------------------------------------------------------
# TODO: add in anything that saves to file
# Get characteristics of R/S line features --------------------------------
# Extracts the length (i.e., number of cells that make up a feature) and
# assigns a unique ID to each ridge and swale line feature as a matrix
#
# Global environment variables:
# R.features[[1]] = Ridge length characteristics
# R.features[[2]] = Ridge unique IDS
# S.features[[1]] = Swale length characteristics
# S.features[[2]] = Swale unique IDS

get.RS.Features()

# Save the ridge and swale features
out <- c(rast(R.features[[1]], ext=ext.Uni, crs=crs.Uni),
         rast(R.features[[2]], ext=ext.Uni, crs=crs.Uni),
         rast(S.features[[1]], ext=ext.Uni, crs=crs.Uni),
         rast(S.features[[2]], ext=ext.Uni, crs=crs.Uni))

names(out) <- c('R.Length', 'R.IDs', 'S.Length', 'S.IDs')
out[out == 0] <- NA
writeRaster(out,paste0(dir,"\\",fn.Uni, "_Preprocess\\", fn.Uni, "_RS_Feats.tif"), overwrite=TRUE)

RS.feat.Ras <- out
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Make model input data ---------------------------------------------------
# Mean values of Elev., RR, nX, nY, nZ, TPI, and TRI for cells along individual
# ridge and swale features
#
# Stores backup as .csvs in the Preprocess folder
#
# Global environment variables:
# R.Means = Ridge features
# S.Means = Swale features

create.Model.Data()
# -------------------------------------------------------------------------

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 2A: Raw                        ###
###                      MODEL DATA AND PROCESSING                      ###
###                                                                     ###
###########################################################################
###########################################################################


# Get raw results from pretrained model -----------------------------------
# push data through the model (R.Means and S.Means) and get predictions
# for each ridge and swale feature
#
# Stores results in ModelResults_Raw folder
#
# Global environment variables:
# p1.Ras = matrix of probabilities that a feature is a toe
# p2.Ras = matrix of probabilities that a feature is a crest
# p3.Ras = matrix of probabilities that a feature is neither toe or crest
# tc.Results = matrix of toe and crest features predicted by the model
# outData = data frame of model results
# S.stack = raster stack of model results for swale features
# R.stack = raster stack of model results for ridge features
#
# pVal = threshold for likelihoods to keep
# modNew = is it the original or new model
# iteration = what iteration model are you using
#             "1_Basemodel": original model, first pass

modelProcessing(pVal=0, modNew=F, iteration = "1_BaseModel")

# ------------------------------------------------------------------------------
############################################################################
############################################################################
###                                                                      ###
###                              SECTION 2B: Intermediate                ###
###            TRAINING DATA AND SITE SPECIFIC MODEL TRAINING            ###
###                                                                      ###
############################################################################
############################################################################

# Initialize, ONLY RUN ONCE AT BEGINNING ----------------------------------
#
# Stores results in ModelResults_Raw folder
#
# c/tP1 = percentile of likelihoods to keep as a 2/1
# c/tP3 = percentile of likelihoods to keep as a 3
#
# Global environment variables:
# trainClass = matrix of training data classes
siteTrainingData(tP1 = 0.90, tP3 = 0.15,
                 cP1 = 0.50, cP3 = 0.10)

# -------------------------------------------------------------------------


# Create subset DEMs to better visualize ridges and swales for
# creating training data. In the pop up window choose how many
# subsamples you would like to make.
# Then click on the map to select the upper left and lower right
# corners of a rectangle that defines a subsample.
# For large datasets choose representative sample locations.
# For smaller datasets it is good practice to include overlap between subsamples.
# ~300m-500m alongshore extent seems to work well, you can always rerun if needed
# Prioritize areas with results that don't seem reasonable
# More subsets with smaller extent is a better practice than fewer subsets
# with a broad extent.
TC.Raster <- rast(tc.Results, ext=ext.Uni, crs=crs.Uni)
c(subs.DEM, polys) %<-% subsetDEM(inel = MT, TC.Raster = TC.Raster,
                                  n = 7, dynamicN = T)

# User defined selection of training samples using subsetted DEM
# In the pop up window make your choice of training samples you want to select.
# Best practice is to choose an n larger than required to account for misclicking,
# NAs and duplicates are removed so a greater n is better, especially for class 3.
# Having the console visible below the window is helpful for seeing if a click
# is missed (returns NA).
# It is important to remember which class you are identifying (1, 2, 3) so that
# your training data is accurate.
# This step may need to be run more than once as you train the model.
#
# The goal is to correct the model results with better samples. So when
# choosing samples it is best to select misclassifications.
# The model does not need a substantial amount of samples, therefore:
# - prioritize misclassifications
# - for classes 1/2 (S1/R2) your number of samples should be a few more than
#   than the real features present in the subsample
# - for class 3 (S3/R3) your number of samples should be a few more than the
#   misclassifications in the subsample, to select these plus additional samples
# - it is always better to have a larger n, if you have succesfully selected
#   the desired points, simply click anywhere that is not a feature to set NA
# - where you click on a line does not matter, it is best to choose a perfectly
#   vertical or horizontal section
# - misclassifications are not colored by ID however line segments may be
#   included under the same ID if separated by only a couple cells, redundancy
#   is best, so assume the individual line segments are separate, but there
#   is no need to select every single segment.
# - model results are colored by ID so it can be advantageous to select these
#   lines at more optimal locations along the line, the colored section is
#   not necessarily where you have to click
# - if a mistake is made (ie not enough samples selected) it is best to choose
#   1 sample for remaining subsample DEMs and click anywhere to reset. This
#   has caused crashes in testing when closing out. OR exit the sample input 
#   using the X button, then do the same for the DEM sample.

# SELECT SWALES THAT ARE REAL TOES
S1 <- getSelection(inel = subs.DEM, feat = RS.feat.Ras$S.IDs, subs.DEM = subs.DEM,
                      TC.Raster = TC.Raster, type = "Swale", goal = 1,
                      n = 2, dynamicN = T)
# SELECT SWALES THAT ARE NOT TOES
S3 <- getSelection(inel = subs.DEM, feat = RS.feat.Ras$S.IDs, subs.DEM = subs.DEM,
                      TC.Raster = TC.Raster, type = "Swale", goal = 3,
                      n = 2, dynamicN = T)
# SELECT RIDGES THAT ARE CRESTS
R2 <- getSelection(inel = subs.DEM, feat = RS.feat.Ras$R.IDs, subs.DEM = subs.DEM,
                      TC.Raster = TC.Raster, type = "Ridge", goal = 2,
                      n = 2, dynamicN = T)
# SELECT RIDGES THAT ARE NOT CRESTS
R3 <- getSelection(inel = subs.DEM, feat = RS.feat.Ras$R.IDs, subs.DEM = subs.DEM,
                      TC.Raster = TC.Raster, type = "Ridge", goal = 3,
                      n = 2, dynamicN = T)

# ONLY if this is the first round of training run this
Swale1 <- S1
Swale3 <- S3
Ridge2 <- R2
Ridge3 <- R3

# ONLY if this is not the first round of training AND you want to include
# previously sample locations run this
Swale1 <- unique(c(Swale1, S1))
Swale3 <- unique(c(Swale3, S3))
Ridge2 <- unique(c(Ridge2, R2))
Ridge3 <- unique(c(Ridge3, R3))


# -------------------------------------------------------------------------

# Create new training data based on selected samples ----------------------
newTraining <- new.Training(inR=R.features[[2]], inS=S.features[[2]])

# Save as raster to intermediate folder
dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Intermediate", sep="")
setwd(dirPre)
out = raster(newTraining,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
projection(out) = crs(M)
writeRaster(out,"NewTraining.tif", overwrite=TRUE)
setwd(dir)

# Train new model the new model is save in the intermediate folder
trainClass <- newTraining
siteModel(trainClass=trainClass, LL=6, ridgePlot=T)

# Push data through the new model and get results saved in the intermediate folder
# it is best to adust the iteration value to account for the adjustments
# (ie iteration = 2: trained model first run, second results)
modelProcessing(pVal=0, modNew=T, iteration = 3)
TC.Raster <- rast(tc.Results, ext=ext.Uni, crs=crs.Uni)

plot(subMT, col=turbo(200), legend = F)
plot(TC.Raster, col=c('darkorchid', 'black'), add=T)

# IF YOU ARE HAPPY WITH RESULTS MOVE ON, IF NOT SELECT MORE SAMPLES AND
# TRAIN AGAIN UNTIL HAPPY (usually 2-3 iterations is enough) IF DIMINISHING
# RETURNS MOVE ON. REMEMBER YOU ARE SELECTING NEW SAMPLES SO THERE IS NO
# NEED TO BE AS THOROUGH. YOU CAN SELECT NEW SUBSAMPLE DEMS AS WELL IF WANTED.
#
# IMPORTANT: the nature of the approach relies on rows/cols so some results
# that seem incorrect may be impossible to fix with training and are later
# corrected so a 'good as I can get it' mentality should be used for training.
# =========================================================================

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 3:                             ###
###                    DATA FILTERING AND FINALIZING                    ###
###                                                                     ###
###########################################################################
###########################################################################


# Select samples to remove ------------------------------------------------
# The removed samples are ones that are not crest or toes and may be the
# result of the approach or the model or some other outside influence.
# It is best to have the results (from the intermediate folder) open in a
# GIS that allows for drawing profiles (Arc, Q, etc.) so that if needed you
# can assess if a line is accurate or not.
# -------------------------------------------------------------------------
# Separate crests and toes ------------------------------------------------
C.ras <- tc.Results
T.ras <- tc.Results

C.ras[C.ras==1] <- NA
C.ras[!is.na(C.ras)] <- 1
T.ras[T.ras==2] <- NA

# Remove single cells
C.ras <- single.RM(C.ras, C.ras*NA)
C.ras <- C.ras[[1]]
T.ras <- single.RM(T.ras, T.ras*NA)
T.ras <- T.ras[[1]]

# Get IDs of features
C.ID <- feat.Extract(el=D, ft=as.matrix(C.ras, wide=T), in.IT=T, inXY=inXY2)[[2]]
T.ID <- feat.Extract(el=D, ft=as.matrix(T.ras, wide=T), in.IT=T, inXY=inXY2)[[2]]
C.ID[C.ID == 0] <- NA
T.ID[T.ID == 0] <- NA

C.ID.ras <- rast(C.ID, ext=ext.Uni, crs=crs.Uni)
T.ID.ras <- rast(T.ID, ext=ext.Uni, crs=crs.Uni)
# -------------------------------------------------------------------------
# Subset DEM for selection of features to remove --------------------------
# Only run this if you need to select more subsample DEMs or different ones
TC.Raster <- rast(tc.Results, ext=ext.Uni, crs=crs.Uni) # MAKE SURE YOU RUN THIS
c(subs.DEM, polys) %<-% subsetDEM(inel = MT, TC.Raster = TC.Raster,
                                  n = 7, dynamicN = T)
# -------------------------------------------------------------------------

# Select results to remove ------------------------------------------------
# If a subsample has nothing to remove set n = 1 and click anything that
# is not a feature for NA. If there are samples to remove, remember it is
# best to have a larger n. 
#
# Sometimes the subsample DEM is too large and clicking on features is tedious,
# it is best to set n to 1 and click an NA value to move on. Then you can redo
# the subsample DEM data to make selecting the features easier in a subsequent
# iteration.

# TOES TO REMOVE
T.RMV <- getSelection(inel = subs.DEM, feat = T.ID.ras, subs.DEM = subs.DEM,
                      TC.Raster = TC.Raster, type = "Swale", goal = 3,
                      n = 2, dynamicN = T, hlShd = T)
# CRESTS TO REMOVE
C.RMV <- getSelection(inel = subs.DEM, feat = C.ID.ras, subs.DEM = subs.DEM,
                      TC.Raster = TC.Raster, type = "Ridge", goal = 3,
                      n = 2, dynamicN = T, hlShd = T)

# Remove samples from results
TC.ras <- tc.Results
TC.ras[which(C.ID %in% C.RMV)] <- NA
TC.ras[which(T.ID %in% T.RMV)] <- NA
tc.Results <- TC.ras
TC.ras <- rast(TC.ras, ext=ext.Uni, crs=crs.Uni)

# IF YOU KNOW THERE ARE STILL SOME THAT NEED TO BE REMOVED START AGAIN,
# THE PREVIOUSLY SELECTED CHOICES STAY REMOVED AND FOCUS SUBSET DEMS

# Visualize
plot(MT, col=grey(1:200/200), legend=F)
TC.p <- (as.points(TC.ras))
plot(TC.p, col=c('seagreen1', 'red')[TC.p$lyr.1], add = T, pch=19, cex = 0.25)


# -------------------------------------------------------------------------

# FROM HERE RESULTS SAVED IN FINAL FOLDER ---------------------------------

# Fill in gaps between features -------------------------------------------
tc.Results.Filtered <- tc.Results
TC.Filled <- suppressWarnings(TC.Filled.fun(RFill=T, TFill=T,
                                            R.RNG=5, S.RNG=0,
                                            doT=T, doC=T))

# Convert to raster
TC.Filled.ras <- rast(TC.Filled, ext=ext.Uni, crs=crs.Uni)

# Make points and visualize to test filled values
TC.ras.pt <- as.points(TC.ras)
TC.Filled.pt <- as.points(TC.Filled.ras)
names(TC.Filled.pt) <- "Feat"

plot(MT, col=grey(1:200/200), legend = F)
plot(TC.Filled.pt, add=T, col=c('white', 'orange')[TC.Filled.pt$Feat], cex=0.25, pch=19)
plot(TC.ras.pt, add=T, cex=0.25, pch=19)

# ADJUST THE R.RNG AND S.RNG VALUES UNTIL HAPPY WITH RESULTS
# (ORANGE AND WHITE LINES SHOULD BE REASONABLE AND NOT JERKY)
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# VERY OPTIONAL: remove sections toe and crest lines based on site specifics
#                like washovers or specific sections you are not interested 
#                having results for. Most datasets shouldn't require this, 
#                only in specific cases.

# Get subsamples of DEM ---------------------------------------------------
# These samples should be slightly larger than the section you want to 
# remove. If you make too many selections (ie n is too large), select a
# region that falls within another section you are removing.
# TC.Raster <- TC.Filled.ras 
# c(subs.DEM, polys) %<-% subsetDEM(inel = MT, TC.Raster = TC.Raster,
#                                   n = 7, dynamicN = T)
# 
# # Select bounds of areas to remove
# # Simply click an location at the alongshore extent of regions you wish
# # to remove, no where in particular just along an imaginary cross-shore line.
# B.RMV <- getSelectionRMV(inel = subs.DEM, 
#                          subs.DEM = subs.DEM, 
#                          TC.Raster = TC.Raster,
#                          sea = sea.Uni)
# 
# # Visualize the bounds
# plot(HS, col=grey(1:200/200), legend=F)
# plot(TC.Filled.pt, add=T, col=c('seagreen1', 'navy')[TC.Filled.pt$Feat], cex=0.25, pch=19)
# abline(v = B.RMV$x, col = rainbow(max(B.RMV$section))[B.RMV$section], lwd=2)
# 
# # Remove data within bounds
# for(k in unique(B.RMV$section)){
#   foo <- B.RMV[B.RMV$section == k, ]
#   if(sea.Uni == "N" || sea.Uni == "S"){
#     TC.Filled.ras[, min(foo$bnd):max(foo$bnd)] <- NA
#   }else{
#     TC.Filled.ras[min(foo$bnd):max(foo$bnd),] <- NA
#   }
#   rm(foo)
# }
# 
# # Make points and visualize 
# TC.Filled.pt <- as.points(TC.Filled.ras)
# names(TC.Filled.pt) <- "Feat"
# 
# plot(HS, col=grey(1:200/200), legend = F)
# plot(TC.Filled.pt, add=T, col=c('red', 'navy')[TC.Filled.pt$Feat], cex=0.25, pch=19)
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Outputs to final folder -------------------------------------------------
# Make crest and toe rasters
C.ras <- TC.Filled.ras
C.ras[C.ras != 2] <- NA
C.ras[!is.na(C.ras)] <- 1
T.ras <- TC.Filled.ras
T.ras[T.ras != 1] <- NA

plot(subMT, col=grey(1:200/200), legend = F)
plot(T.ras, add=T)

# save rasters
dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Final", sep="")
setwd(dirPre)
writeRaster(TC.Filled.ras,paste0(fn.Uni, "_TC_Final.tif"), overwrite=TRUE)
writeRaster(T.ras,paste0(fn.Uni, "_toe_Final.tif"), overwrite=TRUE)
writeRaster(C.ras,paste0(fn.Uni, "_crest_Final.tif"), overwrite=TRUE)
setwd(dir)

# Convert to points and save as shapefiles --------------------------------
# Final points
names(TC.Filled.pt)[1] <- 'Feat'

# Assign variable values as attributes
TC.Filled.pt[["Z"]] <- terra::extract(MT, TC.Filled.pt, ID=F, cells=F, xy=F)
TC.Filled.pt[["RR"]] <- terra::extract(rast(rr, ext=ext.Uni, crs=crs.Uni), TC.Filled.pt, ID=F, cells=F, xy=F)
TC.Filled.pt[["TPI"]] <- terra::extract(rast(tpi, ext=ext.Uni, crs=crs.Uni), TC.Filled.pt, ID=F, cells=F, xy=F)
TC.Filled.pt[["TRI"]] <- terra::extract(rast(tri, ext=ext.Uni, crs=crs.Uni), TC.Filled.pt, ID=F, cells=F, xy=F)

C.pt <- TC.Filled.pt[TC.Filled.pt$Feat == 2,]
T.pt <- TC.Filled.pt[TC.Filled.pt$Feat == 1,]

writeVector(TC.Filled.pt, paste0(dir,"\\",fn.Uni, "_ModelResults_Final\\Final_TC_Points.shp"), overwrite=T)
writeVector(C.pt, paste0(dir,"\\",fn.Uni, "_ModelResults_Final\\Final_Crest_Points.shp"), overwrite=T)
writeVector(T.pt, paste0(dir,"\\",fn.Uni, "_ModelResults_Final\\Final_Toe_Points.shp"), overwrite=T)
# ------------------------------------------------------------------------------
# Plot random samples of "profiles" to PDF -------------------------------------
# maxX = constraint for x axis, keep NA unless you know how for to draw in
# Best practice is to run with minX = 1, maxX = NA, open and then run again
# with ideal x range.
# no.Smp = number of samples to draw
rand.plot.fin(no.Smp = 20, minX=600, maxX=800)
# ------------------------------------------------------------------------------
# Save results as CSVs ---------------------------------------------------------
dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Final", sep="")
setwd(dirPre)
TC.Filled <- as.matrix(TC.Filled.ras, wide = T)
feat.dat <- data.frame(x=c(XE),
                       y=c(YN),
                       z=c(D),
                       feat=c(TC.Filled))
feat.dat <- na.omit(feat.dat)
t.Dat <- feat.dat[feat.dat$feat==1,]
c.Dat <- feat.dat[feat.dat$feat==2,]

write.csv(t.Dat, "Toe_Final.csv",row.names=F)
write.csv(c.Dat, "Crest_Final.csv",row.names=F)
write.csv(feat.dat, "TC_Final.csv",row.names=F)
setwd(dir)
# ------------------------------------------------------------------------------
# Alongshore plot of Crest and Toe ---------------------------------------------
# ww and hh are width and height of the ouput pdf, respectively, in inches
Alongshore.Plot(ww=15, hh=8)
# -------------------------------------------------------------------------

























