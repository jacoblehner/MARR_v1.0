# Find shoreline and draw profiles perpendicular to it
# Load Libraries
#automatic install of packages if they are not installed already
list.of.packages <- c(
  "raster",
  "terra",
  "viridis",
  "segmented",
  "tcltk",
  "scales",
  "ggplot2",
  "ggridges",
  "tidyr",
  "sf",
  "dplyr"
)
# list.of.packages2 <- c("terra")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
    )
  )
}

################################################################################

# Supp --------------------------------------------------------------------
# Dynamic window for selection count
inputs <- function(n = 1){

  xvar <- tclVar(n)

  tt <- tktoplevel()
  tkwm.title(tt,"Sample Count")
  x.entry <- tkentry(tt, textvariable=xvar)

  # reset <- function()
  # {
  #   tclvalue(xvar)<-1
  # }

  # reset.but <- tkbutton(tt, text="Reset", command=reset)

  submit <- function() {
    x <- as.numeric(tclvalue(xvar))
    e <- parent.env(environment())
    e$x <- x
    tkdestroy(tt)
  }
  submit.but <- tkbutton(tt, text="submit", command=submit)

  tkgrid(tklabel(tt,text="How many samples, n?"),columnspan=1)
  tkgrid(tklabel(tt,text="n"), x.entry, pady = 10, padx =10)
  tkgrid(submit.but)#, reset.but)

  tkwait.window(tt)
  return(x)
}
# Color maps
colFn <- function(n=255){

  red <- c(175,255,0,252,128, 120,105,171,255)
  green <- c(240,255,128,186,0,0,48,171,255)
  blue <- c(233,178, 64, 0,0,0,13,171,255)
  relpos <- c(0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 1.0)
  relpos1 <- relpos * n

  cols <- gdalraster::createColorRamp(start_index = relpos1[1],
                                      start_color = c(red[1], green[1], blue[1]),
                                      end_index = relpos1[2],
                                      end_color = c(red[2], green[2], blue[2]))

  for(i in 2:(length(red)-1)){
    tmp <- gdalraster::createColorRamp(start_index = relpos1[i],
                                       start_color = c(red[i], green[i], blue[i]),
                                       end_index = relpos1[i+1],
                                       end_color = c(red[i+1], green[i+1], blue[i+1]))
    cols <- rbind(cols, tmp)
  }

  return(rgb(cols[c(1:n),2:4]/255))

}
gg_color_hue <- function(nn) {
  hues = seq(15, 375, length = nn + 1)
  hcl(h = hues, l = 65, c = 100)[1:nn]
}
# Subsetting DEM
subsetDEM <- function(inel = MT, TC.Raster = TC.Raster,
                      n = ceiling(nrow(MT) / 300), dynamicN = T){

  windows.options(xpos = 0, ypos = 0, width = 15, height = 9)
  # Initialize
  subList <- list()

  names(TC.Raster) <- "Feat"
  TC.p <- (as.points(TC.Raster))
  # TC.p <- extract(TC.Raster, TC.p, xy=T, id=F, method="exact")
  # TC.p <- TC.p %>%
  #   st_as_sf( coords = c("x", "y")) %>%
  #   group_by(Feat) %>%
  #   # summarize() %>%
  #   st_cast("MULTIPOINT") %>%
  #   st_set_crs(crs(inel))

  X11()
  plot(inel, col= grey(1:255/255), legend=F)
  plot(TC.p, col=c('seagreen1', 'red')[TC.p$Feat], add = T, pch=19, cex = 0.25)
  # plot(TC.Raster, col= c('seagreen1', 'darkorchid'), add = T)


  if(dynamicN){
    n <- inputs()
  }
  subTmp <- sel(inel, col='yellow', lwd = 1)
  sel
  p <- as(extent(raster(subTmp)), 'SpatialPolygons')
  p$id <- 1
  subList[[1]] <- subTmp
  # dev.off()

  for(i in 2:n){
    # X11()
    # plot(inel, col= grey(1:255/255), legend=F)
    # plot(TC.Raster, col= c('seagreen1', 'darkorchid'), add = T)
    # plot(p, add=T, border='yellow', lwd=1)

    subTmp <- sel(inel, col='yellow', lwd = 1)
    ptmp <- as(extent(raster(subTmp)), 'SpatialPolygons')
    ptmp$id <- i
    p <- rbind(p, ptmp)
    subList[[i]] <- subTmp
  }
  dev.off()

  plot(inel, col= grey(1:255/255), legend=F)
  plot(TC.Raster, col= c('seagreen1', 'darkorchid'), add = T)
  plot(p, border='yellow', lwd=1, add=T)
  return(list(subList, p))

}
# Feature selection by clicking
getSelection <- function(inel = subs.DEM, feat = RS.feat.Ras, subs.DEM = subs.DEM,
                         TC.Raster = TC.Raster, type = c('Swale', 'Ridge'), goal = c(1, 2, 3),
                         n = 10, dynamicN = T, hlShd = F){
  windows.options(xpos = 0, ypos = 0, width = 15, height = 9)

  if(goal == 3){
    subtit <- paste0("SELECT FALSE TOES/CRESTS: Class 3")
  }else if(goal == 1){
    subtit <- paste0("SELECT TRUE TOE: Class 1")
  }else if(goal == 2){
    subtit <- paste0("SELECT TRUE CREST: Class 2")
  }else{
    return("Goal must be 1, 2, or 3")
  }

  if(type == 'Swale'){
    # feat <- feat$S.IDs
    tc.Foc <- TC.Raster
    tc.Foc[tc.Foc == 2] <- NA
    tc.Foc[!is.na(tc.Foc)] <- 1
    tc.Foc <- tc.Foc * feat
  }else if(type == 'Ridge'){
    # feat <- feat$R.IDs
    tc.Foc <- TC.Raster
    tc.Foc[tc.Foc == 1] <- NA
    tc.Foc[!is.na(tc.Foc)] <- 1
    tc.Foc <- tc.Foc * feat
  }else{
    return("Type must be either 'Ridge' or 'Swale'")
  }

  unqsTOT <- c(na.omit(unique(c(as.matrix(tc.Foc, wide=T)))))

  i <- 1
  HS1 <- shade(terrain(inel[[i]], "slope", unit="radians"),
               terrain(inel[[i]], "aspect", unit="radians"),
               45, 45)
  maintit <- paste0("subDEM ", i, " - ", type, 's')

  x11()
  # cmap <- gg_color_hue(nn = length(unqsTOT))
  plot(HS1, col=grey(1:175/200), legend = F,
       main=paste0(maintit, " || ", subtit), cex.main=1, mar=c(2, 2, 2, 5))
  if(!hlShd){
    plot(inel[[i]], add=T, col=alpha(colFn(255), 0.3), legend = F)
  }
  plot(feat, add = T, legend = F, col= 'black')
  plot(as.factor(tc.Foc), col= rainbow(length(unqsTOT)), add = T)
  # Sys.sleep(5)
  if(dynamicN){
    # n <- as.numeric(dlgInput("How many samples, n?", 1)$res)
    n <- inputs(n)
  }
  clicks <- click(feat, n = n, id=F, xy=F, cell=F, type='p', show=T, col='white')

  clicksDF <- as.data.frame(clicks)

  for(i in 2:length(inel)){
    HS1 <- shade(terrain(inel[[i]], "slope", unit="radians"),
                 terrain(inel[[i]], "aspect", unit="radians"),
                 45, 45)
    maintit <- paste0("subDEM ", i, " - ", type, 's')

    # x11()
    plot(HS1, col=grey(1:175/200), legend = F,
         main=paste0(maintit, " || ", subtit), cex.main=1, mar=c(2, 2, 2, 5))
    if(!hlShd){
      plot(inel[[i]], add=T, col=alpha(colFn(255), 0.3), legend = F)
    }
    plot(feat, col= "black", add = T, legend = F)
    plot(tc.Foc, col= rainbow(length(unqsTOT)), add = T)
    Sys.sleep(1)
    if(dynamicN){
      # n <- as.numeric(dlgInput("How many samples, n?", 1)$res)
      n <- inputs()
    }
    clicks <- click(feat, n = n, id=F, xy=F, cell=F, type='p', show=T, col='white')

    clicksDF <- rbind(clicksDF, clicks)
  }
  dev.off()

  out <- unique(clicksDF[,1])
  out <- out[!is.na(out)]

  return(out)

}
# Feature selection you will remove completely
getSelectionRMV <- function(inel = subs.DEM, subs.DEM = subs.DEM,
                         TC.Raster = TC.Raster, n = 2, sea = sea.Uni){
  windows.options(xpos = 0, ypos = 0, width = 15, height = 9)
  
  
  tc.Foc <- TC.Raster
  unqsTOT <- c(na.omit(unique(c(as.matrix(tc.Foc, wide=T)))))
  
  i <- 1
  HS1 <- shade(terrain(inel[[i]], "slope", unit="radians"),
               terrain(inel[[i]], "aspect", unit="radians"),
               45, 45)
  x11()
  # cmap <- gg_color_hue(nn = length(unqsTOT))
  plot(HS1, col=grey(1:175/200), legend = F, mar=c(2, 2, 2, 5))
  plot(as.factor(tc.Foc), col= c('seagreen1', 'red'), add = T)
  clicks <- click(HS1, n = n, id=T, xy=T, cell=F, type='p', show=T, col='white')
  clicks$section <- i
  clicksDF <- as.data.frame(clicks)
  
  for(i in 2:length(inel)){
    HS1 <- shade(terrain(inel[[i]], "slope", unit="radians"),
                 terrain(inel[[i]], "aspect", unit="radians"),
                 45, 45)
    
    # x11()
    plot(HS1, col=grey(1:175/200), legend = F, cex.main=1, mar=c(2, 2, 2, 5))
    plot(as.factor(tc.Foc), col= c('seagreen1', 'red'), add = T)
    clicks <- click(HS1, n = n, id=T, xy=T, cell=F, type='p', show=T, col='white')
    clicks$section <- i
    clicksDF <- rbind(clicksDF, clicks)
  }
  dev.off()
  
  
  # out <- clicksDF[complete.cases(clicksDF),]
  # out <- clicksDF[,1][unique(clicksDF[,1]),]
  if(sea == 'E'){
    clicksDF$bnd <- rowFromY(TC.Raster, clicksDF$y)
  }
  if(sea == 'W'){
    clicksDF$bnd <- rowFromY(TC.Raster, clicksDF$y)
  }
  if(sea == 'N'){
    clicksDF$bnd <- colFromX(TC.Raster, clicksDF$x)
  }
  if(sea == 'S'){
    clicksDF$bnd <- colFromX(TC.Raster, clicksDF$x)
  }
  
  return(clicksDF)
  
}

# -------------------------------------------------------------------------
# 1A ----------------------------------------------------------------------
# Get scales for RR KEEP
scale.Find <- function(perc=0.1){
  # library(segmented)
  shore <- shore.Mat.Extract(as.matrix(M), sh.Z=mhw.Uni, sea=sea.Uni)
  shoreMask <- shorelineMask(inshore=shore, dist=d2s, sea=sea.Uni, inOM = F)
  # Need the largest window to reach the toe and crest to get a local min at the
  # toe closer to 0.
  # For a random col/row get elevation data based on shore position and distance
  # from shore
  cat("\n", "Random sampling...")
  cat("\n")
  seg.Prof <- function(inProf=NA, inCSD=NA, splits=2, strR=25){

    prof <- inProf
    CS.dis <- inCSD

    # prof <- 2374
    # CS.dis <- 150
    # splits=3
    # thresh=5
    # dim(D)
    DD <- D*shoreMask
    if(sea.Uni=="E"){
      S.p <- which(!is.na(shore[prof,]))
      if(length(S.p) == 0){
        return(c(NA,NA))
      }
      # tmp <- c(D[prof,(S.p):(S.p-CS.dis)])
      tmp <- c(DD[prof,(S.p):1])
    }else if(sea.Uni=="W"){
      S.p <- which(!is.na(shore[prof,]))
      if(length(S.p) == 0){
        return(c(NA,NA))
      }
      # tmp <- c(D[prof,(S.p):(S.p+CS.dis)])
      tmp <- c(DD[prof,(S.p):(ncol(D))])
    }else if(sea.Uni=="N"){
      S.p <- which(!is.na(shore[,prof]))
      if(length(S.p) == 0){
        return(c(NA,NA))
      }
      # tmp <- c(D[(S.p):(S.p+CS.dis), prof])
      tmp <- c(DD[(S.p):(nrow(D)), prof])
    }else if(sea.Uni=="S"){
      S.p <- which(!is.na(shore[,prof]))
      if(length(S.p) == 0){
        return(c(NA,NA))
      }
      # tmp <- c(D[prof,(S.P):(S.p-CS.dis)])
      tmp <- c(DD[prof,(S.P):1])
    }

    # plot(tmp)
    topsZ <- lapply(1:7, function(x) inflect(tmp, threshold = x)$maxima)
    topsZ[[7]] <- topsZ[[7]][which(tmp[topsZ[[7]]] > mean(tmp,na.rm=T))]
    # points(topsZ[[7]],tmp[topsZ[[7]]],pch=19, col='red')
    # abline(h=mean(tmp,na.rm=T))

    if(length(topsZ[[7]]) == 0){
      return(c(NA,NA))
    }

    tmp <- tmp[1:topsZ[[7]][1]]

    df <- data.frame(x=c(1:length(tmp)),
                     y=tmp)

    #fit simple linear regression model
    fit <- lm(y ~ x, data=df)

    #fit piecewise regression model to original model, estimating a breakpoint at x=9
    segmented.fit <- segmented(fit, seg.Z = ~x)

    psi.ind <- c(floor(segmented.fit$psi[,2]))

    # #plot original data
    # if(Viz){
    #   plot(df$x, df$y, pch=16, cex=0.5, xlab="x", ylab="y", main=paste0("Profile: ", inProf))
    #   abline(h=mhw.Uni, lty="dashed")
    #   points(psi.ind, df$y[psi.ind], pch=16, cex=1.5, col='blue')
    #   #add segmented regression model
    #   plot(segmented.fit, add=T, col="blue", lwd=2)
    # }

    # if(length(psi.ind) > 1){
    #   p.foc <- psi.ind[2]
    # }else{
    #
    # }
    p.foc <- psi.ind#[2]
    return(c(inProf, nrow(df) - p.foc))#psi.ind[2])
  }

  # smp.P <- sample(1900:(nrow(D)-800), 250, replace=F)
  if(sea.Uni == "E" || sea.Uni == "W"){
    smp.P <- sample(100:(nrow(D)-100), floor(nrow(D)*perc), replace=F)
  }else{
    smp.P <- sample(100:(ncol(D)-100), floor(ncol(D)*perc), replace=F)
  }

  ind <- data.frame(prof=NA,
                    Dis=NA)

  for(i in 1:length(smp.P)){
    cat("\r", smp.P[i])#, "/", length(smp.P))
    flush.console()
    if(length(min(shore)))
    ind[i,] <- seg.Prof(inProf=smp.P[i], inCSD=500, splits=1, strR=25)
  }

  cat("\n", "Selecting scales...")
  cat("\n")
  # ind <- ind[ind > 10]
  plot(ind$prof, ind$Dis, pch=19, cex=0.5)
  abline(h=c(median(ind$Dis, na.rm=T), mean(ind$Dis, na.rm=T)), col=c("red","blue"), lwd=1.5)
  abline(h=c(mean(ind$Dis, na.rm=T)-sd(ind$Dis, na.rm=T), mean(ind$Dis, na.rm=T)+sd(ind$Dis, na.rm=T)), col=c("forestgreen"), lwd=1.5)
  quantile(ind$Dis, na.rm=T)

  range(ind$Dis, na.rm=T)

  rr.Scales <- round(c(mean(ind$Dis, na.rm=T)-sd(ind$Dis, na.rm=T), mean(ind$Dis, na.rm=T)+sd(ind$Dis, na.rm=T)))
  if((rr.Scales[1] %% 2) == 0) {
    rr.Scales[1] <- rr.Scales[1] - 1
  }
  if((rr.Scales[2] %% 2) == 0) {
    rr.Scales[2] <- rr.Scales[2] + 1
  }
  cat("\n",rr.Scales)
  cat("\n")
  Id.Sc <<- min(65, max(21,rr.Scales[2]))

  mn.Sc <<- max(7, rr.Scales[1])
  shore <<- shore
}
shore.Mat.Extract <- function(inEl = D, sh.Z = 0.8, sea = NA){
  print("Extracting Shoreline Matrix")
  tmp.Shore <- inEl
  tmp.Shore[inEl < sh.Z] <- NA # Set elevation values greater than 0.8 to NA
  tmp.Shore[!is.na(tmp.Shore)] <- 1 # Set elevation values greater than 0.8 to NA
  bndrys <- as.matrix(boundaries(raster(tmp.Shore), type='inner', asNA = T))
  # plot(rast(tmp.Shore), col='cadetblue4')
  # plot(rast(bndrys), add=T, col='black')

  out <- inEl * NA               # Output matrix
  if(sea == "E"){ # Go Right to Left
    for(i in 1:nrow(tmp.Shore)){
      tst <- 1
      tmp <- c(1:ncol(out)) * ((tmp.Shore[i,] * 0)+1)
      if(!is.infinite(suppressWarnings(max(tmp,na.rm=T)))){
        tmp <- max(tmp,na.rm=T)
        out[i,tmp] <- tst
      }
    }
  } else if(sea == "W"){ # Go Left to Right
    for(i in 1:nrow(tmp.Shore)){
      tst <- 1
      tmp <- c(1:ncol(out)) * ((tmp.Shore[i,] * 0)+1)
      if(!is.infinite(suppressWarnings(max(tmp,na.rm=T)))){
        tmp <- min(tmp,na.rm=T)
        out[i,tmp] <- tst
      }
    }
  } else if(sea == "N"){ # Go Top to Bottom
    for(i in 1:ncol(tmp.Shore)){
      tst <- 1
      tmp <- c(1:nrow(out)) * ((tmp.Shore[,i] * 0)+1)
      if(!is.infinite(suppressWarnings(max(tmp,na.rm=T)))){
        tmp <- min(tmp,na.rm=T)
        out[tmp, i] <- tst
      }
    }
  } else if(sea == "S"){ # Go Bottom to Top
    for(i in 1:ncol(tmp.Shore)){
      tst <- 1
      tmp <- c(1:nrow(out)) * ((tmp.Shore[,i] * 0)+1)
      if(!is.infinite(suppressWarnings(max(tmp,na.rm=T)))){
        tmp <- max(tmp,na.rm=T)
        out[tmp, i] <- tst
      }
    }
  } else{
    print("Invalid sea value: sea must be either N, S, E, or W")
  }

  out[out != 1] <- NA

  if(sea=="E"){
    for(i in 1:nrow(out)){
      if(is.infinite(suppressWarnings(max(out[i,], na.rm=T)))){
        tmp <- c(1:ncol(out)) * ((inEl[i,] * 0)+1)
        tmp <- max(tmp,na.rm=T)
        out[i,tmp] <- 1
      }
    }
  }
  if(sea=="W"){
    for(i in 1:nrow(out)){
      if(is.infinite(suppressWarnings(max(out[i,], na.rm=T)))){
        tmp <- c(1:ncol(out)) * ((inEl[i,] * 0)+1)
        tmp <- min(tmp,na.rm=T)
        out[i,tmp] <- 1
      }
    }
  }
  if(sea=="N"){
    for(i in 1:ncol(out)){
      if(is.infinite(suppressWarnings(max(out[,i], na.rm=T)))){
        tmp <- c(1:nrow(out)) * ((inEl[,i] * 0)+1)
        tmp <- min(tmp,na.rm=T)
        out[tmp,i] <- 1
      }
    }
  }
  if(sea=="S"){
    for(i in 1:ncol(out)){
      if(is.infinite(suppressWarnings(max(out[,i], na.rm=T)))){
        tmp <- c(1:nrow(out)) * ((inEl[,i] * 0)+1)
        tmp <- max(tmp,na.rm=T)
        out[tmp,i] <- 1
      }
    }
  }
  print("Done")
  return(out)
}
# Create a mask from the shoreline
shorelineMask <- function(inshore=shore, dist=200, sea=sea.Uni, inOM = T){
  # Find average shoreline angle and set ridge and swales ----------------------
  if(inOM){
    OM <- om.Finder(inXYZ2=shore.XYZ, uniq.Dist = seq(100,1000,50), sea = sea.Uni)
    OM <- mean(abs(OM$Angle), na.rm=T)
    if(OM < 45){
      R <<- RY #* shoreMask
      S <<- SY #* shoreMask
    }else{
      R <<- RX #* shoreMask
      S <<- SX #* shoreMask
    }
  }


  out <- inshore*NA
  inshore[is.na(inshore)] <- 0
  if(sea=="N"){
    for(i in 1:ncol(inshore)){
      if(max(inshore[,i],na.rm=T) != 0){
        tmp <- which((inshore[,i])==1)
        out[tmp:(tmp+dist),i] <- 1
      }
    }
  }
  if(sea=="S"){
    for(i in 1:ncol(inshore)){
      if(max(inshore[,i],na.rm=T) != 0){
        tmp <- which((inshore[,i])==1)
        out[(tmp-dist):tmp,i] <- 1
      }
    }
  }
  if(sea=="E"){
    for(i in 1:nrow(inshore)){
      if(max(inshore[i,],na.rm=T) != 0){
        tmp <- which((inshore[i,])==1)
        out[i, (tmp-dist):tmp] <- 1
      }
    }
  }
  if(sea=="W"){
    for(i in 1:nrow(inshore)){
      if(max(inshore[i,],na.rm=T) != 0){
        tmp <- which((inshore[i,])==1)
        out[i, tmp:(tmp+dist)] <- 1
      }
    }
  }
  out
}
# Find ridge and swales based on 2nd derivative of rr profile
inflect <- function(x, threshold = 1){
  up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
  down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
  a    <- cbind(x,up,down)
  list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
}
# -------------------------------------------------------------------------
# 1B ----------------------------------------------------------------------
preProcess.Morph <- function(makeRR=T, makeSN=T, inMn=7, snWn=5, shore.Ext=150){

  # Create Easting/Northing Matrices For Reference -------------------------------
  cat("\r", "Creating Northing and Easting Rasters...")
  XE <- D * NA
  YN <- D * NA

  ymx <- ymax(M)-(res(M)[1]/2)
  ymn <- ymin(M)+(res(M)[1]/2)

  xmx <- xmax(M)-(res(M)[1]/2)
  xmn <- xmin(M)+(res(M)[1]/2)
  YN.M <- seq(ymx, ymn, by=-1*l)
  XE.M <- seq(xmn, xmx, by=l)
  # length(XE.M)
  #
  # YN.M <- c((ymax(M) - (l/2)):(ymin(M) + (l/2)))
  # XE.M <- c((xmin(M) + (l/2)):(xmax(M) - (l/2)))
  # length()
  # Fill E/N Matrices
  for(i in 1:length(XE.M)){
    XE[,i] <- XE.M[i]
  }

  for(i in 1:length(YN.M)){
    YN[i,] <- YN.M[i]
  }

  # ------------------------------------------------------------------------------

  # CALCULATE RELATIVE RELIEF ----------------------------------------------------
  # Only run this one time, if you already have an RR skip this step
  # TODO: Probably need to make it a function
  shore.XYZ <- shore.Pos(inShore = shore,
                          sea = sea.Uni,
                          inX = XE,
                          inY = YN,
                          inZ = D)
  shore.XYZ <<- shore.XYZ
  shoreMask <- shorelineMask(inshore=shore, dist=shore.Ext, sea=sea.Uni, inOM = F)
  if(makeRR){
    cat("\n", "Creating Relative Relief Raster...")
    # rr.Make <- rr.Calc(inMat=D, inRas=M, mn.S=7, mx.S=Id.Sc, stp=2, save=T, fn=fn.Uni)
    # 495.66-Padre2009 Id.SC=39

    # rr <<- (rr.Make[[1]])
    rr <<- rr.Average(inMn=inMn, inMx=Id.Sc, x=D*shoreMask, inRas=M, save=T)
    # rr[is.na(D)] <<- NA
  }

  # SURFACE NORMAL RAW DATA ------------------------------------------------------
  if(makeSN){
    cat("\n", "Creating DEM Surface Normal Vectors...")
    features.Vec.DEM <<- RS.extract(inData = D, inRas = M, win.sz = snWn,
                                    store.SN = T, store.RS = T,
                                    type.Data = paste(fn.Uni,"_DEM_Raw"),
                                    rm.S=T, prefix = "DEM")
    cat("\n", "Creating RR Surface Normal Vectors...")
    features.Vec.RR <<- RS.extract(inData = rr, inRas = M, win.sz = snWn,
                                  store.SN = T, store.RS = T,
                                  type.Data = paste(fn.Uni,"_RR_Raw"),
                                  rm.S=T, prefix = "RR")
  }

  # 410.00-Padre2009?
  # ----------------------------------------------------------------------------
  # RX <<- features.Vec.DEM[[2]]$ridgeX
  # RY <<- features.Vec.DEM[[2]]$ridgeY
  # # RX <<- features.Vec.RR[[2]]$ridgeX
  # # RY <<- features.Vec.RR[[2]]$ridgeY
  # SX <<- features.Vec.RR[[2]]$swaleX
  # SY <<- features.Vec.RR[[2]]$swaleY
  # UX <- features.Vec.DEM[[1]]$unitX
  # UY <- features.Vec.DEM[[1]]$unitY
  # UZ <<- features.Vec.DEM[[1]]$unitZ
  RX <<- features.Vec.DEM[[2]][[1]]
  RY <<- features.Vec.DEM[[2]][[2]]
  # RX <<- features.Vec.RR[[2]]$ridgeX
  # RY <<- features.Vec.RR[[2]]$ridgeY
  SX <<- features.Vec.RR[[2]][[3]]
  SY <<- features.Vec.RR[[2]][[4]]
  UX <- features.Vec.DEM[[1]][[1]]
  UY <- features.Vec.DEM[[1]][[2]]
  UZ <<- features.Vec.DEM[[1]][[3]]

  if(sea.Uni == "N"){
    UX.o <- UY
    UY.o <- UX*-1
  }else if(sea.Uni=="E"){
    UX.o <- UX
    UY.o <- UY
  }else if(sea.Uni=="S"){
    UX.o <- UY * -1
    UY.o <- UX
  }else{
    UX.o <- UX * -1
    UY.o <- UY * -1
  }
  UX <<- UX.o
  UY <<- UY.o
  # ============================================================================
  # PROFILE EXTRACTION ---------------------------------------------------------
  # ============================================================================
  # Shoreline ------------------------------------------------------------------
  cat("\n", "Creating Shoreline...")
  shore <<- shore.Mat.Extract(D, mhw.Uni, sea.Uni)
  out = raster(shore,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)
  writeRaster(out,paste(fn.Uni,"_Full_shoreline.tif"), overwrite=TRUE)
  # 5.42-Padre2009
  # ----------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # shore.XYZ <<- shore.Pos(inShore = shore,
  #                         sea = sea.Uni,
  #                         inX = XE,
  #                         inY = YN,
  #                         inZ = D)
  XE <<- XE
  YN <<- YN
  shoreMask <- shorelineMask(inshore=shore, dist=shore.Ext, sea=sea.Uni, inOM = T)

  cat("\n", "Creating TPI...")
  tpi <- tpi.Average(inMn=mn.Sc, inMx=Id.Sc, x=D*shoreMask, inRas=M, save=T)
  cat("\n", "Creating TRI...")
  tri <- tri.Average(inMn=mn.Sc, inMx=Id.Sc, x=D*shoreMask, inRas=M, save=T)
  tpi <<- tpi
  tri <<- tri
  shoreMask <<- shoreMask
  # ------------------------------------------------------------------------------
  # cat("\n", "Extracting Profiles...")
  # prof.Extract <<- prof.Extract.init(inXYZ = shore.XYZ,
  #                                    sea = sea.Uni,
  #                                    win.sz = seq(100,1000, 50),
  #                                    prof.L = 150,
  #                                    edge.dist = 150,
  #                                    inRR = rr,
  #                                    inX = XE,
  #                                    inY = YN,
  #                                    inZ = D,
  #                                    inRX = RX,
  #                                    inRY = RY,
  #                                    inSX = SX,
  #                                    inSY = SY,
  #                                    inUX = UX,
  #                                    inUY = UY,
  #                                    inUZ = UZ)
  # # 82.69-Padre2009
  # # ------------------------------------------------------------------------------
  # # Find if the input feature for a profile is a ridge or swale from x/y ---------
  # # prof.Extract.copy <- prof.Extract
  # cat("\n", "Finding input feature directionality...")
  # prof.Extract$R <<- NA
  # prof.Extract$S <<- NA
  #
  # for(i in unique(prof.Extract$ID)){
  #   # cat('\r', "Extracting profile: ", i, "/", length(unique(prof.Extract$ID)))
  #   # flush.console()
  #   data.test <- prof.Extract[prof.Extract$ID == i,]
  #
  #   om.foc <- min(data.test$OM, na.rm=T)
  #
  #   if(abs(om.foc) < 45){
  #     data.test$R <- data.test$RY
  #     data.test$S <- data.test$SY
  #   }else {
  #     data.test$R <- data.test$RX
  #     data.test$S <- data.test$SX
  #   }
  #
  #   prof.Extract[prof.Extract$ID == i,] <<- data.test
  # }
  #
  # prof.R <- prof.Extract
  # prof.R$x <- prof.R$x * prof.R$R
  # prof.R$y <- prof.R$y * prof.R$R
  #
  #
  # prof.R <- prof.R[complete.cases(prof.R$x),]
  #
  # prof.S <- prof.Extract
  # prof.S$x <- prof.S$x * prof.S$S
  # prof.S$y <- prof.S$y * prof.S$S
  # # prof.S$S[prof.S$S==0] <- NA
  # prof.S <- prof.S[complete.cases(prof.S$x),]
  #
  # write.csv(prof.R, "Profile_Ridges.csv")
  # write.csv(prof.S, "Profile_Swales.csv")
  #
  # write.csv(prof.Extract, "Profile_Data.csv")
  #
  # prof.Extract$R[prof.Extract$R!=1] <<- NA
  # prof.Extract$S[prof.Extract$S!=1] <<- NA

  cat("\n", "...Finished preprocessing.")
}
# Create a dataframe that contains the X/Y/Z for each shoreline position
shore.Pos <- function(inShore = shore, sea = NA, inX = XE, inY = YN, inZ = D){
  print("Finding Shore Position")
  # Initialize vectors to use to fill data frame
  tmp.X <- NA
  tmp.Y <- NA
  tmp.Z <- NA
  tmp.i <- NA
  tmp.j <- NA

  ct <- 1
  # Find the shore positions
  for(i in 1:nrow(inShore)){
    for(j in 1:ncol(inShore)){
      if(!is.na(inShore[i,j])){
        if(inShore[i,j] == 1){
          tmp.X[ct] <- inX[i,j]
          tmp.Y[ct] <- inY[i,j]
          tmp.Z[ct] <- inZ[i,j]
          tmp.i[ct] <- i
          tmp.j[ct] <- j
          ct <- ct + 1
        }
      }

    }
  }

  # create dataframe containing positions
  tmp.XYZ <- data.frame(x=tmp.X,
                        y=tmp.Y,
                        z=tmp.Z,
                        i=tmp.i,
                        j=tmp.j)

  # Based on seaward direction sort the data frame and add ID no.
  if(sea == "E" | sea == "W"){ # East or west sort by Northing from large to small
    tmp.XYZ <- tmp.XYZ[order(tmp.XYZ$y, decreasing = T, na.last = F),]
    tmp.XYZ$ID <- c(1:length(tmp.XYZ$y))
    # tmp.XYZ <- 1
  } else if(sea == "N" | sea == "S"){ # North or South sort by Easting small to large
    tmp.XYZ <- tmp.XYZ[order(tmp.XYZ$x, decreasing = F, na.last = F),]
    tmp.XYZ$ID <- c(1:length(tmp.XYZ$y))
  } else if(is.na(sea)){
    print("Invalid sea value: sea must be either N, S, E, or W")
    tmp.XYZ <- NA
  } else {
    print("Invalid sea value: sea must be either N, S, E, or W")
    tmp.XYZ <- NA
  }
  print("Done")
  return(tmp.XYZ)

}
# Average RR function
rr.Average <- function(inMn=7, inMx=Id.Sc, x=D, inRas=M, save=T){
  seq.RR <- seq(inMn, inMx, 2)
  m = (seq.RR - 1) / 2
  mid = ((inMx-1) / 2)
  nr = nrow(x) - mid
  nc = ncol(x) - mid
  out <- x*NA
  for(j in mid:(nc-mid)){
    # cat("\r", j, "/", nc,sep="")
    for(i in mid:(nr-mid)){

      if(!is.na(x[i,j]) & x[i,j] > 0){
        rr.here <- NA
        for(k in 1:length(m)){
          tmp <- c(x[(i-m[k]):(i + m[k]),(j-m[k]):(j + m[k])])
          rr.here[k] <- (x[i,j] - min(tmp, na.rm=T))/(max(tmp, na.rm=T) - min(tmp, na.rm=T))
        }
        out[i,j] = mean(rr.here, na.rm=T)

      }
    }
  }
  # Save surface normal vectors if condition met --------------------------------------------
  if(save){
    # Set name of new output folder
    stat.Dirs <- c(paste(inMn, "_", inMx,"_", fn.Uni, sep=""))
    # Check if folder already exists, if not create it
    for(d in 1:length(stat.Dirs)){
      ifelse(!dir.exists(stat.Dirs[d]), dir.create(stat.Dirs[d]), "Folder exists already")
    }
    # Set new output folder as the working directory
    currDir <- getwd()
    dirRR = paste(currDir, "\\", stat.Dirs,sep="")
    setwd(dirRR)

    # Save average RR as tif
    outRR = raster(out,xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
    projection(outRR) = crs(inRas)
    writeRaster(outRR, filename=paste(stat.Dirs, "_Mean.tif"), overwrite = T)

    # Set the working directory back to the original main directory
    setwd(currDir)
  }

  return(out)
}
# Calculate both the surface normal vectors and ridge and swale components of a surface
RS.extract <- function(inData, inRas, win.sz, store.SN = F, store.RS = F,
                       type.Data = "NA", rm.S = F, prefix=NA){
  # Metadata of the raster used for SurfNorm
  NAvalue(inRas) <- -9999
  nil <- NAvalue(inRas)

  res <- res(inRas)
  L = res[1]

  inData[is.na(inData)] = nil
  # Calculate the components of the surface normal vectors ----------------------------------
  vv = win.sz                            # Window size
  vv2 = vv - 1                           # distance between corners
  surf.x = surfNormX(inData, vv2*L, vv, vv)      # x-component
  surf.y = surfNormY(inData, vv2*L, vv, vv)      # y-component
  surf.El = surfNormEl(inData, 2*L, vv, vv)      # average elevation of corner cells

  surf.z = matrix((vv2*L)^2,nrow(inData),ncol(inData))     # z-component
  surf.z[surf.x==0] <- 0

  den = sqrt((surf.x^2) + (surf.y^2) + (surf.z^2))      # magnitude of the surf. norm. vector

  unitX = surf.x / den                        # x-component of the unit vector
  unitY = surf.y / den                        # y-component of the unit vector
  unitZ = surf.z / den                        # z-component of the unit vector
  NunitZ <- norm.data(unitZ, method = "minmax", mn = NA, mx = NA) # normalized z-component of the unit vec.
  # -----------------------------------------------------------------------------------------
  # Calculate ridge and swale features ------------------------------------------------------
  rs.features <- ridge.swale(inM = inData, inX = unitX, inY = unitY)
  # -----------------------------------------------------------------------------------------
  if(rm.S){
    RS1 <- single.RM(rs.features[[1]], rs.features[[2]])
    RS2 <- single.RM(rs.features[[3]], rs.features[[4]])

    rs.features[[1]] <- RS1[[1]]
    rs.features[[2]] <- RS1[[2]]
    rs.features[[3]] <- RS2[[1]]
    rs.features[[4]] <- RS2[[2]]

  }

  # Save surface normal vectors if condition met --------------------------------------------
  if(store.SN){

    stat.Dirs <- c(paste(vv, "x", vv,"_Surface_Normal_Vectors_", type.Data, sep=""))

    for(d in 1:length(stat.Dirs)){
      ifelse(!dir.exists(stat.Dirs[d]), dir.create(stat.Dirs[d]), "Folder exists already")
    }
    currDir <- getwd()
    dirSN = paste(currDir, "\\", stat.Dirs,sep="")
    setwd(dirSN)

    outEl = raster(surf.El,xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
    projection(outEl) = crs(inRas)

    outX = raster(unitX,xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
    projection(outX) = crs(inRas)

    outY = raster(unitY,xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
    projection(outY) = crs(inRas)

    outZ = raster(unitZ,xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
    projection(outZ) = crs(inRas)

    outNZ = raster(NunitZ,xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
    projection(outNZ) = crs(inRas)


    s.Vec <- stack(c(outX, outY, outZ, outNZ, outEl))
    names(s.Vec) <- paste0(prefix, "_", c("unitX", "unitY", "unitZ", "norm_unitZ", "unitEl"))
    writeRaster(s.Vec, filename=names(s.Vec), bylayer=T, format="GTiff", overwrite = T)

    comp <- c(rast(outX), rast(outY), rast(outZ))
    writeRaster(s.Vec, filename=paste0(prefix, "_", "CompositeXYZ"), format="GTiff", overwrite = T)
    setwd(currDir)
  }
  # -----------------------------------------------------------------------------------------
  # Save Ridge and Swale features if condition met ------------------------------------------
  if(store.RS){

    stat.Dirs <- c(paste(vv, "x", vv,"_Ridge_and_Swale_", type.Data, sep=""))

    for(d in 1:length(stat.Dirs)){
      ifelse(!dir.exists(stat.Dirs[d]), dir.create(stat.Dirs[d]), "Folder exists already")
    }

    currDir <- getwd()
    dirSN = paste(currDir, "\\", stat.Dirs,sep="")
    setwd(dirSN)

    nms <- paste0(prefix, "_", c("ridgeX", "ridgeY", "swaleX", "swaleY",
                                 "ridgeXY", "swaleXY", "ridge.swaleX", "ridge.swaleY"))

    for(i in 1:length(rs.features)){
      out = raster(rs.features[[i]],xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
      projection(out) = crs(inRas)

      fn <- paste(vv, "x", vv, "_", nms[i], ".tif", sep = "")
      writeRaster(out, fn, overwrite = T)
    }

    setwd(currDir)
  }
  # -----------------------------------------------------------------------------------------
  SN.out <- list(unitX, unitY, unitZ, NunitZ)
  names(SN.out) <- paste0(prefix, "_", c("unitX", "unitY", "unitZ", "NunitZ"))

  return(list(SN.out, rs.features))

}
# Using a DEM and window size calculate the surface normal vector components
surfNormX <- function(inx, len, rR=3, cC=3){

  wr = rR - 1
  wc = cC - 1

  nrow = nrow(inx)
  ncol = ncol(inx)

  nrow1 = nrow - wr
  ncol1 = ncol - wc

  outX = matrix(0,nrow,ncol)

  for (i in 1:nrow1){
    for (j in 1:ncol1){

      a = inx[i + wr, j]
      b = inx[i + wr, j + wc]
      c = inx[i,j]
      d = inx[i, j + wc]

      if((a!=nil)&(b!=nil)&(c!=nil)&(d!=nil)){

        va = b - a
        vb = c - a
        vc = b - d
        vd = c - d

        outX[i+(wr/2),j+(wc/2)] = .5*len*(vd - va)
      }
    }
  }
  outX
}
surfNormY <- function(inx, len, rR=3, cC=3){

  wr = rR - 1
  wc = cC - 1

  nrow = nrow(inx)
  ncol = ncol(inx)

  nrow1 = nrow - wr
  ncol1 = ncol - wc

  outY = matrix(0,nrow,ncol)

  for (i in 1:nrow1){
    for (j in 1:ncol1){

      a = inx[i + wr, j]
      b = inx[i + wr, j + wc]
      c = inx[i,j]
      d = inx[i, j + wc]

      if((a!=nil)&(b!=nil)&(c!=nil)&(d!=nil)){
        va = b - a
        vb = c - a
        vc = b - d
        vd = c - d

        outY[i+(wr/2),j+(wc/2)] = .5*len*(vc - vb)
      }
    }
  }
  outY
}
surfNormEl <- function(inx, len, rR=3, cC=3){

  wr = (rR/2) + 0.5
  wc = (cC/2) + 0.5

  nrow = nrow(inx)
  ncol = ncol(inx)

  nrow1 = nrow - wr
  ncol1 = ncol - wc

  outEl = matrix(0,nrow,ncol)

  for (i in 1:nrow1){
    for (j in 1:ncol1){

      a = inx[i + wr, j]
      b = inx[i + wr, j + wc]
      c = inx[i,j]
      d = inx[i, j + wc]

      if((a!=nil)&(b!=nil)&(c!=nil)&(d!=nil)){

        outEl[i+(wr/2),j+(wc/2)] = (a+b+c+d)/4

      }
    }
  }
  outEl
}
# Calculate where ridge and swale features are located
ridge.swale <- function(inM = D, inX = unitX, inY = unitY){
  # variables for processing extent
  nrow = nrow(inM)
  ncol = ncol(inM)

  nrow1 = nrow - 1
  ncol1 = ncol - 1

  nrow2 = nrow - 2
  ncol2 = ncol - 2

  ncol3 = ncol - 3
  nrow3 = nrow - 3

  # # Initialize matrices for storing the first pass data
  # Mult_x = matrix(0,nrow,ncol)
  # Mult_y = matrix(0,nrow,ncol)
  #
  # # For a 3 cell window set the center value to the two end values mutliplied together
  # for(i in 1:nrow){ # X-direction: horizontal window within the x-component of the surf.vec
  #   for(j in 1:ncol2){
  #     Mult_x[i,j+1] = x[i,j] * x[i,j+2]
  #   }
  # }
  #
  # for(i in 1:nrow2){ # Y-direction: vertical window within the y-component of the surf.vec
  #   for(j in 1:ncol){
  #     Mult_y[i+1,j] = inY[i,j] * inY[i+2,j]
  #   }
  # }

  # Set values of the surface normal vectors to +/- 1 based on their sign
  px = inX
  px[px > 0] = 1
  px[px < 0] = -1

  py = inY
  py[py > 0] = 1
  py[py < 0] = -1

  # Initialize ridge and swale matrices for both directional components
  Mult_xx = matrix(0,nrow,ncol)
  Mult_yy = matrix(0,nrow,ncol)

  px[is.na(px)] = 0
  py[is.na(py)] = 0

  # Find where neighboring cells are pos. and neg. and set the value to a ridge or swale
  # based on orientation and elevation value
  # x-component
  for(i in 1:nrow){
    for(j in 1:ncol3){
      # Ridge
      if((px[i,j] + px[i,j+1] < 0) & (px[i,j+2] + px[i,j+3] > 0)){
        if(inM[i,j+1 ] > inM[i,j+2]){
          Mult_xx[i,j+1] = 2
        }else if(inM[i,j+1] < inM[i,j+2]){
          Mult_xx[i,j+2] = 2
        }
      }
      # Swale
      else if((px[i,j] + px[i,j+1] > 0) & (px[i,j+2] + px[i,j+3]<0)){
        if(inM[i,j+1] < inM[i,j+2]){
          Mult_xx[i,j+1] = 1
        }else if(inM[i,j+1] > inM[i,j+2]){
          Mult_xx[i,j+2] = 1
        }
      }
    }
  }
  # y-component
  for(i in 1:nrow3){
    for(j in 1:ncol){
      # Ridge
      if((py[i,j] + py[i+1,j] > 0) & (py[i+2,j] + py[i+3,j] < 0)){
        if(inM[i+1,j] > inM[i+2,j]){
          Mult_yy[i+1,j] = 2
        }else if(inM[i+1,j] < inM[i+2,j]){
          Mult_yy[i+2,j] = 2
        }
      }
      # Swale
      else if((py[i,j] + py[i+1,j]<0) & (py[i+2,j] + py[i+3,j] > 0)){
        if(inM[i+1,j] < inM[i+2,j]){
          Mult_yy[i+1,j] = 1
        }else if(inM[i+1,j] > inM[i+2,j]){
          Mult_yy[i+2,j] = 1
        }
      }
    }
  }


  ##############################################################################

  cMultxx = Mult_xx
  cMultyy = Mult_yy

  tMultxx = Mult_xx
  tMultyy = Mult_yy

  cMultxx[cMultxx<2] = 0
  cMultyy[cMultyy<2] = 0

  cMultxx[cMultxx>0] = 1
  cMultyy[cMultyy>0] = 1

  tMultxx[tMultxx>1] = 0
  tMultyy[tMultyy>1] = 0

  cMultxy = cMultxx + cMultyy
  cMultxy[cMultxy>0] = 1

  tMultxy = tMultxx + tMultyy
  tMultxy[tMultxy>0] = 1

  ##############################################################################

  cMultxx[cMultxx == 0] = NA
  cMultyy[cMultyy == 0] = NA

  tMultxx[tMultxx == 0] = NA
  tMultyy[tMultyy == 0] = NA

  cMultxy[cMultxy == 0] = NA
  tMultxy[tMultxy == 0] = NA

  Mult_xx[Mult_xx == 0] = NA
  Mult_yy[Mult_yy == 0] = NA

  out.List <- list(cMultxx, cMultyy, tMultxx, tMultyy, cMultxy, tMultxy, Mult_xx, Mult_yy)
  names(out.List) <- c("ridgeX", "ridgeY", "swaleX", "swaleY",
                       "ridgeXY", "swaleXY", "ridge.swaleX", "ridge.swaleY")
  return(out.List)
}
# Function for normalizing data
norm.data <- function(inp, method = "minmax", mn = NA, mx = NA){
  if(method == "minmax"){
    out <- (inp - min(inp, na.rm=T)) / (max(inp, na.rm=T) - min(inp, na.rm=T))
  }
  else if(method == "minmaxcus"){
    if(is.na(mx) | is.na(mn)){
      print("ERROR: mn or mx not defined")
      out <- NA
    }else {
      out <- (inp - min(inp, na.rm=T)) / (max(inp, na.rm=T) - min(inp, na.rm=T))
      out <- (mx - mn) * out
      out <- out + mn
    }
  }
  else if(method == "standardize"){
    out <- (inp - mean(inp, na.rm=T)) / sd(inp, na.rm=T)
  }

}
# Remove single cell locations
single.RM <- function(Ridge, Swale){

  Ridge[is.na(Ridge)] <- 0
  Swale[is.na(Swale)] <- 0
  # For ridge and swale features
  for(i in 2:(nrow(Ridge) - 1)){

    # Track progress
    # progress(i, (nrow(Ridge) - 1))#, progress.bar = TRUE)
    # Sys.sleep(0.01)
    # if (i == (nrow(Ridge) - 1)) message("Ridge Done!")

    for(j in 2:(ncol(Ridge) - 1)){
      if(Ridge[i,j] == 1){

        tmp <- sum(c(Ridge[(i-1):(i+1),(j-1):(j+1)]))

        if(tmp == 1){
          Ridge[i,j] = 0
        }

      }
    }
  }
  Ridge[Ridge==0] <- NA

  # For swale features
  for(i in 2:(nrow(Swale) - 1)){

    # Track progress
    # progress(i, (nrow(Swale) - 1))#, progress.bar = TRUE)
    # Sys.sleep(0.01)
    # if (i == (nrow(Swale) - 1)) message("Swale Done!")


    for(j in 2:(ncol(Swale) - 1)){

      if(Swale[i,j] == 1){

        tmp <- sum(c(Swale[(i-1):(i+1),(j-1):(j+1)]))

        if(tmp == 1){
          Swale[i,j] = 0
        }

      }
    }
  }
  Swale[Swale==0] <- NA
  return(list(Ridge,Swale))
}
# For a set distance find the angle tangent to the sequence of shoreline positions
cons.Ang <- function(inXYZ2 = shore.XYZ, uniq.Dist = 100, sea = "E"){

  # u.S <- seq(1, nrow(inXYZ2), uniq.Dist)
  #
  # tmp.ang <- inXYZ2$x*NA
  # ct <- 1
  # for(i in 2:length(u.S)){
  #   section.S <- u.S[i-1]
  #   section.E <- u.S[i] - 1
  #
  #   sec.X <- inXYZ2$x[section.E] - inXYZ2$x[section.S]
  #   sec.Y <- inXYZ2$y[section.E] - inXYZ2$y[section.S]
  #   om <- atan(sec.Y / sec.X) * 180/pi
  #
  #   for(j in section.S:section.E){
  #     tmp.ang[j] <- om
  #   }
  # }
  # tmp.ang[is.na(tmp.ang)] = -90
  # inXYZ2$om <- tmp.ang

  u.S <- seq(1, nrow(inXYZ2), uniq.Dist)

  tmp.ang <- inXYZ2$x*NA
  tmp.r2 <- inXYZ2$x*NA
  ct <- 1
  for(i in 2:length(u.S)){
    # i = 2
    section.S <- u.S[i-1]
    section.E <- u.S[i] - 1

    sec.X <- inXYZ2$x[section.S:section.E]
    sec.Y <- inXYZ2$y[section.S:section.E]

    lm.fit <- lm(sec.Y~sec.X)
    slope <- lm.fit$coefficients[[2]]
    r2 <- summary(lm.fit)$adj.r.squared
    # plot(sec.X, sec.Y, pch=19)
    # abline(a=lm.fit$coefficients[[1]], b=lm.fit$coefficients[[2]])
    # sec.X <- inXYZ2$x[section.E] - inXYZ2$x[section.S]
    # sec.Y <- inXYZ2$y[section.E] - inXYZ2$y[section.S]
    # om <- atan(sec.Y / sec.X) * 180/pi
    om <- atan(slope) * 180/pi
    for(j in section.S:section.E){
      tmp.ang[j] <- om
      tmp.r2[j] <- r2
    }
  }
  tmp.ang[is.na(tmp.ang)] = -90
  inXYZ2$om <- tmp.ang
  inXYZ2$r2 <- tmp.r2

  return(inXYZ2)
}
# For the current section get the segments for all scales
om.Finder <- function(inXYZ2 = shore.XYZ, uniq.Dist = seq(100,1000, 50), sea = "N"){
  #Initialize
  dis.List <- list()
  df.info <- data.frame()
  tmp <- c(NA, NA,1,NA,NA,NA)
  # rm(df.info)

  ct <- 1
  while(ct < 200){
    # cat("\r", ct, "/", 200)

    # Check if you are at the end of the dataset
    if(tmp[3] + 1 < nrow(inXYZ2)){
      strt <- tmp[3] + 1
    }else {
      break
    }
    # strt = 1
    for(i in 1:length(uniq.Dist)){
      dis.List[[i]] <- c(strt,(strt + uniq.Dist[i] - 1))
    }


    tmp.ang <- NA
    tmp.r2 <- NA
    tmp.mae <- NA
    tmp.rmse <- NA

    # For each scale find the angle and R2 values for the given segments

    # par(mfrow=c(4,5))
    for(i in 1:length(dis.List)){
      # i = 8
      u.S <- dis.List[[i]]

      sec.X <- inXYZ2$x[u.S[1]:u.S[2]]
      sec.Y <- inXYZ2$y[u.S[1]:u.S[2]]
      sec.X <- na.omit(sec.X)
      sec.Y <- na.omit(sec.Y)
      lm.fit <- lm(sec.Y~sec.X)
      slope <- lm.fit$coefficients[[2]]
      r2 <- summary(lm.fit)$adj.r.squared

      mae <- sum(abs(sec.Y - lm.fit$fitted.values))/length(sec.Y)
      rmse <- mean((sec.Y - lm.fit$fitted.values)^2, na.rm=T) %>% sqrt()
      # plot(sec.X, sec.Y, pch=19, cex=0.4, xaxt="n",yaxt="n")
      # abline(a=lm.fit$coefficients[[1]], b=lm.fit$coefficients[[2]],lwd=2,col="blue")
      # title(paste(paste("R2: ", round(r2, digits=2), " | MAE: ", round(mae, digits=2), " | RMSE: ", round(rmse, digits=2))),
      #       cex.main=0.8)

      om <- atan(slope) * 180/pi
      tmp.ang[i] <- om
      tmp.r2[i] <- r2
      tmp.mae[i] <- mae
      tmp.rmse[i] <- rmse

    }

    # Which scale gives max R2 value
    ind <- min(which(tmp.r2 == max(tmp.r2, na.rm=T)))

    # Output information about this data point
    tmp <- c(ind, dis.List[[ind]], tmp.ang[ind], tmp.r2[i], uniq.Dist[ind])

    # Initialize/Fill dataframe
    if(nrow(df.info)>=1){
      df.info <- rbind(df.info, tmp)
    }else{
      tmp <- matrix(tmp, nr=1)
      df.info <- as.data.frame(tmp)
      names(df.info) <- c("Indice", "Start", "End", "Angle", "R2", "Scale")
    }
    # Check if you have exceeded the bounds of the input dataset
    df.info$End[df.info$End > nrow(inXYZ2)] <- nrow(inXYZ2)
    ct <- ct + 1
  }
  return(df.info)

}
# Average TPI
tpi.Average <- function(inMn=7, inMx=Id.Sc, x=D, inRas=M, save=T){
  seq.RR <- seq(inMn, inMx, 2)
  m = (seq.RR - 1) / 2
  mid = ((inMx-1) / 2)
  nr = nrow(x) - mid
  nc = ncol(x) - mid
  out <- x*NA
  for(j in mid:(nc-mid)){
    # cat("\r", j, "/", nc,sep="")
    for(i in mid:(nr-mid)){

      if(!is.na(x[i,j]) & x[i,j] > 0){
        rr.here <- NA
        for(k in 1:length(m)){
          tmp <- c(x[(i-m[k]):(i + m[k]),(j-m[k]):(j + m[k])])
          rr.here[k] <- (x[i,j] - mean(tmp, na.rm=T))
        }
        out[i,j] = mean(rr.here, na.rm=T)

      }
    }
  }
  # Save surface normal vectors if condition met --------------------------------------------
  if(save){
    # Set name of new output folder
    stat.Dirs <- c(paste(inMn, "_", inMx,"_", fn.Uni, sep=""))
    # Check if folder already exists, if not create it
    for(d in 1:length(stat.Dirs)){
      ifelse(!dir.exists(stat.Dirs[d]), dir.create(stat.Dirs[d]), "Folder exists already")
    }
    # Set new output folder as the working directory
    currDir <- getwd()
    dirRR = paste(currDir, "\\", stat.Dirs,sep="")
    setwd(dirRR)

    # Save average RR as tif
    outRR = raster(out,xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
    projection(outRR) = crs(inRas)
    writeRaster(outRR, filename=paste(stat.Dirs, "_Mean_TPI.tif"), overwrite = T)

    # Set the working directory back to the original main directory
    setwd(currDir)
  }

  return(out)
}
# Average TRI
tri.Average <- function(inMn=7, inMx=Id.Sc, x=D, inRas=M, save=T){
  seq.RR <- seq(inMn, inMx, 2)
  m = (seq.RR - 1) / 2
  mid = ((inMx-1) / 2)
  nr = nrow(x) - mid
  nc = ncol(x) - mid
  out <- x*NA
  for(j in mid:(nc-mid)){
    # cat("\r", j, "/", nc,sep="")
    for(i in mid:(nr-mid)){

      if(!is.na(x[i,j]) & x[i,j] > 0){
        rr.here <- NA
        for(k in 1:length(m)){
          tmp <- c(x[(i-m[k]):(i + m[k]),(j-m[k]):(j + m[k])])
          rr.here[k] <- sum(abs(x[i,j] - tmp), na.rm=T)/(length(tmp)-1)
        }
        out[i,j] = mean(rr.here, na.rm=T)

      }
    }
  }
  # Save surface normal vectors if condition met --------------------------------------------
  if(save){
    # Set name of new output folder
    stat.Dirs <- c(paste(inMn, "_", inMx,"_", fn.Uni, sep=""))
    # Check if folder already exists, if not create it
    for(d in 1:length(stat.Dirs)){
      ifelse(!dir.exists(stat.Dirs[d]), dir.create(stat.Dirs[d]), "Folder exists already")
    }
    # Set new output folder as the working directory
    currDir <- getwd()
    dirRR = paste(currDir, "\\", stat.Dirs,sep="")
    setwd(dirRR)

    # Save average RR as tif
    outRR = raster(out,xmn=xmin(inRas),xmx=xmax(inRas),ymn=ymin(inRas),ymax(inRas))
    projection(outRR) = crs(inRas)
    writeRaster(outRR, filename=paste(stat.Dirs, "_Mean_TRI.tif"), overwrite = T)

    # Set the working directory back to the original main directory
    setwd(currDir)
  }

  return(out)
}
# Plot input features
plot.Inputs <- function(){
  pdf(paste0(fn.Uni, "_Inputs.pdf"))
  plot(raster(D), col=grey(1:100/100), main="Elevation")
  plot(raster(rr), col=grey(1:100/100), main="RR")
  plot(raster(UX), col=colorspace::divergex_hcl(255), main="UX")
  plot(raster(UY), col=colorspace::divergex_hcl(100), main="UY")
  plot(raster(UZ), col=viridis::plasma(100, direction=-1), main="UZ")
  plot(raster(tpi), col=viridis::plasma(100, direction=1), main="TPI")
  plot(raster(tri), col=viridis::plasma(100, direction=1), main="TRI")
  dev.off()
}
# Create i and j mats
createIJMat <- function(){
  iMat <- D*NA
  jMat <- D*NA
  for(i in 1:(nrow(D))){
    iMat[i,] <- i
  }
  for(j in 1:(ncol(D))){
    jMat[,j] <- j
  }
  iMat <<- iMat
  jMat <<- jMat
  dirPre <- paste(dir,"\\",fn.Uni, "_Preprocess", sep="")
  setwd(dirPre)

  out = raster(iMat,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)

  writeRaster(out,paste0(fn.Uni, "_iMat.tif"), overwrite=TRUE)

  out = raster(jMat,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)

  writeRaster(out,paste0(fn.Uni, "_jMat.tif"), overwrite=TRUE)
  setwd(dir)
}
# Get line features
get.RS.Features <- function(){
  if(sea.Uni=="E" || sea.Uni=="W"){
    inXY2 <<- "X"
  }else{
    inXY2 <<- "Y"
  }

  shoreMask <- shorelineMask(inshore=shore, dist=d2s, sea=sea.Uni)
  R.features <- feat.Extract(el=D, ft=R*shoreMask, in.IT=T, inXY=inXY2)
  S.features <- feat.Extract(el=D, ft=S*shoreMask, in.IT=T, inXY=inXY2)

  dirPre <- paste(dir,"\\",fn.Uni, "_Preprocess", sep="")
  setwd(dirPre)

  R.id.out <- R.features[[2]]
  R.id.out[R.id.out==0] <- NA
  out = raster(R.id.out,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)

  writeRaster(out,paste0(fn.Uni, "_R_IDs.tif"), overwrite=TRUE)

  S.id.out <- S.features[[2]]
  S.id.out[S.id.out==0] <- NA
  out = raster(S.id.out,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)

  writeRaster(out,paste0(fn.Uni, "_S_IDs.tif"), overwrite=TRUE)
  setwd(dir)

  R.features <<- R.features
  S.features <<- S.features

}
# Follow the path of the feature segment to get feature metrics
path.Follow <- function(ii, jj, inel, inft, ink, it, inid, inadj=T, inID, xy = "X"){

  x.vec <- NA # i indices of segment locations
  y.vec <- NA # j indices of segment locations

  k = 1 # starting length

  # while there is adjacency continue to loop
  while(inadj){

    # For MultXX Dataset
    if(xy == "X"){
      # Adjacent cells in the correct direction of path
      # A = inft[ii, (jj - 1)]
      B1 = inft[(ii + 1), (jj - 2)]
      B = inft[(ii + 1), (jj - 1)]
      C = inft[(ii + 1), jj]
      D = inft[(ii + 1), (jj + 1)]
      D1 = inft[(ii + 1), (jj + 2)]

      B2 = inft[(ii + 2), (jj - 1)]
      C2 = inft[(ii + 2), jj]
      D2 = inft[(ii + 2), (jj + 1)]
      # E = inft[ii, (jj + 1)]

      # adj.Cells = c(A, B, C, D, E)
      # tmp.I = c(ii, (ii+1), (ii+1), (ii+1), ii)
      # tmp.J = c((jj-1), (jj-1), jj, (jj+1), (jj+1))
      adj.Cells = c(B, C, D, B1, D1)
      adj.Cells2 = c(B2, C2, D2)
      tmp.I = c((ii+1), (ii+1), (ii+1), (ii+1), (ii+1))
      tmp.J = c((jj-1), jj, (jj+1), (jj-2), (jj+2))
      tmp.I2 = c((ii+2), (ii+2), (ii+2))
      tmp.J2 = c((jj-1), jj, (jj+1))

      sum.Adj = sum(adj.Cells)
      sum.Adj2 = sum(adj.Cells2)
    }

    # For MultYY Dataset
    if(xy == "Y"){
      # Adjacent cells in the correct direction of path
      # A = inft[(ii + 1), jj]
      B1 = inft[(ii + 2), (jj + 1)]
      B = inft[(ii + 1), (jj + 1)]
      C = inft[ii, (jj + 1)]
      D = inft[(ii - 1), (jj + 1)]
      D1 = inft[(ii - 2), (jj + 1)]

      B2 = inft[(ii + 1), (jj + 2)]
      C2 = inft[ii, (jj + 2)]
      D2 = inft[(ii - 1), (jj + 2)]
      # E = inft[(ii - 1), jj]

      # adj.Cells = c(A, B, C, D, E)
      adj.Cells = c(B, C, D, B1, D1)
      adj.Cells2 = c(B2, C2, D2)
      # tmp.I = c((ii+1), (ii+1), ii, (ii-1), (ii-1))
      # tmp.J = c(jj, (jj+1), (jj+1), (jj+1), jj)
      tmp.I = c((ii+1), ii, (ii-1), (ii+2), (ii-2))
      tmp.J = c((jj+1), (jj+1), (jj+1), (jj+1), (jj+1))
      tmp.I2 = c((ii+1), ii, (ii-1))
      tmp.J2 = c((jj+2), (jj+2), (jj+2))
      sum.Adj = sum(adj.Cells)
      sum.Adj2 = sum(adj.Cells2)
    }

    if(is.na(sum.Adj2)){sum.Adj2<-0}
    if(is.na(sum.Adj)){sum.Adj<-0}
    # There is adjacency
    if(sum(sum.Adj) != 0){
      sum.Adj2 <- 0
      # single adjacent cell
      if(sum(sum.Adj) == 1){
        # index of adjacent cell
        cell.ID = match(1, adj.Cells)

        # fill x and y vectors with the i and j indices of the current feature segment location
        x.vec[k] = ii
        y.vec[k] = jj

        # update length value, k
        k = k + 1

        # update the i and j indices to reflect the next segment cell
        ii = tmp.I[cell.ID]
        jj = tmp.J[cell.ID]

      }

      # more than one adjacent cell
      else if(sum(sum.Adj) > 1){
        # indices of adjacent cells
        cell.ID = which(adj.Cells %in% 1)

        # calculate difference in adjacent elev to start
        diff = NA
        for(i in 1:length(cell.ID)){
          diff[i] = abs(inel[ii,jj] - inel[tmp.I[cell.ID[i]], tmp.J[cell.ID[i]]])
        }

        cell.ID = match(min(diff, na.rm=T), diff)

        # fill x and y vectors with the i and j indices of the current feature segment location
        x.vec[k] = ii
        y.vec[k] = jj

        # update length value, k
        k = k + 1

        # update the i and j indices to reflect the next segment cell
        ii = tmp.I[cell.ID]
        jj = tmp.J[cell.ID]
      }
    }
    else if(sum(sum.Adj2) != 0 && sum(sum.Adj) == 0){

      # single adjacent cell
      if(sum(sum.Adj2) == 1){
        # index of adjacent cell
        cell.ID = match(1, adj.Cells2)

        # fill x and y vectors with the i and j indices of the current feature segment location
        x.vec[k] = ii
        y.vec[k] = jj

        # update length value, k
        k = k + 1

        # update the i and j indices to reflect the next segment cell
        ii = tmp.I2[cell.ID]
        jj = tmp.J2[cell.ID]

      }

      # more than one adjacent cell
      else if(sum(sum.Adj2) > 1){
        # indices of adjacent cells
        cell.ID = which(adj.Cells2 %in% 1)

        # calculate difference in adjacent elev to start
        diff = NA
        for(i in 1:length(cell.ID)){
          diff[i] = abs(inel[ii,jj] - inel[tmp.I2[cell.ID[i]], tmp.J2[cell.ID[i]]])
        }

        cell.ID = match(min(diff, na.rm=T), diff)

        # fill x and y vectors with the i and j indices of the current feature segment location
        x.vec[k] = ii
        y.vec[k] = jj

        # update length value, k
        k = k + 1

        # update the i and j indices to reflect the next segment cell
        ii = tmp.I2[cell.ID]
        jj = tmp.J2[cell.ID]
      }
    }
    # There is not adjacency
    # else if(k > 50){
    #   # fill x and y vectors with the i and j indices of the feature segment endpoint location
    #   x.vec[k] = ii
    #   y.vec[k] = jj
    #
    #   for(i in 1:length(x.vec)){
    #
    #     # fill length and id matrices for the feature segment
    #     ink[x.vec[i], y.vec[i]] = k
    #
    #     if(inID){
    #       it[x.vec[i], y.vec[i]] = inid
    #     }
    #
    #
    #     # remove the feature segment from the original input matrix
    #     inft[x.vec[i], y.vec[i]] = 0
    #
    #   }
    #
    #   return(list(ink, it, inft))
    # }
    else if((sum(sum.Adj) == 0 && sum(sum.Adj2) == 0)){

      # fill x and y vectors with the i and j indices of the feature segment endpoint location
      x.vec[k] = ii
      y.vec[k] = jj
      # if(length(x.vec) > 50){
      #   x.vec <- x.vec[1:50]
      #   y.vec <- y.vec[1:50]
      #   k <- 50
      # }

      for(i in 1:length(x.vec)){

        # fill length and id matrices for the feature segment
        ink[x.vec[i], y.vec[i]] = k

        if(inID){
          it[x.vec[i], y.vec[i]] = inid
        }


        # remove the feature segment from the original input matrix
        inft[x.vec[i], y.vec[i]] = 0

      }

      return(list(ink, it, inft))
    }

  }

}
# Extract feature metrics
feat.Extract <- function(el, ft, in.IT = F, inXY = "X"){
  ft[is.na(ft)] <- 0
  ft[1,] <- 0
  ft[nrow(ft),] <- 0
  ft[,1] <- 0
  ft[,ncol(ft)] <- 0

  ft[1,] <- NA
  ft[nrow(ft),] <- NA
  ft[,1] <- NA
  ft[,ncol(ft)] <- NA

  # For MultXX dataset------------------------------------------------------
  # From the top left cell go through the data by row (L->R and T->B)
  # to find segment start locations
  # ------------------------------------------------------------------------
  if(inXY == "X"){
    iter <- 0 # iterator for feature ID
    k <- 0 # iterator for feature length
    tmp.ft <- ft # temporary feature matrix for removal of segments
    cnt <- ft * 0
    id <- ft * 0

    # for(i in 3:(nrow(ft) - 2)){
    for(i in 4:(nrow(ft) - 3)){
      cat("\r", i)
      # # Track progress
      # progress(i, (nrow(ft) - 1))#, progress.bar = TRUE)
      # Sys.sleep(0.01)
      # if (i == (nrow(ft) - 1)) message("Done!")

      # for(j in 3:(ncol(ft) - 2)){
      for(j in 4:(ncol(ft) - 3)){
        # If the cell  is a start point
        if(tmp.ft[i,j] == 1){

          iter = iter + 1 # iterate the id by one
          adj = T         # declare there is adjacency

          # follow the path of the line segment and extract/assign length and id, also remove features in original matrix
          feat.list <- path.Follow(ii = i,
                                   jj = j,
                                   inel = el,
                                   inft = tmp.ft,
                                   ink = cnt,
                                   it = id,
                                   inid = iter,
                                   inadj = adj,
                                   inID = in.IT,
                                   xy = inXY)



          cnt = feat.list[[1]]
          id = feat.list[[2]]
          tmp.ft <- feat.list[[3]]

        }

      }
    }
  }

  # For MultYY dataset------------------------------------------------------
  # From the top left cell go through the data by col (T->B and L->R)
  # to find segment start locations
  # ------------------------------------------------------------------------
  if(inXY == "Y"){
    iter <- 0 # iterator for feature ID
    k <- 0 # iterator for feature length
    tmp.ft <- ft # temporary feature matrix for removal of segments
    cnt <- ft * 0
    id <- ft * 0

    # for(j in 3:(ncol(ft) - 2)){
    for(j in 4:(ncol(ft) - 3)){
      cat("\r", j)

      # # Track progress
      # progress(j, (ncol(ft) - 1))#, progress.bar = TRUE)
      # Sys.sleep(0.01)
      # if (j == (ncol(ft) - 1)) message("Done!")

      # for(i in 3:(nrow(ft) - 2)){
      for(i in 4:(nrow(ft) - 3)){
        # If the cell  is a start point

        if(tmp.ft[i,j] == 1){

          iter = iter + 1 # iterate the id by one
          adj = T         # declare there is adjacency

          # follow the path of the line segment and extract/assign length and id, also remove features in original matrix
          feat.list <- path.Follow(ii = i,
                                   jj = j,
                                   inel = el,
                                   inft = tmp.ft,
                                   ink = cnt,
                                   it = id,
                                   inid = iter,
                                   inadj = adj,
                                   inID = in.IT,
                                   xy = inXY)



          cnt = feat.list[[1]]
          id = feat.list[[2]]
          tmp.ft <- feat.list[[3]]

        }

      }
    }
  }

  return(list(cnt, id))

}
# Create model input data
ID.means <- function(inFeat, inR=T){
  R.IDs <- inFeat[[2]]
  # R.IDs <- R.features[[2]]
  R.IDs[R.IDs==0] <- NA
  R.ID.vec <- sort(na.omit(unique(c(R.IDs))))

  R.D <- NA
  R.RR <- NA
  R.UX <- NA
  R.UY <- NA
  R.UZ <- NA
  R.tpi <- NA
  R.tri <- NA

  # if(inR){
  #   pos.R <- posR
  # }else{
  #   pos.R <- posS
  # }
  df <- data.frame(ID=c(R.IDs),
                   D=c(norm10(D)),
                   RR=c(rr),
                   UX=c(UX),
                   UY=c(UY),
                   UZ=c(UZ),
                   tpi=c(norm10(tpi)),
                   tri=c(norm10(tri)))


  df <- df[complete.cases(df$ID),]

  for(i in 1:length(R.ID.vec)){
    cat("\r", i,"/", length(R.ID.vec))
    # i = 1
    # profvis::profvis({
    tmp <- df[df$ID==R.ID.vec[i],]
    R.D[i] <- mean(tmp$D, na.rm=T)
    R.RR[i] <- mean(tmp$RR, na.rm=T)
    R.UX[i] <- mean(tmp$UX, na.rm=T)
    R.UY[i] <- mean(tmp$UY, na.rm=T)
    R.UZ[i] <- mean(tmp$UZ, na.rm=T)
    R.tpi[i] <- mean(tmp$tpi, na.rm=T)
    R.tri[i] <- mean(tmp$tri, na.rm=T)
  }
  if(inR){
    df.Feat <- data.frame(ID=R.ID.vec,
                          R=1,
                          S=0,
                          D=R.D,
                          RR=R.RR,
                          UX=R.UX,
                          UY=R.UY,
                          UZ=R.UZ,
                          tpi=R.tpi,
                          tri=R.tri)
    # df.Feat <- data.frame(ID=R.ID.vec,
    #                       R=1,
    #                       S=0,
    #                       D=R.D,
    #                       RR=R.RR,
    #                       UX=R.UX,
    #                       UY=R.UY,
    #                       UZ=R.UZ,
    #                       pos=R.Pos)
  }else{
    df.Feat <- data.frame(ID=R.ID.vec,
                          R=0,
                          S=1,
                          D=R.D,
                          RR=R.RR,
                          UX=R.UX,
                          UY=R.UY,
                          UZ=R.UZ,
                          tpi=R.tpi,
                          tri=R.tri)
    # df.Feat <- data.frame(ID=R.ID.vec,
    #                       R=0,
    #                       S=1,
    #                       D=R.D,
    #                       RR=R.RR,
    #                       UX=R.UX,
    #                       UY=R.UY,
    #                       UZ=R.UZ,
    #                       pos=R.Pos)
  }

  return(df.Feat)
}
create.Model.Data <- function(){
  dirPre <- paste(dir,"\\",fn.Uni, "_Preprocess", sep="")
  setwd(dirPre)

  cat("\n", "Creating model input data...")
  cat("\n")
  # These are the data used within a trained model to get predictions
  R.Means <- ID.means(inFeat=R.features, inR=T)
  cat("\n")
  S.Means <- ID.means(inFeat=S.features, inR=F)
  cat("\n")

  # cat("\n", "Making Ridge input data rasters...")
  # cat("\n")
  # # Write the values to raster stacks for optional investigation
  # R.List <- means2ras(inFeat=R.features, inMeans=R.Means)
  # R.stack <- list()
  # for(i in 1:length(R.List)){
  #   out = raster(R.List[[i]],xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  #   projection(out) = crs(M)
  #   R.stack[[i]] <- out
  # }
  # R.stack <- stack(R.stack)
  #
  # writeRaster(R.stack,paste0(fn.Uni, "_R_Means_Rasters.tif"), overwrite=TRUE)
  #
  # cat("\n", "Making Swale input data rasters...")
  # cat("\n")
  # S.List <- means2ras(inFeat=S.features, inMeans=S.Means)
  # S.stack <- list()
  # for(i in 1:length(S.List)){
  #   out = raster(S.List[[i]],xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  #   projection(out) = crs(M)
  #   S.stack[[i]] <- out
  # }
  # S.stack <- stack(S.stack)
  #
  # writeRaster(S.stack,paste0(fn.Uni, "_S_Means_Rasters.tif"), overwrite=TRUE)

  write.csv(R.Means, paste0(fn.Uni, "_R.Means.csv"), row.names=F)
  write.csv(S.Means, paste0(fn.Uni, "_S.Means.csv"), row.names=F)
  R.Means <<- R.Means
  S.Means <<- S.Means

  setwd(dir)
}
norm10 <- function(x, y, v=T){
  if(v){
    return((x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x,na.rm=T)))
  }else {
    return((y - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x,na.rm=T)))
  }

}

# -------------------------------------------------------------------------
# 2A ----------------------------------------------------------------------
# Send data through the model
modelProcessing <- function(pVal=0, modNew=F, iteration){
  if(!modNew){
    dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Raw", sep="")
    setwd(dirPre)
  }else{
    dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Intermediate", sep="")
    setwd(dirPre)
  }

  # First push the data through the master model -------------------------------
  test.Data <- rbind(R.Means, S.Means)
  test.Data <- test.Data[,4:ncol(test.Data)]
  test.Data <- as.matrix(test.Data)

  # Scaling
  # m <- colMeans(test.Data)
  # s <- apply(test.Data, 2, sd)
  # test.Data <- scale(test.Data, center = m, scale = s)
  cat("\n", "Predicting and getting classes...")
  predictions <- model %>% predict(test.Data, verbose=0)
  class_pred <- model %>% predict(test.Data, verbose=0) %>% k_argmax()
  class_pred <- k_get_value(class_pred) + 1
  # ----------------------------------------------------------------------------
  # Next write to raster -------------------------------------------------------
  test.Data <- rbind(R.Means, S.Means)
  outData <- as.data.frame(cbind(test.Data[,1:3], predictions, class_pred))
  names(outData) <- c("ID", "R", "S", "p1", "p2", "p3", "class")

  cat("\n", "Making Ridge Results...")
  cat("\n")
  R.Ras.out <- means2ras.Fin(inFeat=R.features, inMeans = outData, inR=T)
  R.stack <- list()
  for(i in 1:length(R.Ras.out)){
    out = raster(R.Ras.out[[i]],xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
    projection(out) = crs(M)
    R.stack[[i]] <- out
  }
  R.stack <- stack(R.stack)

  writeRaster(R.stack,paste0(fn.Uni, "_R_Results.tif"), overwrite=TRUE)

  cat("\n", "Making Swale Results...")
  cat("\n")
  S.Ras.out <- means2ras.Fin(inFeat=S.features, inMeans = outData, inR=F)

  S.stack <- list()
  for(i in 1:length(S.Ras.out)){
    out = raster(S.Ras.out[[i]],xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
    projection(out) = crs(M)
    S.stack[[i]] <- out
  }
  S.stack <- stack(S.stack)

  # S.Ras.out <<- S.Ras.out
  # R.Ras.out <<- R.Ras.out

  writeRaster(S.stack,paste0(fn.Uni, "_S_Results.tif"), overwrite=TRUE)
  # ----------------------------------------------------------------------------
  write.csv(outData, paste0(fn.Uni, "_outData.csv"), row.names=F)
  # Make Raw Results Raster ----------------------------------------------------
  cat("\n", "Making Raw Class Results...")
  cat("\n")
  p1.Ras <- S.Ras.out[[2]]
  p2.Ras <- R.Ras.out[[2]]
  p1.Ras <<- p1.Ras
  p2.Ras <<- p2.Ras
  tc.Results <- raw.Results(pVal=pVal, iteration = iteration)
  # ----------------------------------------------------------------------------

  tc.Results <<- tc.Results
  outData <<- outData
  S.stack <<- S.stack
  R.stack <<- R.stack
  setwd(dir)
}
# Function to put results into rasters
means2ras.Fin <- function(inFeat, inMeans, inR=T){
  # inFeat <- R.features; inMeans <- outData; inR<-T
  if(inR){
    inMeans <- inMeans[inMeans$R==1,]
  }else{
    inMeans <- inMeans[inMeans$R==0,]
  }

  tmp.C <- D*NA
  tmp.p <- tmp.C

  inMeans <- as.matrix(inMeans)
  # head(inMeans)


  tmpID <- inFeat[[2]]

  # plot(raster(tmp.C))
  for(i in 1:nrow(inMeans)){
    cat("\r",i,"/",nrow(inMeans))
    if(inR){
      tmp.ID <- inMeans[i,1]
      tmp.C[tmpID == tmp.ID] <- inMeans[i,7]
      tmp.p[tmpID == tmp.ID] <- inMeans[i,5]
    }else{
      tmp.ID <- inMeans[i,1]
      tmp.C[tmpID == tmp.ID] <- inMeans[i,7]
      tmp.p[tmpID == tmp.ID] <- inMeans[i,4]
    }

  }
  return(list(tmp.C, tmp.p))
}
# Makes raw results raster, pVal is a threshold of values to keep
raw.Results <- function(pVal=0, iteration = iteration){
  p1.Ras[p1.Ras<pVal] <- NA
  p2.Ras[p2.Ras<pVal] <- NA
  tc.Results <- D*NA
  if(sea.Uni=="E" || sea.Uni=="W"){
    for(i in 1:nrow(D)){
      cat("\r", i)
      # i = 1
      if(!is.infinite(suppressWarnings(max(p2.Ras[i,], na.rm = T)))){
        tmpR <- which(p2.Ras[i,] == max(p2.Ras[i,], na.rm = T))
        tc.Results[i,tmpR] <- 2
        # c.vec[i] <- tmpR
      }
      if(!is.infinite(suppressWarnings(max(p1.Ras[i,], na.rm = T)))){
        tmpS <- which(p1.Ras[i,] == max(p1.Ras[i,], na.rm = T))
        tc.Results[i,tmpS] <- 1
        # t.vec[i] <- tmpS
      }
    }
  }else{
    for(i in 1:ncol(D)){
      cat("\r", i)
      # i = 1
      if(!is.infinite(suppressWarnings(max(p2.Ras[,i], na.rm = T)))){
        tmpR <- which(p2.Ras[,i] == max(p2.Ras[,i], na.rm = T))
        tc.Results[tmpR,i] <- 2
        # c.vec[i] <- tmpR
      }
      if(!is.infinite(suppressWarnings(max(p1.Ras[,i], na.rm = T)))){
        tmpS <- which(p1.Ras[,i] == max(p1.Ras[,i], na.rm = T))
        tc.Results[tmpS,i] <- 1
        # t.vec[i] <- tmpS
      }
    }
  }
  tc.Chk <- tc.Results
  tc.Chk[!is.na(tc.Chk)] <- 1
  tc.Chk <- single.RM(Ridge=tc.Chk, Swale=tc.Chk*NA)
  tc.Chk <- tc.Chk[[1]]
  tc.Results <- tc.Results * tc.Chk

  out = raster(tc.Results,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)
  writeRaster(out,paste0(fn.Uni,"_tc_Results_",iteration,".tif"), overwrite=TRUE)

  return(tc.Results)
}
# -------------------------------------------------------------------------
# 2B ----------------------------------------------------------------------
# Creating Training Data
siteTrainingData <- function(tP1 = 0.95, tP3 = 0.01, cP1 = 0.95, cP3 = 0.01){
  dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Raw", sep="")
  setwd(dirPre)
  TC.Train <- D*NA

  TC.Train[p2.Ras > max(p2.Ras, na.rm=T)*cP1] <- 2
  TC.Train[p2.Ras < max(p2.Ras, na.rm=T)*cP3] <- 3

  TC.Train[p1.Ras > max(p1.Ras, na.rm=T)*tP1] <- 1
  TC.Train[p1.Ras < max(p1.Ras, na.rm=T)*tP3] <- 3

  out = raster(TC.Train,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)

  writeRaster(out, paste0(fn.Uni, "_tc_New_Training.tif"), overwrite=TRUE)
  trainClass <<- TC.Train
  setwd(dir)
}
new.Training <- function(inR, inS){
  trainNew <- trainClass

  for(i in 1:length(Swale1)){
    trainNew[inS==Swale1[i]] <- 1
  }
  for(i in 1:length(Swale3)){
    trainNew[inS==Swale3[i]] <- 3
  }
  for(i in 1:length(Ridge2)){
    trainNew[inR==Ridge2[i]] <- 2
  }
  for(i in 1:length(Ridge3)){
    trainNew[inR==Ridge3[i]] <- 3
  }

  return(trainNew)

}
siteModel <- function(trainClass=trainClass, LL=6, ridgePlot=T){
  dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Intermediate", sep="")
  setwd(dirPre)
  cat("\n", "Creating training data for the new model...")
  cat("\n")
  trainClass[trainClass>3] <- NA
  trainClass[trainClass==0]<- NA
  trainClass_bin <- trainClass
  trainClass_bin[!is.na(trainClass_bin)] <- 1


  # shoreMask <- shorelineMask(inshore=shore, dist=d2s, sea=sea.Uni)
  system.time(T.features <- feat.Extract(el=D, ft=trainClass_bin*shoreMask, in.IT=T, inXY=inXY2))

  t.L <- T.features[[1]]
  trainClass[t.L < LL & trainClass==3] <- NA

  train.Means <- train.ID.means(inFeat=T.features, inTrain=trainClass)
  o.tmp <- train.Means[train.Means$feat==3,]
  train.Data <- rbind(train.Means[train.Means$feat==1,],
                      train.Means[train.Means$feat==2,],
                      o.tmp)
  train.Data <- na.omit(train.Data)
  cat("\n")
  print(table(train.Data$feat))

  write.csv(train.Data, paste0(fn.Uni, "_Training.csv"), row.names=F)
  # # GGRIDGE Visualizations for density plotting
  if(ridgePlot){
    ridge.Plot.new(RS.train=train.Data, scale=F)
  }
  t.L <<- t.L
  # ----------------------------------------------------------------------------
  trainData <- train.Data[,3:ncol(train.Data)]
  trainData <- cbind(trainData, train.Data[,2])
  trainData <- as.matrix(trainData)

  training <-  trainData[,1:(ncol(trainData)-1)]
  trainingtarget <- trainData[, ncol(trainData)]-1
  trainingtarget <- to_categorical(trainingtarget, 3)
  # ----------------------------------------------------------------------------
  cat("\n", "Training new model...")
  # model <- load_model_tf("G:\\RR_Mins\\MorphometricExtraction\\Master\\saved_model/Morph_Model_Master_Lines")
  # model <- load_model_hdf5("D:\\Master\\Morph_Model_Master_Lines2.h5")
  # model <- load_model_hdf5("D:\\Master\\19\\B19_ModelResults_Intermediate\\saved_model\\B19_Morph_Model_Master_Lines")
  history <- model %>% fit(
    training, trainingtarget,
    batch_size = 128,
    epochs = 250,
    verbose = 0,
    validation_split=0,
    # validation_data=list(valid, validtarget),
    shuffle=T
  )
  training1 <-  trainData[,1:(ncol(trainData)-1)]
  trainingtarget1 <- trainData[, ncol(trainData)]-1

  predictions <- model %>% predict(training1, verbose=0)
  class_pred <- model %>% predict(training1, verbose=0) %>% k_argmax()
  class_pred <- k_get_value(class_pred) + 1

  df.pred <- data.frame(act=trainingtarget1+1,
                        pred=class_pred)
  cat("\n")
  print(table(df.pred))
  # ----------------------------------------------------------------------------
  # save_model_tf(model, paste0("saved_model/", fn.Uni, "_Morph_Model_Master_Lines"))
  save_model_hdf5(model, paste0(fn.Uni, "_Morph_Model_Master_Lines.h5"))

  model <<- model
  setwd(dir)
}
train.ID.means <- function(inFeat, inTrain){
  R.IDs <- inFeat[[2]]

  R.L <- inFeat[[1]]
  R.IDs[R.L < 4 & inTrain==3] <- NA
  # R.L[R.L <3] <- NA
  # R.L[!is.na(R.L)] <- 1
  # R.IDs <- R.IDs * R.L
  # R.IDs <- R.features[[2]]
  R.IDs[R.IDs==0] <- NA
  R.ID.vec <- sort(na.omit(unique(c(R.IDs))))

  R.D <- NA
  R.RR <- NA
  R.UX <- NA
  R.UY <- NA
  R.UZ <- NA
  R.tpi <- NA
  R.tri <- NA
  R.ft <- NA

  df <- data.frame(ID=c(R.IDs),
                   D=c(norm10(D)),
                   RR=c(rr),
                   UX=c(UX),
                   UY=c(UY),
                   UZ=c(UZ),
                   tpi=c(norm10(tpi)),
                   tri=c(norm10(tri)),
                   feat=c(inTrain))

  df <- df[complete.cases(df$ID),]

  for(i in 1:length(R.ID.vec)){
    cat("\r", i,"/", length(R.ID.vec))

    tmp <- df[df$ID==R.ID.vec[i],]
    R.D[i] <- mean(tmp$D, na.rm=T)
    R.RR[i] <- mean(tmp$RR, na.rm=T)
    R.UX[i] <- mean(tmp$UX, na.rm=T)
    R.UY[i] <- mean(tmp$UY, na.rm=T)
    R.UZ[i] <- mean(tmp$UZ, na.rm=T)
    R.tpi[i] <- mean(tmp$tpi, na.rm=T)
    R.tri[i] <- mean(tmp$tri, na.rm=T)
    R.ft[i] <- max(tmp$feat, na.rm=T)


  }
  df.Feat <- data.frame(ID=R.ID.vec,
                        feat=R.ft,
                        D=R.D,
                        RR=R.RR,
                        UX=R.UX,
                        UY=R.UY,
                        UZ=R.UZ,
                        tpi=R.tpi,
                        tri=R.tri)
  # df.Feat <- data.frame(ID=R.ID.vec,
  #                       feat=R.ft,
  #                       D=R.D,
  #                       RR=R.RR,
  #                       UX=R.UX,
  #                       UY=R.UY,
  #                       UZ=R.UZ,
  #                       pos=R.pos)
  return(df.Feat)
}
ridge.Plot.new <- function(RS.train=train.Data, scale=T){
  df.RX <- as.data.frame(RS.train[,3:(ncol(RS.train))])
  # for(i in 1:ncol(df.RX)){
  #   df.RX[,i] <- norm10(df.RX[,i])
  # }
  # dim(df.RX)
  if(scale){
    m <- colMeans(df.RX)
    s <- apply(df.RX, 2, sd)

    df.RX <- scale(df.RX, center = m, scale = s)
  }

  df.RX <- as.data.frame(df.RX)

  df.RX$classes <- RS.train[,2]
  # dev.new()
  df.RX %>%
    gather(variable, value, -classes) %>%
    ggplot(aes(y = as.factor(variable),
               fill = as.factor(classes),
               x = percent_rank(value))) + geom_density_ridges(aes(alpha=0.5), scale = 1.2)
}
# -------------------------------------------------------------------------
# 3 -----------------------------------------------------------------------
TC.Filled.fun <- function(RFill=T, TFill=T, R.RNG=5, S.RNG=0, doT=T, doC=T){
  D2 <- D
  D2[is.na(D)] <- 0
  dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Final", sep="")
  setwd(dirPre)
  S.id.out <- S.features[[2]]
  R.id.out <- R.features[[2]]

  # Fill gaps in the toe -------------------------------------------------------
  Toe <- tc.Results.Filtered
  Toe[Toe != 1] <- NA
  chk.Val <- mean(D2*Toe, na.rm=T)
  tmp.Toe <- 0
  if(doT){
    if(sea.Uni=="E" || sea.Uni=="W"){
      for(i in 1:nrow(D2)){
        # cat("\r","Toe: ", i)
        if(!is.infinite(suppressWarnings(max(Toe[i,], na.rm = T)))){
          tmp.Toe[i] <- i
        }else{
          tmp.Toe[i] <- 0
        }
      }
      toe.Fill <- Toe
      for(i in min(tmp.Toe[tmp.Toe>0]):max(tmp.Toe[tmp.Toe>0])){
        cat("\r","Toe: ", i)
        # i = 150
        if(tmp.Toe[i]==0){
          tmp.D <- c((which(toe.Fill[i-1,]==1)-1):(which(toe.Fill[i-1,]==1)+1))
          if(suppressWarnings(max(S[i,tmp.D], na.rm=T) == 1)){
            # dim(D)
            # tmp.Diff <- abs(D[i,tmp.D] - D[i-1, min(which(toe.Fill[i-1,]==1))])
            tmp.Diff <- abs(D2[i,tmp.D] - chk.Val)
            toe.Fill[i, tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))]] <- 1

          }else{
            if(TFill){
              tmp.R <- abs((c(1:ncol(D2))*S[i,]) - which(toe.Fill[i-1,]==1))
              # plot(tmp.R)
              if(min(tmp.R,na.rm=T) < S.RNG && !is.infinite(min(tmp.R,na.rm=T))){
                toe.Fill[i,which(tmp.R == min(tmp.R,na.rm=T))] <- 1
              }else{
                # tmp.Diff <- abs(D[i,tmp.D] - D[i-1,which(toe.Fill[i-1,]==1)])
                tmp.Diff <- abs(D2[i,tmp.D] - chk.Val)
                toe.Fill[i,tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))]] <- 1
              }

            }else{
              # tmp.Diff <- abs(D[i,tmp.D] - D[i-1,which(toe.Fill[i-1,]==1)])
              tmp.Diff <- abs(D2[i,tmp.D] - chk.Val)
              toe.Fill[i,tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))]] <- 1
            }
          }
        }
      }
    }else{
      for(i in 1:ncol(D)){
        # cat("\r","Toe: ", i)
        if(!is.infinite(suppressWarnings(max(Toe[,i], na.rm = T)))){
          tmp.Toe[i] <- i
        }else{
          tmp.Toe[i] <- 0
        }
      }
      # plot(tmp.Toe, pch=19, cex=0.5)
      toe.Fill <- Toe
      for(i in min(tmp.Toe[tmp.Toe>0]):max(tmp.Toe[tmp.Toe>0])){
        cat("\r","Toe: ", i)
        # i = 86
        # i = 150
        if(tmp.Toe[i]==0){
          # print(i)
          # plot(toe.Fill[, i-1])
          tmp.D <- c((which(toe.Fill[,i-1]==1)-1):(which(toe.Fill[,i-1]==1)+1))
          if(suppressWarnings(max(R[tmp.D,i], na.rm=T) == 1)){
            # dim(D)
            # plot(raster(D), col=grey(1:100/100))
            # tmp.Diff <- abs(D[tmp.D,i] - D[(which(toe.Fill[,i-1]==1)),i-1])
            tmp.Diff <- abs(D[tmp.D,i] - chk.Val)
            toe.Fill[tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))],i] <- 1
            # crest.Fill[tmp.D[which(R[tmp.D,i] == max(R[tmp.D,i], na.rm=T))],i] <- 2
          }else{
            if(RFill){
              tmp.R <- abs((c(1:nrow(D))*S[,i]) - which(toe.Fill[,i-1]==1))
              # plot(tmp.R)
              if(min(tmp.R,na.rm=T) < S.RNG && !is.infinite(min(tmp.R,na.rm=T))){
                toe.Fill[which(tmp.R == min(tmp.R,na.rm=T)),i] <- 1
              }else{
                # tmp.Diff <- abs(D[tmp.D,i] - D[which(toe.Fill[,i-1]==1),i-1])
                tmp.Diff <- abs(D[tmp.D,i] - chk.Val)
                toe.Fill[tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))],i] <- 1
              }

            }else{
              # tmp.Diff <- abs(D[tmp.D,i] - D[which(toe.Fill[,i-1]==1),i-1])
              tmp.Diff <- abs(D[tmp.D,i] - chk.Val)
              toe.Fill[tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))],i] <- 1
            }
          }
        }
      }
    }
  }else{
    toe.Fill <- Toe
  }

  # plot(tmp.Crest)
  out = raster(toe.Fill,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)
  writeRaster(out,paste0(fn.Uni, "_toe_Final.tif"), overwrite=TRUE)
  # ----------------------------------------------------------------------------
  # Fill gaps in the crest -----------------------------------------------------
  Crest <- tc.Results.Filtered
  Crest[Crest != 2] <- NA
  chk.Val <- mean(D*Crest, na.rm=T)
  tmp.Crest <- 0

  if(doC){
    if(sea.Uni=="E" || sea.Uni=="W"){
      for(i in 1:nrow(D)){
        # cat("\r","Crest: ", i)
        if(!is.infinite(suppressWarnings(max(Crest[i,], na.rm = T)))){
          tmp.Crest[i] <- i
        }else{
          tmp.Crest[i] <- 0
        }
      }
      # plot(tmp.Crest, pch=19, cex=0.5)
      crest.Fill <- Crest
      for(i in min(tmp.Crest[tmp.Crest>0]):max(tmp.Crest[tmp.Crest>0])){
        # i=2672
        cat("\r","Crest: ", i)
        # if(i+1 < ncol(D)){
        # i = 254
        # print(tmp.Crest[i])
        if(tmp.Crest[i]==0){
          tmp.D <- c((which(crest.Fill[i-1,]==2)-1):(which(crest.Fill[i-1,]==2)+1))
          if(suppressWarnings(max(R[i,tmp.D], na.rm=T) == 1)){
            # dim(D)
            tmp.Diff <- abs(D[i,tmp.D] - D[i-1, min(which(crest.Fill[i-1,]==2))])
            crest.Fill[i, tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))]] <- 2
            # crest.Fill[tmp.D[which(R[tmp.D,i] == max(R[tmp.D,i], na.rm=T))],i] <- 2
          }else{
            if(RFill){
              tmp.R <- abs((c(1:ncol(D))*R[i,]) - which(crest.Fill[i-1,]==2))
              # plot(tmp.R)
              if(min(tmp.R,na.rm=T) < R.RNG && !is.infinite(min(tmp.R,na.rm=T))){
                crest.Fill[i,which(tmp.R == min(tmp.R,na.rm=T))] <- 2
              }else{
                tmp.Diff <- abs(D2[i,tmp.D] - D2[i-1,which(crest.Fill[i-1,]==2)])
                crest.Fill[i,tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))]] <- 2
              }

            }else{
              tmp.Diff <- abs(D2[i,tmp.D] - D2[i-1,which(crest.Fill[i-1,]==2)])
              crest.Fill[i,tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))]] <- 2
            }
          }

        }
        # }

      }
    }else{
      for(i in 1:ncol(D)){
        # cat("\r","Crest: ", i)
        if(!is.infinite(suppressWarnings(max(Crest[,i], na.rm = T)))){
          tmp.Crest[i] <- i
        }else{
          tmp.Crest[i] <- 0
        }
      }
      # plot(tmp.Crest, pch=19, cex=0.5)
      crest.Fill <- Crest
      for(i in min(tmp.Crest[tmp.Crest>0]):max(tmp.Crest[tmp.Crest>0])){
        cat("\r","Crest: ", i)
        # i = 950
        # print(i)
        # RFill=T
        # dim(crest.Fill)
        # plot(crest.Fill[,i-1])
        # if(i+1 < ncol(D)){
        if(tmp.Crest[i]==0){
          # print(i)
          tmp.D <- c((which(crest.Fill[,i-1]==2)-1):(which(crest.Fill[,i-1]==2)+1))
          if(suppressWarnings(max(R[tmp.D,i], na.rm=T) == 1)){
            # dim(D)
            tmp.Diff <- abs(D[tmp.D,i] - D[(which(crest.Fill[,i-1]==2)),i-1])
            crest.Fill[tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))],i] <- 2
            # crest.Fill[tmp.D[which(R[tmp.D,i] == max(R[tmp.D,i], na.rm=T))],i] <- 2
          }else{
            if(RFill){
              tmp.R <- abs((c(1:nrow(D))*R[,i]) - which(crest.Fill[,i-1]==2))
              # plot(tmp.R)
              if(min(tmp.R,na.rm=T) < R.RNG && !is.infinite(min(tmp.R,na.rm=T))){
                crest.Fill[which(tmp.R == min(tmp.R,na.rm=T)),i] <- 2
              }else{
                tmp.Diff <- abs(D[tmp.D,i] - D[which(crest.Fill[,i-1]==2),i-1])
                crest.Fill[tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))],i] <- 2
              }

            }else{
              tmp.Diff <- abs(D[tmp.D,i] - D[which(crest.Fill[,i-1]==2),i-1])
              crest.Fill[tmp.D[which(tmp.Diff == min(tmp.Diff,na.rm=T))],i] <- 2
            }
          }
        }
        # }
      }
    }
  }else{
    crest.Fill <- Crest
  }


  out = raster(crest.Fill,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)
  writeRaster(out,paste0(fn.Uni, "_crest_Final.tif"), overwrite=TRUE)
  # ----------------------------------------------------------------------------
  # Combine Toe and Crest ------------------------------------------------------
  # mnB <- max(min(tmp.Crest[tmp.Crest>0]), min(tmp.Toe[tmp.Toe>0]))
  # mxB <-
  TC.Final <- D*NA
  TC.Final[!is.na(toe.Fill)] <- 1
  TC.Final[!is.na(crest.Fill)] <- 2

  out = raster(TC.Final,xmn=xmin(M),xmx=xmax(M),ymn=ymin(M),ymax(M))
  projection(out) = crs(M)
  writeRaster(out,paste0(fn.Uni, "_TC_Final.tif"), overwrite=TRUE)

  Crest <<- crest.Fill
  Toe <<- toe.Fill
  setwd(dir)
  return(TC.Final)

}
rand.plot.fin <- function(no.Smp = 100, minX=1, maxX=NA){
  dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Final", sep="")
  setwd(dirPre)
  pdf(paste0(fn.Uni,"_Profiles_Final.pdf"), width=8, height=4)

  minSmp <- minX
  maxSmp <- maxX

  # D.plot <- norm10(D)
  D.plot <- D
  # rr.plot <- norm10(rr)
  if(sea.Uni=="E"||sea.Uni=="W"){
    r.foc <- sample(minSmp:nrow(D.plot), no.Smp)
  }else{
    r.foc <- sample(minSmp:ncol(D.plot), no.Smp)
  }

  for(i in 1:length(r.foc)){
    # i = 1
    if(sea.Uni=="E"||sea.Uni=="W"){
      r.foc <- sample(minSmp:nrow(D.plot), no.Smp)
      if(!is.infinite(suppressWarnings(min(which(!is.na(D.plot[r.foc[i],])))))){
        if(is.na(maxSmp)){
          st <- c(min(which(!is.na(D.plot[r.foc[i],]))):max(which(!is.na(D.plot[r.foc[i],]))))
        }else{
          st <- c(maxSmp:max(which(!is.na(D.plot[r.foc[i],]))))
        }

        plot((D.plot[r.foc[i],st]), pch=19,cex=0.5, xlab='x', ylab='y',type='l', main=paste0("Profile: ", r.foc[i]))#, ylim=c(-0.1,1))
        # lines(D.plot[r.foc[i],500:800])
        points(D.plot[r.foc[i],st]*Toe[r.foc[i],st], pch=19,cex=1, col="red")
        points(D.plot[r.foc[i],st]*(Crest[r.foc[i],st]/2), pch=19,cex=1, col="blue")
        # lines(rr.plot[r.foc[i],st], col='forestgreen')
        # abline(h=5, lty="dashed", col="blue")
      }

    }else{

      if(!is.infinite(suppressWarnings(min(which(!is.na(D.plot[,r.foc[i]])))))){
        st <- c(min(which(!is.na(D.plot[,r.foc[i]]))):max(which(!is.na(D.plot[,r.foc[i]]))))
        plot((D.plot[st,r.foc[i]]), pch=19,cex=0.5, xlab='x', ylab='y',type='l', main=paste0("Profile: ", r.foc[i]))#, ylim=c(-0.1,1))
        # lines(D.plot[r.foc[i],500:800])
        points(D.plot[st,r.foc[i]]*Toe[st,r.foc[i]], pch=19,cex=1, col="red")
        points(D.plot[st,r.foc[i]]*(Crest[st,r.foc[i]]/2), pch=19,cex=1, col="blue")
        # lines(rr.plot[r.foc[i],st], col='forestgreen')
        # abline(h=5, lty="dashed", col="blue")
      }
    }
  }
  dev.off()
  setwd(dir)
}
Alongshore.Plot <- function(ww=15, hh=5){
  dirPre <- paste(dir,"\\",fn.Uni, "_ModelResults_Final", sep="")
  setwd(dirPre)
  feat.dat <- data.frame(x=c(XE),
                         y=c(YN),
                         z=c(D),
                         feat=c(TC.Filled))
  feat.dat <- na.omit(feat.dat)
  t.Dat <- feat.dat[feat.dat$feat==1,]
  c.Dat <- feat.dat[feat.dat$feat==2,]

  pdf(paste0(fn.Uni, "_Planar_AS_Lines.pdf"), width=ww, height=hh)
  # dev.new()
  layout(matrix(c(1,2), 1, 2, byrow = TRUE),
         widths=c(1,3), heights=c(1.5,1))

  if(sea.Uni=="E"||sea.Uni=="W"){
    t.Dat <- t.Dat[order(t.Dat$y),]
    c.Dat <- c.Dat[order(c.Dat$y),]

    plot(c.Dat$x, c.Dat$y, type='l',
         xlab="Easting (m)", ylab="Northing (m)", main="Planar")
    lines(t.Dat$x, t.Dat$y, col='blue')


    plot(c.Dat$y, c.Dat$z, type='l', ylim=c(0,max(c.Dat$z, na.rm=T)),
         xlab="Northing (m)", ylab="Elevation (m)", main="Alongshore")
    lines(t.Dat$y, t.Dat$z, col='blue')
  }else{
    t.Dat <- t.Dat[order(t.Dat$x),]
    c.Dat <- c.Dat[order(c.Dat$x),]

    plot(c.Dat$x, c.Dat$y, type='l',
         xlab="Easting (m)", ylab="Northing (m)", main="Planar")
    lines(t.Dat$x, t.Dat$y, col='blue')

    plot(c.Dat$x, c.Dat$z, type='l', ylim=c(0,max(c.Dat$z, na.rm=T)),
         xlab="Easting (m)", ylab="Elevation (m)", main="Alongshore")
    lines(t.Dat$x, t.Dat$z, col='blue')
  }
  dev.off()
  setwd(dir)
}
# -------------------------------------------------------------------------
