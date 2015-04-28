# Overlay all of the POP and MCS clusters, and produce a list (for each cluster ID and timestep) of the overlapping IDs

library(raster)
library(rasterVis)
source("getMyData.R")
source("functions/adjCoords.R")
source("functions/getCorrectDates.R")
source("functions/mycrosstab.R")
source("functions/makeLines.R")

overwrite <- FALSE
timestep <- '5min'
mydata <- getMyData(timestep=timestep, var="lsrain", overwrite=F)
rb5216.4km.std <- adjCoords(mydata)
models <- "rb5216.4km.std"
id <- "s"

# Output results file
current <- "20112014" # format(Sys.time(), "%d%m%Y")
all_overlaps_rdata <- paste("/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/all_overlaps_",current,".RData", sep="")

# Get Dates
gcd <- getCorrectDates(rb5216.4km.std)
alldates <- gcd[[1]]
inbr <- gcd[[2]]

# ID, x, y, timestep, overlaps_in_t-1
all_overlaps <- data.frame(ID=integer(), x=numeric(), y=numeric(), timestep=as.POSIXct, count=integer(), tminus1_overlaps=character())

for (dt in 2:(length(alldates)-1)){
    
    lastdt <- format(alldates[dt-1], "%d.%H%M")
    thisdt <- format(alldates[dt], "%d.%H%M")
    
    print(thisdt)
    
    # Output files
    this_all.file <- paste('/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/AllClusters/allIDs_',thisdt,'.tif',sep="")
    last_all.file <- paste('/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/AllClusters/allIDs_',lastdt,'.tif',sep="")
    this_allpts.file <- paste('/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/AllClusters/allIDs_',thisdt,'_points.RData',sep="")
    
    # Load the data
    last_mcs <- raster(paste('/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_1000km_5min/mcs_tracking_1000km_',lastdt,'.tif', sep=''))
    last_pop <- raster(paste('/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_1000km_5min/pop_backtracked_',lastdt,'.tif', sep=''))
    
    this_mcs <- raster(paste('/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_1000km_5min/mcs_tracking_1000km_',thisdt,'.tif', sep=''))
    this_pop <- raster(paste('/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_1000km_5min/pop_backtracked_',thisdt,'.tif', sep=''))
    
    if (file.exists(last_all.file) & !overwrite){
        last_all <- raster(last_all.file)
    } else {
        last_all <- overlay(last_mcs, last_pop, fun=function(x,y){
            x[is.na(x)] <- 0
            i <- which(x == 0 & y > 0)
            x[i] <- y[i]
            x[x==0] <- NA
            return(x) }, filename=last_all.file, overwrite=T)
    }
    
    if (file.exists(this_all.file) & !overwrite){
        this_all <- raster(this_all.file)
    } else {
        this_all <- overlay(this_mcs, this_pop, fun=function(x,y){
            x[is.na(x)] <- 0
            i <- which(x == 0 & y > 0)
            x[i] <- y[i]
            x[x==0] <- NA
            return(x) }, filename=this_all.file, overwrite=T)        
    }
    
    # Ideal output:
    ## ID, x, y, timestep, overlaps_in_t-1
    # Cross-tabulate
    myxtab <- mycrosstab(last_all, this_all) # rows, columns
    # Get the frequency
    myfreq <- as.data.frame(freq(this_all, useNA="no"))
    # Rasterize clusters to get the centroid
    if (file.exists(this_allpts.file) & !overwrite){
        load(this_allpts.file)
        #mypts <- readOGR(dsn=this_allpts.file, layer=strsplit(basename(this_allpts.file), extension(this_allpts.file))[[1]])
    } else {
        file.remove(this_allpts.file)
        this_all.poly <- rasterToPolygons(this_all, dissolve=T)
        # Match the order of points in myxtab
        i <- match(as.integer(colnames(myxtab)), this_all.poly@data[,1])
        mypts <- coordinates(this_all.poly)[i,]
        # Just in case there's only 1 polygon ...
        if (inherits(mypts, "numeric")){
            mypts <- matrix(mypts, nrow = length(i), ncol = 2)
        }
        save(mypts, ascii=TRUE, file=this_allpts.file)
    }
    
    for (tclust in 1:dim(myxtab)[2]){
        # Get this ID
        id <- as.integer(colnames(myxtab)[tclust])
        # Get overlaps from myxtab
        ri <- which(myxtab[,tclust] > 0)
        tclust_overlaps <- as.numeric(rownames(myxtab))[ri]
        all_overlaps <- rbind(all_overlaps, data.frame(ID=id, x=mypts[tclust,1], y=mypts[tclust,2], timestep=alldates[dt], count=myfreq[tclust,"count"], tminus1_overlaps=paste(tclust_overlaps, collapse=",")))
    }
    #file.remove(all_overlaps_rdata)
    save(all_overlaps, file=all_overlaps_rdata)
    
}
save(all_overlaps, file=all_overlaps_rdata)

# Make lines from the dataframe
load(all_overlaps_rdata)
id_gt50 <- as.integer(attributes(which(idtab[order(idtab)] > 50))$names)
mylines <- makeLines_v2(subset(all_overlaps, ID %in% id_gt50))

levelplot(mycl.1bf, att="Landcover", maxpixels=600000, main=paste("MCS Tracks Over Vegetation Classes", sep=""), xlim=c(-12,10), ylim=c(4,18), xlab=NULL, ylab=NULL, col.regions=as.character(ftab$Colours), colorkey=list(space="right", labels=list(rot=0))) + 
    latticeExtra::layer(sp.polygons(land_simple, lty=2)) +
    latticeExtra::layer(sp.lines(mylines[["ll"]]))
