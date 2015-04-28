# Load and process data from all 4 model variants
# Need to re-work many of old functions so that they read and write to/from the correct locations.

library(raster)
library(rasterVis)
library(rgdal)
source("functions/consec_drywet.R")
source("functions/diurnal_cycle_v2.R")
source("functions/loadVeg.R") # Returns myveg = fractional veg cover for each pft tile
source("functions/loadOtherAncils.R")
source("functions/makeBoxes.R") 
source("functions/vegPrep.R") # Returns allveg = vegetation classess (1 to 6) AND veg classes intersected with zones (i.e. boxes)
source("functions/patches.R")
source("functions/movies.R")
# source("functions/tracker.R")
source("functions/tracker_v2.R")
source("functions/mcsStats.R")
source("functions/popStats.R")
source("functions/initiations.R")
source("functions/makeLines.R")
source("getMyData.R")
source("trackCheck.R")
source("functions/adjCoords.R")
source("functions/getLUT.R")
source("functions/trackInitiationsBackwards.R")
source("functions/loadAllAncils.R")
# source("mcsIntensity.R")
# hello

if (Sys.info()[["sysname"]] == "Darwin"){
	indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
	dlresultsdir <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/"
	resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
	scratchdir <- "/Users/ajh235/Work/Scratch/"
} else {
	indatadir <- "/data/local/hadhy/ModelData/WAFR/"
	dlresultsdir <- "/data/local/hadhy/Projects/InternalSabbatical/Results/"
	resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
	scratchdir <- "/data/local/hadhy/Scratch/"
	require(PP,lib.loc="/project/ukmo/rhel6/R")
}

rasterOptions(tmpdir=scratchdir, todisk=F)

timestep <- "5min" # "avg" "5min"
threshold <- 1000 # 1000 # 16 # Size threshold above which a patch of rain is classified as MCS
myproj <- "rp" # "ll" #

if (myproj == "rp"){
    land_simple <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer="land_rp") # Rotated Pole
} else {
    land_simple <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer="land_ll") # Lat Long
#     land_simple <- get(data(wrld_simpl)) # LatLong
#     land_country <- getData("countries") # LatLong
}

myLUT <- getLUT(nBndClass=3) # or 1

mydata <- getMyData(timestep=timestep, var="lsrain", overwrite=F)
if (timestep == "10min"){
    rb5216.1km.std <- adjCoords(mydata[[1]])
    rb5216.4km.std <- adjCoords(mydata[[2]])
    rb5216.4km.50k <- adjCoords(mydata[[3]])
    rb5216.4km.300k <- adjCoords(mydata[[4]])    
    models <- c("rb5216.4km.std", "rb5216.4km.50k", "rb5216.4km.300k") # [2:3]
    id <- "s" # c("s","w","y")#[2:3]
}
if (timestep == "5min"){
    rb5216.4km.std <- adjCoords(mydata)
    models <- "rb5216.4km.std"
    id <- "s"
}

# u = 1km std veg
# s = 4km std veg
# w = 4km smoothed 50km veg
# y = 4km smoothed 300km veg

# models <- c("avg5216.4km.std", "avg5216.1km.std", "avg5216.4km.300k", "avg5216.4km.50k")
# models <- c("rb5216.4km.std", "rb5216.1km.std", "rb5216.4km.300k", "rb5216.4km.50k")
# id <- c("s","u","y","w")

ancils <- loadAllAncils(myproj=myproj, nBndClass=1, model="rb5216.4km.std", overwrite=F)
mycl <- ancils$mycl
myveg <- ancils$myveg
blfrac <- myveg[[1]]; blfrac[is.na(blfrac)] <- 0
grassfrac <- myveg[[3]] + myveg[[4]]; grassfrac[is.na(grassfrac)] <- 0
spp  <- ancils[["spp"]]

# Recode "boundary tree" and "boundary grass" to "boundary"
mycl[mycl==5] <- 4
mycl[mycl==6] <- 4

for (x in 1:length(id)){
	print(models[x])
	# load data
	inbr <- get(models[x])
    
	allpatch <- mypatch(threshold=threshold, inbr=inbr, id=id[x], land_simple, timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, spp=spp, sppa=spp, overwrite=F)
    
# 	mcs <- allpatch[[1]]
# 	pop <- allpatch[[2]]
    mcs.infile <- paste(indatadir,"djzx",id[x],"/patches/",threshold,"km_",timestep,"/allmcs.",threshold,"km.vrt", sep="")
    pop.infile <- paste(indatadir,"djzx",id[x],"/patches/",threshold,"km_",timestep,"/allpop.",threshold,"km.vrt", sep="")
    if (myproj == "rp"){
        mcs <- adjCoords(brick(mcs.infile))
        pop <- adjCoords(brick(pop.infile))
    } 
    if (myproj == "ll"){
        # Following is not tested ...
        # Need to reproject rp to ll ...
        mcs <- system(paste('source $HOME/scitools/bin/activate ; python /Users/ajh235/Scripts/Python/getRotPoleDetails.py "',mcs.infile, '"', sep=""), intern=T)
        pop <- system(paste('source $HOME/scitools/bin/activate ; python /Users/ajh235/Scripts/Python/getRotPoleDetails.py "',pop.infile, '"', sep=""), intern=T)
    }
	
# 	pop.stats(pop=pop, id=id[x], threshold=threshold, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)
    
# 	mymcsstats <- mcs.stats(mcs=mcs, precip=inbr, id=id[x], threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)
	
#     initiations(inpop=pop, id=id[x], threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T)
    
	results <- tracker2(threshold=threshold, mcs, inbr, id=id[x], timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T)

    realInitiations <- trackInitiationsBackwards(results=results, ancils=ancils, threshold=threshold, mcs, pop, inbr, id=id[x], timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T)

    # Check MCS tracks match MCS blocks
#     trackCheck(results, mcs, inbr, id=id[x], threshold=threshold, timestep=timestep)
    
# 	mcsfollow <- mcs.stats(mcs=mcs, precip=inbr, id=id[x], threshold=1000, timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)
    
	# Create plots of diurnal cycle
# 	diurnalcycle2(inbr, type="all", patch=F, model.nm=models[x], id=id[x], spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=F) # Creates pdf plots
#     diurnalcycle2(inbr, type="intense", patch=F, model.nm=models[x], id=id[x], spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=F) # Creates pdf plots


# 	diurnalcycle(inbr, type="pop", patch=pop, model.nm=models[x], id=id[x], spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=T) # Creates pdf plots
	
    # Consecutive dry and wet days ...
# 	consec_drywet(inbr, id=id[x], indatadir=indatadir, scratchdir=scratchdir, overwrite=T)
    
    # Within MCS intensity
    # source("mcsIntensity.R")

    # Make stats on MCS Initiations 
#     source("initiation_analysis.R")

}

# Plot time series of avg.5216 for each model ...
# makeMovie(rb5216.4km.std, rb5216.4km.std, rb5216.4km.50k, rb5216.4km.300k, land_simple, sppa, dlresultsdir=dlresultsdir, models=models, boxes=F, type="precip", overwrite=T)
