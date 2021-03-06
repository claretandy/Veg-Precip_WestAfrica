# Calculate consecutive dry or wet minutes ...
# input is time series raster brick of hourly accumulated precipitation

library(raster)
library(rasterVis)
#require(PP,lib.loc="/project/ukmo/rhel6/R")
library(rgdal)
library(sm)
library(plotrix)
library(RColorBrewer)
library(classInt)
# source("~/Rscripts/defaultFunctions.R")

consec_drywet <- function(inbr, id="s", indatadir="/data/local/hadhy/ModelData/WAFR/", scratchdir="/data/local/hadhy/Scratch/", overwrite=T){

	setOptions(tmpdir=scratchdir, todisk=T)

	# Set indir and outdir
	datadir <- paste(indatadir,"djzx",id,"/",sep="") 
	outdir <- paste(indatadir,"djzx",id,"/derived/",sep="")
	if (!file.exists(outdir)){ dir.create(outdir) }
	
	# Load some saved objects ...
	#load("~/Projects/InternalSabbatical/Scripts/myWorkspace_06082012.Rdata")
	
	
	# Loop through all layers
	for (f in 1:(nlayers(inbr))){
		# Set current date and time
		mydate <- as.POSIXct(gsub("X","",layerNames(inbr)[f]), format="%Y.%m.%d.%H.%M.%S")
		print(layerNames(inbr)[f])
		
		# Set output filenames
		out.time <- as.character(mydate, format="%d%H.%M")
		acc.p.f <- paste(outdir, "acc_precip.", out.time, ".tif", sep="")
		acc.dry.f <- paste(outdir, "acc_dry.", out.time, ".tif", sep="")
		acc.wet.f <- paste(outdir, "acc_wet.", out.time, ".tif", sep="")
		consec.dry.f <- paste(outdir, "consec_dry.", out.time, ".tif", sep="")
		consec.wet.f <- paste(outdir, "consec_wet.", out.time, ".tif", sep="")
		
		# Load data for current time step
		r <- inbr[[f]]
		
		# Convert precip from kg/m2/s to total mm in last 60 minutes
		r <- calc(r, fun=function(x){x*60*60}, filename=paste(scratchdir,"myr.tif",sep=""), overwrite=T, format="GTiff")
		
		# Do the business ...
		if (f == 1){
			# Create new raster layers for the results
			acc.p      <- r ; layerNames(acc.p) <- "Accumulated Precipitation"
			acc.dry <- setValues(acc.p, 0); layerNames(acc.dry) <- "Accumulated dry minutes"
			acc.wet <- setValues(acc.p, 0); layerNames(acc.wet) <- "Accumulated wet minutes"
			consec.dry <- setValues(acc.p, 0); layerNames(consec.dry) <- "Consecutive dry minutes"
			consec.wet <- setValues(acc.p, 0); layerNames(consec.wet) <- "Consecutive wet minutes"
		} else {
			# NB: No rain ~ 1mm per day (i.e. 1/24 )
			
			# Accumulated rain (corrected above to include accumulated rain over last 5 minutes)
			print("Accumulated rain ...")
			if (!file.exists(acc.p.f) | overwrite==T){
				acc.p <- overlay(acc.p, r, fun=function(x,y){return(x+z)}, na.rm=T, filename=acc.p.f, overwrite=T, format="GTiff")
			}
			# NB: For all of the acc.* and consec.* below, need to multiply result by 5 to get total minutes dry
			# Correct by /60/24 to get fractions of dry days
			
			# Accumulated dry minutes
			#if (!file.exists(acc.dry.f)){
			print("Accumulated dry ...")
			acc.dry <- overlay(acc.dry, r, fun=function(x,y){
					z <- y<1/24 ; # Returns 1 if precip < 1/24mm per hour
					z <- z*60 ; # Number of minutes in time step
					return(x+z)}, 
					na.rm=T, filename=acc.dry.f, overwrite=T, format="GTiff") #}
			# Consecutive dry minutes
			#if (!file.exists(consec.dry.f)){
			print("Consecutive dry ...")
			consec.dry <- overlay(consec.dry, r, fun=function(x,y){
					z <- y<1/24 ; # Returns 1 if precip < 1/24mm per hour
					z[z==0] <- NA ;
					z <- z*60 ; # Number of minutes in time step
					newx <- x+z ;
					newx[is.na(newx)] <- 0 ; # Grid cells currently wet revert to 0
					return(newx)}, 
					na.rm=T, filename=consec.dry.f, overwrite=T, format="GTiff") #}
			
			# Accumulated wet minutes
			#if (!file.exists(acc.wet.f)){
			print("Accumulated wet ...")
			acc.wet <- overlay(acc.wet, r, fun=function(x,y){
					z <- y>=1/24 ; # Returns 1 if precip >= 1/24mm per hour
					z <- z*60 ; # Number of minutes in time step
					return(x+z)}, 
					na.rm=T, filename=acc.wet.f, overwrite=T, format="GTiff") #}
			# Consecutive wet minutes
			#if (!file.exists(consec.wet.f)){
			print("Consecutive wet ...")
			consec.wet <- overlay(consec.wet, r, fun=function(x,y){
					z <- y>=1/24 ; # Returns 1 if precip >= 1/24mm per hour
					z[z==0] <- NA ;
					z <- z*60 ; # Number of minutes in time step
					newx <- x+z ;
					newx[is.na(newx)] <- 0 ; # Grid cells currently dry revert to 0 
					return(newx)}, 
					na.rm=T, filename=consec.wet.f, overwrite=T, format="GTiff") #}
			# browser()
		}
	}
}





