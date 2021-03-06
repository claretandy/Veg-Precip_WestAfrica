library(raster)
library(rasterVis)
library(rgdal)
library(sm)
library(plotrix)
library(RColorBrewer)
library(classInt)


# Produce line plots of Precip vs time for different veg types 
# Time: Both 16/08 0Z to 19/08 6Z and averaged per hour
# Precip: Both total precip and precip split by MCS/POP # This has not been done yet!!!
diurnalcycle <- function(inbr, type="all", patch=F, model.nm="avg5216.4km.300k", id="s", spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=F){
	print("Plotting diurnal cycle ...")
	
	if (Sys.info()[["sysname"]] == "Darwin"){
		indatadir <- "/Volumes/MYBOOK/WAFR/"
		resultsdir <- "/Volumes/MYBOOK/Projects/InternalSabbatical/Results/"
		scratchdir <- "/Volumes/MYBOOK/Scratch/"
	} else {
		indatadir <- "/data/local/hadhy/ModelData/WAFR/"
		resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
		scratchdir <- "/data/local/hadhy/Scratch/"
		require(PP,lib.loc="/project/ukmo/rhel6/R")
	}
	setOptions(tmpdir=scratchdir, todisk=F)
	# Set indir and outpdf
# 	datadir <- paste("/data/local/hadhy/ModelData/WAFR/djzx",id,"/",sep="") 
	outpdf <- paste(resultsdir,"timeseries_byzone3-stdveg_",id,".pdf", sep="")
	
	# Load masks etc ...
# 	load("~/Projects/InternalSabbatical/Scripts/myWorkspace_29082012.Rdata")
# 	# Load land boundaries
# 	land_simple <- readOGR(dsn="/data/local/hadhy/WAFR3/ancils", layer="land_rp")
	# List box IDs
	myboxes <- 1:15
	
	# Open PDF 
	pdf(outpdf, width=12, height=6)
	for (bx in myboxes){

		print(bx)
#         browser()
		mymsk <- spp.r == bx
		mymsk[mymsk == 0] <- NA
# 		browser()
		e <- extent(sppa[bx,])
		# Plot close-up of grid box vegetation ...
        print("Plot close-up of grid box vegetation ...")
#         browser()
		plot(shift(mask(mycl, mymsk), x=-360), col=c("green", "yellow", "grey", "brown","red","blue"), main=paste("Vegetation boundaries masked by orography > 500m; Zone:",bx), legend=F, xlim=c(e@xmin,e@xmax), ylim=c(e@ymin,e@ymax))
		# text(centroids, labels=1:nrow(centroids), cex=3)
		plot(sppa, add=T)
		plot(land_simple, add=T)
		legend(x=e@xmax+1, y=e@ymin+((e@ymax - e@ymin)*0.75), c("tree","grass","mix","boundary","boundary tree","boundary grass")[c(1,5,3,4,6,2)], fill=c("green", "yellow", "grey", "brown","red","blue")[c(1,5,3,4,6,2)], xpd=T, cex=0.8)
		
		
		# 1) From 16/08 0Z to 19/08 6Z by total precip
        print("pr.bx ...")
		pr.bx <- crop(inbr, spp[bx,], filename=paste(scratchdir,"pr.bx.tif",sep=""), overwrite=T)
		mycl.bx <- crop(mycl, spp[bx,])
		mystats <- zonal(pr.bx, mycl.bx)
	# 	for (x in 1:nlayers(pr.ss.5min)){
	# 		print(x)
	# 		tmp <- zonal(pr.ss.5min[[x]], mask(mycl, mymsk))[,2]
	# 		if (x == 1){
	# 			mystats <- tmp
	# 		} else {
	# 			mystats <- cbind(mystats, tmp)
	# 		}
	# 	}
# 		if (bx==3 & model.nm=="avg5216.4km.300k"){ browser() }
		# Plot the full time series for this box ...
		mycols <- c("green", "yellow", "grey", "brown","red","blue")
		plot(c(1,(ncol(mystats)-1)), c(0,max(mystats[,-1])*60*60), type="n", xaxt="n", ylab="Precipitation (mm/hr)", xlab="Time", main=paste("Time series of mean precipitation by vegetation type in zone",bx)) #, ylim=c(min(mystats), max(mystats)))
		for (x in 1:nrow(mystats)){
			lines(1:(ncol(mystats)-1), mystats[x,-1]*60*60, col=mycols[x], type="l")
		}

		daylab <- which(format(getZ(inbr), "%H:%M") == "00:00")
		h6lab <- which(format(getZ(inbr), "%H:%M") == "06:00")
		h12lab <- which(format(getZ(inbr), "%H:%M") == "12:00")
		h18lab <- which(format(getZ(inbr), "%H:%M") == "18:00")
		axis(1, format(getZ(inbr)[daylab], "%b-%d"), at=daylab, tcl=-0.7)
		axis(1, format(getZ(inbr)[h6lab], "%H"), at=h6lab, tcl=-0.4, cex.axis=0.8)
		axis(1, format(getZ(inbr)[h12lab], "%H"), at=h12lab, tcl=-0.4, cex.axis=0.8)
		axis(1, format(getZ(inbr)[h18lab], "%H"), at=h18lab, tcl=-0.4, cex.axis=0.8)

		legend(x="topright", c("tree","grass","mix","boundary","boundary tree","boundary grass")[c(1,5,3,4,6,2)], fill=c("green", "yellow", "grey", "brown","red","blue")[c(1,5,3,4,6,2)], xpd=T, cex=0.8) # x=ncol(mystats)-100, y=max(mystats*60*60)

		# Plot the average for each hour ...
		hrs <- format(getZ(inbr), "%H")
		mystats.hrs <- matrix(NA, nrow=nrow(mystats), ncol=24)
		for (x in 1:nrow(mystats)){
			mystats.hrs[x,] <- tapply(mystats[x,-1], hrs, FUN=mean)
		}
		
		plot(c(0,23), c(0,max(mystats.hrs)*60*60), type="n", xaxt="n", ylab="Precipitation (mm/hr)", xlab="Time", main=paste("Mean diurnal cycle of precipitation by vegetation type in zone",bx))
		for (x in 1:nrow(mystats)){
			lines(0:23, mystats.hrs[x,]*60*60, col=mycols[x], type="l")
		}
		axis(1, c(0,6,12,18,24), at=c(0,6,12,18,24), tcl=-0.7, cex.axis=1.1)
		axis(1, c(1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23), at=c(1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23), tcl=-0.2, cex.axis=0.6, mgp=c(0,0.1,0))
		legend(x="topright", c("tree","grass","mix","boundary","boundary tree","boundary grass")[c(1,5,3,4,6,2)], fill=c("green", "yellow", "grey", "brown","red","blue")[c(1,5,3,4,6,2)], xpd=T, cex=0.8) # x=ncol(mystats)-100, y=max(mystats*60*60)
	}
	dev.off()

}