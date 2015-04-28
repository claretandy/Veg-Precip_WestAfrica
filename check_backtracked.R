source("functions/makeLines.R")

# Check backtracked initiation points

# Select all POP initiations that are linked to MCS

outresults_file <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/backtracked_results_1000km_5min.Rdata"
load(outresults_file)

bt_results$Hour <- as.numeric(format(bt_results$timestep, "%H"))
bt_results$Day <- as.numeric(format(bt_results$timestep, "%d"))

#pts <- realInitiations[realInitiations$type == "POP" & realInitiations$class == "Generation",c("x","y")]

bt_results_rp <- SpatialPointsDataFrame(coords=bt_results[,c("x","y")], data=bt_results, proj4string=CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere"))
# Reproject all points
bt_results_ll <- spTransform(bt_results_rp, CRSobj=CRS("+init=epsg:4326"), use_ob_tran=T)
# Write new lat long points to data.frame
bt_results_ll@data[,c("x","y")] <- coordinates(bt_results_ll)
lines_ll <- makeLines(bt_results_ll@data)[[1]]
# Write out points
writeOGR(bt_results_ll, dsn='/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/LatLongMapData/', layer='bt_results_ll', driver='ESRI Shapefile', morphToESRI=T, check_exists=T, overwrite_layer=T, delete_dsn=T)
# writeOGR(bt_results_rp, dsn='/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/LatLongMapData/', layer='bt_results_rp', driver='ESRI Shapefile', morphToESRI=T, check_exists=T, overwrite_layer=T, delete_dsn=T)
# Write out lines
writeOGR(lines_ll, dsn='/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/LatLongMapData/', layer='bt_results_lines_ll', driver='ESRI Shapefile', check_exists=T, overwrite_layer=T, delete_dsn=T, morphToESRI=T)



# Get ancils
ancils <- loadAllAncils(myproj="ll", nBndClass=1, model="rb5216.4km.std", overwrite=F)
land_rp <- ancils[[9]]
# Get all MCS IDs
mcsids <- unique(bt_results$ID[bt_results$type == "MCS"])
# Get all generations
mcs_start <- bt_results[grep("MCS start",bt_results$class),]
pop_start <- bt_results[bt_results$class == "Generation",]
pop_start_pm <- bt_results_ll[bt_results_ll$class == "Generation" & bt_results_ll$Hour >= 16 & bt_results_ll$Hour <= 18,]

# Plot POP Generations ove land cover ...
myLUT <- data.frame(ID=c(1,2,3,4,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "orography"), levels=c("tree", "grass", "sparse", "boundary", "orography")[c(3,2,4,1,5)]), Colours=c("dark green", "yellow", "orange", "sienna", "dark grey"), plotOrder=c(4,2,1,3,5))

mycl.1b <- ancils$mycl
mycl.1b[mycl.1b == 5] <- 4
mycl.1b[mycl.1b == 6] <- 4
mycl.1bf <- as.factor(mycl.1b)
ftab <- myLUT[myLUT$ID %in% levels(ancils$mycl.f)[[1]]$ID, ]
levels(mycl.1bf) <- ftab

print(
    levelplot(mycl.1bf, maxpixels=600000, par.settings=rasterTheme(region=myLUT$Colours), xlab=NULL, ylab=NULL, xlim=c(-12,10), ylim=c(4,18), main="Vegetation classes and afternoon MCS initiation points") + # , scales=list(draw=FALSE),  xlim=c(-24,15), ylim=c(-1,26), 
        latticeExtra::layer(sp.polygons(ancils$land_simple, col="black", lty=1)) + 
        latticeExtra::layer(sp.points(pop_start_pm, pch=3, cex=0.8, col="black"))
#         latticeExtra::layer(sp.polygons(spp)) +
#         latticeExtra::layer(sp.text(loc=coordinates(spp), txt=1:nrow(spp@data), cex=3)) +
#         latticeExtra::layer(sp.polygons(land12k.pol)) +
#         latticeExtra::layer(sp.polygons(land4k.pol))
)


# Plot examples of POP generations
t_start <- as.POSIXct("2006-08-17 17:00:00")
for (dn in 600*(0:-6)){
    thisdate <- t_start + dn
    outfile <- paste("/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_10min/", "pop_backtracked_",format(thisdate, "%d.%H%M"), ".tif", sep="")
    r <- raster(outfile)
    rss <- crop(r, extent(11,16,3.5,8.5))
    rss[rss == 0] <- NA
    # What are the MCS initiations in this timestep?
    if (t_start == thisdate){
        tracking_ids <- mcs_start[mcs_start$timestep == t_start,"ID"]
        pts_now <- mcs_start[mcs_start$timestep == thisdate,c("x","y")] 
        init_pts <- pts_now
    } else {
        pts_now <- bt_results[bt_results$timestep == thisdate,c("x","y")]
    }
    print(levelplot(rss, col.regions=rev(terrain.colors(20)), margin=F, main=list(thisdate)) + layer(sp.points(init_pts, col="black", pch=17)) + layer(sp.points(pts_now, col="black")))
}


# follow each MCS backwards, and ID the earliest initiation point
for (x in rev(mcs_start$ID)[1:5]){
    repeat{
        print(x)
        # Statements
        myrelated <- bt_results[bt_results$ID %in% x,"related"]
        myrelated <- myrelated[!is.na(myrelated)]
        myrelated <- as.numeric(unlist(strsplit(paste(myrelated, collapse=","), ",")))
        print(paste(myrelated, collapse=","))
        class_list <- bt_results[bt_results$ID %in% myrelated,"class"]
        if (any(class_list %in% "Generation")){
            print(bt_results[which(class_list %in% "Generation"),])
            browser()
            break
        }
    }
}



bt_results[bt_results$ID %in% c(1059,1239,1325,1241),]

# Reproject points
# coordinates(pop_start) <- ~ x+y
allgen_rp <- SpatialPointsDataFrame(coords=pop_start_pm[,c("x","y")], data=pop_start_pm, proj4string=CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere"))
# Reproject afternoon initiation points
initpts_ll <- spTransform(allgen_rp, CRSobj=CRS("+init=epsg:4326"), use_ob_tran=T)
initpts_ll <- initpts_ll[!is.na(extract(mycl.f, initpts_ll)),]
initpts_ll_lc <- extract(mycl.1bf, initpts_ll, method='simple', sp=T)

initpts_ll_lc$Landcover <- myLUT[match(initpts_ll_lc$mycl_masked_ll_0.3, myLUT$ID), "Landcover"]
initpts_ll_lc$Colours <- myLUT[match(initpts_ll_lc$mycl_masked_ll_0.3, myLUT$ID), "Colours"]
initpts_ll_lc$plotOrder <- myLUT[match(initpts_ll_lc$mycl_masked_ll_0.3, myLUT$ID), "plotOrder"]

initplot <- melt(initpts_ll_lc, id.vars=c("ID","Landcover","Colours","plotOrder"))
ggplot(initpts_ll_lc@data, aes(x="Landcover")) + geom_bar(stat="identity")


ptdens <- rasterize(initpts_ll, mycl.1bf, fun=function(x, ...)length(x))
ptdens_agg <- aggregate(ptdens, fact=6, fun=sum, na.rm=T)
