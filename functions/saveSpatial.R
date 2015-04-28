saveSpatial <- function(bt_results, filename, doLines=F){
    
    source("functions/makeLines.R")
    require(sp)
    require(rgdal)
    
    bt_results$Hour <- as.numeric(format(bt_results$timestep, "%H"))
    bt_results$Day <- as.numeric(format(bt_results$timestep, "%d"))
    #bt_results$class <- sub("minor parts from ID [0-9]{1,6}", "minor parts", bt_results$class)
    
    bt_results_rp <- SpatialPointsDataFrame(coords=bt_results[,c("x","y")], data=bt_results, proj4string=CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere"))
    # Reproject all points
    bt_results_ll <- spTransform(bt_results_rp, CRSobj=CRS("+init=epsg:4326"), use_ob_tran=T)
    # Write new lat long points to data.frame
    bt_results_ll@data[,c("x","y")] <- coordinates(bt_results_ll)
    
    # Write out points
    writeOGR(bt_results_ll, dsn='/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/LatLongMapData/', layer=paste(filename,'_ll',sep=""), driver='ESRI Shapefile', morphToESRI=T, check_exists=T, overwrite_layer=T, delete_dsn=T)
    proj4string(bt_results_rp) <- proj4string(bt_results_ll)
    writeOGR(bt_results_rp, dsn='/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/RPMapData/', layer=paste(filename,'_rp',sep=""), driver='ESRI Shapefile', morphToESRI=T, check_exists=T, overwrite_layer=T, delete_dsn=T)
    
    
    if (doLines){
        lines_ll <- makeLines(bt_results_ll@data)[[1]]
        # Write out lines
        writeOGR(lines_ll, dsn='/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/LatLongMapData/', layer=paste(filename,'_lines_ll',sep=""), driver='ESRI Shapefile', check_exists=T, overwrite_layer=T, delete_dsn=T, morphToESRI=T)
        
        return(list(bt_results_rp=bt_results_rp, bt_results_ll=bt_results_ll, lines_ll=lines_ll))
    } else {
        return(list(bt_results_rp=bt_results_rp, bt_results_ll=bt_results_ll))
    }
    
    
}

