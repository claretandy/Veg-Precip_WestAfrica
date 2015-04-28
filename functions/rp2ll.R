rp2ll <- function(inraster, ancils, outfile, overwrite=F){
    require(rgdal)
    require(raster)
    require(sp)
    if (!file.exists(outfile) | overwrite){
        rp_proj <- CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere")
        ll_proj <- CRS("+init=epsg:4326")
        inraster.spdf <- SpatialPointsDataFrame(coordinates(inraster), data=data.frame(values(inraster)), proj4string=CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere"))
        inraster.spdf_ll <- spTransform(inraster.spdf, CRSobj=CRS("+init=epsg:4326"), use_ob_tran=T)
        inraster.ll_gaps <- rasterize(inraster.spdf_ll, y = ancils$mycl, field=colnames(inraster.spdf_ll@data)[1],fun=mean)
        inraster.ll <- resample(inraster.ll_gaps, ancils$mycl, method='bilinear', filename=outfile, overwrite=T)        
    } else {
        inraster.ll <- raster(outfile)
    }
    return(inraster.ll)
}


# Shell method...
# pr.sum.mcs.f_ll <- paste(indatadir, "/djzxs/derived/pr.sum.mcs.",timestep,".",threshold,"km.",hr,"_",plen,"_ll3.tif",sep="")
# pr.sum.pop.f_ll <- paste(indatadir, "/djzxs/derived/pr.sum.pop.",timestep,".",threshold,"km.",hr,"_",plen,"_ll.tif",sep="")
# system('echo GEOGCS[\"GCS_WGS_1984\", DATUM[\"WGS_1984\", SPHEROID[\"WGS_1984\",6378137.0,298.257223563]], PRIMEM[\"Greenwich\",0.0], EXTENSION[\"PROJ4\",\"+proj=ob_tran +o_proj=latlon +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +a=1 +to_meter=1\"], UNIT[\"Degree\",0.0174532925199433]] > in_proj_file.txt')
# system(paste('gdalwarp -s_srs in_proj_file.txt -t_srs "epsg:4326" -tr 0.036 0.036 -r bilinear ',pr.sum.mcs.f,' ',pr.sum.mcs.f_ll, sep=''))
# pr.sum.mcs.ll <- raster(pr.sum.mcs.f_ll)

