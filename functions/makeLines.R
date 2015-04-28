makeLines <- function(mcstracks){

    require(sp)
    
    jjj <- unique(mcstracks$ID)
    li <- 1
    mylines <- list()
#     if(exists("mytext")){rm(mytext)}
    for (jj in jjj){
        if(!exists("this.mytext")){
            this.mytext <- data.frame(mcstracks[mcstracks$ID == jj,c("x","y","ID")][1,])
        } else {
            this.mytext <- rbind(this.mytext, mcstracks[mcstracks$ID == jj,c("x","y","ID")][1,])
        }
        myline <- Line(mcstracks[mcstracks$ID == jj,c("x","y")])
        mylines[[li]] <- Lines(list(myline), ID=li) # or ID=li
        #         assign(x=paste("myline",li,sep=""), value=mylines)
        li <- li+1
    }
    myspl <- SpatialLines(mylines, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    mylines <- SpatialLinesDataFrame(sl=myspl, data=data.frame(lineid=1:length(jjj), MCSid=jjj))
    
    return(list(mylines, this.mytext))
}

makeLines_v2 <- function(all_overlaps){
    
    require(sp)
    require(raster)
    require(rgdal)
    
    li <- 1
    lineids <- integer()
    mylines <- list()
    for (id in unique(all_overlaps$ID)){
        ss <- all_overlaps[unlist(lapply(as.character(all_overlaps$tminus1_overlaps), FUN = function(x){
            xvec <- as.numeric(strsplit(x, split=",")[[1]])
            i <- which(xvec %in% id)
            if(length(i)==0){
                return(FALSE)
            } else {
                return(TRUE)
            }})),]
        if (nrow(ss) > 1){
            ss <- checkCoords(ss)
            # If checkCoords returns more than one ID, then create new line for each ID
            for (ssid in unique(ss$ID)){
                print(ssid)
                myline <- Line(ss[ss$ID == ssid,c("x","y")])
                mylines[[li]] <- Lines(list(myline), ID=li) # or ID=li
                lineids <- c(lineids, ssid)
                li <- li+1                                
            } 
                
        }
    }
    
    myspl <- SpatialLines(mylines, proj4string=CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere"))
    mylines_rp <- SpatialLinesDataFrame(sl=myspl, data=data.frame(lineid=1:(li-1), ClusterID=lineids))
    mylines_ll <- spTransform(mylines_rp, CRSobj=CRS("+init=epsg:4326"), use_ob_tran=T)
    return(list("ll"=mylines_ll, "rp"=mylines_rp))
}

checkCoords <- function(ss){
    
    # Tidy up lines
    
    require(sp)
    require(raster)
    
    # Check for duplicates in the timestep
    dupes  <- duplicated(ss$timestep)
    ss_nodupes <- ss[!dupes,]
    ss_nodupes <- ss_nodupes[order(ss_nodupes$timestep),]
    times      <- ss_nodupes$timestep
    times_m1   <- c(NA,times[-length(times)])
    timediff   <- as.integer(format(times - times_m1, "%M"))
    
    # Test if all time differences are 5 minutes
    if (any(timediff != 5, na.rm = T)){
        tdgt5 <- which(timediff > 5)
        start <- 1
        for(x in tdgt5){
            end <- x-1
            ss_nodupes[start:end,"ID"] <- as.integer(paste("888",ss_nodupes[start:end,"ID"],sep=""))
            start <- x
        }
    }

    # Test if the same Id is being used for far away clusters, so need to check the distance ...
    xy    <- ss_nodupes[,c("x","y")]
    xy_m1 <- rbind(c(NA,NA), xy[-nrow(xy),])
    xy$dist  <- pointDistance(xy, xy_m1, lonlat=TRUE) / 1000
    if (any(xy$dist > 500, na.rm=T)){
        # Do something here to separate lines
        dgt500 <- which(xy$dist > 500)
        start <- 1
        for (x in dgt500){
            if (ss_nodupes[x-1,"ID"] == ss_nodupes[x,"ID"]){
                end <- x-1
                ss_nodupes[start:end,"ID"] <- as.integer(paste("999",ss_nodupes[start:end,"ID"],sep=""))
                start <- x
            }
        }
    }
    
    # Check the timesteps are consecutive ...
    if (!all(sort(ss_nodupes$timestep) == ss_nodupes$timestep)){
        browser()
    }
    
    return(ss_nodupes)
}

