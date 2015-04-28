getMyData <- function(timestep="avg", var="lsrain", indatadir="/Users/ajh235/Work/DataLocal/ModelData/WAFR/", overwrite=T){
	
	if (timestep=="avg"){
				
		# avg.5216
		avg5216.1km.std <- brick(paste(indatadir,"djzxu/layers_nc/avg.5216.nc", sep=""))
		avg5216.4km.std <- brick(paste(indatadir,"djzxs/layers_nc/avg.5216.nc", sep=""))
		avg5216.4km.50k <- brick(paste(indatadir,"djzxw/layers_nc/avg.5216.nc", sep=""))
		avg5216.4km.300k <- brick(paste(indatadir,"djzxy/layers_nc/avg.5216.nc", sep=""))
		
		names(avg5216.1km.std) <- as.POSIXct("2006-08-16 00:00")+0:(nlayers(avg5216.1km.std)-1)*60*60
		names(avg5216.4km.std) <- as.POSIXct("2006-08-16 00:00")+0:(nlayers(avg5216.4km.std)-1)*60*60
		names(avg5216.4km.50k) <- as.POSIXct("2006-08-16 00:00")+0:(nlayers(avg5216.4km.50k)-1)*60*60
		names(avg5216.4km.300k) <- as.POSIXct("2006-08-16 00:00")+0:(nlayers(avg5216.4km.300k)-1)*60*60
		
		avg5216.1km.std <- setZ(avg5216.1km.std, as.POSIXct("2006-08-16 00:00")+0:(nlayers(avg5216.1km.std)-1)*60*60)
		avg5216.4km.std <- setZ(avg5216.4km.std, as.POSIXct("2006-08-16 00:00")+0:(nlayers(avg5216.4km.std)-1)*60*60)
		avg5216.4km.50k <- setZ(avg5216.4km.50k, as.POSIXct("2006-08-16 00:00")+0:(nlayers(avg5216.4km.50k)-1)*60*60)
		avg5216.4km.300k <- setZ(avg5216.4km.300k, as.POSIXct("2006-08-16 00:00")+0:(nlayers(avg5216.4km.300k)-1)*60*60)
		
		return(list(avg5216.1km.std, avg5216.4km.std, avg5216.4km.50k, avg5216.4km.300k))
	}
	if (timestep=="10min"){
		
		# 5216 10min
		m5216.1km.std <- brick(paste(indatadir,"djzxu/layers_nc/5216.nc", sep=""))
		m5216.4km.std <- brick(paste(indatadir,"djzxs/layers_nc/5216.nc", sep=""))
		m5216.4km.50k <- brick(paste(indatadir,"djzxw/layers_nc/5216.nc", sep=""))
		m5216.4km.300k <- brick(paste(indatadir,"djzxy/layers_nc/5216.nc", sep=""))
# 		browser()
		names(m5216.1km.std) <- as.POSIXct("2006-08-16 00:00")+0:(nlayers(m5216.1km.std)-1)*60*10
		names(m5216.4km.std) <- as.POSIXct("2006-08-16 00:00")+0:(nlayers(m5216.4km.std)-1)*60*10
		names(m5216.4km.50k) <- as.POSIXct("2006-08-16 00:00")+0:(nlayers(m5216.4km.50k)-1)*60*10
		names(m5216.4km.300k) <- as.POSIXct("2006-08-16 00:00")+0:(nlayers(m5216.4km.300k)-1)*60*10
		
		m5216.1km.std <- setZ(m5216.1km.std, as.POSIXct("2006-08-16 00:00")+0:(nlayers(m5216.1km.std)-1)*60*10)
		m5216.4km.std <- setZ(m5216.4km.std, as.POSIXct("2006-08-16 00:00")+0:(nlayers(m5216.4km.std)-1)*60*10)
		m5216.4km.50k <- setZ(m5216.4km.50k, as.POSIXct("2006-08-16 00:00")+0:(nlayers(m5216.4km.50k)-1)*60*10)
		m5216.4km.300k <- setZ(m5216.4km.300k, as.POSIXct("2006-08-16 00:00")+0:(nlayers(m5216.4km.300k)-1)*60*10)
	
		return(list(m5216.1km.std, m5216.4km.std, m5216.4km.50k, m5216.4km.300k))

	}

    if (timestep=="5min"){
        m5216.4km.std <- brick(paste(indatadir,"djzxs/5min/4203.nc",sep=""))
        indates <- as.POSIXct(getZ(m5216.4km.std))
        i <- format(indates, "%S") == "59"
        if(any(i)){ 
            indates[i] <- indates[i]+1
            m5216.4km.std <- setZ(m5216.4km.std, indates)
            print("Correcting dates ...")
        }
        names(m5216.4km.std) <- indates
        return(m5216.4km.std)
    }
	
}