getCorrectDates <- function(inbr){
    require(raster)
    indates <- as.POSIXct(getZ(inbr))
    i <- format(indates, "%S") == "59"
    if(any(i)){ 
        indates[i] <- indates[i]+1
        inbr <- setZ(inbr, indates)
        names(inbr) <- indates
        print("Correcting dates ...")
    }
    return(list(indates,inbr))
}
