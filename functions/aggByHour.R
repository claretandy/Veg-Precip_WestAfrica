aggByHour <- function(smalldata, xx, infun){
    try(
        treepr <- aggregate(Precip ~ time, data=smalldata[smalldata$ID == xx & smalldata$Class==1,], FUN=infun) # smalldata$Hour > 14 & smalldata$Hour < 21 & 
    , silent=T)
    try(
        grasspr <- aggregate(Precip ~ time, data=smalldata[smalldata$ID == xx & smalldata$Class==2,], FUN=infun) # smalldata$Hour > 14 & smalldata$Hour < 21 & 
        , silent=T)
    try(
        bndrypr <- aggregate(Precip ~ time, data=smalldata[smalldata$ID == xx & smalldata$Class==4,], FUN=infun) # smalldata$Hour > 14 & smalldata$Hour < 21 & 
        , silent=T)
    
    tmp <- merge(treepr, grasspr, by="time", all=T, sort=T)
    tmp <- merge(tmp, bndrypr, by="time", all=T, sort=T)
    colnames(tmp) <- c("time", "PrecipTree", "PrecipGrass", "PrecipBndry")
    
    tmp$greatest <- apply(tmp[,2:4], 1, FUN=which.max)
    tmp$Hour <- as.numeric(format(tmp$time, "%H"))
    
    return(tmp)
}