# When an MCS occurs over forest, boundary and grass, where does it rain the most intensely?

source("functions/makeLines.R")
source("functions/multiplot.R")
source("functions/getLUT.R")
source("functions/loadAllAncils.R")
source("getMyData.R")
source("functions/tracker_v2.R")
source("functions/aggByHour.R")
source("functions/adjCoords.R")
library(ggplot2)
library(reshape)
library(Hmisc)
  
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

inthr <- 10
timestep <- "10min"
# pdf(paste("/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_LC_comparison_gt",inthr,"mm_v1.pdf",sep=""), width=14, height=9)
mcsrst.path <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_10min/"
gntheme <- rasterTheme(pch=19, region=brewer.pal(9, 'YlGn'))
myLUT <- getLUT(nBndClass=1)

# Load veg fractions
# return(list(myveg, mylandfrac, myorog, mycl, mycl.f, spp, spp.r, myboxes, land_simple, mycl.z))
ancils <- loadAllAncils(myproj="rp", nBndClass=1, model="rb5216.4km.std", overwrite=F)
mycl <- ancils[[4]]
myveg <- ancils[[1]]

blfrac <- myveg[[1]]; blfrac[is.na(blfrac)] <- 0
grassfrac <- myveg[[3]] + myveg[[4]]; grassfrac[is.na(grassfrac)] <- 0

# Recode "boundary tree" and "boundary grass" to "boundary"
mycl[mycl==5] <- 4
mycl[mycl==6] <- 4

# Get unique list of MCS ids
mcs.infile <- paste(indatadir,"djzxs/patches/1000km_10min/allmcs.1000km.vrt", sep="")
mcs <- adjCoords(brick(mcs.infile))
mydata <- getMyData(timestep=timestep, var="lsrain", overwrite=F)
rb5216.4km.std <- adjCoords(mydata[[2]])
results <- tracker2(threshold=1000, mcs, rb5216.4km.std, id="s", timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)

mcs.ids <- unique(results$ID)
alllines <- makeLines(results)[[1]]

results$TOD <- format(results$timestep, "%H:%M")
mcslist <- sort(table(results[,"ID"])) # All times of day # results$TOD > 15 & results$TOD <= 21
iii <- as.numeric(attributes(mcslist[mcslist > 6])$dimnames[[1]]) # These are all the MCS that occur # between 15:00 and 21:00

spresults <- makeLines(results[results$ID %in% iii,])
myspdf <- spresults[[1]] # Spatial Lines
mytext <- spresults[[2]] # Spatial points of starting point

# Get coordinates for starting points
coordinates(mytext) <- ~x+y
# Get the longest tracks for plotting starting points
iiilongest <- as.numeric(attributes(mcslist[mcslist > 50])$dimnames[[1]])

print(
    levelplot(blfrac, main="MCS tracks longer than 1 hour (16/08 to 18/08)", margin=F, par.settings=gntheme) + 
          latticeExtra::layer(sp.polygons(land_simple, col="grey")) + 
          latticeExtra::layer(sp.lines(myspdf)) + 
          latticeExtra::layer(sp.text(loc=coordinates(mytext[mytext$ID %in% iiilongest,]), txt=mytext[mytext$ID %in% iiilongest,]$ID)) 
    )

if (exists("bigdata")){ rm(bigdata) }

# Loop through each track
for (ii in iii){
    
    mydates <- results[results$ID == ii,"timestep"]
    if (exists("mydf")){ rm(mydf) }
    
    print(paste("MCS: ",ii, "; ",length(mydates)," timesteps",sep=" "))
    
    for (x in 1:length(mydates)){
#         print(mydates[x])
        
        mcs.now <- shift(raster(paste(mcsrst.path,"mcs_tracking_1000km_",format(mydates[x], "%d.%H%M"), ".tif",sep="")), x=-360)
        mcs.now <- adjCoords(mcs.now)
        
        mcsiinow <- mcs.now == ii
        mcsiinow[mcsiinow == 0] <- NA
        precip.now <- shift(rb5216.4km.std[[which(getZ(rb5216.4km.std) == mydates[x])]], x=-360)
        
        # Mask precip that's not in the MCS
        mcsprecip   <- raster::mask(precip.now, mcsiinow)*3600

        # Where the rainfall is most intense, what's the rate over grass and the rate over tree?
        # 1. Classify intense rain (>10mm)
        mcsprecip.intense <- calc(mcsprecip, fun=function(x){x[x<inthr]<-NA; return(x)})
        
        # 2. Mask veg classes to intense precip patches
        mycl.mask <- raster::mask(mycl, mcsprecip.intense)
        
        # Check that some cells have data in ...
        # If MCS occurs in Sea, this will be FALSE
        haveData <- length(which(is.na(getValues(mycl.mask)))) != ncell(mycl.mask)
        
        if (haveData){
            
            # Write all data to a big dataframe
            loc <- which(!is.na(getValues(mcsprecip.intense)))
            prvals <- getValues(mcsprecip.intense)[loc]
            clvals <- getValues(mycl.mask)[loc]
            if (!exists("bigdata")){
                bigdata <- data.frame("ID"=ii, "time"=mydates[x], "TOD"=format(mydates[x], "%H:%M"), "Precip"=prvals, "Class"=clvals)
            } else {
                bigdata <- rbind(bigdata, data.frame("ID"=ii, "time"=mydates[x], "TOD"=format(mydates[x], "%H:%M"), "Precip"=prvals, "Class"=clvals))
            }
            
        } else {
#             print("We don\'t have any data for this timestep")
        }
        
    } # End of mydates
        
} # End of ii


# For all MCS that occur for longer than 1 hour, how many intense rainfall points occur over each cover type?
bigdata$Hour <- as.numeric(format(bigdata$time, "%H"))
try(rm(out))
for (i in 0:23){
    print(i)
    myfreq <- table(bigdata[bigdata$Hour >= i & bigdata$Hour < i+1,"Class"])
    myfreq <- data.frame(myfreq)
    colnames(myfreq)[2] <- paste("Hr", i, sep="_")
    if (!exists("out")){
        out <- cbind(myLUT[which(myLUT$ID %in% myfreq$Var1),], myfreq)
        out$Var1 <- NULL
    } else {
        out$smthg <- NA
        out[which(out$ID %in% myfreq$Var1),"smthg"] <- myfreq[,2]
        z <- length(colnames(out))
        colnames(out)[z] <- paste("Hr", i, sep="_")
    }
    out$Var1 <- NULL
}
outmelt <- melt(out, id.vars=c("ID","Landcover","Colours","plotOrder"))
ggplot(outmelt, aes(x=variable, y=value, fill=Landcover)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(values=as.character(myLUT$Colours[order(myLUT$plotOrder)]),
                      name="Land cover", 
                      labels=capwords(as.character(myLUT$Landcover[order(myLUT$plotOrder)]))) +
    labs(title="Frequency of intense precipitation over land cover types", y="Frequency", x="Hour") +
    scale_x_discrete(labels=0:23) 

# barplot(as.matrix(out[,-c(1:4)]), names.arg=0:23, col=as.character(out$Colours))

# Find MCS that occur over tree/boundary/grass at the same time, and compare mean precip rate
try(rm(smalldata))
print("Smalldata ...")
try(load("smalldata_gt1hour_3cl_gt10pc.Rdata"))
if(!exists("smalldata")){
    
    for (x in unique(bigdata$ID)){
    print(x)
    alltime <- unique(bigdata[bigdata$ID == x, "time"])
    for (y in 1:length(alltime)){
        indata <- bigdata[bigdata$ID == x & bigdata$time == alltime[y], c("Class", "Precip")]
        # Test if >10% of each class
        infreq <- data.frame(table(indata$Class))
        infreq$Freq <- infreq$Freq/sum(infreq$Freq)
        this.lc <- myLUT[which(myLUT$ID %in% infreq$Var1), "Landcover"]
        len.lc <- length(which(c("tree", "grass", "boundary") %in% this.lc))
        if (len.lc >= 3){
            gt10pc <- length(which(infreq[infreq$Var1 %in% c(1,2,4),"Freq"] > 0.1))
            if (gt10pc >= 3){
                # Here we are!!!
                if (!exists("smalldata")){
                    smalldata <- bigdata[bigdata$ID == x & bigdata$time == alltime[y], ]
                } else {
                    smalldata <- rbind(smalldata, bigdata[bigdata$ID == x & bigdata$time == alltime[y], ])
                }
            }
        }
    }
}

} # End of if
browser()

# Idea: Plot scatterplot of tree precip vs grass precip for only MCS where >10% of each cover type exists
ids <- as.numeric(dimnames(table(smalldata$ID[!is.na(smalldata$Class)]))[[1]])
q90 <- function(x, na.rm=T){quantile(x, probs=0.9, na.rm=T)}
try(rm(list=c("compall", "compall2")))
for (xx in ids){
    print(xx)
    tmp <- aggByHour(smalldata, xx, q90)
    tmp2 <- aggByHour(smalldata, xx, mean)
    
    if (!exists("compall")){
        compall <- data.frame(ID=xx, tmp)
    } else {
        compall <- rbind(compall, data.frame(ID=xx, tmp))
    }

    if (!exists("compall2")){
        compall2 <- data.frame(ID=xx, tmp2)
    } else {
        compall2 <- rbind(compall2, data.frame(ID=xx, tmp2))
    }

}

ggplot(compall, aes(factor(Hour), fill=factor(greatest))) + 
    geom_bar() +
    labs(title="Frequency of preference of intense precipitation for land cover types by time of day when MCS covers all 3 cover types", x="Hour", y="Count") +
    scale_fill_manual(values = as.character(myLUT[c(1,2,4),"Colours"]),
                      name = "Land Cover",
                      labels = c("Forest", "Grass","Boundary"))
    
ggplot(compall2, aes(factor(Hour), fill=factor(greatest))) + 
    geom_bar() +
    labs(title="Frequency of preference of intense precipitation for land cover types by time of day when MCS covers all 3 cover types", x="Hour", y="Count") +
    scale_fill_manual(values = as.character(myLUT[c(1,2,4),"Colours"]),
                      name = "Land Cover",
                      labels = c("Forest", "Grass","Boundary"))

print("90th percentile")
round(table(compall$PrecipTree > compall$PrecipGrass & compall$PrecipTree > compall$PrecipBndry) / len, 2)
round(table(compall$PrecipGrass > compall$PrecipBndry & compall$PrecipGrass > compall$PrecipTree) / len, 2)
round(table(compall$PrecipBndry > compall$PrecipGrass & compall$PrecipBndry > compall$PrecipTree) / len, 2)

print("Mean")
len <- length(compall$PrecipTree)
binconf(x=table(compall2$PrecipTree > compall2$PrecipGrass & compall2$PrecipTree > compall2$PrecipBndry)[2], n=len, alpha=0.05, method="wilson")
binconf(table(compall2$PrecipGrass > compall2$PrecipBndry & compall2$PrecipGrass > compall2$PrecipTree)[2], n=len, alpha=0.05, method="wilson")
binconf(table(compall2$PrecipBndry > compall2$PrecipGrass & compall2$PrecipBndry > compall2$PrecipTree)[2], n=len, alpha=0.05, method="wilson")
# round(table(compall2$PrecipTree > compall2$PrecipGrass & compall2$PrecipTree > compall2$PrecipBndry) / len, 2)
# round(table(compall2$PrecipGrass > compall2$PrecipBndry & compall2$PrecipGrass > compall2$PrecipTree) / len, 2)
# round(table(compall2$PrecipBndry > compall2$PrecipGrass & compall2$PrecipBndry > compall2$PrecipTree) / len, 2)

print("Mean; Afternoon")
j <- which(compall2$Hour >= 13 & compall$Hour < 20)
lenj <- length(compall2$PrecipTree[j])
round(binconf(table(compall2$PrecipTree[j] > compall2$PrecipGrass[j] & compall2$PrecipTree[j] > compall2$PrecipBndry[j])[2], n=lenj, alpha=0.05, method="wilson"), 3)
round(binconf(table(compall2$PrecipGrass[j] > compall2$PrecipBndry[j] & compall2$PrecipGrass[j] > compall2$PrecipTree[j])[2], n=lenj, alpha=0.05, method="wilson"), 3)
round(binconf(table(compall2$PrecipBndry[j] > compall2$PrecipGrass[j] & compall2$PrecipBndry[j] > compall2$PrecipTree[j])[2], n=lenj, alpha=0.05, method="wilson"), 3)
# round(table(compall2$PrecipTree[j] > compall2$PrecipGrass[j] & compall2$PrecipTree[j] > compall2$PrecipBndry[j]) / lenj, 2)
# round(table(compall2$PrecipGrass[j] > compall2$PrecipBndry[j] & compall2$PrecipGrass[j] > compall2$PrecipTree[j]) / lenj, 2)
# round(table(compall2$PrecipBndry[j] > compall2$PrecipGrass[j] & compall2$PrecipBndry[j] > compall2$PrecipTree[j]) / lenj, 2)

print("Mean; Early morning")
j <- which(compall2$Hour >= 1 & compall$Hour < 8)
lenj <- length(compall2$PrecipTree[j])
round(binconf(table(compall2$PrecipTree[j] > compall2$PrecipGrass[j] & compall2$PrecipTree[j] > compall2$PrecipBndry[j])[2], n=lenj, alpha=0.05, method="wilson"), 3)
round(binconf(table(compall2$PrecipGrass[j] > compall2$PrecipBndry[j] & compall2$PrecipGrass[j] > compall2$PrecipTree[j])[2], n=lenj, alpha=0.05, method="wilson"), 3)
round(binconf(table(compall2$PrecipBndry[j] > compall2$PrecipGrass[j] & compall2$PrecipBndry[j] > compall2$PrecipTree[j])[2], n=lenj, alpha=0.05, method="wilson"), 3)
# round(table(compall2$PrecipTree[j] > compall2$PrecipGrass[j] & compall2$PrecipTree[j] > compall2$PrecipBndry[j]) / lenj, 2)
# round(table(compall2$PrecipGrass[j] > compall2$PrecipBndry[j] & compall2$PrecipGrass[j] > compall2$PrecipTree[j]) / lenj, 2)
# round(table(compall2$PrecipBndry[j] > compall2$PrecipGrass[j] & compall2$PrecipBndry[j] > compall2$PrecipTree[j]) / lenj, 2)


dev.off()

