library(raster)
library(ggplot2)
library(plyr)
source('functions/adjCoords.R')
source("getMyData.R")


# Plot fractional cover of BL + Shrub compared to mean total surface evapotranspiration
in_data_dir <- '/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/'
timestep <- "5min" # "10min" # "avg"
threshold <- 1000 # Threshold for minimum size of an MCS
myproj <- "rp" # Should be "ll" for paper
models <- "rb5216.4km.std"
id <- "s"
vegThreshold <- 0.3 
bndDef       <- 2
nBuf         <- 3

pdf("../../Results/treeFrac_vs_evap.pdf", width=12, height=9)

# Get precip and evap data
mydata <- getMyData(timestep=timestep, var="lsrain", overwrite=F)
pre <- adjCoords(mydata)
evap <- adjCoords(brick('/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/3223_gt20060815.nc'))

# Get tree fraction
ancils_rp <- loadAllAncils(myproj="rp", nBndClass=1, model="rb5216.4km.std", vegThreshold=vegThreshold, bndDef=bndDef, nBuf=nBuf, overwrite=F)
myveg      <- ancils_rp$myveg
tree_frac <- myveg[[1]] + myveg[[2]] + myveg[[5]]
tree_bins <- cut(tree_frac, breaks=seq(0,1,0.1)) / 10

print(
    levelplot(crop(tree_frac, extent(ancils_rp$spp)), main='Tree fractional cover') + layer(sp.polygons(ancils_rp$spp))    
    )

# Get P - E
p_minus_e_file <- paste(in_data_dir,'p_minus_e.nc', sep='')
if (!file.exists(p_minus_e_file)){
    system(paste('cdo sub ',in_data_dir,'4203.nc ',in_data_dir,'3223_gt20060815.nc ',p_minus_e_file, sep=''))
}
p_minus_e <- adjCoords(brick(p_minus_e_file))
p_minus_e_timmean <- adjCoords(brick('/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/p_minus_e_timmean.nc'))

print(
    levelplot(crop(p_minus_e_timmean, extent(ancils_rp$spp)) * 60*60, main='Mean p - e (mm/hour)') + layer(sp.polygons(ancils_rp$spp))    
    )


# Plot some stats for P-E over tree bins and zones ...

mydf <- data.frame(tree_bins=getValues(crop(tree_bins, extent(ancils_rp$spp))), 
                   zones=getValues(crop(ancils_rp$spp.r, extent(ancils_rp$spp))), 
                   pme=as.vector(getValues(crop(p_minus_e_timmean, extent(ancils_rp$spp)))*60*60))

mydata <- ddply(mydf, .(zones, tree_bins), summarize,
                N    = sum(!is.na(pme)),
                mean = mean(pme, na.rm=T),
                sd   = sd(pme, na.rm=T),
                se   = sd / sqrt(N)
                )

print(
    ggplot(mydata, aes(x=tree_bins, y=mean)) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) + geom_line() + geom_vline(xintercept=0.3, colour='red') + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="P-E (Precipitation minus Total Surface Moisture Flux) by zone and tree fraction", x="Tree Fraction", y="Mean P-E (mm/hour)")
)

alldata <- cbind(mydata, Var='P-E')

# Plot the same for mean evapotranspiration ....
evap_mean <- adjCoords(brick('/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/3223_timmean.nc'))

mydf <- data.frame(tree_bins=getValues(crop(tree_bins, extent(ancils_rp$spp))), 
                   zones=getValues(crop(ancils_rp$spp.r, extent(ancils_rp$spp))), 
                   evap=as.vector(getValues(crop(evap_mean, extent(ancils_rp$spp)))*60*60))

mydata <- ddply(mydf, .(zones, tree_bins), summarize,
                N    = sum(!is.na(evap)),
                mean = mean(evap, na.rm=T),
                sd   = sd(evap, na.rm=T),
                se   = sd / sqrt(N)
)
print(
    ggplot(mydata, aes(x=tree_bins, y=mean)) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) + geom_line() + geom_vline(xintercept=0.3, colour='red') + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="Mean Surface Moisture Flux by Zone and Tree Fraction", x="Tree Fraction", y="Mean Surface Moisture Flux (mm/hour)")    
    )

alldata <- rbind(alldata, cbind(mydata, Var='Evapotranspiration'))

# Plot the same for mean precipitation ....
precip_mean <- adjCoords(brick('/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/4203_timmean.nc'))

mydf <- data.frame(tree_bins=getValues(crop(tree_bins, extent(ancils_rp$spp))), 
                   zones=getValues(crop(ancils_rp$spp.r, extent(ancils_rp$spp))), 
                   pre=as.vector(getValues(crop(precip_mean, extent(ancils_rp$spp)))*60*60))

mydata <- ddply(mydf, .(zones, tree_bins), summarize,
                N    = sum(!is.na(pre)),
                mean = mean(pre, na.rm=T),
                sd   = sd(pre, na.rm=T),
                se   = sd / sqrt(N)
)

print(
    ggplot(mydata, aes(x=tree_bins, y=mean)) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) + geom_line() + geom_vline(xintercept=0.3, colour='red') + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="Mean Precipitation by Zone and Tree Fraction", x="Tree Fraction", y="Mean Precipitation Rate (mm/hour)")
    )

alldata <- rbind(alldata, cbind(mydata, Var='Precipitation'))

print(
    ggplot(alldata, aes(x=tree_bins, y=mean, colour=Var)) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) + geom_line() + geom_vline(xintercept=0.3, colour='black', linetype="longdash") + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="Moisture Fluxes by Zone and Tree Fraction", x="Tree Fraction", y="Moisture Flux (mm/hour)")
)


# What if we exclude grid cells where it has rained in the last 3/6/9/12 hours?
# 1. Time since last precipitation (>1mm/hour) ...
timeSinceLastRain <- function(pre){
#     lrd <- brick(pre); lrd[] <- 0 # as.numeric(getZ(pre)[1])
#     tslr <- brick(pre); tslr[] <- 0
    tslr_file <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/timeSinceLastRain.tif"
    
    if (file.exists(tslr_file)){
        tslr <- brick(tslr_file)
    } else {
        lrd <- vector('list')
        lrd[[1]] <- raster(pre, layer=0)
        lrd[[1]][] <- 0
        
        tslr <- vector('list')
        tslr[[1]] <- raster(pre, layer=0)
        tslr[[1]][] <- 0
        
        for(i in 2:nlayers(pre)){
            print(i)
            #thisdate <- as.numeric(getZ(pre)[i])
            thisdate <- i
            lrd[[i]] <- overlay(pre[[i]], lrd[[i-1]], fun=function(x, y){ii <- which(x>(1/3600)); x[ii] <- thisdate; x[-ii] <- y[-ii]; return(x)})
            
            tslr[[i]] <- calc(lrd[[i]], fun=function(x){ z <- thisdate - x; return(z) })
        }
        tslr <- stack(tslr)
        tslr <- writeRaster(tslr, filename=tslr_file)
    }
    
    return(tslr)
}
tslr <- timeSinceLastRain(pre)

# NB: All af these three alternatives were replaced by a python script that was MUCH faster 
evap_norain_alt1 <- function(tslr, evap, hrs=3, overwrite=F){
    evap_norainXhrs_file <- paste("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/evap_norain_",hrs,"hrs.tif", sep="")
    if (file.exists(evap_norainXhrs_file) & !overwrite){
        evap_norainXhrs <- brick(evap_norainXhrs_file)
    } else {
        evap_norainXhrs <- vector("list")
        for (i in 1:nlayers(tslr)){
            print(i)
            tslr_xhrs <- tslr[[i]] > 5*60*hrs
            evap_norainXhrs[[i]] <- tslr_xhrs * evap[[i]]
        }
        print("Creating a raster stack ...")
        evap_norainXhrs <- stack(evap_norainXhrs)
        print("Writing data to file ...")
        evap_norainXhrs <- writeRaster(evap_norainXhrs, filename = evap_norainXhrs_file, overwrite=T)
    }
    return(evap_norainXhrs)
}

evap_norain_alt2 <- function(tslr, evap, hrs=3, overwrite=F){
    evap_norainXhrs_file <- paste("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/evap_norain_",hrs,"hrs.tif", sep="")
    if (file.exists(evap_norainXhrs_file) & !overwrite){
        evap_norainXhrs <- brick(evap_norainXhrs_file)
    } else {
        tslr_Xhrs <- tslr < 5*60*hrs
        evap_norainXhrs <- tslr_Xhrs * evap
        print("Writing data to file ...")
        evap_norainXhrs <- writeRaster(evap_norainXhrs, filename = evap_norainXhrs_file, overwrite=T)
    }
    return(evap_norainXhrs)
}

evap_norain_alt3 <- function(tslr, evap, hrs=3, overwrite=F){
    evap_norainXhrs_file <- paste("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/evap_norain_",hrs,"hrs.tif", sep="")
    if (file.exists(evap_norainXhrs_file) & !overwrite){
        evap_norainXhrs <- brick(evap_norainXhrs_file)
    } else {
        evap_norainXhrs <- overlay(evap, tslr, fun=function(x, y){
            ii <- which(y < 5*60*3)
            x[ii] <- NA
            return(x)},
            filename=evap_norainXhrs_file,
            overwrite=T)
    }
    return(evap_norainXhrs)
}

# Speed tests ...
# system.time(evap_norain_alt1(subset(tslr, 1:5), subset(evap, 1:5), hrs=6, overwrite=T))
# system.time(evap_norain_alt2(subset(tslr, 1:5), subset(evap, 1:5), hrs=6, overwrite=T))
# system.time(evap_norain_alt3(subset(tslr, 1:5), subset(evap, 1:5), hrs=6, overwrite=T))
# 
# evap_norainXhrs <- evap_norain_alt1(subset(tslr, 1:5), subset(evap, 1:5), hrs=6, overwrite=T)
# system.time(calc(evap_norainXhrs, fun=mean, na.rm=T, filename=paste("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/evap_norain_",hrs,"hrs_MEAN.tif", sep="")))

for (x in c(3,6,12,24,36,48)){
    print(paste("Hours since last rain:",x))
    #evap_norainXhrs <- evap_norain_alt1(tslr, evap, hrs=x, overwrite=F)
    running_tot <- paste("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/evap_norain_",x,"hrs_gt20060815_RUNNING_TOT.nc", sep="")
    running_sum <- paste("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/evap_norain_",x,"hrs_gt20060815_RUNNING_SUM.nc", sep="")
    outmean <- paste("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/evap_norain_",x,"hrs_gt20060815_TIMEMEAN_BYDAY.nc", sep="")
    if (file.exists(outmean)){
        evap_norainXhrs_mean <- adjCoords(raster(outmean))
        evap_norainXhrs_running_tot <- adjCoords(raster(running_tot))
        evap_norainXhrs_running_sum <- adjCoords(raster(running_sum))
        evap_norainXhrs_running_tot[evap_norainXhrs_running_tot < 144] <- NA
        evap_norainXhrs_running_sum <- mask(evap_norainXhrs_running_sum, evap_norainXhrs_running_tot)
        evap_norainXhrs_mean <- mask(evap_norainXhrs_mean, evap_norainXhrs_running_tot)
    } else {
        print("Can\'t find time mean from iris")
        browser()
        #evap_norainXhrs_mean <- calc(evap_norainXhrs, fun=mean, na.rm=T, filename=outmean, overwrite=T)        
    }
    assign(paste("evap_norain",x,"hrs_mean",sep=""), evap_norainXhrs_mean)
    assign(paste("evap_norain",x,"hrs_running_tot",sep=""), evap_norainXhrs_running_tot)
    assign(paste("evap_norain",x,"hrs_running_sum",sep=""), evap_norainXhrs_running_sum)
}

tree_bins_cropped <- crop(tree_bins, extent(ancils_rp$spp))
zones_cropped     <- crop(ancils_rp$spp.r, extent(ancils_rp$spp))

mydf <- data.frame(tree_bins=getValues(tree_bins_cropped), 
                zones=getValues(zones_cropped), 
                evap=as.vector(getValues(crop(evap_norain3hrs_mean, extent(ancils_rp$spp)))/24),
                tot=as.vector(getValues(crop(evap_norain3hrs_running_tot, extent(ancils_rp$spp)))/24),
                evap_sum=as.vector(getValues(crop(evap_norain3hrs_running_sum, extent(ancils_rp$spp)))/24),
                tslr="3 Hours")

mydf <- rbind(mydf, 
              data.frame(tree_bins=getValues(tree_bins_cropped), 
                zones=getValues(zones_cropped), 
                evap=as.vector(getValues(crop(evap_norain6hrs_mean, extent(ancils_rp$spp)))/24),
                tot=as.vector(getValues(crop(evap_norain6hrs_running_tot, extent(ancils_rp$spp)))/24),
                evap_sum=as.vector(getValues(crop(evap_norain6hrs_running_sum, extent(ancils_rp$spp)))/24),
                tslr="6 Hours") 
)

mydf <- rbind(mydf, 
              data.frame(tree_bins=getValues(tree_bins_cropped), 
                zones=getValues(zones_cropped), 
                evap=as.vector(getValues(crop(evap_norain12hrs_mean, extent(ancils_rp$spp)))/24),
                tot=as.vector(getValues(crop(evap_norain12hrs_running_tot, extent(ancils_rp$spp)))/24),
                evap_sum=as.vector(getValues(crop(evap_norain12hrs_running_sum, extent(ancils_rp$spp)))/24),
                tslr="12 Hours"
              )
)

mydf <- rbind(mydf, 
              data.frame(tree_bins=getValues(tree_bins_cropped), 
                         zones=getValues(zones_cropped), 
                         evap=as.vector(getValues(crop(evap_norain24hrs_mean, extent(ancils_rp$spp)))/24),
                         tot=as.vector(getValues(crop(evap_norain24hrs_running_tot, extent(ancils_rp$spp)))/24),
                         evap_sum=as.vector(getValues(crop(evap_norain24hrs_running_sum, extent(ancils_rp$spp)))/24),
                         tslr="24 Hours")
)

mydf <- rbind(mydf, 
              data.frame(tree_bins=getValues(tree_bins_cropped), 
                         zones=getValues(zones_cropped), 
                         evap=as.vector(getValues(crop(evap_norain36hrs_mean, extent(ancils_rp$spp)))/24),
                         tot=as.vector(getValues(crop(evap_norain36hrs_running_tot, extent(ancils_rp$spp)))/24),
                         evap_sum=as.vector(getValues(crop(evap_norain36hrs_running_sum, extent(ancils_rp$spp)))/24),
                         tslr="36 Hours")
)

mydf <- rbind(mydf, 
              data.frame(tree_bins=getValues(tree_bins_cropped), 
                         zones=getValues(zones_cropped), 
                         evap=as.vector(getValues(crop(evap_norain48hrs_mean, extent(ancils_rp$spp)))/24),
                         tot=as.vector(getValues(crop(evap_norain48hrs_running_tot, extent(ancils_rp$spp)))/24),
                         evap_sum=as.vector(getValues(crop(evap_norain48hrs_running_sum, extent(ancils_rp$spp)))/24),
                         tslr="48 Hours")
)

mydf <- subset(mydf, !is.na(evap) & !is.na(zones))

mydata <- ddply(mydf, .(zones, tree_bins, tslr), summarize,
                Ngridcells = length(tot),
                Nmeancount = mean(tot),
                N    = length(tot) * mean(tot),
                median_evap = median(evap_sum, na.rm=T),
                mean_evap = mean(evap_sum),
                mean = mean(evap_sum / tot, na.rm=T),
                sd   = sd(evap_sum / tot, na.rm=T),
                se   = sd / sqrt(Ngridcells)
)

print(
    ggplot(mydata, aes(x=tree_bins, y=mean*3600, colour=tslr)) + geom_errorbar(aes(ymin=3600*(mean-se), ymax=3600*(mean+se))) + geom_line() + geom_vline(xintercept=0.3, colour='black', linetype='dashed') + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="Mean Evapotranspiration by Zone and Tree Fraction", x="Tree Fraction", y="Mean Evaporation Rate (mm/hour)") + scale_colour_discrete(name="Time Since\nLast Rain")
)

print(
    ggplot(mydata, aes(x=tree_bins, y=mean_evap*3600, colour=tslr)) + geom_line() + geom_vline(xintercept=0.3, colour='black', linetype='dashed') + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="Mean Accumulated Evapotranspiration by Zone and Tree Fraction", x="Tree Fraction", y="Mean Evaporation Rate (mm/hour)") + scale_colour_discrete(name="Time Since\nLast Rain")
)

print(
    ggplot(mydata, aes(x=tree_bins, y=N*5, colour=tslr)) + geom_line() + geom_vline(xintercept=0.3, colour='black', linetype='dashed') + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="Total number of minutes contributing towards evaporation totals by zone and tree fraction", x="Tree Fraction", y="Minutes") + scale_colour_discrete(name="Time Since\nLast Rain")
)

print(
    ggplot(subset(mydata, tslr=='24 Hours'), aes(x=tree_bins, y=mean*3600, colour=tslr)) + geom_errorbar(aes(ymin=3600*(mean-se), ymax=3600*(mean+se))) + geom_line() + geom_vline(xintercept=0.3, colour='black', linetype='dashed') + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="Mean Evapotranspiration by Zone and Tree Fraction", x="Tree Fraction", y="Mean Evaporation Rate (mm/hour)") + scale_colour_discrete(name="Time Since\nLast Rain")
)

print(
    ggplot(subset(mydata, tslr=='24 Hours'), aes(x=tree_bins, y=mean*3600)) + geom_errorbar(aes(ymin=3600*(mean-se), ymax=3600*(mean+se))) + geom_line() + geom_vline(xintercept=0.3, colour='black', linetype='dashed') + scale_x_continuous(breaks=seq(0,1,0.2)) + facet_wrap(~zones, ncol=5) + labs(title="Mean Evapotranspiration by Zone and Tree Fraction", x="Tree Fraction", y="Mean Evaporation Rate (mm/hour)") + scale_colour_discrete(name="Time Since\nLast Rain")
)

dev.off()

