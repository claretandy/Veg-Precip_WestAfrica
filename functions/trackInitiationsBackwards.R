# Using MCS initiation points, track backwards to find the grid cell at which convection started.
source("functions/tracker_v2.R")
source("functions/adjCoords.R")
source("functions/mcsBackTrack.R")
source("functions/getCorrectDates.R")

# Runs with the following: 
# rm(list=ls()); source('~/Work/Projects/InternalSabbatical/Scripts/allruns/loadData_diurnal.R')

# for (x in list.files(path=paste("/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_10min"), pattern="mcs_tracking_1000km_*", full.names=T)){
#     print(basename(x))
#     r <- adjCoords(raster(x))
#     writeRaster(r, filename=extension(x, "rp.tif"), overwrite=T)
#     browser()
# }

getCluster <- function(thisdate, type, threshold, timestep){
    if (type == "mcs"){
        infile <- paste('/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_',threshold,'km_',timestep,'/mcs_tracking_',threshold,'km_', format(thisdate, "%d.%H%M"), '.tif', sep="")
        clust <- raster(infile)
    }
    if (type == "pop"){
        infile <- paste('/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/patches/',threshold,'km_',timestep,'/pop.',threshold,'km_',format(thisdate, "%d%H.%M"),'.tif',sep="")
        pop_binary <- adjCoords(raster(infile))
        clust <- clump(pop_binary, directions=8, gaps=F)
    }
    
    return(clust)
}


trackInitiationsBackwards <- function(results=results, ancils=ancils, threshold=threshold, mcs=mcs, pop=pop, inbr=inbr, id=id[x], timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F){
    
    resultsdir <- paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,"/", sep="")
    genpts <- results[results$class == "MCS Generation" | results$class == "MCS Generation / dissipation",]
    coordinates(genpts) <- ~x+y
            
    # print(levelplot(ancils$mycl.f, par.settings=rasterTheme(region=myLUT$Colours), xlab=NULL, ylab=NULL, margin=F) + layer(sp.polygons(land_simple)) + layer(sp.points(genpts, col="black")) + layer(sp.polygons(ancils$spp)))
    
    mydates  <- sort(unique(genpts$timestep), decreasing=T)
    gcd <- getCorrectDates(inbr)
    alldates <- gcd[[1]]
    inbr <- gcd[[2]]

    # Using the existing results dataframe, we need to set all MCS 'initiations' to 'regular tracking' because we will track all MCS backwards to find the exact starting point of all MCS
    bt_results <- results
    bt_results$related <- NA
    bt_results$type <- "MCS"
    bt_results[bt_results$class == "generation","class"] <- "regularTracking (MCS start)"
    bt_results[bt_results$class == "generation / dissipation","class"] <- "dissipation (MCS start)"
    
#     for (x in (length(alldates)-1):(length(alldates)-30)){
    for (x in (length(alldates)-1):1){ 
        
        # Get the timestep
        thisdate <- alldates[x]
        thisdate.p1 <- alldates[x+1]
        print(paste("thisdate",thisdate))
        print(paste("thisdate.p1",thisdate.p1))
        
        # Get the MCS initiations at this timestep
        thisgen <- genpts[genpts$timestep == thisdate,]
        thisgen.p1 <- genpts[genpts$timestep == thisdate.p1,]
        
        # Get pop and mcs clusters for this timestep
        pop.t <- getCluster(thisdate, "pop", threshold, timestep)*10000
        pop.t[is.na(pop.t)] <- 0
        
        # If we're on the last time step (i.e. first in the backtracking loop), the previous clusters (pop.new) need to come from the MCS in the last time step
        # if (alldates[x+1] == alldates[length(alldates)]){
        #     mcs.t <- getCluster(alldates[length(alldates)], "mcs", threshold)
        #     pop.new <- mcs.t
        # }
        
        # Get clusters from t+1 clumps
        # pop.fut = small scale convection + MCS initiations in future timestep
        if (exists("pop.new")){
            pop.fut <- pop.new
            pop.fut[is.na(pop.fut)] <- 0            
        } else {
            pop.fut <- raster(pop.t)
            pop.fut <- setValues(x = pop.fut, values = 0)
        }
        
        # If any MCS generations exist, get the MCS cluster image for this timestep 
        if (length(thisgen.p1) > 0){
            mcs.fut <- getCluster(thisdate.p1, "mcs", threshold, timestep)
            # Extract only the clusters that are initiations
            mcs.futinit <- mcs.fut * mcs.fut %in% thisgen.p1$ID
            mcs.futinit[is.na(mcs.futinit)] <- 0
            # browser()
            # Add MCS initiations to pop for current time step
            pop.fut <- pop.fut + mcs.futinit
        }
        
        ### Include MCS dissipation in pop.t so that we catch cases where an patch is > 1000km, then <1000km, then >1000km again. 
        ### In this case, they would be coded as 2 separate MCS, but we want only 1 ignition point.
        # Get all MCS in this timestep
        i_mcs_now <- which(bt_results$timestep == thisdate)
        i_dissipation <- grep("dissipation", bt_results[,"class"])
        diss_ids <- bt_results[i_dissipation[which(i_dissipation %in% i_mcs_now)],"ID"]
        if (length(diss_ids) > 0){
            mcs.t <- getCluster(thisdate, "mcs", threshold, "5min")
            starters <- cellStats(pop.t, "max")+10000
            mcs.t <- calc(mcs.t, fun = function(x){
                x[is.na(x)] <- 0; 
                x[!(x %in% diss_ids)] <- 0; 
                x[(x %in% diss_ids)] <- x[(x %in% diss_ids)] + starters; 
                return(x)})
            pop.t <- mcs.t + pop.t
        }
        
        
        # If there were no MCS initiations in t+1, and pop.fut doesn't exist, then move on to next timestep ...
        if (!exists("pop.fut")){
            if (length(thisgen) > 0){
                pop.fut <- mcsinit
                # bt_results <- cbind(thisgen@data, coordinates(thisgen)) # Starts the results data.frame for backtracking
            }
            print("Moving on to next timestep ...")
            next
        } else {
                                        
            # Do the overlay
            # mycrosstab comes from tracker_v2.R
            myxtab <- mycrosstab(pop.t, pop.fut)
            
            # Something like this, but will probably need to add in splitting, merging, etc
            # pop.clump : unique small scale features with unique IDs
            # pop.fut  : the small scale that we had from the future (previous) time step
            # pop.new   : a copy of pop.fut for writing the result to
            # myxtab    : crosstabulation between pop.fut (cols) and pop.clump (rows)
            # bt_results: backtracking results data.frame to be filled in
            # uniqID    : 
            
            myxtab$recodeTo <- NA
            myxtab$class    <- NA
            
            trk <- mcsBackTrack(myxtab, bt_results, thisdate, thisdate.p1)
            
            myxtab     <- trk[[1]]
            bt_results <- trk[[2]]
            out_rctab  <- trk[[3]]
            
            # For the last time step, check that myxtab doesn't have lots of "dissipation" rows in it
            print(myxtab[myxtab$recodeTo != 0,])
            print(tail(bt_results))
            print(out_rctab)
            
#             browser()
            
            if(nrow(out_rctab)>0){
                # Generations in future timestep ...
                ids <- out_rctab$fromID
                for (i in ids){
                    # Recode according to out_rctab. 
                    # Includes Generations and Splits (minor parts)
                    ri <- bt_results$ID == i & bt_results$timestep == thisdate.p1
                    oi <- out_rctab$fromID == i
                    if (length(out_rctab[oi,"toID"])>1){browser()}
                    bt_results[ri,"ID"] <- out_rctab[oi,"toID"]
                    bt_results[ri,"class"] <- as.character(out_rctab[oi,"class"])
                    bt_results[ri,"related"] <- as.character(out_rctab[oi,"related"])
                }
            }
            
            ##############################################################################
            # Now, do something with myxtab ...
            
            # Recode current timestep image to "toID"
            # Set up recoding data.frame
            rctab <- data.frame("from"=as.numeric(rownames(myxtab)), "to"=as.numeric(myxtab$recodeTo))
            
            # Substitutes from with to, and sets all other cells to NA. 
            print("Recoding pop.t using rctab ...")
            pop.rc <- subs(pop.t, rctab, by="from", which="to", subWithNA=T)
            
            # Write to bt_results
            # NB: For generations, records in the future timestep need to have their class changed to "Generation"
            # A few more fields for the results table
            print("Creating polygons from pop.rc to get centroids ...")
            pop4poly <- pop.rc
            pop4poly[pop4poly == 0] <- NA
            myfreq <- freq(pop4poly, useNA="no")
            trc.poly <- rasterToPolygons(pop4poly, dissolve=T)
            mypts <- coordinates(trc.poly)
            
            trc.data <- data.frame(ID=trc.poly@data$to, x=mypts[,1], y=mypts[,2])
            
            # Append to results data.frame
            if (all(trc.data$ID == myfreq[,"value"])){
                pixcount <- myfreq[,"count"]
            } else {
                pixcount <- as.numeric(myfreq[which(myfreq[,"value"] %in% trc.data$ID), "count"])               
            }
            myorder <- match(trc.data$ID, myxtab$recodeTo)
            newresults <- data.frame("ID"=trc.data$ID, "pixcount"=pixcount, "timestep"=rep(thisdate, length(trc.data$ID)),"class"=myxtab[myorder,"class"], "x"=trc.data$x, "y"=trc.data$y, "related"=myxtab[myorder,"related"], "type"=rep("POP",length(trc.data$ID)))
            
            ####################
            # Could argue here that we don't need to include backtracked POP that don't turn into MCS, but let's leave it in for the moment ....
            ##################
            bt_results <- rbind(bt_results, newresults)
            
            # Set pop.new for next timestep and write raster to file
            outfile <- paste("/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_",threshold,"km_",timestep,"/pop_backtracked_",format(thisdate, "%d.%H%M"), ".tif", sep="")
            pop.new <- writeRaster(pop.rc, filename=outfile, overwrite=T)
            
        } # end of actions if statement
        
        
        outresults_file <- paste("/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/backtracked_results_",threshold,"km_",timestep,".Rdata", sep="")
        save(bt_results, file=outresults_file)
        # Create backup with timestamp
        bkup_file <- sub(".Rdata",paste("_",format(Sys.time(), "%d%m%Y"), ".Rdata", sep=""), outresults_file)
        file.copy(from=outresults_file, to=bkup_file, overwrite=T)

    } # end of date loop
    return(bt_results)
} # end of function

