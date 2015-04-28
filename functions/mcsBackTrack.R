mcsBackTrack <- function(myxtab, bt_results, thisdate, thisdate.p1){
    
    # Go through myxtab, and recode each cluster to either a new ID or the overlapping ID
    
    # Cases may include (NB: going backwards in time, previous time step = 1 step in future):
    #   - Regular tracking  : one overlap between current time step and previous
    #   - Dissipation_bk    : Occurred in current time step, but not 1 time step in the future
    #   - Generation        : Occurred 1 time step in the future, but not in current time step
    #   - Splitting         : One cluster in current time step overlaps with 2 or more in future time step
    #   - Merging           : Two or more clusters in current time step overlap with only one in the future time step
    #   - Anything else     : May include combinations of splits and merges

    
    # Add extra column to myxtab
    myxtab$recodeTo <- NA
    myxtab$class    <- NA
    myxtab$related  <- NA
    myxtab$rowSums  <- rowSums(myxtab, na.rm=T)
    myxtab_colSums  <- colSums(myxtab, na.rm=T)
    
    # Create empty vectors for recoding clusters in t+1 (where splits, merges or initiations occur)
    recode_p1.fromID  <- vector("numeric")
    recode_p1.toID    <- vector("numeric")
    recode_p1.class   <- vector("character")
    recode_p1.related <- vector("character")    
    
    # Get row and col IDs that equal 0
    col0 <- which(colnames(myxtab)==0)
    row0 <- which(rownames(myxtab)==0)
    
    # Get next available unique ID
    maxID <- max(bt_results$ID)
    newid <- maxID + 1
    
    ###########################################################################################
    # Always recode 0 in current time step to 0
    myxtab[row0,"recodeTo"] <- 0
    
    ###########################################################################################
    #   - Regular tracking  : one overlap between current time step and previous
    
    # Get index of rows and columns with data in
    i_col <- which(!colnames(myxtab) %in% c("0","recodeTo","rowSums","class","related"))
    i_row <- which(!rownames(myxtab) %in% "0")
    
    # Loop through columns
    for (i in i_col){
        # If one row contains > 95% of non-zero pixels, then it is regular tracking
        pc_non0 <- myxtab[i_row,i] / sum(myxtab[i_row,i]) # % of non zero changes
        x <- which(pc_non0 > 0.95) # Only 1 or empty
        if (length(x) == 1){ 
            toID <- as.numeric(colnames(myxtab)[i])
            myxtab[i_row,][x,"recodeTo"] <- toID
            myxtab[i_row,][x,"class"] <- "regularTracking"
            # If any small fractions, set to 0
            ri <- which(pc_non0 <= 0.05 & pc_non0 != 0)
            if (any(ri)){
                for (r in ri){
                    # Make sure that this column holds > 95% of total pixels in the row
                    r_pc <- myxtab[i_row,i][r] / myxtab[i_row,"rowSums"][r]
                    if (r_pc > 0.95){
                        myxtab[i_row,"recodeTo"][r] <- toID
                        myxtab[i_row,"class"][r] <- "regularTracking (tiny split)"
                    }
                }                
            }
        } 
    }
    
    
    ###########################################################################################
    #   - Dissipation_bk    : Occurred in current time step, but not 1 time step in the future
    i_diss      <- myxtab[,col0] == myxtab$rowSums
    i_disslen   <- length(which(i_diss==TRUE))
    if (i_disslen > 0){
        newids <- newid + 0:(i_disslen-1)
        # myxtab[i_diss, "recodeTo"] <- newids
        # myxtab[i_diss, "class"] <- "dissipation"
        myxtab[i_diss, "recodeTo"] <- 0
        myxtab[i_diss, "class"] <- "POP not related to an MCS"
        
        newid <- max(newids) + 1
    }        

    
    ###########################################################################################
    #   - Generation        : Occurred 1 time step in the future, but not in current time step
    # THIS IS THE ACTUAL POINT OF INITIATION (AT FUTURE TIME STEP)
    # NEED TO CHANGE DESCRIPTION OF FUTURE TIMESTEP TO "GENERATION"
    # If Future timestep cell count == current timestep 
    i_inits <- which(myxtab[row0,i_col] == myxtab_colSums[i_col])
    inits   <- as.numeric(attributes(myxtab[,i_col][i_inits])$name)
    recode_p1.fromID  <- c(recode_p1.fromID, inits)
    recode_p1.toID    <- c(recode_p1.toID, inits)
    recode_p1.class   <- c(recode_p1.class, rep("Generation", length(inits)))
    recode_p1.related <- c(recode_p1.related, rep(NA, length(inits)))
    #recode_p1.related <- c(recode_p1.related, rep(I(list(NA)), length(inits)))
    
    ###########################################################################################
    #   - Splitting         : One cluster in current time step overlaps with 2 or more in future time step
        
    # Get the number of splits per ID @ current timestep
    if(class(myxtab[i_row, i_col]) == "data.frame"){
        splitcount <- apply(myxtab[i_row, i_col], 1, FUN=function(x){
            y<-100*(x/sum(x)); 
            z <- which(y>0 & as.numeric(attributes(x)$names))!=0; 
            return(length(z))})
    } else {
        # If there's only 1 column, the object class is "integer", and there would be no future splits, because a split needs at least 2 columns (i.e. 1 present cluster splitting into 2 or more future clusters)
        splitcount <- 0
    }
    # Loop through split count, pick out all rows > 1, and check if the cols qualify for a merge
    if(any(splitcount>1)){
        # Which rows have > 1 corresponding column?
        rowsplit <- attributes(splitcount[splitcount>1])$names
        for (rs in rowsplit){
            # For this row, which columns overlap?
            alloverlaps <- as.numeric(colnames(myxtab[rs,i_col][which(myxtab[rs,i_col]>0)]))
            # For this row, which col(s) have the largest overlap
            max_i <- myxtab[rs,i_col] %in% max(myxtab[rs,i_col])
            largest_overlap <- as.numeric(attributes(myxtab[rs,i_col][max_i])$names)
            if (length(largest_overlap) > 1){
                print(paste("More than 1 columns are tied for the greatest overlap with row",rs))
                
                # Which columns are tied?
                i_tied <- which(myxtab[rs,i_col] == max(myxtab[rs,i_col]))
                tied <- attributes(myxtab[rs,i_col][i_tied])$name
                
                # Which column has the greatest % overlap with this row?
                # print(100 * myxtab[rs,tied] / colSums(myxtab[,tied]))
                # Following returns the column name with the maximum value, or if it's a tie, the column name with the smallest value
                tied_max <- attributes(which.max(100 * myxtab[rs,tied] / colSums(myxtab[,tied])))$name
                # Set the largest_overlap to this column
                largest_overlap <- tied_max
                
            }
            # Is the row's largest_overlap in agreement with the column's largest overlap?
            # i.e. does the future timestep also merge?
            ip1 <- which.max(myxtab[i_row,as.character(largest_overlap)])
            maxrow_forthiscol <- rownames(myxtab[i_row,])[ip1]
            # Get the remaining overlaps (i.e. those that are not the largest)
            remaining_overlaps <- alloverlaps[-which(alloverlaps %in% largest_overlap)]
            
            if (rs != maxrow_forthiscol){
                print("A merge occurs at the same time as a split")
                # This is doing the job of the merge section below ...
                myxtab[maxrow_forthiscol,"recodeTo"] <- largest_overlap
                myxtab[maxrow_forthiscol,"class"] <- "Merges and splits in next timestep"
                # No related clusters here
                                
                if (any(remaining_overlaps)){
                    myxtab[rs,"recodeTo"] <- newid
                    myxtab[rs,"class"] <- paste("Smaller part of a merge to ID ",largest_overlap,", and a split to ID(s) ",paste(remaining_overlaps, collapse=","),sep="")
                    myxtab[rs,"related"] <- paste(alloverlaps, collapse=",")
                    newid <- newid + 1
                }
            } 
            
            if (rs == maxrow_forthiscol) {
                # Change ID of this row to column with largest overlap
                myxtab[rs,"recodeTo"] <- largest_overlap
                # Set class of this row to "regular tracking / splits in next timestep"
                myxtab[rs,"class"]    <- "regular tracking / splits in next timestep (largest overlap)"
                # Get minor overlaps (i.e. splits in next timestep that do not retain same ID)
                # i_overlaps_minor <- which(myxtab[rs,i_col] > 0 & as.numeric(colnames(myxtab[rs,i_col])) != largest_overlap)
                # overlaps_minor   <- as.numeric(colnames(myxtab[rs,i_col][i_overlaps_minor]))
                myxtab[rs,"related"]  <- paste(remaining_overlaps, collapse=",")
            }

            # if(length(remaining_overlaps)>1){ browser() }
            dupes <- which(recode_p1.fromID %in% remaining_overlaps)
            if(any(dupes)){
                for (d in dupes){
                    recode_p1.class[d] <- "Merges, but also includes one or more minor split"
                    # Related can be either a newid or the largest_overlap
                    recode_p1.related[d] <- paste(recode_p1.related[d], myxtab[rs,"recodeTo"],sep=",")
                }
                
            } else {
                # If this ID doesn't already exist in the recode table ...
                recode_p1.fromID <- c(recode_p1.fromID, remaining_overlaps)
                recode_p1.toID <- c(recode_p1.toID, remaining_overlaps)
                recode_p1.class <- c(recode_p1.class, rep(paste("regular tracking / splits this timestep (minor parts from ID ",myxtab[rs,"recodeTo"],")",sep=""), length(remaining_overlaps)))
                recode_p1.related <- c(recode_p1.related, rep(myxtab[rs,"recodeTo"], length(remaining_overlaps)))
            }
        }
    } else {
        print("No splits in this timestep")
    }
    
    ###########################################################################################
    #   - Merging           : Two or more clusters in current time step overlap with only one in the future time step

    # Get the number of merges per ID @ future timestep
    if(class(myxtab[i_row, i_col]) == "data.frame"){
        mergecount <- apply(myxtab[i_row, i_col], 2, FUN=function(x){y<-100*(x/sum(x)); z <- which(y>0 & as.numeric(attributes(x)$names))!=0; length(z)})
    } else {
        # If there's only 1 column, the object class is "integer"
        # There may still be multiple rows in the current timestep that overlap with 1 cluster in the future timestep
        mergecount <- length(which(myxtab[i_row, i_col]>0))
    }
    
    if (any(mergecount > 1)){
        print("We have some clusters in t that merge in t+1")
        
        # All merges need a new ID in the current timestep (since we are tracking backwards)
        # Also, need to recode t+1 entries to "Merges in this timestep"
        merges <- attributes(mergecount[which(mergecount > 1)])$name
        for (m in merges){
            # Which rows are the merge rows for this column?
            mr <- rownames(myxtab[i_row,][which(myxtab[i_row,m]>0),])
            # Which of these rows contributes the most to the merge?
            mr.grtst <- rownames(myxtab[i_row,])[which.max(myxtab[i_row,m])]
            # Write info out to myxtab
            myxtab[mr.grtst,"class"] <- "Merges in next timestep (largest part)"
            myxtab[mr.grtst,"recodeTo"] <- as.numeric(m)
            # Remove mr.grtst from mr
            mr.smlst <- mr[-which(mr == mr.grtst)]
            # What are the new IDs?
            newids <- newid + 0:(length(mr.smlst)-1)
#             browser()
            # Write related 
            for (mri in mr.smlst){
                # Set the ID for the largest cluster in this merge to 
                #if (thisdate == as.POSIXct("2006-08-19 02:50:00")){browser()}
                # Write info out to myxtab
                myxtab[mri,"class"] <- "Merges in next timestep (smaller part)"
                myxtab[mri,"recodeTo"] <- newids[which(mr.smlst %in% mri)]
                
#                 if (!any(relatedvec %in% m)){
                if(is.na(myxtab[mri,"related"])){ 
                    myxtab[mri,"related"] <- m 
                } else {
                    relatedvec <- unlist(strsplit(myxtab[mri,"related"], split=","))
                    myxtab[mri,"related"] <- paste(c(relatedvec, m), collapse=",")
                }
#                 }
            }
            
            # Get the next available new ID
            newid <- max(newids) + 1
            
        }
    }
    
    ###########################################################################################
    #   - Anything else     : May include combinations of splits and merges
    out.rctab <- data.frame(fromID=recode_p1.fromID, toID=recode_p1.toID, class=recode_p1.class, related=recode_p1.related)
#     for (x in 1:length(recode_p1.related)){
#         out.rctab[x,"related"] <- I(list(recode_p1.related[[x]]))
#     }
    
    if(length(which(is.na(myxtab$recodeTo))) > 0 ){
        print(which(is.na(myxtab$recodeTo)))
        print("This means there are still some cases that I have missed")
        browser()        
    } else {
        print("All cases have been successfully coded")
        return(list(myxtab, bt_results, out.rctab))
    }
        
}