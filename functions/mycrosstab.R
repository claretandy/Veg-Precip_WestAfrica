mycrosstab <- function(r1, r2){
    
    myxtab <- crosstab(r1, r2)
    myxtab2 <- reshape(myxtab, v.names="Freq", idvar="Var1", timevar="Var2", direction="wide")
    colnames(myxtab2) <- sub("Freq.","",colnames(myxtab2))
    myxtab2 <- myxtab2[,colnames(myxtab2)!="NA"]
    myxtab2 <- myxtab2[!is.na(myxtab2$Var1),] # Drop the row with NA
    rownames(myxtab2) <- myxtab2$Var1
    myxtab2 <- myxtab2[,-1]
    
    return(myxtab2)
}
