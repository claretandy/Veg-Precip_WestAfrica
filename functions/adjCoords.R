adjCoords <- function(r){
    require(raster)
    
    this.res <- res(r)[2]
    
    if ( this.res <= 0.0136){
        
        r <- raster::shift(r, x=0.0068, y=0.00675)
        
    } else {
        # Check that x coordinates are located correctly
        if (round(xmin(r) - round(xmin(r), 1), 5) != 0){
            r <- raster::shift(r, x=0.01798)
        }
        
        # Check that y coordinates are located correctly
        if (round(ymin(r) - round(ymin(r), 1), 5) != 0){
            r <- raster::shift(r, y=0.018)
        }
        
        # Check that pixel resolution is correct
        if (any(res(r) - c(0.036, 0.036) != c(0,0))){
            res(r) <- 0.036
        }
        
    }
    
    # Finally, check that the extent is centred around 0, not 360
    if (extent(r)@xmax > 360){
        r <- raster::shift(r, x=-360)
    }
    
    return(r)
}
