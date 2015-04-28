densify <- function(xy,n=5){
    ## densify a 2-col matrix
    cbind(dens(xy[,1],n=n),dens(xy[,2],n=n))
}

dens <- function(x,n=5){
    ## densify a vector
    out = rep(NA,1+(length(x)-1)*(n+1))
    ss = seq(1,length(out),by=(n+1))
    out[ss]=x
    for(s in 1:(length(x)-1)){
        out[(1+ss[s]):(ss[s+1]-1)]=seq(x[s],x[s+1],len=(n+2))[-c(1,n+2)]
    }
    out
}

simplecentre <- function(xyP,dense){
    require(deldir)
    require(splancs)
    require(igraph)
    require(rgeos)
    
    ### optionally add extra points
    if(!missing(dense)){
        xy = densify(xyP,dense)
    } else {
        xy = xyP
    }
    
    ### compute triangulation
    d=deldir(xy[,1],xy[,2])
    
    ### find midpoints of triangle sides
    mids=cbind((d$delsgs[,'x1']+d$delsgs[,'x2'])/2,
               (d$delsgs[,'y1']+d$delsgs[,'y2'])/2)
    
    ### get points that are inside the polygon 
    sr = SpatialPolygons(list(Polygons(list(Polygon(xyP)),ID=1)))
    ins = over(SpatialPoints(mids),sr)
    
    ### select the points
    pts = mids[!is.na(ins),]
    
    dPoly = gDistance(as(sr,"SpatialLines"),SpatialPoints(pts),byid=TRUE)
    pts = pts[dPoly > max(dPoly/1.5),]
    
    ### now build a minimum spanning tree weighted on the distance
    G = graph.adjacency(as.matrix(dist(pts)),weighted=TRUE,mode="upper")
    T = minimum.spanning.tree(G,weighted=TRUE)
    
    ### get a diameter
    path = get.diameter(T)
    
    if(length(path)!=vcount(T)){
        stop("Path not linear - try increasing dens parameter")
    }
#     browser()
    ### path should be the sequence of points in order
    list(pts=pts[path,],tree=T)
    
}

onering=function(p, i){p@polygons[[1]]@Polygons[[i]]@coords}

capture = function(){p=locator(type="l")
                     SpatialLines(list(Lines(list(Line(cbind(p$x,p$y))),ID=1)))}

s = capture()
p = gBuffer(s,width=0.2)
plot(p,col="#cdeaff")
plot(s,add=TRUE,lwd=3,col="red")
# scp = simplecentre(onering(p))
scp = simplecentre(onering(rbnd.pol, 1))
lines(scp$pts,col="white")

# Get polygons with number of vertices > 45
maxid <- 0
maxval<- 0
for (i in 1:length(r2p@polygons[[1]]@Polygons)){
    nr <- nrow(onering(r2p, i))
    if (nr > 45){ print(paste(i,nr, sep=" : "))}
    if (nr > maxval){maxval <- nr ; maxid <- i}
}
print(paste("Max ID : ",maxid, "(",maxval,")", sep=""))

# Make a new feature for each polygon ...
ancils <- loadAllAncils(myproj=myproj, nBndClass=1, model="rb5216.4km.std", vegThreshold=vegThreshold, bndDef=bndDef, nBuf=1, overwrite=F)
r2p <- rasterToPolygons(ancils$mycl.f, fun = function(x){x==4}, n=8, dissolve=T)
plist <- vector("list")
for (i in 1:length(r2p@polygons[[1]]@Polygons)){
    plist[[i]] <- Polygons(list(r2p@polygons[[1]]@Polygons[[i]]), paste("s",i,sep=""))
}
plist.sp <- SpatialPolygons(plist, 1:length(plist))
