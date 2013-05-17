#' Tools and Analysis Functions for the Bering-Arctic Seal Surveys
#'
#' bass is a central package for creating basic/global functions and data
#' related to the Bering-Arctic Seal Surveys. The 2012 BOSS surveys are the
#' first surveys to rely on this package
#'
#' @docType package
#' @name bass
#' @aliases bass bass-package
#' @author Paul Conn & Josh London
NULL

#' Add a column to the SpatialPolygonsDataFrame object giving the proportion of area
#' of each cell that is covered by land
#'
#' This function returns an object of type \code{SpatialPolygonsDataFrame} that takes
#' an original SpatialPolygonsDataFrame and computes the proportion of area for each 
#' subgeometry that is covered by land.  "Land" must be represented as a list of SpatialPolygonsDataFrames.
#' 
#' @param Grid	master SpatialPolygonsDataFrame object for which land calculations are desired (separate calculations occur for each subgeometry)
#' @param Land  a list of SpatialPolygonsDataFrames that represent land cover 
#' @param rel   if TRUE (default), returns proportional area; if FALSE, returns absolute area
#' @author Paul Conn \email{paul.conn@@noaa.gov}
#' @return SpatialPolygonsDataFrame object, which includes proportion land cover as an additional column
#' @export
add.prop.land<-function(Grid,Land,rel=TRUE){
	n_cells=dim(Grid)[1]
	Land.area=rep(0,n_cells)
	for(icell in 1:n_cells){
		for(iobj in 1:length(Land)){ #loop over list
			I.intersect=gIntersects(Grid[icell,],Land[[iobj]],byid=TRUE)
			if(sum(I.intersect)>0){
				Int_id=which(I.intersect==TRUE)
				for(iind in 1:length(Int_id))Land.area[icell]=Land.area[icell]+gArea(gIntersection(Grid[icell,],Land[[iobj]][Int_id[iind],]))
			}
		}
	}
	if(rel==TRUE)Land.area=Land.area/gArea(Grid,byid=TRUE)
	Grid[["land_cover"]]=Land.area
	Grid
}

#' Add line to SpatialPolygonsDataFrame that gives distance from "southern ice edge" for the centroid of each grid cell.
#' This version function is depracated - see function get.edge.Bering for the current version
#'
#' @param Grid  master SpatialPolygonsDataFrame object to add to
#' @param Grid_points centroids of Grid (SpatialPoints object)
#' @param IceExtent  a SpatialPolygonsDataFrames giving ice extent (e.g. as imported from NIC shapefile)
#' @param proj  current projection
#' @param mean_adjust if TRUE, standardize distance covariate by dividing by its mean
#' @author Paul Conn \email{paul.conn@@noaa.gov}
#' @return SpatialPolygonsDataFrame object, which includes distance from southern ice edge as an additional column
#' @export
add.dist.s.ice.edge<-function(Grid,Grid_points,IceExtent,proj,mean_adjust=TRUE){
  bb=bbox(Grid)
  bbox.p=Polygon(cbind(c(bb[1,1],bb[1,1],bb[1,2],bb[1,2],bb[1,1]),c(bb[2,1],bb[2,2],bb[2,2],bb[2,1],bb[2,1])))
  bbox.ps=Polygons(list(bbox.p),1)
  bbox.sps = SpatialPolygons(list(bbox.ps)) 
  proj4string(bbox.sps)=CRS(proj) 
  ice_edge<-gIntersection(gBuffer(IceExtent,byid=TRUE,width=0), gBuffer(bbox.sps,width=0))  #buffer can fix issues with self-intersection (use gIsValid to figure out if imported geometries are valid)
  bb=bbox(ice_edge)
  southern_line=Line(cbind(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,1])))
  southern_line=Lines(list(southern_line),1)
  southern_line=SpatialLines(list(southern_line))
  proj4string(southern_line)=CRS(proj)
  dist_ice_edge=as.vector(gDistance(Grid_points,southern_line,byid=TRUE))
  if(mean_adjust)dist_ice_edge=dist_ice_edge/mean(dist_ice_edge)
  Grid=spCbind(Grid,dist_ice_edge)
  Grid
}

#' calculate southern ice edge from spatial polygons object giving sea ice concentration
#' @param SpDF  SpatialPolygonsDataFrame object giving sea ice concentration
#' @param proj projection associated with the object
#' @author Paul Conn \email{paul.conn@@noaa.gov}
#' @return SpatialLines object giving approx southern ice edge
#' @export
get.edge.Bering<-function(SpDF,proj){
  coords=coordinates(SpDF)
  x=sort(unique(coords[,1]))
  y=unique(coordinates(SpDF)[,2])
  y.min=min(y)
  y.max=max(y)
  n.points=length(x)
  Point.vec=matrix(0,n.points,2)
  for(i in 1:n.points){
    which.row=which(coords[,1]==x[i])
    cur.y=coords[which.row,2]
    cur.mat=data.frame(matrix(0,length(which.row),2))
    cur.mat[,1]=as.numeric(cur.y)
    cur.mat[,2]=as.numeric(SpDF[["ice_conc"]][which.row])
    cur.mat=cur.mat[order(cur.mat[,1]),]
    cur.mat[,2]=(cur.mat[,2]<=0.1)
    cur.mat[1,2]=1  #southernmost cell gets a 1 by default
    if(nrow(cur.mat)==1)cur.y=cur.mat[1,1]
    else{
      flag=0
      j=1
      while(flag==0){
        j=j+1
        if(j==nrow(cur.mat)){
          cur.y=cur.mat[nrow(cur.mat),1]
          flag=1          
        }
        else{
          if(cur.mat[j,2]==0){
            cur.y=cur.mat[j-1,1]
            flag=1
          }
        }
      }
    }
    Point.vec[i,]=c(x[i],cur.y)
  }
  L=SpatialLines(list(Lines(list(Line(Point.vec)),ID="a")),proj4string=CRS(proj))
  L
}


#' Produce an RW1 adjacency matrix for a rectangular grid for use with areal spatial models (queens move)
#' @param x number of cells on horizontal side of grid
#' @param y number of cells on vertical side of grid
#' @param byrow If TRUE, cell indices are filled along rows (default is FALSE)
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn \email{paul.conn@@noaa.gov}
rect_adj <- function(x,y,byrow=FALSE){
  Ind=matrix(c(1:(x*y)),y,x,byrow=byrow)
  if(byrow==TRUE)Ind=t(Ind)
  n.row=nrow(Ind)
  n.col=ncol(Ind)
  Adj=matrix(0,x*y,x*y)
  for(i in 1:n.row){
    for(j in 1:n.col){
      if(i==1 & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i==1 & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i==1 & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i>1 & i<n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i>1 & i<n.row & j>1 & j<n.col){
        cur.nums=c(Ind[i,j]-n.row-1,Ind[i,j]-n.row,Ind[i,j]-n.row+1,Ind[i,j]-1,Ind[i,j]+1,Ind[i,j]+n.row-1,Ind[i,j]+n.row,Ind[i,j]+n.row+1)
        Adj[Ind[i,j],cur.nums]=1
      }
      if(i>1 & i<n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
        
      }
      if(i==n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
      }
      if(i==n.row & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
      if(i==n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
    }
  }
  if(byrow==TRUE)Adj=t(Adj)
  return(Adj)
}

#' Produce an RW2 Adjacency matrix for a rectangular grid for use with areal spatial models.
#' This formulation uses cofficients inspired by a thin plate spline, as described in Rue & Held, section 3.4.2
#' Here I'm outputting an adjacency matrix of 'neighbor weights' which makes Q construction for regular latices
#' easy to do when not trying to make inference about all cells (i.e., one can just cut off 
#' eliminate rows and columns associated with cells one isn't interested in and set Q=-Adj+Diag(sum(Adj)) 
#' @param x number of cells on horizontal side of grid
#' @param y number of cells on vertical side of grid
#' @param byrow If TRUE, cell indices are filled along rows (default is FALSE)
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn \email{paul.conn@@noaa.gov}
rect_adj_RW2 <- function(x,y,byrow=FALSE){
  cur.x=x+4  #make calculations on a larger grid and then cut off rows/columns at end
  cur.y=y+4
  Ind=matrix(c(1:(cur.x*cur.y)),cur.y,cur.x,byrow=byrow)
  if(byrow==TRUE)Ind=t(Ind)
  n.row=nrow(Ind)
  n.col=ncol(Ind)
  Adj=matrix(0,cur.x*cur.y,cur.x*cur.y)
  for(i in 3:(n.row-2)){
    for(j in 3:(n.col-2)){
        #kings move
        Adj[Ind[i,j],Ind[i,j]+1]=8
        Adj[Ind[i,j],Ind[i,j]+n.row]=8
        Adj[Ind[i,j],Ind[i,j]-n.row]=8
        Adj[Ind[i,j],Ind[i,j]-1]=8
        #bishops move        
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=-2
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=-2
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=-2
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=-2
        #kings move + 1
        Adj[Ind[i,j],Ind[i,j]+2]=-1
        Adj[Ind[i,j],Ind[i,j]+2*n.row]=-1  
        Adj[Ind[i,j],Ind[i,j]-2]=-1
        Adj[Ind[i,j],Ind[i,j]-2*n.row]=-1  
    }
  }
  #compile list of cells that need to be removed
  I.rem=matrix(0,n.row,n.col)
  I.rem[c(1,2,n.row-1,n.row),]=1
  I.rem[,c(1,2,n.col-1,n.col)]=1
  Adj=Adj[-which(I.rem==1),-which(I.rem==1)]
  if(byrow==TRUE)Adj=t(Adj)
  return(Adj)
}

