#' Tools and Analysis Functions for the Bering-Arctic Seal Surveys
#'
#' bass is a central package for creating basic/global functions and data
#' related to the Bering-Arctic Seal Surveys. The 2012 BOSS surveys are the
#' first surveys to rely on this package
#'
#' @docType package
#' @name bass
#' @aliases bass bass-package
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
#' @author Paul Conn (paul.conn@noaa.gov)
#' @return SpatialPolygonsDataFrame object, which includes proportion land cover as an additional column
#' @examples
#' New_grid <- add.prop.land(Grid,Land,rel=TRUE)
#' New_grid
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

#' Add line to SpatialPolygonsDataFrame that gives distance from "southern ice edge" for the centroid of each grid cell
#'
#' @param Grid  master SpatialPolygonsDataFrame object to add to
#' @param Grid_points centroids of Grid (SpatialPoints object)
#' @param IceExtent  a SpatialPolygonsDataFrames giving ice extent (e.g. as imported from NIC shapefile)
#' @param proj  current projection
#' @author Paul Conn (paul.conn@noaa.gov)
#' @return SpatialPolygonsDataFrame object, which includes distance from southern ice edge as an additional column
#' @examples
#' New_grid <- add.dist.s.ice.edge(Grid,IceExtent)
#' New_grid
#' @export
add.dist.s.ice.edge<-function(Grid,Grid_points,IceExtent,proj){
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
  dist_ice_edge=t(gDistance(Grid_points,southern_line,byid=TRUE))
  Grid=spCbind(Grid,dist_ice_edge)
  Grid
}

#' Produce an adjacency matrix for a rectangular grid for use with areal spatial models (queens move)
#' @param x number of cells on horizontal side of grid
#' @param y number of cells on vertical side of grid
#' @param byrow If TRUE, cell indices are filled along rows (default is FALSE)
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn
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
