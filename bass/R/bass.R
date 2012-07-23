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
				for(iind in 1:length(Int_id))Land.area[icell]=Land.area[icell]+gArea(gIntersection(Grid.poly[icell,],Land[[iobj]][Int_id[iind],]))
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
  Grid[["dist_ice_southern"]]=t(gDistance(Grid_points,southern_line,byid=TRUE))
  Grid
}
