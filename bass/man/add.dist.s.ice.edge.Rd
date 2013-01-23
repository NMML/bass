\name{add.dist.s.ice.edge}
\alias{add.dist.s.ice.edge}
\title{Add line to SpatialPolygonsDataFrame that gives distance from "southern ice edge" for the centroid of each grid cell.
This version function is depracated - see function get.edge.Bering for the current version}
\usage{
  add.dist.s.ice.edge(Grid, Grid_points, IceExtent, proj,
    mean_adjust = TRUE)
}
\arguments{
  \item{Grid}{master SpatialPolygonsDataFrame object to add
  to}

  \item{Grid_points}{centroids of Grid (SpatialPoints
  object)}

  \item{IceExtent}{a SpatialPolygonsDataFrames giving ice
  extent (e.g. as imported from NIC shapefile)}

  \item{proj}{current projection}

  \item{mean_adjust}{if TRUE, standardize distance
  covariate by dividing by its mean}
}
\value{
  SpatialPolygonsDataFrame object, which includes distance
  from southern ice edge as an additional column
}
\description{
  Add line to SpatialPolygonsDataFrame that gives distance
  from "southern ice edge" for the centroid of each grid
  cell. This version function is depracated - see function
  get.edge.Bering for the current version
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}

