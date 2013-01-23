\name{add.prop.land}
\alias{add.prop.land}
\title{Add a column to the SpatialPolygonsDataFrame object giving the proportion of area
of each cell that is covered by land}
\usage{
  add.prop.land(Grid, Land, rel = TRUE)
}
\arguments{
  \item{Grid}{master SpatialPolygonsDataFrame object for
  which land calculations are desired (separate
  calculations occur for each subgeometry)}

  \item{Land}{a list of SpatialPolygonsDataFrames that
  represent land cover}

  \item{rel}{if TRUE (default), returns proportional area;
  if FALSE, returns absolute area}
}
\value{
  SpatialPolygonsDataFrame object, which includes
  proportion land cover as an additional column
}
\description{
  This function returns an object of type
  \code{SpatialPolygonsDataFrame} that takes an original
  SpatialPolygonsDataFrame and computes the proportion of
  area for each subgeometry that is covered by land.
  "Land" must be represented as a list of
  SpatialPolygonsDataFrames.
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}

