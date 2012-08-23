#script file for windoze

library(sp)
#library(proj4)
library(rgdal)
library(rgeos)
library(raster)
library(nPacMaps)  #from Josh London
library(maptools)


#r <- raster("~/Dropbox/sea_ice/ssmi_ease/nt_20120401_f17_nrt_n.bin.reproj.tif")
r<-raster("//afsc/akc-nmml/Polar/Data/Environ/SeaIce/SSMI_SIC/2012/nt_20120401_f17_nrt_n.bin.reproj.tif")

#extent(r)<-c(-5922823,6109589,-5874646,6157766) 
#projection(r) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
		
fun <- function(x) {
	ifelse(x < 251, x/2.5, NA)
}
sic_raster <- calc(r, fun)

plot(sic_raster, main = "Northern Hemisphere SSM/I Sea-Ice Concentrationnn01 April 2012")



laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
		"+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
sic_raster <- projectRaster(sic_raster, res = c(25067.53, 25067.53),
		crs = laea_180_proj)
data(alaska_dcw)
data(russia_dcw)

alaska_dcw <- spTransform(alaska_dcw, CRS(laea_180_proj))
russia_dcw <- spTransform(russia_dcw, CRS(laea_180_proj))
plot(sic_raster, col = colorRampPalette(c("dark blue", "deepskyblue","skyblue", "lightskyblue", "white"))(20), main = "Northern Hempisphere SSM/I Sea-Ice Concentration\n01 April 2012")
plot(alaska_dcw, add = TRUE)
plot(russia_dcw, add = TRUE)

x_min <- 0
x_max <- 1400000
y_min <- -3800000
y_max <- -1600000

pep_ext <- extent(x_min, x_max, y_min, y_max)
sic_raster <- crop(sic_raster, pep_ext, snap = "near")
pep_ext <- extent(sic_raster)

plot(sic_raster, col = colorRampPalette(c("dark blue", "deepskyblue","skyblue", "lightskyblue", "white"))(20), main = "SSM/I Sea-Ice Concentration and PEP BOSS Extent\n01 April 2012")
plot(alaska_dcw, col = "black", add = TRUE)
plot(russia_dcw, col = "black", add = TRUE)
plot(pep_ext, col = "red", add = TRUE)

#define separate SpPolyDFs for mainland
Area_alaska=gArea(alaska_dcw,byid=TRUE)
Area_russia=gArea(alaska_dcw,byid=TRUE)
Alaska_mainland=alaska_dcw[which(Area_alaska==max(Area_alaska)),]
Russia_mainland=russia_dcw[which(Area_russia==max(Area_russia)),]

#Attach proportion land for each cell
Grid_poly<-rasterToPolygons(sic_raster,na.rm=FALSE) #convert to SpPolyDF for compatibility with rgeos
Land=list(alaska=alaska_dcw,russia=russia_dcw)
#the following takes awhile.  Instead, consider loading cur.Grid.Rdat
Grid_poly=add.prop.land(Grid=Grid_poly,Land=Land)
#remove initial layer
Grid_poly=Grid_poly[,-1]
save.image("cur.Grid.Rdat")
#load("cur.Grid.Rdat")

#Attach distance to mainland for each cell (distance from grid cell centroid to landscape polygon)
Grid_points=gCentroid(Grid_poly,byid=TRUE)
Dist_AK=gDistance(Grid_points,Alaska_mainland,byid=TRUE)
Dist_Rus=gDistance(Grid_points,Russia_mainland,byid=TRUE)
Dist_mainland=apply(cbind(as.vector(Dist_AK),as.vector(Dist_Rus)),1,'min')
Grid_poly[["dist_mainland"]]=Dist_mainland

#Attach distance to land (including islands) for each cell
Dist_AK=apply(gDistance(Grid_points,alaska_dcw,byid=TRUE),2,'min')
Dist_Rus=apply(gDistance(Grid_points,russia_dcw,byid=TRUE),2,'min')
Dist_land=apply(cbind(as.vector(Dist_AK),as.vector(Dist_Rus)),1,'min')
Grid_poly[["dist_land"]]=Dist_land
save.image("cur.Grid.Rdat")

#Attach distance to 10% Sea-ice Concentration Contour for each cell
sic_contour <- rasterToContour(sic_raster,levels=c(10))
Dist_contour <- gDistance(Grid_points,sic_contour,byid=TRUE)
Grid_poly[["dist_contour"]]=as.vector(Dist_contour)

#input and attach distance to "southern ice edge latitude" (determined by creating a line corresponding to southernmost ice edge and determining distance to this line)
#note to do this for real, we'd have to loop over filenames and attach values for each date
IceExtent=readOGR(dsn="c:/users/paul.conn/git/bass/ice",layer="nic_autoc2012136n_pl_a")
IceExtent=spTransform(IceExtent, CRS(laea_180_proj))
Grid_poly=add.dist.s.ice.edge(Grid=Grid_poly,Grid_points=Grid_points,IceExtent=IceExtent,proj=laea_180_proj)

#convert Grid_poly to SpatialPixelsDataFrame and then to a RasterStack
Grid_pix <- SpatialPixelsDataFrame(coordinates(Grid_poly),data=Grid_poly@data)
Grid_stack<-stack(Grid_pix)
