#script file to produce BOSS data for 2012 power analysis
#authors: Paul Conn & Josh London
#output: A list, "Data", comprised of the following objects:
# 1) "Grid" - this is a list of SpatialPolygonDataFrames, with one SPDF for each date of the survey.
#     Each polygon corresponds to a 25x25k grid cell. Presently, the following covariates are included in each data frame:
#     land_cover - proportion land covered by grid cell
#     dist_mainland - minimum distance from grid cell centroid to russia or alaska coast (excluding islands)
#     dist_land - minimum distance to any land (including islands)
#     dist_shelf - distance from 1000 meter depth contour
#     ice_conc - Percent of the grid cell covered by sea ice (from NSIDC)
#     dist_contour - distance from grid cell centroid to the closest 10% sea ice contour (a surrogate for distance to open water)
#     dist_ice_edge - distance from the latitude of the southernmost ice edge
#  The final three variables vary over time.
# 2) "Adj" - this is an adjacency matrix calculated with a queen's move for use with areal spatial models
# 3) "Meta" - this holds a few metadata entries
# 4) "Photos" - not currently implemented, but will eventually hold a spatial points dataframe giving locations,
#     and seal species counts, etc. needed for abundance estimation
library(sp)
#library(proj4)
library(rgdal)
library(rgeos)
library(raster)
library(nPacMaps)  #from Josh London
library(maptools)

DEBUG=TRUE
MEAN_ADJUST=TRUE  #if TRUE, standardizes "distance to" covariates to their means

if(DEBUG)source("c:/users/paul.conn/git/bass/bass/R/bass.R")

#tmp <- tempfile(fileext = ".tif",tmpdir="d:/temp")
#download.file(url=paste("file:////afsc/akc-nmml/Polar/Data/Environ/SeaIce/SSMI_SIC/2012/","nt_20120401_f17_nrt_n.bin.reproj.tif", sep = ""), destfile=tmp,mode="wb")
#r <- raster(tmp)
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
#y_max <- -2630000

pep_ext <- extent(x_min, x_max, y_min, y_max)
sic_raster <- crop(sic_raster, pep_ext, snap = "near")
dim_raster=dim(sic_raster)
pep_ext <- extent(sic_raster)

plot(sic_raster, col = colorRampPalette(c("dark blue", "deepskyblue","skyblue", "lightskyblue", "white"))(20), main = "SSM/I Sea-Ice Concentration and PEP BOSS Extent\n01 April 2012")
plot(alaska_dcw, col = "black", add = TRUE)
plot(russia_dcw, col = "black", add = TRUE)
plot(pep_ext, col = "red", add = TRUE)

#define separate SpPolyDFs for mainland
Area_alaska=gArea(alaska_dcw,byid=TRUE)
Area_russia=gArea(russia_dcw,byid=TRUE)
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
if(MEAN_ADJUST)Dist_mainland=Dist_mainland/mean(Dist_mainland)
Grid_poly[["dist_mainland"]]=Dist_mainland

#Attach distance to land (including islands) for each cell
Dist_AK=apply(gDistance(Grid_points,alaska_dcw,byid=TRUE),2,'min')
Dist_Rus=apply(gDistance(Grid_points,russia_dcw,byid=TRUE),2,'min')
Dist_land=apply(cbind(as.vector(Dist_AK),as.vector(Dist_Rus)),1,'min')
if(MEAN_ADJUST)Dist_land=Dist_land/mean(Dist_land)
Grid_poly[["dist_land"]]=Dist_land

#Attach distance to shelf break (1000m depth contour) for each cell
Shelf_break=readOGR("c:/users/paul.conn/git/bass/shapefiles/1000m_depth_contour.shp",layer="1000m_depth_contour")
Shelf_break=spTransform(Shelf_break, CRS(laea_180_proj))
Shelf_break=gBoundary(Shelf_break[1,])
#remove part of Shelf break that intersect with land or the EEZ line
Shelf_break=gDifference(Shelf_break,gBuffer(alaska_dcw,width=10000,byid=TRUE))
Shelf_break=gDifference(Shelf_break,gBuffer(russia_dcw,width=10000,byid=TRUE))
EEZs=readOGR("c:/users/paul.conn/git/bass/shapefiles/EEZ_Albers83_BS.shp",layer="EEZ_Albers83_BS")
EEZs=spTransform(EEZs, CRS(laea_180_proj))
EEZ_Alaska=EEZs[1,]  #limit to Alaska EEZ
Shelf_break=gDifference(Shelf_break,gBuffer(gBoundary(EEZ_Alaska),width=30000))
Dist_shelf=apply(gDistance(Grid_points,Shelf_break,byid=TRUE),2,'min')
if(MEAN_ADJUST)Dist_shelf=Dist_shelf/mean(Dist_shelf)
Grid_poly[["dist_shelf"]]=Dist_shelf

#output grid - this outputs spatial polygons dataframe for "large" grid with all time-constant
#covariates (e.g. distance to land).  We'll define time-dependent covariates on a smaller
#grid since this will have to be a list of spatial polygon DFs and will take up considerably more
#disk space
save(Grid_poly,file="BOSS_Grid_static.Rdat")

#create spatial connectivity matrix for CAR/ICAR modeling
Adj=rect_adj(x=dim_raster[2],y=dim_raster[1],byrow=TRUE)

#define indicator vector for which cells are to be included in analysis (exclude all cells north of xx latitude (Bering Straight), those that are 100% land, and those that have 100% water over course of study)
Include=(coordinates(Grid_points)[,2]<(-2630000)) #dont include cells above Bering Strait
Include[which(Grid_poly[["land_cover"]]>0.99)]=0  #dont include points with > 99% land
#scroll through sea ice data to determine which cells always occur in open water
#Date=c(paste("040",c(4:9),sep=''),paste("04",c(10:30),sep=''),paste("050",c(1:9),sep=''),paste("05",c(10:22),sep='')) #limit dates to 4/4-5/22 (first-last BOSS flights in 2012)
Date=c(paste("04",c(20:27),sep='')) #limit dates to 4/20-4/27 (first-last BOSS flights used in power analysis)
str1="//afsc/akc-nmml/Polar/Data/Environ/SeaIce/SSMI_SIC/2012/nt_2012"
str2="_f17_nrt_n.bin.reproj.tif"
Ice1=rep(0,length(Include))
for(idate in 1:length(Date)){
  cat(paste("\n date",idate,"of",length(Date),"\n"))
  filename=paste(str1,Date[idate],str2,sep='')
  r<-raster(filename)
  tmp_raster <- calc(r, fun)
  tmp_raster <- projectRaster(tmp_raster, res = c(25067.53, 25067.53),
                            crs = laea_180_proj)
  tmp_raster <- crop(tmp_raster, pep_ext, snap = "near")
  tmp_raster[is.na(tmp_raster)]=0
  Ice1=Ice1+values(tmp_raster>0)
}
Include[which(Ice1==0)]=0
#import shapefile specifying russian/alaskan EEZs
EEZs=readOGR("c:/users/paul.conn/git/bass/shapefiles/EEZ_Albers83_BS.shp",layer="EEZ_Albers83_BS")
EEZs=spTransform(EEZs, CRS(laea_180_proj))
EEZ_Alaska=EEZs[1,]  #limit to Alaska EEZ
I.Alaska=gIntersects(Grid_poly,EEZ_Alaska,byid=TRUE)
Include[which(I.Alaska==0)]=0

#reduce "grid" to those cells of interest for analysis, reformatting spatial connectivity matrix and declaring list of spatial polygon dataframe objects
Grid_reduced=Grid_poly[Include==1,]
Adj_reduced=Adj[Include==1,Include==1]
plot(Grid_reduced,add=TRUE)
Data=list(Adj=Adj_reduced)
Data$Grid=vector("list",length(Date))  #allocate space for one SpatialPolygonsDataframe for each date in study
for(idate in 1:length(Date))Data$Grid[[idate]]=Grid_reduced

#attach daily sea ice concentration & distance to 10% sea-ice concentration contour for each cell
for(idate in 1:length(Date)){
  cat(paste("\n date",idate,"of",length(Date),"\n"))
  filename=paste(str1,Date[idate],str2,sep='')
  r<-raster(filename)
  tmp_raster <- calc(r, fun)
  tmp_raster <- projectRaster(tmp_raster, res = c(25067.53, 25067.53),
                              crs = laea_180_proj)
  tmp_raster <- crop(tmp_raster, pep_ext, snap = "near")
  sic_contour=rasterToContour(tmp_raster,levels=c(10))
  dist_contour=as.vector(gDistance(Grid_points[Include==1,],sic_contour,byid=TRUE))
  ice_conc=values(tmp_raster)[Include==1]
  if(MEAN_ADJUST)dist_contour=dist_contour/mean(dist_contour)
  Data$Grid[[idate]]=spCbind(spCbind(Data$Grid[[idate]],ice_conc),dist_contour)
}

#input and attach distance to "southern ice edge latitude" (determined by creating a line corresponding to southernmost ice edge and determining distance to this line)
#note to do this for real, we'd have to loop over filenames and attach values for each date
str1="nic_autoc2012"
str2="n_pl_a"
Date=c(111:118)  #April 20-27 in Julian
for(idate in 1:length(Date)){
  filename=paste(str1,Date[idate],str2,sep='')
  IceExtent=readOGR(dsn="x:/polar/data/environ/seaice/nic/2012",layer=filename)
  IceExtent=spTransform(IceExtent, CRS(laea_180_proj))
  Data$Grid[[idate]]=add.dist.s.ice.edge(Grid=Data$Grid[[idate]],Grid_points=Grid_points[Include==1,],IceExtent=IceExtent,proj=laea_180_proj,mean_adjust=MEAN_ADJUST)
}  

#add some metadata
Data$Meta=list(date.start="20Apr2012",date.end="27Apr2012")
Data$Meta$info="2012 PEP BOSS survey data for power analysis"

#save dataset
save(Data,file="Power_Data_noPhotos.rdat")

#output grid data from one year as a shapefile for Mike & other Arc users
writeOGR(Grid_reduced,dsn=".",layer="Grid2012",driver="ESRI Shapefile")

#convert Grid_poly to SpatialPixelsDataFrame and then to a RasterStack ... this might be a better object class as the foundation for remaining analysis
Grid_pix <- SpatialPixelsDataFrame(coordinates(Grid_poly),data=Grid_poly@data)
Grid_stack<-stack(Grid_pix)


