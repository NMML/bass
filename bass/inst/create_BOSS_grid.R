#script file for windoze

library(sp)
library(proj4)
library(rgdal)
library(rgeos)
library(raster)
library(nPacMaps)  #from Josh London


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

pep_ext <- extent(x_min, x_max, y_min, y_max)
sic_raster <- crop(sic_raster, pep_ext, snap = "near")
pep_ext <- extent(sic_raster)

plot(sic_raster, col = colorRampPalette(c("dark blue", "deepskyblue","skyblue", "lightskyblue", "white"))(20), main = "SSM/I Sea-Ice Concentration and PEP BOSS Extent\n01 April 2012")
plot(alaska_dcw, col = "black", add = TRUE)
plot(russia_dcw, col = "black", add = TRUE)
plot(pep_ext, col = "red", add = TRUE)

