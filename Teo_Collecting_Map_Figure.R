library(geosphere)
library(raster)
library(ggmap)
library(rgdal)
library(dismo)

range01=function(x){
newnums <- sapply(1:length(x), function(z) ifelse(!is.na(x[z]), (x[z]-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)), NA ))
return(newnums)
}
opar <- par()

#get bioclim data
sites <- read.csv("~/Dropbox/AO-1/population locations.csv",header=T)
#sort site file in meaningful way
sites$popabbr <- c("TZ","FP","MT","DA","MC","M","TC","TX","CL","CU")
sitesO <- sites[order(sites$popabbr),]
Ptpops <- data.frame(lon=sitesO$lon,lat=sitesO$lat)
SPtpops <- SpatialPoints(Ptpops,proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
#read in ENV data
files <- list.files(path ="~/bio_22_tif/", pattern=".tif", full.names = T) #Redistributing worldclim files is not a permitted use, but the lead author (Anna O'Brien) will likely retain the downloaded files for some time, please reach out if there are questions.
bioclimvals <- matrix(nrow=nrow(Ptpops),ncol=length(files))
for(i in 1:length(files)){
	bioclimvals[,i]=extract(raster(files[i]),SPtpops) #bioclim crs is essentially crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
}
popbiodat <- data.frame(bioclimvals) #sorted by how are in sitesO, so alphabetical.
#columns are what the bioclim levels are(1 through 19): TAnn TmeanWarmQ TmeanCQ Pann PwetM PDM Pseas PwetQ PDQ PWarmQ PCQ TDiRange TIsother  TSD TmaxWarmM TminCM TannR TmeanWetQ TmeanDQ

MxCity <- SpatialPoints(data.frame(-99.1333,19.4328),proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

mxalt <- crop(raster("~/maketeofigure/alt_22.tif"),extent(-112,-97,16,26)) #check if projection crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
mxsub <- crop(mxalt,extent(-100,-98.5,18.65,19.75))

bw <- colorRampPalette(c(rgb(0,0,0),rgb(1,1,1)))

par(opar)
pdf("~/Dropbox/AO-1/PopLocs.pdf", width=6,height=4.9)
	plot(mxalt,col=bw(1000))
	mtext("Latitude",side=2,line=2)
	mtext("Longitude",side=1,line=2)
	points(MxCity, cex=2,pch='*',col=rgb(0,0,0))
	points(SPtpops, col = rgb(range01(popbiodat[,1]),0,1-(range01(popbiodat[,1]))), pch = 16,cex=0.5)
	par(new = TRUE)
	par(fig = c(0.16, 0.55, 0.25, 0.65)) 
	plot(mxsub,col=bw(1000),legend=FALSE,zlim=c(0,5000))
	par(fig = c(0.16, 0.55, 0.25, 0.65))
	points(MxCity, cex=3,pch='*',col=rgb(0,0,0))
	par(fig = c(0.16, 0.55, 0.25, 0.65))
	points(SPtpops, col = rgb(range01(popbiodat[,1]),0,1-(range01(popbiodat[,1]))), pch = 16)
	points(SPtpops[c(5,8:10),], col = rgb(range01(popbiodat[,1])[c(5,8:10)],0,1-range01(popbiodat[,1])[c(5,8:10)]), pch = 1, cex=2)
dev.off()

sitetann <- sort( popbiodat[,1] / 10)
rb<-colorRampPalette(c(rgb(1,0,0),rgb(0,0,1)))
pdf("~/Dropbox/AO-1/baronly.pdf",width=5.5,height=1)
par(mar=c(3.5,1,0,1))
image(as.matrix(sitetann), col=rev(rb(100)),zlim=c(min(sitetann),max(sitetann)),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side=1,at=seq(from=0, to=1,length.out=10),labels=sitetann)
	mtext("Site mean annual temperature", side=1,line=2,adj=0.5)
	text(seq(from=0, to=1,length.out=10)[c(2,8:10)],rep(0.125,times=4),rep("*",times=4),col=rgb(1,1,1),cex=3)
dev.off()
#combined in keynote, which is simpler

