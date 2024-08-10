#######################
###Libraries and functions
####################
###

library(MCMCglmm) #for linear models
library(coda) #for HPDinterval, dependency

#a remove missing data function:
getfull <- function(dat){
	whichfull <- which(sapply(1:nrow(dat), function(z) any(is.na(dat[z,]) ) )==FALSE)
	return(whichfull)
}

#add a shortcut for coda HPDinterval 
HPDi <- function(vect,prob){
	int <- HPDinterval(as.mcmc(vect),prob=prob)
	return(int)
}

bufferX <- function(x,p) { #provide a buffer around a range, e.g. for axes limits when plotting
	r<- range(x,na.rm=T)
	add <- c(-1,1)*p*(r[2]-r[1])
	return(r+add)
	}	

range01=function(x){ # rescale a vector to fall within 0 and 1, useful for plotting with color assignments based on a variable
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

std.error <- function(dat, na.rm=TRUE) {sd(dat,na.rm=na.rm)/sqrt(length(dat))}
# basic standard error function, defaults to na.rm=T

setwd("~/Variable_teosinte-biota-interactions/") 

##Growth Curves function - second function in this section, relies on first (getcoeffs) to extract coefficients
#for each individual, makes a growth curve between days alive and height,
#then pulls coefficients from that curve and stores them in a dataframe, rows accding to plant number.
getcoeffs <- function(x1, y1,type,force0=F){
	y <- y1[!is.na(y1)]
	x <- x1[!is.na(y1)]
	nobs<-length(y)
	if(type == "parabola"){
		if(nobs>=3){
			if(force0==T ){ 
		 	coeffs<-as.vector(lm(y~ 0 + x + I(x^2))$coefficients)		
			} else if(force0==F){
		 	coeffs<-as.vector(lm(y~ x + I(x^2))$coefficients)
			} else {print("force0 must be T or F")}
		} else if(nobs<3 & force0==T) {
			coeffs <- as.numeric(c(NA,NA))
		} else if(nobs<3 & force0==F) {
		 	coeffs <- as.numeric(c(NA,NA,NA))
		} else {print("force0 must be T or F; or can't work with the y vector")}		
	} else if(type == "line"){
		if(nobs>=2){		
			if(force0==T){
			coeffs<-as.vector(lm(y~ 0 + x)$coefficients)
			} else if(force0==F){
		 	coeffs<-as.vector(lm(y~ x)$coefficients)
			} else {print("force0 must be T or F")}
		} else if(nobs<2 & force0==T) {
			coeffs <- as.numeric(c(NA))
		} else if(nobs<2 & force0==F) {
		 	coeffs <- as.numeric(c(NA,NA))
		} else {print("force0 must be T or F; or can't work with the y vector")}	
	} else {print("type options are line or parabola")}
}
GrowthCurves<-function(size,timepoints,type,force0=F){	
	Growth.Curves <- sapply(1:nrow(size), function(z) getcoeffs(as.numeric(t(timepoints[z,])),as.numeric(t(size[z,])),type=type,force0=force0))			
	return(Growth.Curves)}



#######################
###Read in data, numbers reported in methods
####################
###

#data on collecting sites
PopEnvDat <- read.csv("PopEnvDat.csv",header=T) 
#bioclim data originally extracted from worldclim. bioclim rasters are not authorized for redistribution, extracted numbers are included
#the data also include columns for data provided by soil testing facility.
#below is the same data, but sampled (copied rows) according to the experimental design for plants OR inocula sources
PlantSourceEnvDat <- read.csv("Plant_EnvDat.csv",header=T,nrows=720)# pop env dat sampled for plant ID
InocSourceEnvDat <- read.csv( "Inoc_EnvDat.csv",header=T,nrows=720)# same but for biota

#treatment design for experimental plants in greenhouse
treatments <- read.csv("treatments.csv",header=T,nrows=720)
colnames(treatments)[1]<-"Individual"

#measured data from experimental plants
final_measurements <- read.csv("final_measurements.csv",header=T,nrows=720)
#leaf data in 3 matching matrices. Leaves were measured at 2 timepoints. "L" is length "W" is width, and "whichlf" documents which leaf was measured (1 = 1st true leaf, 2 = 2nd true leaf, etc)
lfsizeL_points <- read.csv("lfsizeL_points.csv",header=T,nrows=720)#only first measure in leaf dataframes is always the same (2nd) leaf. Not all plants made a 4th leaf
lfsizeW_points <- read.csv("lfsizeW_points.csv",header=T,nrows=720)
lfsize_whichlf <- read.csv("lfsize_whichlf.csv",header=T) 
#height data. "days" is age of the plant, "points" is matching matrix giving the height value
height_points <- read.csv("heightpoints.csv",header=T)
height_days <- read.csv("height_days.csv",header=T)
#ionomics data
ao1ionomics <- read.csv("ionomics.csv",header=T,nrows=720)

#harvest timing variation reported in methods
sum(final_measurements$harvestage>= 48 & final_measurements$harvestage <= 52,na.rm=T)
#654
sum(is.na(final_measurements$harvestage))
#52
654/(720-52) # 
# %
#reporting co-correlated environmental variables, methods.
cor(PopEnvDat$TAnn, PopEnvDat[,c(4,21,22,27,33)])

#######################
###Setting up vectors and dataframes
####################
###

#generating growth rate curves, Figure S2
Height.Curves <- t(GrowthCurves(height_points,height_days,type="parabola",force0=F))
#height curves plot for supplement
ispos <- as.numeric(Height.Curves[,3] > 0)
#average age at harvest for plants harvested out of sync. add 2 days for time between planting and first germinant
mean(c(height_days[which(height_days[,6]>52),6], height_days[which(height_days[,6]<49),6]))
range(c(height_days[which(height_days[,6]>52),6], height_days[which(height_days[,6]<49),6]))
#reported in MS.
#Figure S3
pdf("HeightCurves.pdf",width=3,height=6)
par(mfrow=c(2,1))
par(mar=c(2,4,3,1))
par(oma=c(2,0,0,0))
plot(1~c(1), pch = NA, ylim=range(height_points,na.rm=T),xlim=range(height_days,na.rm=T),ylab="",xlab="")
	invisible(sapply((1:nrow(height_points))[ispos==1], function(z) lines(unlist(height_points[z,])~unlist(height_days[z,]), col=rgb(0,0,0,alpha=0.25)) ))
	mtext("Height (cm)",side=2,line=2.5,adj=-0.9)
	mtext("Late growth pattern",side=3,line=0.5)
plot(1~c(1), pch = NA, ylim=range(height_points,na.rm=T),xlim=range(height_days,na.rm=T),ylab="",xlab="")
	invisible(sapply((1:nrow(height_points))[ispos==0], function(z) lines(unlist(height_points[z,])~unlist(height_days[z,]), col=rgb(0,0,0,alpha=0.25)) ))
	mtext("Days since germination",side=1,line=2.5)
	mtext("Early growth pattern",side=3,line=0.5)
dev.off()

#useful vectors describing treatments
SAN <- ifelse(treatments$soil.full != "N",as.character(treatments$SA),NA)#sympatric, allopatric, or no inocula
isinoc <- ifelse(treatments$soil.full=="N",0,1)
tann <- sort( PopEnvDat$TAnn / 10)

###Tissue elements data. Includes figure S4
#checking for normality, taking the log where this improves normality, and expanding to match final measurement data dimensions with NAs where samples were lost during ionomics analysis
subions1 <- ao1ionomics[ao1ionomics$sampleid%in%(1:720) , 14:33] #remove information columns, just ions.
subionsid <- ao1ionomics$sampleid[ao1ionomics$sampleid%in%(1:720)]
subions <- subions1 #duplication to retain subions1 for later, has the untransformed data
for(i in 1:20){if(any(subions1[,i]<0,na.rm=T)){subions[,i] <- subions[,i] + abs(min(subions[,i],na.rm=T)) + min(subions[subions[,i]>0,i],na.rm=T)  }} 
#using the minimum transformed dataframes to check for normality
logimprovesionA <- sapply( 1:ncol(subions) , function(colmn) shapiro.test(subions[,colmn])$statistic) < sapply( 1:ncol(subions) , function(colmn) shapiro.test(log(subions[,colmn]))$statistic)
wstatsA <- cbind(sapply( 1:ncol(subions) , function(colmn) shapiro.test(subions[,colmn])$statistic), sapply( 1:ncol(subions) , function(colmn) shapiro.test(log(subions[,colmn]))$statistic)) #minimum are 0.49 As and 0.31 Se, rest are above 0.85, so pretty normal or log-normal
#most are log normal, As is closer to normal than log-normal
#add back sample id information
ao1subions <- cbind(data.frame(sampleid=subionsid),subions) #add to unscaled data
ao1mixions <- cbind(data.frame(sampleid=subionsid),log(subions)) #logged element concentrations
ao1mixions[,which(!logimprovesionA)+1] <- ao1subions[,which(!logimprovesionA)+1] #put back unlogged when log does not make more normal; plus 1 for samp id column
ionnames <- c("B","Na","Mg","Al","P","S","K","Ca","Fe","Mn","Co","Ni","Cu","Zn","As","Se","Rb","Sr","Mo","Cd")
#not all 720 plants have samples, and this data is not perfectly sorted, so need to fill in new matrices with sorted data.
subionsfull <- as.data.frame(matrix(NA,ncol=ncol(subions1) + 1,nrow=720)) #not logged, for means, SEs and raw data supplemental figs.
colnames(subionsfull) <- colnames(ao1subions)
for(i in 1:720){ 
		if( length(which(ao1mixions$sampleid == i))!= 0 ){
			subionsfull[i,2:ncol(subionsfull)] <- subions1[which(ao1mixions$sampleid == i),] 
			subionsfull$sampleid[i] <- i
		} else{subionsfull$sampleid[i] <- i}
}
ionsM <- as.data.frame(matrix(NA,ncol=ncol(ao1mixions),nrow=720)) #logged (only when makes more normal) 
colnames(ionsM) <- colnames(ao1mixions)
for(i in 1:720){ 
		if( length(which(ao1mixions$sampleid == i))!= 0 ){
			ionsM[i,] <- ao1mixions[which(ao1mixions$sampleid == i),] 
		} else{ionsM$sampleid[i] <- i}
}
#Figure S4: normality tests, histograms, for supplement
pdf("TissueElementNormality.pdf",height=10,width=8)
par(mfrow=(c(5,4)))
par(mar=(c(7,4,3,1)))
for(i in 1:20){
	hist(subions1[,i],breaks=20,main=ionnames[i],xlab="") #colnames(ao1ionomics[,14:33])[i]
	mtext(paste("W=",round(wstatsA[i,1],digits = 2),ifelse(!logimprovesionA[i],"*",""), sep=""),	 side=1, line = 3)
	mtext(paste("ln W=",round(wstatsA[i,2],digits = 2),ifelse(logimprovesionA[i],"*",""), sep=""),	 side=1, line = 5)
}
dev.off()#

###combining plant trait data, includes Supplemental Figure S3
##extracts variables, puts in single dataframe, logs as appropriate for data normality
#get second leaf measurements, if the second leaf is measured twice, use second measure.
ll2 <- ifelse(lfsize_whichlf[1:720,2]==2,lfsizeL_points[,2],lfsizeL_points[,1]) # the first measurement is always the second leaf; the second measurement is sometimes also the first leaf
lw2 <- ifelse(lfsize_whichlf[1:720,2]==2,lfsizeW_points[,2],lfsizeW_points[,1])
#plant traits across all plants
Plant.Traits <- data.frame(biomass = final_measurements$biomass, SR = final_measurements$shoot.root,
							Ht = final_measurements$final.height, Wt = final_measurements$final.stem, 
							Lfnum = final_measurements$final.lvs, LL = ll2,
							LW =lw2, GrthDelay = Height.Curves[,3],
							germday = final_measurements$germday, hairs = final_measurements$hairs,
							Lgermday = log(final_measurements$germday), Lhairs = log(final_measurements$hairs+1,) )#log germday, and log hairs if you want them normally distr. otherwise ok.
p1vec4hairsGD <- c(0,0,0,0,0,0,0, abs(min(Plant.Traits$GrthDelay,na.rm=T)) + min(Plant.Traits$GrthDelay[Plant.Traits$GrthDelay>0],na.rm=T) ,0,1) #adding minimum values in order to take log and not lose datapoints falling at or below 0
wstatsplanttrts <- cbind(sapply( c(1:10) , function(colmn) shapiro.test(Plant.Traits[,colmn])$statistic),sapply( c(1:10) , function(colmn) shapiro.test(log(Plant.Traits[,colmn] + p1vec4hairsGD[colmn]))$statistic))
# For hairs and germination day is it logged much better
	#other plant traits generally normal, sometimes log is only slightly higher (stem width and S:R) but then both are normal,
traitnames <- c("Biomass","S:R","Height","Width","Leaf#","LeafL","LeafW","GrwTime","GermDay","Hairs")
traitnamesNL <- c("Biomass","Shoot:Root","Stem Height","Stem Width","Leaf Number",
				"Leaf Length","Leaf Width","Growth Timing","Germination Day","Stem Hairs")
logORnoTrati <- c(0,0,0,0,0,0,0,0,1,1 )#rule for traits is that only log transform if originally W <.9 AND log improves normality
pdf("BiomassandTraitNormality.pdf",height=6,width=8)
par(mfrow=(c(3,4)))
par(mar=(c(7,4,3,1)))
for(i in 1:10){
	hist(Plant.Traits[,i],breaks=20,main=traitnamesNL[i],xlab="") #colnames(ao1ionomics[,14:33])[i]
	mtext(paste("W=",round(wstatsplanttrts[i,1],digits = 2),ifelse(logORnoTrati[i]==0,"*",""), sep=""),	 side=1, line = 3)
	mtext(paste("ln W=",round(wstatsplanttrts[i,2],digits = 2),ifelse(logORnoTrati[i]==1,"*",""), sep=""),	 side=1, line = 5)
}
dev.off()#


#split into inoculated and uninoculated data
#for some analyses, we use inoculated only
Plant.Traits.live <- Plant.Traits[!is.na(SAN),]
Plant.Traits.UNinoc <- Plant.Traits[is.na(SAN),]
inoc_trts <- treatments[which(treatments$soil.abv != "N"),]
ioninoc <- ionsM[treatments$soil.abv !="N",]
ionUNinoc <- ionsM[treatments$soil.abv =="N",]
inoc_Penv <- PlantSourceEnvDat[which(treatments$soil.abv != "N"),]
inoc_Ienv <- InocSourceEnvDat[which(treatments$soil.abv != "N"),]
inoc_Ienv_I <- inoc_Ienv
colnames(inoc_Ienv_I) <- paste(colnames(inoc_Ienv_I),"_I",sep="")
inoc_Penv_P <- inoc_Penv
colnames(inoc_Penv_P) <- paste(colnames(inoc_Penv_P),"_P",sep="")

shapiro.test(PopEnvDat$InorganicN.ppm); shapiro.test(log(PopEnvDat$InorganicN.ppm))
shapiro.test(PopEnvDat$P.Bray.ppm); shapiro.test(log(PopEnvDat$P.Bray.ppm))
shapiro.test(PopEnvDat$K.ppm); shapiro.test(log(PopEnvDat$K.ppm))
#in Shapiro test, P.Bray.ppm and K.ppm in from the soil data are less than 0.9 and log increases W;. so these are logged below

#get a few combined treatment and data dataframes. 
TraitIon1 <- data.frame(cbind(ioninoc[,2:21], 
	Plant.Traits.live, inoc_Ienv_I[,c(1,4)],inoc_Penv_P[,c(1,4)]), #Tann and Pann
	lPp = log(inoc_Penv$P.Bray.ppm) ,  lPi = log(inoc_Ienv$P.Bray.ppm), lKp = log(inoc_Penv$K.ppm), lKi =  log(inoc_Ienv$K.ppm),
	mom=factor(inoc_trts$Mom),SA= as.numeric(as.factor(SAN[!is.na(SAN)]) ) ) 
TraitIonInocCom <- TraitIon1[getfull(TraitIon1[,c(1:20,22:28,31,32)]),] #complete cases for element and trait variables (other columns are environmental data columns)
TI_trts <- TraitIonInocCom
TI_trts$Plant <- inoc_trts$Plant[getfull(TraitIon1[,c(1:20,22:28,31,32)])]#get treatment data for the complete cases of inoculated plants 
TI_trts$Soil <- inoc_trts$soil.full[getfull(TraitIon1[,c(1:20,22:28,31,32)])]#get treatment data for the complete cases of inoculated plants 
TI_trts$SYM <- TI_trts$SA -1

#######################
###Main effect of inoculation
####################
###

###Main effect of inoculation
##biomass, P, RDA, responses to inoculation; using logged but not scaled elements for tests
phosdat <- data.frame(y=ionsM[,5+1],x=isinoc)
phosdat <- phosdat[getfull(phosdat),]
phosmod <- MCMCglmm(y~x,data=phosdat,verbose=F)#defaul nitt, burnin, thin. 13000, 3000, 10.
biodat <- data.frame(y=Plant.Traits[,1],x=isinoc)
biodat <- biodat[getfull(biodat),]
biomod <- MCMCglmm(y~x,data=biodat,verbose=F)#default nitt, burnin, thin. 13000, 3000, 10.; 

#for reporting in text: Figure S5
#Means and SEs
bioUImn <- tapply(Plant.Traits$biomass,isinoc ,mean, na.rm=T)# 3.50 2.59 with SEs 0.036 and 0.094 for inoc and uninoc, respectively
bioUIse <- tapply(Plant.Traits$biomass,isinoc ,std.error )
pUImn <- tapply(subionsfull$P31, isinoc,mean, na.rm=T) # 977 and 569, with SE 11.8 and 12.1
pUIse <- tapply(subionsfull$P31,isinoc ,std.error )
#vectors & getting intervals, T.s used also in other sections
T.s <- seq(from= min(TI_trts$TAnn_P), to=max(TI_trts$TAnn_P),length.out=1000)/10
pdf("InoculationEffects.pdf",height=4,width=5)
layout(matrix(1:2,ncol=2,byrow=F))
par(oma=c(7,0,2,2))
par(mar=c(0,4,0,0))
plot(bioUImn~c(1,2),ylim=bufferX(c(bioUImn+bioUIse,bioUImn-bioUIse),0.1),pch=16,xlim=c(0.5,2.5), xaxt="n",ylab="")
	mtext("Biomass g",side=2,line=2)
	arrows(c(1,2),y0=bioUImn-bioUIse,y1=bioUImn+bioUIse,length=0)
	axis(side=1,at=c(1,2),las=2,labels=c("Uninoculated","Inoculated"))
plot(pUImn~c(1,2),ylim=bufferX(c(pUImn+pUIse,pUImn-pUIse),0.1),pch=16,xaxt="n",xlim=c(0.5,2.5),ylab="")
	mtext(expression("Phosphorus"*mu*"/mg"),side=2,line=2)
	arrows(c(1,2),y0=pUImn-pUIse,y1=pUImn+pUIse,length=0)	
	axis(side=1,at=c(1,2),las=2,labels=c("Uninoculated","Inoculated"))
dev.off()	

#####################
####Sources of variation
### 
####################


###sums of squares
##treatment variables
plantsite <- as.character(TI_trts$TAnn_P)
biotasite <- as.character(TI_trts$TAnn_I)
plantandbiota <- paste(plantsite,biotasite)

ssbyvar <- function(response,category.vec){ #sums of squares function
	means <- tapply(response,category.vec,mean,na.rm=T) #take the means by category
	ssresid <- sum(sapply(sort(unique(category.vec)), function(z) sum( (response[category.vec==z] - means[names(means)==z])^2,na.rm=T ))) #square of difference of each datapoint from its associated treatment mean (residual variation)
	sstot <- sum((response-mean(response,na.rm=T))^2,na.rm=T) #square of difference of each datapoint from the grand mean (total variation)
	sst <- (sstot-ssresid) # total variation - residual variation = treatment variation
	return(sst/sstot) # treatment variance as a fraction of total variation
	}

limcols <- c(21:28,31,32,1:20)
ssP <- sapply(limcols, function(z) ssbyvar(TI_trts[,z],plantsite))
ssB <- sapply(limcols, function(z) ssbyvar(TI_trts[,z],biotasite))
ssPB <- sapply(limcols, function(z) ssbyvar(TI_trts[,z],plantandbiota))
ssM <- sapply(limcols, function(z) ssbyvar(TI_trts[,z],TI_trts$mom))
amntVars <- cbind(ssP,ssB,ssPB)
##testing sig of var explained with anova
trait_sourcevar_Asig <- matrix(NA,ncol=30,nrow=3)
 datatmp <- TI_trts[,c(21:28,31,32,1:20,33:ncol(TI_trts))]
for(i in c(1:30)){
 datatmp$y <- datatmp[,i]
trait_sourcevar_Asig[1,i] <- anova(lm(y~Plant,data=datatmp))$"Pr(>F)"[1]
trait_sourcevar_Asig[2,i] <- anova(lm(y~Soil,data=datatmp))$"Pr(>F)"[1]
trait_sourcevar_Asig[3,i] <- anova(lm(y~Soil+Plant,data=datatmp),lm(y~Plant*Soil,data=datatmp))$"Pr(>F)"[2]
}
colnames(trait_sourcevar_Asig) <- colnames(TI_trts[,c(21:28,31,32,1:20)])
ssPsig <- trait_sourcevar_Asig[1,]<0.05
ssBsig <- trait_sourcevar_Asig[2,]<0.05
ssPBsig <- trait_sourcevar_Asig[3,]<0.05


###Randomization of the above to create null distribution
################may take upwards of 5 hrs#############
set.seed(1) #to produce identical objects as in manuscript
sampAsig <- array(NA, dim=c(length(limcols),3,1000))#1000
sampAamnt <- array(NA, dim=c(length(limcols),3,1000))#1000
toSS <- list(plantsite, biotasite, plantandbiota)
for(j in 1:dim(sampAsig)[3]){
	sampdat_datonly <- TI_trts[sample(1:nrow(TI_trts),replace=F),limcols]
	sampAamnt[,,j] <- sapply(1:length(toSS), function(SofV) sapply(1:ncol(sampdat_datonly), function(z) ssbyvar(sampdat_datonly[,z],toSS[[SofV]])) )
	sampdat <- as.data.frame(cbind(sampdat_datonly, TI_trts[,33:ncol(TI_trts)]))
	for(i in 1:ncol(sampdat_datonly)){
		sampdat$yvar <- sampdat[,i]
		sampAsig[i,1,j] <- anova(lm(yvar~Plant,data=sampdat))$"Pr(>F)"[1]
		sampAsig[i,2,j] <- anova(lm(yvar~Soil,data=sampdat))$"Pr(>F)"[1]
		sampAsig[i,3,j] <- anova(lm(yvar~Soil+Plant,data=sampdat),lm(yvar~Plant*Soil,data=sampdat))$"Pr(>F)"[2]
	}
}
save(sampAamnt,file="sampAamnt.Rdata")
save(sampAsig,file="sampAsig.Rdata")
load("sampAamnt.Rdata")
amtintervals <- sapply(1:nrow(amntVars), function(z) sapply(1:ncol(amntVars), function(p) findInterval(amntVars[z,p],sort(sampAamnt[z,p,]))/(dim(sampAamnt)[3])   )  )
isamntVsig <- amtintervals > 0.95 #amount of variance explained by factor in real data higher than in 95% of randomizations
##alternative evaluation
# load("sampAsig.Rdata")
# sigvarintervals <- sapply(1:ncol(trait_sourcevar_Asig), function(z) sapply(1:nrow(trait_sourcevar_Asig), function(p) findInterval(trait_sourcevar_Asig[p,z],sort(sampAsig[z,p,]))/(dim(sampAsig)[3])   )  )
# issigVsig <- sigvarintervals < 0.05 #p-value of model comparison more significant that in 95% of randomizations


###Fitting the biomass model across all traits (see below for model, results go in figure 2, which is generated here)
set.seed(1) #to produce identical objects as in manuscript
Deps <- scale(TI_trts[,c(21:28,31,32,1:20)]) 
#above scales dependent variables for faster fitting (important in permutations, see below)
#scaling does not affect results, as the data in each variable are simply divided by all the same number (the mean)
Inds <- data.frame(TAnn_I= TI_trts$TAnn_I, TAnn_P = TI_trts$TAnn_P, SYM = TI_trts$SYM, mom = TI_trts$mom)
sigs <- matrix(NA, nrow=ncol(Deps),ncol=5)
slopes <- matrix(NA, nrow=ncol(Deps),ncol=5)
for(i in 1:ncol(Deps)){
	Inds$Dep <- Deps[,i]
	mod <- MCMCglmm(Dep ~I(TAnn_P/10) +I(TAnn_I/10) + SYM + SYM:I(TAnn_P/10), random = ~ mom  , data = Inds,verbose=FALSE,nitt=10000, thin = 10, burnin=1000,pr=T)
	sigs[i,] <- summary(mod)$solutions[,5]
	slopes[i,] <- summary(mod)$solutions[,1]
} #

###Randomization of the above to create null distribution
################may take upwards of 5 hrs#############
sampInds <- data.frame(TAnn_I= TI_trts$TAnn_I, TAnn_P = TI_trts$TAnn_P, SYM = TI_trts$SYM, mom = TI_trts$mom) #NOT permuted
sampsigs <- array(NA, dim=c(ncol(Deps),5,1000))
sampslopes <- array(NA, dim=c(ncol(Deps),5,1000))
for(j in 1:dim(sampsigs)[3]){
	sampDeps <- Deps[sample(1:nrow(Deps),replace=F),] #permuted
	for(i in 1:ncol(sampDeps)){
		sampInds$sampDep <- sampDeps[,i]
		sampmod <- MCMCglmm(sampDep ~I(TAnn_P/10) +I(TAnn_I/10) + SYM + SYM:I(TAnn_P/10), random = ~ mom  , data = sampInds,verbose=FALSE,nitt=10000, thin = 10, burnin=1000,pr=T)
		sampsigs[i,,j] <- summary(sampmod)$solutions[,5]
		sampslopes[i,,j] <- summary(sampmod)$solutions[,1]
	}
}
save(sampslopes,file="sampslopes.Rdata")
save(sampsigs,file="sampsigs.Rdata")
load("sampslopes.Rdata")
slopeintervals <- sapply(1:nrow(slopes), function(z) sapply(1:ncol(slopes), function(p) findInterval(slopes[z,p],sort(sampslopes[z,p,]))/(dim(sampslopes)[3])   )  )
permslopemn <- sapply(1:nrow(slopes), function(z) sapply(1:ncol(slopes), function(p) mean(sampslopes[z,p,]) )  )
isslopesig <- slopeintervals < 0.05 | slopeintervals > 0.95 &  t(abs(slopes)) > abs(permslopemn) 
##alternative evalution:
# load("sampsigs.Rdata")
# sigintervals <- sapply(1:nrow(sigs), function(z) sapply(1:ncol(sigs), function(p) findInterval(sigs[z,p],sort(sampsigs[z,p,]))/(dim(sampsigs)[3])   )  )
# issigsig <- sigintervals < 0.05


###Figure 2, main text
pdf("sumsofsquares_randsig.pdf",height=4,width=6) 
 par(oma=c(0,0,0,0))
par(mar=c(7.5,3,1,0))
dummy <- rbind(ssP,ssB,ssPB-ssP-ssB)
bg <- barplot((dummy*100),ylim=c(0,50),ylab="",xaxt="n",col=c(rgb(0.33,0.33,0.33),rgb(0.66,0.66,0.66),rgb(1,1,1)))
legend(15,50,c("Plant Source","Biota Source","Plant X Biota Source"),fill = c(rgb(0.33,0.33,0.33),rgb(0.66,0.66,0.66),rgb(1,1,1)),bty="n")
points(bg-0.4, ssPB*100+2, pch=ifelse(isamntVsig[1,],16,NA),col=rgb(0.33,0.33,0.33),cex=0.75)
points(bg, ssPB*100+2, pch=ifelse(isamntVsig[2,],16,NA),col=rgb(0.66,0.66,0.66),cex=0.75)
points(bg+0.4, ssPB*100+2, pch=ifelse(isamntVsig[3,],1,NA),cex=0.6)
mtext("% of variation explained",side=2,line=2)
axis(side=1,at=bg,las=2,labels=c(traitnamesNL,ionnames))
dev.off()

##Figure 4, main text
pdf("mods_randsig.pdf",height=2.5,width=7) 
layout(matrix(c(1,2),ncol=2),widths=c(5,1.6))
par(oma=c(0,0,0,0))
par(mar=c(7.5,2.6,0.1,0.25))
image(matrix(c(rep(c(0:2),each=30),rep(2, times=30)),nrow=30,ncol=4,byrow=F),col=c(rgb(0.4,0.4,0.4),rgb(0.8,0.8,0.8),rgb(1,1,1)),zlim=c(0,2),xaxt="n",yaxt="n")
abline(h=seq(from=0,to=1,length.out=7)[c(2,4,6)])
abline(v=seq(from=0,to=1,length.out=59)[seq(from=2,to=58,by=2)])
axis(side=1,at=seq(from=0,to=1,length.out=30),las=2,labels=c(traitnamesNL,ionnames))
text(seq(from=0,to=1,length.out=30),1,sapply(1:30, function(z) if(isslopesig[5,z]){ifelse(slopes[z,5]>0,"+","-")} else {""} ))
text(seq(from=0,to=1,length.out=30),0.666,sapply(1:30, function(z) if(isslopesig[4,z]){ifelse(slopes[z,4]>0,"+","-")} else {""} ))
text(seq(from=0,to=1,length.out=30),0.333,sapply(1:30, function(z) if(isslopesig[3,z]){ifelse(slopes[z,3]>0,"+","-")} else {""} ))
text(seq(from=0,to=1,length.out=30),0,sapply(1:30, function(z) if(isslopesig[2,z]){ifelse(slopes[z,2]>0,"+","-")} else {""} ))
axis(side=2,at=seq(from=0,to=1,length.out=4),las=2,labels=c(expression(beta[E[P]]),expression(beta[E[B]]),expression(beta[S]),expression(beta[EXS])))
par(mar=c(0,0,0,0))
plot(c(1:10)~I(c(1:10)),pch=NA,xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
legend(0,10,c("Plant Source","Biota Source","Plant X Biota Source"),cex=0.9, fill = c(rgb(0.4,0.4,0.4),rgb(0.8,0.8,0.8),rgb(1,1,1)),bty="n")
dev.off()


##REPORTING NUMBERS
range(round(100*(ssB)))
range(round(100*(ssP)))
range(round(100*(ssPB-ssP-ssB)))
range(round(100*(ssPB)))
mean(round(100*(ssB)))
mean(round(100*(ssP)))
mean(round(100*(ssPB-ssP-ssB)))
mean(round(100*(ssPB)))

sort(tapply(TI_trts$biomass,biotasite,mean,na.rm=T))
sort(tapply(TI_trts$biomass,biotasite,std.error,na.rm=T))

## ###including mom as source of variation - not as in manuscript, but of potential interest
# trait_sourcevarMom_Asig <- matrix(NA,ncol=32,nrow=5)
# for(i in 1:32){
#  datatmp <- TI_trts
#  datatmp$y <- TI_trts[,i]
# trait_sourcevarMom_Asig[1,i] <-  anova(lm(y~mom,data=datatmp))$"Pr(>F)"[1]
# trait_sourcevarMom_Asig[2,i] <- anova(lm(y~Soil,data=datatmp))$"Pr(>F)"[1]
# trait_sourcevarMom_Asig[3,i] <- anova(lm(y~Soil,data=datatmp),lm(y~mom + Soil,data=datatmp))$"Pr(>F)"[2]
# trait_sourcevarMom_Asig[4,i] <- anova(lm(y~mom,data=datatmp),lm(y~mom + Soil,data=datatmp))$"Pr(>F)"[2]
# trait_sourcevarMom_Asig[5,i] <- anova(lm(y~Soil+mom,data=datatmp),lm(y~mom*Soil,data=datatmp))$"Pr(>F)"[2]
# }
# colnames(trait_sourcevarMom_Asig) <- colnames(TI_trts[,1:32])
# sum(trait_sourcevarMom_Asig[1,]<0.05)
# sum(trait_sourcevarMom_Asig[2,]<0.05)
# sum(trait_sourcevarMom_Asig[5,]<0.05) #one of these is unloggeg hairs (i.e. violates normality assumption)
# round(100*(ssM-ssP)) #plant family explains an additional percent of the variance.
# round(100 - 100*(ssM + ssB + (ssPB-ssP-ssB))) #while remaining variance in measured variables not ascribed to treatments ranges widely






#####################
###Is variation in biomass, related to sympatry and/or environmental descriptors of plant and biota source site?
####################

####ANALYSIS OF BIOMASS MODELS, 

#################
set.seed(1) 
bio.ttsx <- MCMCglmm(biomass ~I(TAnn_P/10) +I(TAnn_I/10) + SYM + SYM:I(TAnn_P/10), random = ~ mom  , data = TI_trts,verbose=FALSE,nitt=1000000, thin = 100, burnin=100000,pr=T)
bio.pcsx <- MCMCglmm(biomass ~Pann_P +Pann_I + SYM + SYM:Pann_P, random = ~ mom  , data = TI_trts,verbose=FALSE,nitt=1000000, thin = 100, burnin=100000,pr=T)#1082
bio.phsx <- MCMCglmm(biomass ~lPp +lPi + SYM + SYM:lPp, random = ~ mom  , data = TI_trts,verbose=FALSE,nitt=1000000, thin = 100, burnin=100000,pr=T)#1139
bio.kksx <- MCMCglmm(biomass ~lKp +lKi + SYM + SYM:lKp, random = ~ mom  , data = TI_trts,verbose=FALSE,nitt=1000000, thin = 100, burnin=100000,pr=T) #1119
##note that allowing ns plant source to be moved to random effect of plant (repl E_SxS with Inoc source) slightly changes things. Pann becomes sig for SYM & ints.
##options for responding to reviewer
#1 presenting results for a slimmed 4x4 factorial
##  if allow removal of TAnn_P, Sym and Sym:I remain sig and of similar-ish sign and mag. if replace with random ~Plant, they are marginal and likewise similar. Pann likewise marginal, but similarly sig as for full dataset changes to plant source
######lP similarly n.s. for any sympatric term whether allow removal or not. (w/ remove simplifies to plant + inoc). 
######lKi actually becomes marginally significant for lKp, SYM and SYM:lKp. (in the inocK pos, sym pos, inocKxSym neg)
#2  removals of treatments to show that local adaptation estimate varies by MAT
##removing TZ and ML plants and inocula results in only TAnn_P sig (pos), unless allow removals, in which case ends at sig pos eff SYM, and of TAnn_P
##removing TC and TX plants and inocula results in the same model (regardless of removal TAnn_P)


###get objects ready to plot BIOMASS fitted models
rb<-colorRampPalette(c(rgb(1,0,0),rgb(0,0,1)))
#make posteriors for best models
postbiomid <- bio.ttsx$Sol
mommnbio <- rowMeans(postbiomid[,6:122])
#get HPDIs, means, biomass
bio.predtimax.tp.m <- sapply(1:length(T.s), function(z) mean(postbiomid[,1] + postbiomid[,2]*T.s[z] + postbiomid[,3]*max(TI_trts$TAnn_I/10) + postbiomid[,4]*0 + postbiomid[,5]*0*T.s[z] + mommnbio)) 
bio.predtimax.tp.ci <- sapply(1:length(T.s), function(z) HPDi(postbiomid[,1] + postbiomid[,2]*T.s[z] + postbiomid[,3]*max(TI_trts$TAnn_I/10) + postbiomid[,4]*0 + postbiomid[,5]*0*T.s[z] + mommnbio  ,.95))  
bio.predti186.tp.m <- sapply(1:length(T.s), function(z) mean(postbiomid[,1] + postbiomid[,2]*T.s[z] + postbiomid[,3]*18.6 + postbiomid[,4]*0 + postbiomid[,5]*0*T.s[z] + mommnbio)) 
bio.predti186.tp.ci <- sapply(1:length(T.s), function(z) HPDi(postbiomid[,1] + postbiomid[,2]*T.s[z] + postbiomid[,3]*18.6 + postbiomid[,4]*0 + postbiomid[,5]*0*T.s[z] + mommnbio  ,.95))  
bio.predti153.tp.m <- sapply(1:length(T.s), function(z) mean(postbiomid[,1] + postbiomid[,2]*T.s[z] + postbiomid[,3]*15.3 + postbiomid[,4]*0 + postbiomid[,5]*0*T.s[z] + mommnbio)) 
bio.predti153.tp.ci <- sapply(1:length(T.s), function(z) HPDi(postbiomid[,1] + postbiomid[,2]*T.s[z] + postbiomid[,3]*15.3 + postbiomid[,4]*0 + postbiomid[,5]*0*T.s[z] + mommnbio  ,.95))  
bio.predtimin.tp.m <- sapply(1:length(T.s), function(z) mean(postbiomid[,1] + postbiomid[,2]*T.s[z] + postbiomid[,3]*min(TI_trts$TAnn_I/10) + postbiomid[,4]*0 + postbiomid[,5]*0*T.s[z] + mommnbio))
bio.predtimin.tp.ci <- sapply(1:length(T.s), function(z) HPDi(postbiomid[,1] + postbiomid[,2]*T.s[z] + postbiomid[,3]*min(TI_trts$TAnn_I/10) + postbiomid[,4]*0 + postbiomid[,5]*0*T.s[z] + mommnbio ,.95)) 
bio_allopat_equalT_hpdi		<- sapply(1:length(tann), function(z) HPDi(postbiomid[,1] + postbiomid[,2]*tann[z] + postbiomid[,3]*tann[z] + postbiomid[,4]*0 + postbiomid[,5]*0*tann[z] + mommnbio ,.95)) 
allobiomn <-tapply(TI_trts$biomass[TI_trts$SYM==0],paste(TI_trts$TAnn_P[TI_trts$SYM==0]/10,TI_trts$TAnn_I[TI_trts$SYM==0]/10,sep="_"),mean,na.rm=T) 
allobiose <-tapply(TI_trts$biomass[TI_trts$SYM==0],paste(TI_trts$TAnn_P[TI_trts$SYM==0]/10,TI_trts$TAnn_I[TI_trts$SYM==0]/10,sep="_"),std.error,na.rm=T) 
TPallobio <-  as.numeric(lapply(1:36, function(z) strsplit(names(allobiomn),"_")[[z]][1] ))
TIallobio <-  as.numeric(lapply(1:36, function(z) strsplit(names(allobiomn),"_")[[z]][2] ))

#Plot biomass. Figure 3, main text
pdf("Biomass_modelpred_avgpoints.pdf",height=3,width=6) 
layout(matrix(1:3,ncol=3),widths =c(1,1,0.75))
par(mar=c(0,0,0,0))
par(oma=c(4,4,3,0))
plot(TI_trts$biomass ~ I(TI_trts$TAnn_P/10), pch = NA, xlim=bufferX(TI_trts$TAnn_P/10,0.01),ylim=bufferX(allobiomn,0.05),
		ylab="Biomass",xlab="",cex.axis=1.2)
	lines(bio.predtimax.tp.m~T.s,col=rgb(1,0,0))
	lines(bio.predtimin.tp.m~T.s,col=rgb(0,0,1))
	polygon(c((T.s), rev(T.s)), c( (bio.predtimax.tp.ci[1,]),rev(bio.predtimax.tp.ci[2,])), col=rgb(1,0,0,alpha=0.5), border = NA)
	polygon(c((T.s), rev(T.s)), c( (bio.predtimin.tp.ci[1,]),rev(bio.predtimin.tp.ci[2,])), col=rgb(0,0,1,alpha=0.5), border = NA)
	points(allobiomn ~TPallobio,pch=16	,col = rgb(range01(TIallobio),0,1-range01(TIallobio),alpha=1), cex=2 )
 	mtext("Biomass",side=2,line=2.5)
 	mtext("Allopatric",side=3,line=0.5)
plot(tapply(TI_trts$biomass[TI_trts$SYM==1],TI_trts$TAnn_P[TI_trts$SYM==1]/10,mean)~
	 tapply(TI_trts$TAnn_P[TI_trts$SYM==1]/10,TI_trts$TAnn_P[TI_trts$SYM==1]/10,mean),
	 col=rgb(range01(tann),0,1-range01(tann)), pch=16, cex=2,
	 xlab="", ylab="",xlim=bufferX(TI_trts$TAnn_P/10,0.01),ylim=bufferX(allobiomn,0.05),main="",yaxt="n",cex.axis=1.2)
	arrows(tann,bio_allopat_equalT_hpdi[1,], x1=tann,y1=bio_allopat_equalT_hpdi[2,],length=0,	
		col=rgb(range01(tann),0,1-range01(tann),alpha=0.5)) #,lty=2)		
 	mtext("Sympatric",side=3,line=0.5)
 	mtext("Plant source site mean annual temperature",side=1,line=2.5,at =13)
par(mar=c(2,5,2,5))
image(as.matrix(t(tann)), col=rev(rb(100)),zlim=c(min(tann),max(tann)),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side=4,at=seq(from=0, to=1,length.out=10),labels=tann,las=2)
	mtext("Biota Source MAT",side = 1, line=1)
	text(rep(0,times=4),seq(from=0, to=1,length.out=10)[c(2,8:10)],rep("*",times=4),col=rgb(1,1,1),cex=2)
dev.off()
#extract predictions for text
mtp <- mean(tann)
#biomass predictions for inocula from the lowest and highest temp site, for average plant source site
sapply(c(T.s[1],T.s[length(T.s)]), function(z) mean(postbiomid[,1] + postbiomid[,2]*mtp + postbiomid[,3]*z + postbiomid[,4]*0 + postbiomid[,5]*0*mtp + mommnbio)) 
sapply(c(T.s[1],T.s[length(T.s)]), function(z) HPDi(postbiomid[,1] + postbiomid[,2]*mtp + postbiomid[,3]*z + postbiomid[,4]*0 + postbiomid[,5]*0*mtp + mommnbio  ,.95)) 
#predicted means for coldest and highest temp site plants with allopatric biota from the same temp site
sapply(c(T.s[1],T.s[length(T.s)]), function(z) mean(postbiomid[,1] + postbiomid[,2]*z + postbiomid[,3]*z + postbiomid[,4]*0 + postbiomid[,5]*0*z + mommnbio)) 
sapply(c(T.s[1],T.s[length(T.s)]), function(z) HPDi(postbiomid[,1] + postbiomid[,2]*z + postbiomid[,3]*z + postbiomid[,4]*0 + postbiomid[,5]*0*z + mommnbio, .95)) 
#predicted means for the coldest and warmest temp site plants with sympatric biota
sapply(c(T.s[1],T.s[length(T.s)]), function(z) mean(postbiomid[,1] + postbiomid[,2]*z + postbiomid[,3]*z + postbiomid[,4]*1 + postbiomid[,5]*1*z + mommnbio)) 
sapply(c(T.s[1],T.s[length(T.s)]), function(z) HPDi(postbiomid[,1] + postbiomid[,2]*z + postbiomid[,3]*z + postbiomid[,4]*1 + postbiomid[,5]*1*z + mommnbio, .95)) 

symbiomn <- tapply(TI_trts$biomass[TI_trts$SYM==1],TI_trts$TAnn_P[TI_trts$SYM==1]/10,mean)
symbiose <- tapply(TI_trts$biomass[TI_trts$SYM==1],TI_trts$TAnn_P[TI_trts$SYM==1]/10,std.error)

##Expanded version, Figure S6
pdf("bio_Allpoints.pdf",height=4,width=8) 
layout(matrix(c(1,2,3),ncol=3,byrow=TRUE), widths=c(1,1,0.5))
par(oma=c(5,6,3,0))
par(mar=c(1,0,0,0))
plot(TI_trts$biomass ~ I(TI_trts$TAnn_P/10), pch = NA, xlim=bufferX(TI_trts$TAnn_P/10,0.01),ylim=c(2.2,4.6),
		ylab="Biomass",xlab="",cex.axis=1.2)
	polygon(c((T.s), rev(T.s)), c( (bio.predtimax.tp.ci[1,]),rev(bio.predtimax.tp.ci[2,])), col=rgb(1,0,0,alpha=0.5), border = NA)
	polygon(c((T.s), rev(T.s)), c( (bio.predtimin.tp.ci[1,]),rev(bio.predtimin.tp.ci[2,])), col=rgb(0,0,1,alpha=0.5), border = NA)
	polygon(c((T.s), rev(T.s)), c( (bio.predti186.tp.ci[1,]),rev(bio.predti186.tp.ci[2,])), col=rgb(range01(c(12.9,19.8,18.6))[3],0,1-range01(c(12.9,19.8,18.6))[3],alpha=.5), border = NA)
	polygon(c((T.s), rev(T.s)), c( (bio.predti153.tp.ci[1,]),rev(bio.predti153.tp.ci[2,])), col=rgb(range01(c(12.9,19.8,15.3))[3],0,1-range01(c(12.9,19.8,15.3))[3],alpha=.5), border = NA)
	lines(bio.predtimax.tp.m~T.s,col=rgb(1,0,0))
	lines(bio.predtimin.tp.m~T.s,col=rgb(0,0,1))
	lines(bio.predti186.tp.m~T.s,col=rgb(range01(c(12.9,19.8,18.6))[3],0,1-range01(c(12.9,19.8,18.6))[3]))
	lines(bio.predti153.tp.m~T.s,col=rgb(range01(c(12.9,19.8,15.3))[3],0,1-range01(c(12.9,19.8,15.3))[3]))
	points(allobiomn ~TPallobio,pch=16	,col = rgb(range01(TIallobio),0,1-range01(TIallobio),alpha=1), cex=2 )
	arrows(TPallobio,y0=allobiomn-allobiose,y1=allobiomn+allobiose,length=0,col = rgb(range01(TIallobio),0,1-range01(TIallobio),alpha=1))
 	mtext("Biomass",side=2,line=2.5)
 	mtext("Allopatric",side=3,line=0.5)
plot(TI_trts$biomass[TI_trts$SYM==1]~ I(TI_trts$TAnn_P[TI_trts$SYM==1]/10),ylim=c(2.2,4.6),
	 col=rgb(range01(TI_trts$TAnn_I),0,1-range01(TI_trts$TAnn_I),alpha=1)[TI_trts$SYM==1] , pch=NA,
	 xlab="", ylab="",xlim=bufferX(TI_trts$TAnn_P/10,0.01),main="",yaxt="n",cex.axis=1.2)
	sapply(1:10, function(z) polygon( tann[z]+c(-0.05,0.05,0.05,-0.05),rep(c(bio_allopat_equalT_hpdi[1,z],bio_allopat_equalT_hpdi[2,z]),each=2 ),
		col=rgb(range01(tann),0,1-range01(tann),alpha=0.5)[z],border=NA ) )#,lty=2)		
 	points(symbiomn~tann,pch=16,cex=2,col=rgb(range01(tann),0,1-range01(tann)))
 	arrows(tann,y0=symbiomn-symbiose,y1=symbiomn+symbiose,length=0,col=rgb(range01(tann),0,1-range01(tann)))
 	mtext("Sympatric",side=3,line=0.5)
	mtext("Plant Source MAT",side = 1, line=4,at=12.5)
par(mar=c(5,4.5,5,4))
image(as.matrix(t(tann)), col=rev(rb(100)),zlim=c(min(tann),max(tann)),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side=4,at=seq(from=0, to=1,length.out=10),labels=tann,las=2,cex.axis=1.25)
	mtext("Biota Source MAT",side = 1, line=1,at=0)
	text(rep(0,times=4),seq(from=0, to=1,length.out=10)[c(2,8:10)],rep("*",times=4),col=rgb(1,1,1),cex=3)
dev.off()





##################
###univariate phenotypic "selection" gradients
###############

tfittrt <- data.frame(cbind(Plant.Traits[,1:10],subionsfull[,-1]), TAnn_I = InocSourceEnvDat$TAnn)[InocSourceEnvDat$TAnn%in%c(130,153,186,198,NA),]
 #unlogged data only, log transformation not appropriate for selection analysis. Hansen & Houle 2008. 
 #results are not reported for variables that were extremely non-normal (W =<0.5), as these produce results driven heavily by outliers
 	#these are Na, Ni, Zn, As, & Se 
 #uses expanded leaf tissue concentration dataset including all 720 ids, some of which have missing data
 #then subsets to only the inocula treatments that were applied to every plant source
 
 tfittrt$TAnn_I[is.na(tfittrt$TAnn_I)] <- "none" #re-code uninoculated treatments

set.seed(1)
#global stand trait
fitmodsGT <- list()
fitmodsbGT <- list()
scaleddfsGT <- list()
for(trait in c(2:30)){
	datall <- data.frame(uy = tfittrt[,1], ux= tfittrt[,trait],inoc = as.character(tfittrt$TAnn_I))
	#inocula have strong effects on biomass, but this is not of interest here, remove this effect by scaling biomass "locally"
	 datall$y <- sapply(1:nrow(datall), function(z)  datall$uy[z] / mean(datall$uy[datall$inoc==datall$inoc[z]], na.rm=T))
	 datall$x <- sapply(1:nrow(datall), function(z)  datall$ux[z] / mean(datall$ux, na.rm=T) )
	scaleddfsGT[[trait]] <- datall
		datf <- datall[getfull(datall),]
		#main effect of inocula skipped in models, as we removed it above.
		fitmodsGT[[trait]] <- MCMCglmm(y~x:inoc, data = datf, verbose=FALSE,nitt=10000, thin = 10, burnin=100)
		fitmodsbGT[[trait]] <- MCMCglmm(y~x:inoc + I(x^2):inoc, data = datf, verbose=FALSE,nitt=10000, thin = 10, burnin=100)
		#nitt is low because must match permutation to be DIC comparable
}

set.seed(1)
#CHECK that linear and quadratic models make non-overlapping predictions in the range of the observed data
DiffPredGT <- array(NA,dim=c(5,1000,29))
numpointsGT <- matrix(NA, nrow=29,ncol=5)
for(i in 1:29){
	for(TI in 1:5){
	datall <- scaleddfsGT[[i+1]]
	yrange <- range(datall$y,na.rm=T)
	xrange <- range(datall$x,na.rm=T)
	datall <- datall[as.numeric(as.factor(datall$inoc))==TI,]
	xvals <- seq(from=min(datall$x,na.rm=T),to=max(datall$x,na.rm=T),length.out=1000)
	hpdilin <- sapply(1:1000, function(z) HPDi(fitmodsGT[[i+1]]$Sol[,1] + fitmodsGT[[i+1]]$Sol[,TI+1]*xvals[z], 0.95))
	hpdiquad <- sapply(1:1000, function(z) HPDi(fitmodsbGT[[i+1]]$Sol[,1] + fitmodsbGT[[i+1]]$Sol[,TI+1]*xvals[z] + fitmodsbGT[[i+1]]$Sol[,TI+6]*(xvals[z]^2),0.95))
	DiffPredGT[TI,,i] <- sapply(1:1000, function(z) hpdilin[2,z] < hpdiquad[1,z] | hpdilin[1,z] > hpdiquad[2,z] )
	if(sum(DiffPredGT[TI,,i])>0){ 
		numpointsGT[i,TI] <- length(which( sapply(datall$x, function(z) any( abs(z-xvals[which(DiffPredGT[TI,,i])])<0.01)) ))  
	} else{numpointsGT[i,TI] <- 0}
	}
}
QLpreddiffGT <- sapply(1:29, function(z) any(rowSums(DiffPredGT[,,z])>0) )

DICmodtype_matGT <- cbind(unlist(lapply(fitmodsGT, function(x) x$DIC)), unlist(lapply(fitmodsbGT, function(x) x$DIC)))
quadterm_p_matGT <-  sapply(2:30, function(x) any(summary(fitmodsbGT[[x]])$solutions[7:11,5] < 0.05))
traitswquadGT <- which(quadterm_p_matGT & DICmodtype_matGT[,2]<DICmodtype_matGT[,1] & QLpreddiffGT)
#  [1]  1  2  3  4  5  6 15 18 21 22 26 28
#These are the traits that have potentially non-linear [single-optimum/quadratic] selection gradients


fivecols <- c(rgb(range01(tann)[c(2,8:10)],0,1-range01(tann)[c(2,8:10)]),"black")

#####INSPECT linear and quadratic relationships, where quadratic relationship appears to be significant
par(mar=c(0,0,0,0))
par(oma=c(4,4,1,1))
layout(matrix(1:60,nrow=5,byrow=F))
for(i in traitswquadGT){
	for(TI in 1:5){
	datall <- scaleddfsGT[[i+1]]
	yrange <- range(datall$y,na.rm=T)
	xrange <- range(datall$x,na.rm=T)
	datall <- datall[as.numeric(as.factor(datall$inoc))==TI,]
	plot(datall$y~datall$x, pch=16,cex=0.5,ylab="",xlab="",xaxt="n",yaxt="n",ylim=yrange,xlim=xrange,col=fivecols[TI])
		if(i==1){axis(side=2)}
		if(i==1){mtext("Biomass",side=2,line=2)}
		if(TI==5){
			axis(side=1)
			mtext(c(traitnamesNL,ionnames)[i+1],side=1,line=2.5)
			}
		xvals <- seq(from=min(datall$x,na.rm=T),to=max(datall$x,na.rm=T),length.out=1000)
		lines(sapply(1:1000, function(z) mean(fitmodsGT[[i+1]]$Sol[,1] + fitmodsGT[[i+1]]$Sol[,TI+1]*xvals[z])) ~ xvals, 
			col = fivecols[TI])
		hpdilin <- sapply(1:1000, function(z) HPDi(fitmodsGT[[i+1]]$Sol[,1] + fitmodsGT[[i+1]]$Sol[,TI+1]*xvals[z], 0.95))
 		polygon(c(xvals,rev(xvals)),c(hpdilin[1,],rev(hpdilin[2,])), col=rgb(0,0,0,alpha=0.25), border = NA )
		lines(sapply(1:1000, function(z) mean(fitmodsbGT[[i+1]]$Sol[,1] + fitmodsbGT[[i+1]]$Sol[,TI+1]*xvals[z] + fitmodsbGT[[i+1]]$Sol[,TI+6]*(xvals[z]^2))) ~ xvals, 
			col = fivecols[TI],lty=2)
		hpdiquad <- sapply(1:1000, function(z) HPDi(fitmodsbGT[[i+1]]$Sol[,1] + fitmodsbGT[[i+1]]$Sol[,TI+1]*xvals[z] + fitmodsbGT[[i+1]]$Sol[,TI+6]*(xvals[z]^2),0.95))
 		polygon(c(xvals,rev(xvals)),c(hpdiquad[1,],rev(hpdiquad[2,])), col=rgb(0,0,0,alpha=0.25), border = NA )
	}
}
lapply(fitmodsbGT[traitswquadGT+1],summary) #also check model output
## based on plots and model summaries, 
	#leaf number has too few unique x-values (max 5, most points in only 2 different values) and too few points at the extremes to make fitting parabolas reasonable.
	#while Nickel is also problematic, this result is excluded due to non-normality already
	#Fe also appears problematic, but many of the selection parameters are n.s. in model or after permutations

#####INSPECT remaining linear relationships
par(mar=c(0,0,0,0))
par(oma=c(4,4,1,1))
layout(matrix(1:85,nrow=5,byrow=F))
for(i in which(!1:29%in%traitswquadGT)){
	for(TI in 1:5){
	datall <- scaleddfsGT[[i+1]]
	yrange <- range(datall$y,na.rm=T)
	xrange <- range(datall$x,na.rm=T)
	datall <- datall[as.numeric(as.factor(datall$inoc))==TI,]
	plot(datall$y~datall$x, pch=16,cex=0.5,ylab="",xlab="",xaxt="n",yaxt="n",ylim=yrange,xlim=xrange,col=fivecols[TI])
		if(i==1){axis(side=2)}
		if(i==1){mtext("Biomass",side=2,line=2)}
		if(TI==5){
			axis(side=1)
			mtext(c(traitnamesNL,ionnames)[i+1],side=1,line=2.5)
			}
		xvals <- seq(from=min(datall$x,na.rm=T),to=max(datall$x,na.rm=T),length.out=1000)
		lines(sapply(1:1000, function(z) mean(fitmodsGT[[i+1]]$Sol[,1] + fitmodsGT[[i+1]]$Sol[,TI+1]*xvals[z])) ~ xvals, 
			col = fivecols[TI])
		hpdilin <- sapply(1:1000, function(z) HPDi(fitmodsGT[[i+1]]$Sol[,1] + fitmodsGT[[i+1]]$Sol[,TI+1]*xvals[z], 0.95))
 		polygon(c(xvals,rev(xvals)),c(hpdilin[1,],rev(hpdilin[2,])), col=rgb(0,0,0,alpha=0.25), border = NA )
	}
}
#remaining variables that were highly non-normal show problematic results (Na, Zn, As, Se), these are not reported

#given results of inspection, shift leaf number model to linear only, exclude from quadratic
traitswquadGTn <- traitswquadGT[-4]

mn_tfslopesGT <- sapply(c(2:30), function(trait) summary(fitmodsGT[[trait]])$solutions[,1])#table cols are in order of traits
lhpdi_tfslopesGT <- sapply(c(2:30), function(trait) summary(fitmodsGT[[trait]])$solutions[,2])#table cols are in order of traits
uhpdi_tfslopesGT <- sapply(c(2:30), function(trait) summary(fitmodsGT[[trait]])$solutions[,3])#table cols are in order of traits

mn_tfslopesQGT <- sapply(c(2:30), function(trait) summary(fitmodsbGT[[trait]])$solutions[,1])#table cols are in order of traits
lhpdi_tfslopesQGT <- sapply(c(2:30), function(trait) summary(fitmodsbGT[[trait]])$solutions[,2])#table cols are in order of traits
uhpdi_tfslopesQGT <- sapply(c(2:30), function(trait) summary(fitmodsbGT[[trait]])$solutions[,3])#table cols are in order of traits

mn_tfslopesOPT <- sapply(c(2:30), function(trait) sapply(1:5, function(inoc) 
	mean( (-1*fitmodsbGT[[trait]]$Sol[,1+inoc]) / (2*fitmodsbGT[[trait]]$Sol[,6+inoc])  ) ) )#table cols are in order of traits
uhpdi_tfslopesOPT <- sapply(c(2:30), function(trait) sapply(1:5, function(inoc) 
	HPDi( (-1*fitmodsbGT[[trait]]$Sol[,1+inoc]) / (2*fitmodsbGT[[trait]]$Sol[,6+inoc]),0.95  )[2] ) )
lhpdi_tfslopesOPT <- sapply(c(2:30), function(trait) sapply(1:5, function(inoc) 
	HPDi( (-1*fitmodsbGT[[trait]]$Sol[,1+inoc]) / (2*fitmodsbGT[[trait]]$Sol[,6+inoc]),0.95  )[1] ) )
#-b / 2a , b is linear, a is quad


##permutations
# ######EXPECTED TO TAKE ~4 hours on a laptop, depending on specific conditions
# #permute each trait with respect to biomass. do this within treatments, rather than at global level, 
# #so that can compare DIC between intercept model. we want both: DIC real < DIC permuted AND in real data, that DIC with slope < DIC with intercept
# # OR do we want both: slope real > slope permuted AND in real data, that DIC with slope < DIC with intercept
set.seed(1)
tfittrt_130 <- tfittrt[tfittrt$TAnn_I=="130",]   
tfittrt_153 <- tfittrt[tfittrt$TAnn_I=="153",]   
tfittrt_186 <- tfittrt[tfittrt$TAnn_I=="186",]   
tfittrt_198 <- tfittrt[tfittrt$TAnn_I=="198",]  
tfittrt_none <- tfittrt[tfittrt$TAnn_I=="none",]   
permdicsGT <- array(NA,dim=c(ncol(tfittrt)-2,1000))#1000
permslopesGT <- array(NA,dim=c(11,ncol(tfittrt)-2,1000))#1000
for(i in 1:dim(permslopesGT)[3]){
 	newd <- data.frame(biomass = tfittrt[,1],  #add biomass, not permuted
 		rbind(tfittrt_130[sample(1:nrow(tfittrt_130),replace=F),2:(ncol(tfittrt_130)-1)], #permute traits wrt to biomass within each treatment
 		tfittrt_153[sample(1:nrow(tfittrt_153),replace=F),2:(ncol(tfittrt_153)-1)], 
 		tfittrt_186[sample(1:nrow(tfittrt_186),replace=F),2:(ncol(tfittrt_186)-1)], 
 		tfittrt_198[sample(1:nrow(tfittrt_198),replace=F),2:(ncol(tfittrt_198)-1)], 
 		tfittrt_none[sample(1:nrow(tfittrt_none),replace=F),2:(ncol(tfittrt_none)-1)]),  #sampled trait data.	
 		TAnn_I=tfittrt[,ncol(tfittrt)])
	for(trait in c(2:30)){
		samp_trait <- data.frame(uy = newd[,1], ux= newd[,trait],inoc = as.character(newd$TAnn_I))
		#inocula have strong effects on biomass, remove this effect by scaling biomass and traits "locally"
		samp_trait$y <- sapply(1:nrow(samp_trait), function(z)  samp_trait$uy[z] / mean(samp_trait$uy[samp_trait$inoc==samp_trait$inoc[z]], na.rm=T))
		samp_trait$x <- sapply(1:nrow(samp_trait), function(z)  samp_trait$ux[z] / mean(samp_trait$ux[samp_trait$inoc==samp_trait$inoc[z]], na.rm=T))
		samp_traitf <- samp_trait[getfull(samp_trait),]
	if((trait-1) %in% traitswquadGTn){		
		sampfitmod <- MCMCglmm(y~x:inoc + I(x^2):inoc, data = samp_traitf, verbose=FALSE,nitt=10000, thin = 10, burnin=100)
		permslopesGT[,trait-1,i] <- summary(sampfitmod)$solutions[,1]  #intercept included, though in real and permuted data we expect this to be 1 (the average of mean-scaled biomass data should be 1)
	}
	else{
		sampfitmod <- MCMCglmm(y~x:inoc, data = samp_traitf, verbose=FALSE,nitt=10000, thin = 10, burnin=100)
		permslopesGT[,trait-1,i] <- c(summary(sampfitmod)$solutions[,1], rep(NA, times=5))  #intercept included, though in real and permuted data we expect this to be 1 (the average of mean-scaled biomass data should be 1)	
	}
		permdicsGT[trait-1,i] <-  sampfitmod$DIC
	}
print(i)
}
 save(permdicsGT,file="permdicsGTQ.Rdata")
 save(permslopesGT,file="permslopesGTQ.Rdata")


load("permslopesGTQ.Rdata")
slopeintervalSelGT <- matrix(NA, nrow=29,ncol=11)
permslopemnSelGT <- matrix(NA, nrow=29,ncol=11)
isslopesigSelGT <- matrix(NA, nrow=29,ncol=11)
for(trait in 1:29){
	if(trait %in% traitswquadGTn){
		slopeintervalSelGT[trait,] <- sapply(1:nrow(mn_tfslopesQGT), function(z) findInterval(mn_tfslopesQGT[z,trait],sort(permslopesGT[z,trait,]))/(dim(permslopesGT)[3])   )
		permslopemnSelGT[trait,] <- sapply(1:nrow(mn_tfslopesQGT), function(z)  mean(permslopesGT[z,trait,]) )  
		isslopesigSelGT[trait,] <- (slopeintervalSelGT[trait,] < 0.05 | slopeintervalSelGT[trait,] > 0.95) &  abs(mn_tfslopesQGT[,trait]) > abs(permslopemnSelGT[trait,]) #real slope needs to be both outside the permutation window AND further from 0
	} else {
		slopeintervalSelGT[trait,1:6] <- sapply(1:nrow(mn_tfslopesGT), function(z) findInterval(mn_tfslopesGT[z,trait],sort(permslopesGT[z,trait,]))/(dim(permslopesGT)[3])   )
		permslopemnSelGT[trait,1:6] <- sapply(1:nrow(mn_tfslopesGT), function(z) mean(permslopesGT[z,trait,]) )
		isslopesigSelGT[trait,1:6] <- (slopeintervalSelGT[trait,1:6] < 0.05 | slopeintervalSelGT[trait,1:6] > 0.95) &  abs(mn_tfslopesGT[,trait]) > abs(permslopemnSelGT[trait,1:6]) #real slope needs to be both outside the permutation window AND further from 0
	}
}

##Figure S7
fivecols <- c(rgb(range01(tann)[c(2,8:10)],0,1-range01(tann)[c(2,8:10)]),"black")
pdf("Trait-fitness-slopes_wquad.pdf",height=6,width=5.5) 
reportedtraits <- c(1:10,12:20,22,26:29)
layout(matrix(c(1:25),ncol=5,byrow=TRUE))
par(oma=c(1,4,0.25,0))
par(mar=c(0, 0, 0, 0))
for(i in reportedtraits){
	if(i %in% traitswquadGTn){	
	plot(mn_tfslopesQGT[-1,i]~c(1:5,1:5+0.25),xaxt="n", pch = c(rep(16,times=5),rep(17,times=5)),
		xlim=c(0,5.5), ylim=c(-2.25,3),yaxt="n", col=fivecols)
	arrows(c(1:5,1:5+0.25),y0=lhpdi_tfslopesQGT[-1,i],y1=uhpdi_tfslopesQGT[-1,i],length=0, col=fivecols)	
 	text(c(1:5,1:5+0.25),c(rep(3,times=5),rep(-1.4,times=5)),ifelse(isslopesigSelGT[i,-1],"*",""),col= fivecols)
	} else {
	plot(mn_tfslopesGT[-1,i]~c(1:5),pch=16,xaxt="n", xlim=c(0,5.5), ylim=c(-2.25,3),yaxt="n", col=fivecols)
	arrows(c(1:5),y0=lhpdi_tfslopesGT[-1,i],y1=uhpdi_tfslopesGT[-1,i],length=0, col=fivecols)
 	text(c(1:5),3,ifelse(isslopesigSelGT[i,-1],"*",""),col= fivecols)
	}
	abline(h=0,lty=3)
	if(i %in% c(1,6,12,17,26)){axis(side=2)}
	if(i ==12){mtext("Estimated phenotypic selection",side=2,line=2)}else{}
	text(0,-2,c(traitnamesNL,ionnames)[i+1],adj=0)
}
plot(1:10~c(1:10),pch=NA,bty="n",xaxt="n",yaxt="n")
legend(0,10,c("13.0","15.3","18.6","19.8","none"), bty="n",fill=fivecols)
	text(5.25,9.75,"Biota Source MAT")
legend(0.25,3.5,c("linear","quadratic"), bty="n",pch=c(16,17))
	text(5.5,3.25,"Selection term form")
dev.off()


#fewer traits, Figure 5, main text
reduceset <- c(2,3,6,9,14,15,16,26)
pdf("Trait-fitness-slopes_sub_wquad.pdf",height=3.75,width=3.75) 
layout(matrix(c(1:9),ncol=3,byrow=TRUE))
par(oma=c(1,4,0.25,0))
par(mar=c(0, 0, 0, 0))
for(i in reduceset){
	if(i %in% traitswquadGTn){	
	plot(mn_tfslopesQGT[-1,i]~c(1:5,1:5+0.25),xaxt="n", pch = c(rep(16,times=5),rep(17,times=5)),
		xlim=c(0,5.5), ylim=c(-2.25,3),yaxt="n", col=fivecols)
	arrows(c(1:5,1:5+0.25),y0=lhpdi_tfslopesQGT[-1,i],y1=uhpdi_tfslopesQGT[-1,i],length=0, col=fivecols)	
 	text(c(1:5,1:5+0.25),c(rep(3,times=5),rep(-1.4,times=5)),ifelse(isslopesigSelGT[i,-1],"*",""),col= fivecols)
	} else {
	plot(mn_tfslopesGT[-1,i]~c(1:5),pch=16,xaxt="n", xlim=c(0,5.5), ylim=c(-2.25,3),yaxt="n", col=fivecols)
	arrows(c(1:5),y0=lhpdi_tfslopesGT[-1,i],y1=uhpdi_tfslopesGT[-1,i],length=0, col=fivecols)
 	text(c(1:5),3,ifelse(isslopesigSelGT[i,-1],"*",""),col= fivecols)
	}
	abline(h=0,lty=3)
	if(i %in% c(2,9,16)){axis(side=2)} 
	if(i ==9){mtext("Estimated phenotypic selection",side=2,line=2)}else{}
	text(0,-2,c(traitnamesNL,ionnames)[i+1],adj=0)
}
plot(1:10~c(1:10),pch=NA,bty="n",xaxt="n",yaxt="n")
legend(0,10,c("13.0","15.3","18.6","19.8","none"), bty="n",fill=fivecols)
	text(4.9,9.75,"Biota Source MAT")
legend(0.25,3.5,c("linear","quadratic"), bty="n",pch=c(16,17))
	text(5.25,3.3,"Selection term form")
dev.off()

##Figure S8, panel a
pdf("Trait-fitness-optima_quad.pdf",height=3.5,width=3.75) 
par(mfrow=c(3,3))
par(mar=c(0,0,0,0))
par(oma=c(0.5,4,0.5,0.5))
for(trait in traitswquadGTn[-c(7,8)] ) { #one not reported due to data, the other has no optima (all significant quadratic terms suggest diversifying sel.)
plot(1~c(1), pch=NA, ylim=c(-5,7),xlim=c(0,6), ylab="",xlab="",yaxt="n", xaxt="n")
	for(z in 1:5){ 
		if( isslopesigSelGT[trait,z+6] & mn_tfslopesQGT[z+6,trait]<0 ){
		points(z, mn_tfslopesOPT[z,trait], pch=1, cex=1.25, col=fivecols[z]) 
		arrows(z, y0=lhpdi_tfslopesOPT[z,trait], y1=uhpdi_tfslopesOPT[z,trait],length=0,col=fivecols[z],lwd=1.25)} else { 
		text(z,0,"") } 
	}
	if(trait %in% c(1,5,22)){axis(side=2)}
	if(trait==5){mtext("Estimated (relative) optimum trait value",side=2,line=2)}
	if(trait==15){
		legend(-0.5,6.5,c("13.0","15.3","18.6","19.8","none"), bty="n",fill=fivecols,ncol=2)
		text(0,6.25,"Biota Source MAT",adj=c(0,0))		}
	text(0,-4.6,c(traitnamesNL,ionnames)[trait+1],adj=0)
	abline(h=0,lty=3)
	}
dev.off()

##Figure S8, panel b
#stem height and rubidium
pdf("Trait-fitness-optdiff_quad.pdf",height=2,width=5) 
fivecollite <- rgb(c(0,.33,.66,1,0), 0, c(1,0.66,0.33,0,0),alpha=0.1)
par(mar=c(0,0,0,0))
par(oma=c(4,4,1,1))
par(mfrow=c(1,3))
#plot stem height, inoc 13.0, 15.3 and none
	i = 2
	datall <- scaleddfsGT[[i+1]]
	yrange <- range(datall$y,na.rm=T)
	xrange <- range(datall$x,na.rm=T)
	plot(datall$y~datall$x, pch=NA,cex=0.5,ylab="",xlab="",ylim=yrange,xlim=xrange,col=fivecols[TI])
			mtext(c(traitnamesNL,ionnames)[i+1],side=1,line=2.5)
	mtext("Relative fitness", side=2,line=2)
	for(TI in c(1,3,5)){
		datallsub <- datall[as.numeric(as.factor(datall$inoc))==TI,]
		points(datallsub$y~datallsub$x, pch=16,cex=0.5,col=fivecols[TI])
		xvals <- seq(from=min(datall$x,na.rm=T),to=max(datall$x,na.rm=T),length.out=1000)
		lines(sapply(1:1000, function(z) mean(fitmodsbGT[[i+1]]$Sol[,1] + fitmodsbGT[[i+1]]$Sol[,TI+1]*xvals[z] + fitmodsbGT[[i+1]]$Sol[,TI+6]*(xvals[z]^2))) ~ xvals, 
			col = fivecols[TI],lty=1)
		hpdiquad <- sapply(1:1000, function(z) HPDi(fitmodsbGT[[i+1]]$Sol[,1] + fitmodsbGT[[i+1]]$Sol[,TI+1]*xvals[z] + fitmodsbGT[[i+1]]$Sol[,TI+6]*(xvals[z]^2),0.95))
 		polygon(c(xvals,rev(xvals)),c(hpdiquad[1,],rev(hpdiquad[2,])), col=fivecollite[TI], border = NA )
	}
#plot rubidium, inoc 15.3 and 18.6
	i = 26
	datall <- scaleddfsGT[[i+1]]
	yrange <- range(datall$y,na.rm=T)
	xrange <- range(datall$x,na.rm=T)
	plot(datall$y~datall$x, pch=NA,cex=0.5,ylab="",xlab="",yaxt="n",ylim=yrange,xlim=xrange,col=fivecols[TI])
			mtext(c(traitnamesNL,ionnames)[i+1],side=1,line=2.5)
	for(TI in c(2,3)){
		datallsub <- datall[as.numeric(as.factor(datall$inoc))==TI,]
		points(datallsub$y~datallsub$x, pch=16,cex=0.5,col=fivecols[TI])
		xvals <- seq(from=min(datall$x,na.rm=T),to=max(datall$x,na.rm=T),length.out=1000)
		lines(sapply(1:1000, function(z) mean(fitmodsbGT[[i+1]]$Sol[,1] + fitmodsbGT[[i+1]]$Sol[,TI+1]*xvals[z] + fitmodsbGT[[i+1]]$Sol[,TI+6]*(xvals[z]^2))) ~ xvals, 
			col = fivecols[TI],lty=1)
		hpdiquad <- sapply(1:1000, function(z) HPDi(fitmodsbGT[[i+1]]$Sol[,1] + fitmodsbGT[[i+1]]$Sol[,TI+1]*xvals[z] + fitmodsbGT[[i+1]]$Sol[,TI+6]*(xvals[z]^2),0.95))
 		polygon(c(xvals,rev(xvals)),c(hpdiquad[1,],rev(hpdiquad[2,])), col=fivecollite[TI], border = NA )
	}
plot(1~c(1),ylim=c(0,10),xlim=c(0,5),pch=NA,xaxt="n",yaxt="n",bty="n")
		legend(0,8,c("13.0","15.3","18.6","19.8 (not plotted)","none"), bty="n",fill=fivecols)
		text(0,9,"Biota Source MAT",adj=c(0,0))		
dev.off()

############## #Alternative method: multiple traits, true phenotypic selection "gradients"
############fits model to biomass using all traits and reverse elimination until models no longer contain n.s. terms
###########REPORTED IN SUPPLEMENT
#checking assumptions
sign(cor(tfittrt[getfull(tfittrt),-31])>0.8) #Fe and Al, Ca and Sr, 
# not all traits can be passed in, as some are too tightly correlated. Aluminium and strontium not used as explanatory variables.
# Also eliminate Na, Ni, Zn, As, and Se due to normality issues mentioned above
# Include quadratic term in all starting models when it was significant for any single-trait models for any inoculum treatment above

scaledtraitmat <- as.data.frame(sapply(2:30, function(z) scaleddfsGT[[z]]$x))
colnames(scaledtraitmat) <- colnames(tfittrt)[2:30]
scaledtraitmat$biomass <- scaleddfsGT[[2]]$y #can be any scaleddfsGT, all have same y and inoc values
scaledtraitmat$inoc <- scaleddfsGT[[2]]$inoc #can be any scaleddfsGT, all have same y and inoc values

scaledtrait130 <- scaledtraitmat[scaledtraitmat$inoc=="130",]
scaledtrait153 <- scaledtraitmat[scaledtraitmat$inoc=="153",]
scaledtrait186 <- scaledtraitmat[scaledtraitmat$inoc=="186",]
scaledtrait198 <- scaledtraitmat[scaledtraitmat$inoc=="198",]
scaledtraitnone <- scaledtraitmat[scaledtraitmat$inoc=="none",]

nutrterms <- c("B11", "Mg26", "P31", "S34", "I(S34^2)", "K39", "Ca44", "Fe54", "I(Fe54^2)", "Mn55", 
			"Co59", "Cu63", "I(Cu63^2)", "Rb85", "I(Rb85^2)", "Mo98", "I(Mo98^2)", "Cd111")
trtterms <- c("SR", "I(SR^2)", "Ht", "I(Ht^2)", "Wt", "I(Wt^2)", 
			"Lfnum", "LL", "I(LL^2)", "LW", "I(LW^2)", "GrthDelay", "germday", "hairs")


#Split traits into morphology and nutrition (leaf tissue element concentrations) to limit number of traits in any one model
###nutr
summary(MCMCglmm(biomass ~ B11 + Mg26 + P31 + S34 + I(S34^2) + K39 + Ca44 + Fe54 + I(Fe54^2) + Mn55 + 
			Co59 + Cu63 + I(Cu63^2) + Rb85 + I(Rb85^2) + Mo98 + I(Mo98^2) + Cd111, 
			data = scaledtrait130[getfull(scaledtrait130),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: Fe54^2, Fe54, Mo98, Mo98^2, Cd111, B11 (slight worse), Mg26(worse), P31 (worse), Co59 (worse), no further
			#best stopped before B11 (had B11 and all subs rms)
nutr_130 <- (MCMCglmm(biomass ~ B11 + Mg26 + P31 + S34 + I(S34^2) + K39 + Ca44 + Mn55 + 
			Co59 + Cu63 + I(Cu63^2) + Rb85 + I(Rb85^2), data = scaledtrait130[getfull(scaledtrait130),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
summary(MCMCglmm(biomass ~ B11 + Mg26 + P31 + S34 + I(S34^2) + K39 + Ca44 + Fe54 + I(Fe54^2) + Mn55 + 
			Co59 + Cu63 + I(Cu63^2) + Rb85 + I(Rb85^2) + Mo98 + I(Mo98^2) + Cd111, 
			data = scaledtrait153[getfull(scaledtrait153),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop:Mo98 (Mn close), Mn, Fe^2, Co59, Cd111 (Mo^2 close), Mo^2, B11 (Rb^2 close, DIC worse),
				# Mg (Rb^2 close, DIC worse), Rb^2 (worse), Ca44 (worse), stop. best had B11 and all subs rms
nutr_153 <- (MCMCglmm(biomass ~ B11 + Mg26 + P31 + S34 + I(S34^2) + K39 + Ca44 + Fe54 + Cu63 + I(Cu63^2) + Rb85 + I(Rb85^2), 
			data = scaledtrait153[getfull(scaledtrait153),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
summary(MCMCglmm(biomass ~ B11 + Mg26 + P31 + S34 + I(S34^2) + K39 + Ca44 + Fe54 + I(Fe54^2) + Mn55 + 
			Co59 + Cu63 + I(Cu63^2) + Rb85 + I(Rb85^2) + Mo98 + I(Mo98^2) + Cd111, 
			data = scaledtrait186[getfull(scaledtrait186),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: Mo^2, Ca44, Mo, S^2 (Mn, Co), Co, Mg (B, Cd, Rb^2), Mn, Cd (Cu,B), B (Cu, Rb^2), Rb^2, Cu, 
				#P (DIC worse), Cu^2 (worse), Fe (worse),  Fe^2 (recovered to close but not quite best, mod only S, K, Rb)
				#best has P and all subs rms
nutr_186 <- (MCMCglmm(biomass ~ P31 + S34 + K39 + Fe54 + I(Fe54^2) + I(Cu63^2) + Rb85, 
			data = scaledtrait186[getfull(scaledtrait186),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T))
summary(MCMCglmm(biomass ~ B11 + Mg26 + P31 + S34 + I(S34^2) + K39 + Ca44 + Fe54 + I(Fe54^2) + Mn55 + 
			Co59 + Cu63 + I(Cu63^2) + Rb85 + I(Rb85^2) + Mo98 + I(Mo98^2) + Cd111, 
			data = scaledtrait198[getfull(scaledtrait198),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: Mn, S, Cd, Mo^2, Ca (Rb^2), Rb^2, Fe (2/3, P close), Fe^2, P, Rb, Mo, Mg, 
			# B11 (worse), Cu^2 (worse), Cu (DIC BEST YET), S^2 (DIC worse), no more can be rm
nutr_198 <- (MCMCglmm(biomass ~I(S34^2) + K39 + Co59 , data = scaledtrait198[getfull(scaledtrait198),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T))
summary(MCMCglmm(biomass ~ B11 + Mg26 + P31 + S34 + I(S34^2) + K39 + Ca44 + Fe54 + I(Fe54^2) + Mn55 + 
			Co59 + Cu63 + I(Cu63^2) + Rb85 + I(Rb85^2) + Mo98 + I(Mo98^2) + Cd111, 
			data = scaledtraitnone[getfull(scaledtraitnone),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: Mo^2 (Fe^2), Fe^2, S, Mo (B, P), P, B (Co, bit worse), S^2 (DIC recover), Mn (worse), Co (worse)
				#Fe, no more can be rm, best was model immed after rm S^2 (has Mn and all subs rms)
nutr_none <- (MCMCglmm(biomass ~ Mg26 + K39 + Ca44 + Fe54 + Mn55 + Co59 + Cu63 + I(Cu63^2) + Rb85 + I(Rb85^2) + Cd111, 
			data = scaledtraitnone[getfull(scaledtraitnone),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T))  

##morph
summary(MCMCglmm(biomass ~ SR + I(SR^2) + Ht + I(Ht^2) + Wt + I(Wt^2) + 
			Lfnum + LL + I(LL^2) + LW + I(LW^2) + GrthDelay + germday + hairs, 
			data = scaledtrait130[getfull(scaledtrait130),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: LW^2, Wt^2, hairs, LW, LL^2 (DIC slight worse), LL (DIC not better or worse), SR^2 (DIC worse)
morph_130 <- (MCMCglmm(biomass ~ SR + I(SR^2) + Ht + I(Ht^2) + Wt +  Lfnum + LL + I(LL^2) + GrthDelay + germday, 
			data = scaledtrait130[getfull(scaledtrait130),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T))
summary(MCMCglmm(biomass ~ SR + I(SR^2) + Ht + I(Ht^2) + Wt + I(Wt^2) + 
			Lfnum + LL + I(LL^2) + LW + I(LW^2) + GrthDelay + germday + hairs, 
			data = scaledtrait153[getfull(scaledtrait153),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: LW^2 (SR, SR^2, LW), SR (SR^2), Wt (LW), LW, 
				#Lfnum (hairs, slightly worse after rm Lfnum), hairs (worse), LL (worse), LL^2 (worse), stop
morph_153 <- (MCMCglmm(biomass ~ I(SR^2) + Ht + I(Ht^2) + Wt + Lfnum + LL + I(LL^2) + GrthDelay + germday + hairs, 
			data = scaledtrait153[getfull(scaledtrait153),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T))
summary(MCMCglmm(biomass ~ SR + I(SR^2) + Ht + I(Ht^2) + Wt + I(Wt^2) + 
			Lfnum + LL + I(LL^2) + LW + I(LW^2) + GrthDelay + germday + hairs, 
			data = scaledtrait186[getfull(scaledtrait186),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: Wt^2, hairs (lfnum), Lfnum, LL (worse after), LL^2 (DIC now better than before),
				#germday (DIC tiny bit worse), Ht^2 (worse), SR^2 (worse), LW (worse), LW^2 (worse), no more poss
				#best includes germday and all subs rm.
morph_186 <- (MCMCglmm(biomass ~ SR + I(SR^2) + Ht + I(Ht^2) + Wt + LW + I(LW^2) + GrthDelay + germday , 
			data = scaledtrait186[getfull(scaledtrait186),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
summary(MCMCglmm(biomass ~ SR + I(SR^2) + Ht + I(Ht^2) + Wt + I(Wt^2) + 
			Lfnum + LL + I(LL^2) + LW + I(LW^2) + GrthDelay + germday + hairs, 
			data = scaledtrait198[getfull(scaledtrait198),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: SR^2, Ht^2, Lfnum, germday, LL^2 (LL), LL, Wt (worse), GrthDelay (worse) no more can be rm
				#best had Wt and all subs rm'd 
morph_198 <- (MCMCglmm(biomass ~ SR + Ht + Wt + I(Wt^2) + LW + I(LW^2) + GrthDelay + hairs, 
			data = scaledtrait198[getfull(scaledtrait198),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T))
summary(MCMCglmm(biomass ~ SR + I(SR^2) + Ht + I(Ht^2) + Wt + I(Wt^2) + 
			Lfnum + LL + I(LL^2) + LW + I(LW^2) + GrthDelay + germday + hairs, 
			data = scaledtraitnone[getfull(scaledtraitnone),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 
			#drop: Wt^2, Lfnum, germday (GrthDelay), LW^2, GrthDelay (SR^2), SR^2 (DIC smallest bit worse),
				#LL^2 (worse), LL (recovery to almost before dropped SR^2), no more can be rm.
morph_none <- (MCMCglmm(biomass ~ SR + I(SR^2) + Ht + I(Ht^2) + Wt +
			LL + I(LL^2) + LW + hairs, 
			data = scaledtraitnone[getfull(scaledtraitnone),],verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)) 

listnutrmod <- list(nutr_130, nutr_153, nutr_186, nutr_198, nutr_none)
listmorphmod <- list(morph_130, morph_153, morph_186, morph_198, morph_none)

nutrterms <- c("B11", "Mg26", "P31", "S34", "I(S34^2)", "K39", "Ca44", "Fe54", "I(Fe54^2)", "Mn55", 
			"Co59", "Cu63", "I(Cu63^2)", "Rb85", "I(Rb85^2)", "Mo98", "I(Mo98^2)", "Cd111")
morphterms <- c("SR", "I(SR^2)", "Ht", "I(Ht^2)", "Wt", "I(Wt^2)", 
			"Lfnum", "LL", "I(LL^2)", "LW", "I(LW^2)", "GrthDelay", "germday", "hairs")

morph_pres <- sapply(listmorphmod,function(mod) 
	sapply(morphterms, function(z) z%in%rownames(summary(mod)$solutions)) )
nutr_pres <- sapply(listnutrmod,function(mod) 
	sapply(nutrterms, function(z) z%in%rownames(summary(mod)$solutions)) )
morphnames <- traitnamesNL[c(2,2,3,3,4,4,5,6,6,7,7,8,9,10)]
nutrnames <- ionnames[c(1,3,5,6,6,7,8,9,9,10,11,13,13,17,17,19,19,20)]

morph_slo <- sapply(1:length(listmorphmod),function(mod) 
					sapply(1:length(morphterms), function(z) ifelse(morph_pres[z,mod], mean(listmorphmod[[mod]]$Sol[,
						which(rownames(summary(listmorphmod[[mod]])$solutions)==morphterms[z])]), NA ) ) )
morph_lo <- sapply(1:length(listmorphmod),function(mod) 
					sapply(1:length(morphterms), function(z) ifelse(morph_pres[z,mod], HPDi(listmorphmod[[mod]]$Sol[,
						which(rownames(summary(listmorphmod[[mod]])$solutions)==morphterms[z])],prob=0.95)[1], NA ) ) )
morph_up <- sapply(1:length(listmorphmod),function(mod) 
					sapply(1:length(morphterms), function(z) ifelse(morph_pres[z,mod], HPDi(listmorphmod[[mod]]$Sol[,
						which(rownames(summary(listmorphmod[[mod]])$solutions)==morphterms[z])],prob=0.95)[2], NA ) ) )
nutr_slo <- sapply(1:length(listnutrmod),function(mod) 
					sapply(1:length(nutrterms), function(z) ifelse(nutr_pres[z,mod], mean(listnutrmod[[mod]]$Sol[,
						which(rownames(summary(listnutrmod[[mod]])$solutions)==nutrterms[z])]), NA ) ) )
nutr_lo <- sapply(1:length(listnutrmod),function(mod) 
					sapply(1:length(nutrterms), function(z) ifelse(nutr_pres[z,mod], HPDi(listnutrmod[[mod]]$Sol[,
						which(rownames(summary(listnutrmod[[mod]])$solutions)==nutrterms[z])],prob=0.95)[1], NA ) ) )
nutr_up <- sapply(1:length(listnutrmod),function(mod) 
					sapply(1:length(nutrterms), function(z) ifelse(nutr_pres[z,mod], HPDi(listnutrmod[[mod]]$Sol[,
						which(rownames(summary(listnutrmod[[mod]])$solutions)==nutrterms[z])],prob=0.95)[2], NA ) ) )
nutr_slo <- sapply(1:length(listnutrmod),function(mod) 
					sapply(1:length(nutrterms), function(z) ifelse(nutr_pres[z,mod], mean(listnutrmod[[mod]]$Sol[,
						which(rownames(summary(listnutrmod[[mod]])$solutions)==nutrterms[z])]), NA ) ) )

############
##permutations
#######MAY TAKE several hours on a laptop, depending on specific conditions
##permute biomass with respect to all other traits
##refit the best model - one per inocula
set.seed(1)
permslope_nutr <- array(NA,dim=c(19,5,1000))# no. terms in model if all potential traits included, no. inocula, 1000
permslope_morph <- array(NA,dim=c(15,5,1000))#no. terms in model if all potential traits included, no. inocula, 1000
for(inoc in 1:5){
	odat <- list(scaledtrait130,scaledtrait153,scaledtrait186,scaledtrait198,scaledtraitnone)[[inoc]]
	bestmodnutr <- listnutrmod[[inoc]]$Fixed$formula
	bestmodmorph <- listmorphmod[[inoc]]$Fixed$formula
	for(j in 1:1000){#get working before move to 1000
		sampbio <- sample(odat$biomass, size = nrow(odat), repl=F) #permute bio wrt to all others
		sampdat <- odat
		sampdat$biomass <- sampbio
		sampdat <- sampdat[getfull(sampdat),]
		sampmodmorph <- MCMCglmm(bestmodmorph, data = sampdat,verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)
		sampmodnutr <- MCMCglmm(bestmodnutr, data = sampdat,verbose=FALSE,nitt=50000, thin = 10, burnin=100,pr=T)
		permslope_morph[c(1,1+ which(morph_pres[,inoc])),inoc,j] <- summary(sampmodmorph)$solutions[,1] 
		permslope_nutr[c(1,1+ which(nutr_pres[,inoc])),inoc,j] <- summary(sampmodnutr)$solutions[,1] 
	}
}
save(permslope_nutr,file="permslope_nutr.RData")
save(permslope_morph,file="permslope_morph.RData")
############

load("permslope_nutr.RData")
load("permslope_morph.RData")
slopeinterval_nutr <- matrix(NA, nrow=18,ncol=5)
slopeinterval_morph <- matrix(NA, nrow=14,ncol=5)
permslopemn_nutr <- matrix(NA, nrow=18,ncol=5)
permslopemn_morph <- matrix(NA, nrow=14,ncol=5)
isslopesig_nutr <- matrix(NA, nrow=18,ncol=5)
isslopesig_morph <- matrix(NA, nrow=14,ncol=5)
for(trait in 1:18){
		slopeinterval_nutr[trait,] <- sapply(1:ncol(nutr_slo), function(z) findInterval(nutr_slo[trait,z],sort(permslope_nutr[trait+1,z,]))/(dim(permslope_nutr)[3])   )#recall permslope includes intercepts while nutr_slo does not
		permslopemn_nutr[trait,] <- sapply(1:ncol(nutr_slo), function(z)  mean(permslope_nutr[trait+1,z,]) )  
		isslopesig_nutr[trait,] <- (slopeinterval_nutr[trait,] < 0.05 | slopeinterval_nutr[trait,] > 0.95) &  abs(nutr_slo[trait,]) > abs(permslopemn_nutr[trait,]) #real slope needs to be both outside the permutation window AND further from 0
}
for(trait in 1:14){
		slopeinterval_morph[trait,] <- sapply(1:ncol(morph_slo), function(z) findInterval(morph_slo[trait,z],sort(permslope_morph[trait+1,z,]))/(dim(permslope_morph)[3])   )#recall permslope includes intercepts while morph_slo does not
		permslopemn_morph[trait,] <- sapply(1:ncol(morph_slo), function(z)  mean(permslope_morph[trait+1,z,]) )  
		isslopesig_morph[trait,] <- (slopeinterval_morph[trait,] < 0.05 | slopeinterval_morph[trait,] > 0.95) &  abs(morph_slo[trait,]) > abs(permslopemn_morph[trait,]) #real slope needs to be both outside the permutation window AND further from 0
}

###Figure S9
pdf("Trait-fitness-slopes_multitrait_morph.pdf",height=4.5,width=3.75) 
layout(matrix(c(1:12),ncol=3,byrow=TRUE),heights=c(1,1,1,0.75))
par(oma=c(0,4,0.25,0))
par(mar=c(0, 0, 0, 0))
for(i in c(1,3,5)){
	plot(morph_slo[i,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	points(morph_slo[i+1,]~c(1:5 + 0.25), pch = 17, col=fivecols)
	arrows(1:5,y0=morph_lo[i,],y1=morph_up[i,],length=0, col=fivecols)
	arrows(c(1:5 + 0.25),y0=morph_lo[i+1,],y1=morph_up[i+1,],length=0, col=fivecols)
	abline(h=0,lty=3)
	if(i ==1){axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))}
	text(0,-2.8,morphnames[i],adj=0)
	text(1:5,4.1,ifelse(isslopesig_morph[i,], "*",""),col=fivecols)
	text(c(1:5 + 0.25),-2.4,ifelse(isslopesig_morph[i+1,], "*",""),col=fivecols)
}
	plot(morph_slo[7,]~c(1:5),xaxt="n", pch =16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	arrows(c(1:5),y0=morph_lo[7,],y1=morph_up[7,],length=0, col=fivecols)
	axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))
	abline(h=0,lty=3)
	text(0,-2.8,morphnames[7],adj=0)
	text(1:5,4.1,ifelse(isslopesig_morph[7,], "*",""),col=fivecols)
	mtext("Estimated phenotypic selection",side=2,line=2)
for(i in c(8,10)){
	plot(morph_slo[i,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	points(morph_slo[i+1,]~c(1:5 + 0.25), pch = 17, col=fivecols)
	arrows(1:5,y0=morph_lo[i,],y1=morph_up[i,],length=0, col=fivecols)
	arrows(c(1:5 + 0.25),y0=morph_lo[i+1,],y1=morph_up[i+1,],length=0, col=fivecols)
	abline(h=0,lty=3)
	text(1:5,4.1,ifelse(isslopesig_morph[i,], "*",""),col=fivecols)
	text(c(1:5 + 0.25),-2.4,ifelse(isslopesig_morph[i+1,], "*",""),col=fivecols)
	text(0,-2.8,morphnames[i],adj=0)
}
for(i in 12:14){	
	plot(morph_slo[i,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	arrows(1:5,y0=morph_lo[i,],y1=morph_up[i,],length=0, col=fivecols)
	abline(h=0,lty=3)
	if(i ==12){axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))}
	text(0,-2.8,morphnames[i],adj=0)
	text(1:5,4.1,ifelse(isslopesig_morph[i,], "*",""), col = fivecols)
}
plot(1:10~c(1:10),pch=NA,bty="n",xaxt="n",yaxt="n")
legend(0,10,c("13.0","15.3","18.6","19.8","none"), bty="n",fill=fivecols)
	text(4.9,9.75,"Biota Source MAT")
plot(1:10~c(1:10),pch=NA,bty="n",xaxt="n",yaxt="n")
legend(1,9.5,c("linear","quadratic"), bty="n",pch=c(16,17))
	text(5,9.75,"Selection term form")
dev.off()

###Figure S10
pdf("Trait-fitness-slopes_multitrait_nutr.pdf",height=5.75,width=3.75) 
layout(matrix(c(1:15),ncol=3,byrow=TRUE))
par(oma=c(1,4,0.25,0))
par(mar=c(0, 0, 0, 0))
for(i in c(1:3)){
	plot(nutr_slo[i,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	arrows(1:5,y0=nutr_lo[i,],y1=nutr_up[i,],length=0, col=fivecols)
	abline(h=0,lty=3)
	if(i ==1){axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))}
	text(0,-2.8,nutrnames[i],adj=0)
	text(1:5,4.1,ifelse(isslopesig_nutr[i,], "*",""),col=fivecols)
}
	plot(nutr_slo[4,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	points(nutr_slo[4+1,]~c(1:5 + 0.25), pch = 17, col=fivecols)
	arrows(1:5,y0=nutr_lo[4,],y1=nutr_up[4,],length=0, col=fivecols)
	arrows(c(1:5 + 0.25),y0=nutr_lo[4+1,],y1=nutr_up[4+1,],length=0, col=fivecols)
	abline(h=0,lty=3)
	text(0,-2.8,nutrnames[4],adj=0)
	text(1:5,4.1,ifelse(isslopesig_nutr[4,], "*",""),col=fivecols)
	text(c(1:5 + 0.25),-2.4,ifelse(isslopesig_nutr[4+1,], "*",""),col=fivecols)
	axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))
for(i in c(6:7)){
	plot(nutr_slo[i,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	arrows(1:5,y0=nutr_lo[i,],y1=nutr_up[i,],length=0, col=fivecols)
	abline(h=0,lty=3)
	if(i ==1){axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))}
	text(0,-2.8,nutrnames[i],adj=0)
	text(1:5,4.1,ifelse(isslopesig_nutr[i,], "*",""),col=fivecols)
}
	plot(nutr_slo[8,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	points(nutr_slo[8+1,]~c(1:5 + 0.25), pch = 17, col=fivecols)
	arrows(1:5,y0=nutr_lo[8,],y1=nutr_up[8,],length=0, col=fivecols)
	arrows(c(1:5 + 0.25),y0=nutr_lo[8+1,],y1=nutr_up[8+1,],length=0, col=fivecols)
	abline(h=0,lty=3)
	text(0,-2.8,nutrnames[8],adj=0)
	axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))
	text(1:5,4.1,ifelse(isslopesig_nutr[8,], "*",""),col=fivecols)
	text(c(1:5 + 0.25),-2.4,ifelse(isslopesig_nutr[8+1,], "*",""),col=fivecols)
	mtext("Estimated phenotypic selection",side=2,line=2)
for(i in c(10:11)){
	plot(nutr_slo[i,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	arrows(1:5,y0=nutr_lo[i,],y1=nutr_up[i,],length=0, col=fivecols)
	abline(h=0,lty=3)
	if(i ==1){axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))}
	text(0,-2.8,nutrnames[i],adj=0)
	text(1:5,4.1,ifelse(isslopesig_nutr[i,], "*",""),col=fivecols)
}
for(i in c(12,14,16)){
	plot(nutr_slo[i,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	points(nutr_slo[i+1,]~c(1:5 + 0.25), pch = 17, col=fivecols)
	arrows(1:5,y0=nutr_lo[i,],y1=nutr_up[i,],length=0, col=fivecols)
	arrows(c(1:5 + 0.25),y0=nutr_lo[i+1,],y1=nutr_up[i+1,],length=0, col=fivecols)
	abline(h=0,lty=3)
	if(i ==12){axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))}
	text(1:5,4.1,ifelse(isslopesig_nutr[i,], "*",""),col=fivecols)
	text(c(1:5 + 0.25),-2.4,ifelse(isslopesig_nutr[i+1,], "*",""),col=fivecols)
	text(0,-2.8,nutrnames[i],adj=0)
}
	plot(nutr_slo[18,]~c(1:5),xaxt="n", pch = 16, xlim=c(0,5.5), ylim=c(-3,4),yaxt="n", col=fivecols)
	arrows(1:5,y0=nutr_lo[18,],y1=nutr_up[18,],length=0, col=fivecols)
	abline(h=0,lty=3)
	if(i ==1){axis(side=2, at=c(-2,0,2,4), labels=c("-2","0","2","4"))}
	text(1:5,4.1,ifelse(isslopesig_nutr[i,], "*",""),col=fivecols)
	text(0,-2.8,nutrnames[18],adj=0)
plot(1:10~c(1:10),pch=NA,bty="n",xaxt="n",yaxt="n")
legend(0,10,c("13.0","15.3","18.6","19.8","none"), bty="n",fill=fivecols)
	text(4.9,9.75,"Biota Source MAT")
plot(1:10~c(1:10),pch=NA,bty="n",xaxt="n",yaxt="n")
legend(1,9.5,c("linear","quadratic"), bty="n",pch=c(16,17))
	text(5,9.75,"Selection term form")
dev.off()



#####################
####Analyses for SUPPLEMENTAL note
### Similarity of GH plants to field data, nutritional limitation/not of plants in GH, is AMF colonization correlated with plant size in the field  
####################

#FIELD DATA 
#process field ionomics data identically to main experiment data
#note that we treat field data as GH data were treated ---
#but there are some differences bt field and GH over whether logging is better. P: gh y, f medium no. 17th element: gh y f slight no; using GH mixed method, As and Se deviate from normal, P also now under 0.85 at W=0.81
fionomics <- read.csv("~/Variable_teosinte-biota-interactions/August 2013 field leaves ionomics.csv",header=T)
fionomics[which(fionomics$pop==1 & fionomics$individual==20),15:34] <- NA #take out mixed up samples
fsubions1 <- fionomics[,15:34]
fsubions <- fsubions1
#adjustment of the minimum value identical to subions object in greenhouse data processing above
for(i in 1:20){if(any(subions1[,i]<0,na.rm=T)){fsubions[,i] <- fsubions1[,i] + abs(min(subions1[,i],na.rm=T)) + min(subions1[subions1[,i]>0,i],na.rm=T)  }} 
fmixions <- cbind(data.frame(sampleid = fionomics$sample),log(fsubions)) 
fmixions[,which(!logimprovesionA)+1] <- fsubions[,which(!logimprovesionA)] #log identically as in subions (GH data)
#note that this results in NA values because not all negative values in the field samples are brought above zero
##this happens for elements where the GH samples did not have values below 0, but the field samples did, 
###or where the GH samples do have values below 0, but minimum of the field samples is lower than the minimum of the GH samples. 
#

#read in field site trait data things
field <- read.csv("~/Variable_teosinte-biota-interactions/FieldData.csv")
#read in slide data, this is arbusular mycorrhizal colonization data in root intersections. H for hyphae, A arbuscule, V vesicle, AV arbuscule or vesicle, T all/any
slide <- read.csv("~/Variable_teosinte-biota-interactions/root_slides_24Feb14.csv")
slide <- slide[-which(slide$pop==4 & slide$sample==9),]#take out possibly mixed up sample, the others that are duplicated have NA values at the duplicate, and so are not really duplicated.
slide <- slide[!is.na(slide$pop),]
pAV <- (slide$A+slide$V)/(slide$Total*2)
ionpAV <-     sapply(1:nrow(fionomics), function(z) mean(pAV[which(slide$pop==fionomics$pop[z] & slide$sample==fionomics$individual[z])],na.rm=T))
height <-    sapply(1:nrow(fionomics), function(z) field$altura[which(field$pop==fionomics$pop[z] & field$plant==fionomics$individual[z])])
stemwidth <- sapply(1:nrow(fionomics), function(z) field$anchura[which(field$pop==fionomics$pop[z] & field$plant==fionomics$individual[z])])
fionenv <- data.frame(MAT = PopEnvDat$TAnn[fionomics$pop],MAP = PopEnvDat$Pann[fionomics$pop], N = PopEnvDat$InorganicN.ppm[fionomics$pop], P = PopEnvDat$P.Bray.ppm[fionomics$pop], K = PopEnvDat$K.ppm[fionomics$pop])
fpop <- fionomics$pop
##trait differences between field and greenhouse
SW.fg <- data.frame(Ht = c(Plant.Traits$Ht,height),Wt = c(Plant.Traits$Wt,stemwidth), FG = c(rep("g",times=nrow(Plant.Traits)),rep("f",times=length(height))) )
tapply(SW.fg$Ht, SW.fg$FG,mean,na.rm=T); tapply(SW.fg$Ht, SW.fg$FG,std.error)
tapply(SW.fg$Wt, SW.fg$FG,mean,na.rm=T); tapply(SW.fg$Wt, SW.fg$FG,std.error)

#project elemental profiles by multiplying identially scaled (logged and adjusted or not at the same time above, and then using the mean as cc function does automatically) element scores by the CCA ycoefs for each element and summing.
#cc/cancor automatically centers -- so to project any new data, the column means of the original dataset passed to cc have to be subtracted.
# to be projected: fmixions
Iionmeans <- apply(TraitIonInocCom[,1:20],2,mean)

fieldiontrait <- fmixions
fieldiontrait$Wt <- stemwidth
fieldiontrait$Ht <- height
fieldiontrait$pAV <- ionpAV
fieldiontrait$MAT <- fionenv$MAT

###not reported in manuscript. of potential interest
# for(i in c(2:21)){  
# 	tmpdat <- data.frame(x=fieldiontrait[,i],Ht=fieldiontrait$Ht,MAT=fieldiontrait$MAT)
# 	tmpdat <- tmpdat[getfull(tmpdat),]
# 	print(summary(lm(Ht~ x+ MAT+ MAT:x,data=tmpdat)) )
# }
# for(i in c(2:21)){  
# 	tmpdat <- data.frame(y=fieldiontrait[,i],MAT=fieldiontrait$MAT,pAV=fieldiontrait$pAV)
# 	tmpdat <- tmpdat[getfull(tmpdat),]
# 	print(summary(lm(y~  MAT*pAV,data=tmpdat)) )
# }


#note fpop number is in order of sampling: TZ FP MT DA MC ML TC TX CL CU  
tapply(ionpAV,fionenv$MAT/10,mean,na.rm=T)
range(ionpAV,na.rm=T)
ftraitenv <- data.frame(height=height, stemwidth = stemwidth, HW = height/stemwidth, MAT=fionenv$MAT/10)
ftraitenv.col <- data.frame(height=height, stemwidth = stemwidth, HW = height/stemwidth, MAT=fionenv$MAT/10,pAV = ionpAV)
fte <- ftraitenv[getfull(ftraitenv),]
fte.col <- ftraitenv.col[getfull(ftraitenv.col),]
T.sf <- seq(from = min(fionenv$MAT),to=max(fionenv$MAT),length.out=1000)

###Model fitting. Reported in Table S4
#stem width
summary(MCMCglmm(stemwidth~MAT*pAV,data = fte.col,verbose=FALSE,nitt=100000, thin = 100, burnin=10000))#stop, mat-pav interaction is sig and neg, pav main sig pos. mat main nonsig pos.
#height
summary(MCMCglmm(height ~MAT*pAV,data = fte.col,verbose=FALSE,nitt=100000, thin = 100, burnin=10000))#best model, same as width: pav, and interaction terms sig, pav pos, int neg, main eff MAT ns and pos.
pAV.s <- seq(from= min(ionpAV,na.rm=T),to=max(ionpAV,na.rm=T),length.out=1000)
fp.w.AVxmod  <- MCMCglmm(stemwidth ~MAT*pAV, 		data = fte.col,verbose=FALSE,nitt=1000000, thin = 100, burnin=100000)
wpost <- fp.w.AVxmod$Sol
wpred.tmax.ci <- sapply(1:length(pAV.s), 	function(z) HPDi(wpost[,1] + wpost[,2]*max(fte.col$MAT) + wpost[,3]*pAV.s[z] + wpost[,4]*pAV.s[z]*max(fte.col$MAT),.95) ) 
wpred.tmax.mean <- sapply(1:length(pAV.s), function(z) mean(wpost[,1] + wpost[,2]*max(fte.col$MAT) + wpost[,3]*pAV.s[z] + wpost[,4]*pAV.s[z]*max(fte.col$MAT)) ) 
wpred.tmin.ci <- sapply(1:length(pAV.s), 	function(z) HPDi(wpost[,1] + wpost[,2]*min(fte.col$MAT) + wpost[,3]*pAV.s[z] + wpost[,4]*pAV.s[z]*min(fte.col$MAT),.95) ) 
wpred.tmin.mean <- sapply(1:length(pAV.s), function(z) mean(wpost[,1] + wpost[,2]*min(fte.col$MAT) + wpost[,3]*pAV.s[z] + wpost[,4]*pAV.s[z]*min(fte.col$MAT)) ) 
#
fp.h.AVxmod	<- MCMCglmm(height ~MAT*pAV, 		data = fte.col,verbose=FALSE,nitt=1000000, thin = 100, burnin=100000)#best model, same as width: pav, and interaction terms sig, pav pos, int neg, main eff MAT ns and pos.
hpost <- fp.h.AVxmod$Sol
hpred.tmax.ci <- sapply(1:length(pAV.s), 	function(z) HPDi(hpost[,1] + hpost[,2]*max(fte.col$MAT) + hpost[,3]*pAV.s[z] + hpost[,4]*pAV.s[z]*max(fte.col$MAT),.95) ) 
hpred.tmax.mean <- sapply(1:length(pAV.s), function(z) mean(hpost[,1] + hpost[,2]*max(fte.col$MAT) + hpost[,3]*pAV.s[z] + hpost[,4]*pAV.s[z]*max(fte.col$MAT)) ) 
hpred.tmin.ci <- sapply(1:length(pAV.s), 	function(z) HPDi(hpost[,1] + hpost[,2]*min(fte.col$MAT) + hpost[,3]*pAV.s[z] + hpost[,4]*pAV.s[z]*min(fte.col$MAT),.95) ) 
hpred.tmin.mean <- sapply(1:length(pAV.s), function(z) mean(hpost[,1] + hpost[,2]*min(fte.col$MAT) + hpost[,3]*pAV.s[z] + hpost[,4]*pAV.s[z]*min(fte.col$MAT)) ) 

#Figure S13
pdf("Field_Mod_Data.pdf",height=4.5,width=4)
layout(matrix(c(1,2,3,3),ncol=2,byrow=F),widths=c(1,0.75))
par(oma=c(4,3,1,0))
par(mar=c(0,1,0,0.5))
plot(stemwidth~pAV,data=fte.col,pch=NA,ylab= "",xlab="",xaxt="n",ylim=c(-3,38))#,main="Field Macrotraits")
	polygon(c(pAV.s,rev(pAV.s)), c(wpred.tmin.ci[1,],rev(wpred.tmin.ci[2,])), col=rgb(0,0,1,alpha=.5), border = NA)
	polygon(c(pAV.s,rev(pAV.s)), c(wpred.tmax.ci[1,],rev(wpred.tmax.ci[2,])), col=rgb(1,0,0,alpha=.5), border = NA)
	lines(wpred.tmax.mean~pAV.s,,col=rgb(1,0,0))
	lines(wpred.tmin.mean~pAV.s,col=rgb(0,0,1))
	mtext("Stem width (mm)",side=2,line=2)
	points(stemwidth~pAV,col=rgb(range01(MAT),0,1-range01(MAT)),data=fte.col,pch=16,cex=1)#,main="Field Macrotraits")
par(mar=c(0,1,0,0.5))
plot(height~pAV,data=fte.col,pch=NA,ylab="",main="",ylim=c(-35,325),xlab="")
	polygon(c(pAV.s,rev(pAV.s)), c(hpred.tmin.ci[1,],rev(hpred.tmin.ci[2,])), col=rgb(0,0,1,alpha=.5), border = NA)
	polygon(c(pAV.s,rev(pAV.s)), c(hpred.tmax.ci[1,],rev(hpred.tmax.ci[2,])), col=rgb(1,0,0,alpha=.5), border = NA)
	lines(hpred.tmax.mean~pAV.s,,col=rgb(1,0,0))
	lines(hpred.tmin.mean~pAV.s,col=rgb(0,0,1))
	mtext("Height (cm)",side=2,line=2)
	mtext("Proportion AMF Colonized",side=1,line=2.5)
	points(height~pAV,col=rgb(range01(MAT),0,1-range01(MAT)),data=fte.col,pch=16,cex=1)
par(mar=c(5,2.5,2.5,5))
image(as.matrix(t(tann)), col=rev(rb(100)),zlim=c(min(tann),max(tann)),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side=4,at=seq(from=0, to=1,length.out=10),labels=tann,las=2)
	mtext("Site MAT",side = 1, line=1,at=1.5)
dev.off()

#####################
###supplemental figures for nutrients 
####################

#read in deficiency info
deftox <- read.csv("~/Variable_teosinte-biota-interactions/marschnerDefToxConcentrations.csv",stringsAsFactors=F)
# sort same as element data
deftox2 <- matrix(NA,ncol=20,nrow=4)
colnames(deftox2) <- colnames(fsubions1)
whichdft <-  sapply(deftox[,1], function(z) which(colnames(fsubions1)==z) ) 
for(i in 1:length(whichdft)){deftox2[,(whichdft[i])] <- unlist(deftox[i,2:5])}

#make table of treatment level means.
uionmns<- sapply(2:ncol(subionsfull), function(z) tapply(subionsfull[isinoc==0,z],PlantSourceEnvDat$TAnn[isinoc==0],mean,na.rm=T))
uionses<- sapply(2:ncol(subionsfull), function(z) tapply(subionsfull[isinoc==0,z],PlantSourceEnvDat$TAnn[isinoc==0],std.error,na.rm=T))
ghionmns <- sapply(2:21, function(i) 
				sapply(1:10, function(pt) 
					sapply(1:10, function(pi)  
						mean(subionsfull[PlantSourceEnvDat$TAnn==tann[pt]*10 & InocSourceEnvDat$TAnn==tann[pi]*10,i],na.rm=T) )))
tps <- rep(tann, each =10)
tis <- rep(tann, times = 10)

#plot, supplemental figure S11
elmranges <- rbind( ghionmns, uionmns, deftox2[1,])
pdf("~/Variable_teosinte-biota-interactions/Elements_LIM.pdf",height=8,width=8)
layout(matrix(c(1:20,21,21,22,22), ncol=4, byrow=T),heights=c(2,2,2,2,2,1))
par(mar=c(3,1,1,1))
par(oma=c(1,3,1,0))
for(i in 1:20){
	plot(ghionmns[,i]~tps,ylab="",xlab="",#
		ylim=bufferX(range(elmranges[,i],na.rm=T),p=0.1),pch=1,col=rgb(range01(tis),0,1-range01(tis)) )
 	points(uionmns[,i]~tann,pch=1,col=rgb(0.5,0.5,0.5))
	mtext(ionnames[i],side=3)
	abline(h=(deftox2[1,i]),lty=1, lwd=2)
	abline(h=(deftox2[2,i]),lty=2,lwd=2)
if(i==4){legend(14,100,c("Deficient","Possibly deficient"),lty=c(1,2),bty="n")}
if(i==9 ){mtext("micrograms per gram dry weight",side=2,line=2.5)}
if(i==18){mtext("Plant source site mean annual temperature",side=1,line=2.5,at=21)}
}
par(mar=c(0,0,3,0))
plot(1:3~c(1:3),bty="n",pch=NA,xaxt="n",yaxt="n",xlab="",ylab="")
	legend(1,3.5,"Uninoculated",fill=rgb(0.5,0.5,0.5),bty="n",cex=1.5)
par(mar=c(2.5,0,1.5,1))
image(as.matrix(tann), col=rev(rb(100)),zlim=c(min(tann),max(tann)),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side=1,at=seq(from=0, to=1,length.out=10),labels=tann)
	mtext("Biota source site mean annual temperature",side = 1, line=2.25)
dev.off()


#supplemental figure S12
gfmns <- rbind(sapply(2:21, function(i) mean(subionsfull[,i],na.rm=T)),
			sapply(1:20, function(i) mean(fsubions1[,i],na.rm=T) ))			
gfses <- rbind(sapply(2:21, function(i) std.error(subionsfull[,i])),
			sapply(1:20, function(i) std.error(fsubions1[,i]) ))
gfrange <- rbind(gfmns+gfses,gfmns-gfses,deftox2[1,])
pdf("~/Variable_teosinte-biota-interactions/Field_elmt.pdf",height=5,width=5)
layout(matrix(c(1:20,21,21,21,21),ncol=4,byrow=T),)
par(oma=c(0,1,0,0))
par(mar=c(1,1,2,1))
for(i in 1:20){
	plot(c(gfmns[,i])~c(1,2),xlim=c(0.5,2.5),ylim=bufferX(range(gfrange[,i],na.rm=T),p=0.1),
	 ylab="",xlab="",bty="l",xaxt="n",col=c("black",rgb(0.85,0.75,0)),pch=1)
	 arrows(c(1:2),gfmns[,i]-gfses[,i], y1=gfmns[,i]+gfses[,i] ,col=c("black",rgb(0.85,0.75,0)),length=0,lwd=2)
	abline(h=(deftox2[1,i]),lty=1, lwd=2)
	abline(h=(deftox2[2,i]),lty=2,lwd=2)
	mtext(ionnames[i],side=3,line=0)
}
par(mar=c(0,0,0,0))
plot(1:10~c(1:10),xlim=c(1,10),xaxt="n",yaxt="n",bty="n",pch=NA)
legend(2,10,c("Inoculated greenhouse plants","Field plants"),fill=c(rgb(0,0,0,alpha=1),rgb(0.85,0.75,0,alpha=1)),bty="n")
dev.off()
