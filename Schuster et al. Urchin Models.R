#### Schuster et al. Tropicalization of temperate reef fish communities facilitated 
#    by urchin grazing and diversity of thermal affinities
 
library('ggplot2'); library('nlme'); library('mgcv'); library('lmeSplines'); library('maps') 
library('dplyr'); library('pspline')

### Summarizing global RLS data; extracting urchin sites and fish data
dataM1<-read.csv("RLS_Data_ALL_M1.csv",header=TRUE) # Pelagic fish data RLS
dataM2<-read.csv("RLS_Data_ALL_M2.csv",header=TRUE) # Invertebrate + cryptic fish data RLS
all_data<-rbind(dataM1,dataM2)

# Extract global sea urchin observations
Urchin_Sites<-all_data %>% filter(CLASS == "Echinoidea") %>% filter(ORDER != "Clypeasteroida") %>% filter(ORDER != "")

# Extract global fish observations (for thermal ranges of fishes)
ALLSPP<-all_data %>% filter(PHYLUM=="Chordata") %>% droplevels() %>% filter(CLASS!='Ascidiacea',CLASS!='Aves',CLASS!='Mammalia',CLASS!='Reptilia') %>% filter(SPECIES_NAME!="No species found") %>% droplevels() # select fishes only

### EXTRACT TEMPs (SST) AT EACH GLOBAL SURVEY SITE
Site_Metadata<-aggregate(TotalAbundance~SiteCode+ECOregion+SiteLat+SiteLong+SurveyDate+SurveyID,data=all_data,FUN=sum) #aggregate species occurances to survey site level (SurveyID)

####### Weekly SST Temperature ########
# EXTRACT TEMPERATURES AT EACH RLS SURVEY SITE

weekly_sst<-open.nc("sst_wkmean_1981-March2019.nc")
sst<-read.nc(weekly_sst,unpack=T)
# Vector of each variable
long<-data.frame(sst[3])
lat<-data.frame(sst[2])
time<-data.frame(sst[1])

##### MAKE GRID OF ALL POSSIBLE LAT AND LONG COMBINATIONS. REMOVE -180 from Longs 
all.lat.long1<-expand.grid(lat[,1],long[,1])
all.lat.long<-all.lat.long1
all.lat.long[,2]<-all.lat.long1[,2]-180

## Make sst array
sst.array<-abind(sst[4])

SP<-SpatialPoints(cbind(Site_Metadata$SiteLong,Site_Metadata$SiteLat))
crs(SP)<- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

##### DATE CONVERSION FROM DATE TO JULIAN (RLS) (weekly till 2019 data)
write.xlsx(Site_Metadata,file="Site_Metadata_RLS.xlsx",sheetName="RLS All unique Sites",col.names=TRUE)

## Split SurveyDate column in excel to remove time; add Origin column (1981-11-08 = start date of weekly sst)
# note: date_diff calculation does not include end date (if wanted, add 1 day)
Site_Metadata<-read.csv("Site_Metadata_RLS.csv")
Site_Metadata$date_diff <- as.Date(as.character(Site_Metadata$SurveyDate), format="%Y-%m-%d")-
  as.Date(as.character(Site_Metadata$Origin), format="%Y-%m-%d")
Site_Metadata$Julian<-as.numeric(Site_Metadata$date_diff + 2444917)

##### CREATE BLANK MATRIX
temp.matrix<-matrix(0,dim(Site_Metadata)[1],4)
colnames(temp.matrix)<-c("Survey.T.Match","Mean.Survey.T.4wkPrior","Site.T.Min","Site.T.Max")

###### MATCHING RLS SURVEY SITES W/ TEMP DATA BY TIME POINTS
for (i in 1:dim(Site_Metadata)[1]){
  
  diff.times<-abs(Site_Metadata$Julian[i]-time[,1]) #calculates difference in time difference between RLS days and sst days
  time.point<-which(diff.times==min(diff.times)) #chooses time point where time diff is minimal (closest to survey date)
  time.point4wkprior<-time.point-3 #chooses time point for 4 weeks before survey
  time.point1yrprior<-time.point-52 #choose time point 1 year before survey
  #### EUCLIDEAN DISTANCE 
  values<-sqrt(((Site_Metadata[i,]$SiteLat+90)-(all.lat.long[,1]+90))^2+((Site_Metadata[i,]$SiteLong+180)-(all.lat.long[,2]+180))^2)
  lat.long<-all.lat.long[which(values==min(values)),]
  ### second [1] "first element" is added here because there's two identitical spatial points for lat.long[,1]
  lat1<-which(lat[,1]==lat.long[,1][1])
  long1<-which(long[,1]-180==lat.long[,2][1])
  temp.matrix[i,1]<-sst.array[long1,lat1,time.point]
  temp.matrix[i,2]<-mean(sst.array[long1,lat1,time.point4wkprior:time.point],na.omit=T)
  temp.matrix[i,3]<-min(sst.array[long1,lat1,time.point:time.point1yrprior]) # Min temp in 1 year prior survey period
  temp.matrix[i,4]<-max(sst.array[long1,lat1,time.point:time.point1yrprior]) # max temp in 1 year prior 
  
}

Site_Temps<-cbind(Site_Metadata,temp.matrix) # merge temp data with survey site data
#### use this snipped to remove no-longer-needed columns (Origin, days since origin); remove TotalABundance from this, as it's total abundance per SiteCode, not SurveyID
Site_Temps<-Site_Temps[,-7]

# Merge Site_Temps back to fish and urchin community data for next analysis steps
UrchinSiteAbundance<-left_join(Urchin_Sites,Site_Temps,by=c("SurveyID"))
# Calculate urchin site abundance & add barren-former score (barren forming Y/N)
UrchinSiteAbundance<-aggregate(TotalAbundance~SiteCode+SurveyID+SurveyDate+ECOregion+GENUS+SPECIES_NAME+Survey.T.Match+SiteLat,data=UrchinSiteAbundance,FUN=sum)
Urchin_Scores<-read.csv("Urchin_Scores_GE.csv") # barren former scores (Y/N); see supplementary Table S2 for species scores
UrchinSiteAbundance<-left_join(UrchinSiteAbundance,Urchin_Scores,by="SPECIES_NAME")
UrchinSiteAbundance<-UrchinSiteAbundance %>% filter(Barren_forming=="Y") # remove urchin species that don't form barrens
# Survey site-level abundance of sea urchins (i.e. total number of barren-forming urchins, summed across species)
UrchinSiteAbundance<-plyr::ddply(UrchinSiteAbundance,c("SiteCode","SurveyID","SurveyDate",'GENUS'),summarize,Barren_former_AB=sum(TotalAbundance))

# Merge Site temp data with global species (fishes) observations
ALLSPP<-left_join(ALLSPP,Site_Temps,by=c("SurveyID"))

### EXTRACT GBIF OCCURRENCES 
# FishSpecies is a list of fish species names extracted from Fish_Communities (i.e., fishes occurring at included sites)
library(spocc)

dfall <- occ(query = FishSpecies, from = c('gbif'),gbifopts=list(hasCoordinate=TRUE)) #extract GBIF fish occurrences for all fish species recorded in RLS data
SPP_OCCURall<-occ2df(dfall) #makes df with most important info

dfall2 <- fixnames(dfall, how = 'query') # returns original query name (RLS spp name) to each record
SPP_OCCURall2<-occ2df(dfall2) 
SPP_OCCURall2<-SPP_OCCURall2[-c(2:5)]
names(SPP_OCCURall2)[names(SPP_OCCURall2) == "name"] <- "SPECIES_NAME" #rename to match with RLS species column name

SPP_OCCall<-left_join(SPP_OCCURall,SPP_OCCURall2,by="key") # joins GBIF spp names with RLS spp names
write.csv(SPP_OCCall,file="FishGBIF.csv",col.names=TRUE)

FishGBIF <- FishGBIF %>%
  filter(year > 1983) %>% filter(is.na(depth)| depth<20) %>%
  filter(!is.na(longitude)) %>%
  filter(!is.na(latitude)) # remove records from before temp record starts, non coastal records & records without spatial info

####### Weekly SST Temperature ########
# EXTRACT TEMPERATURES AT EACH GBIF OCCURRENCE SITE 
# same procedure as above, but using GBIF fish occurrence locations

weekly_sst<-open.nc("sst_wkmean_1981-March2019.nc")
sst<-read.nc(weekly_sst,unpack=T)
# Vector of each variable
long<-data.frame(sst[3])
lat<-data.frame(sst[2])
time<-data.frame(sst[1])

##### MAKE GRID OF ALL POSSIBLE LAT AND LONG COMBINATIONS. REMOVE -180 from Longs 
all.lat.long1<-expand.grid(lat[,1],long[,1])
all.lat.long<-all.lat.long1
all.lat.long[,2]<-all.lat.long1[,2]-180

## Make sst array
sst.array<-abind(sst[4])

SP<-SpatialPoints(cbind(FishGBIF$longitude,FishGBIF$latitude))
crs(SP)<- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

##### DATE CONVERSION FROM DATE TO JULIAN (RLS) (weekly till 2019 data)
## Split SurveyDate column in excel to remove time; add Origin (1981-11-08 = start date of weekly sst)
# note: date_diff calculation does not include end date (if wanted, add 1 day)
setwd("/Users/jmschuster/Desktop/PhD/Herbivory Paper/Urchin Analysis")
FishGBIF$date_diff <- as.Date(as.character(FishGBIF$eventDate), format="%Y-%m-%d")-
  as.Date(as.character(FishGBIF$Origin), format="%Y-%m-%d")
FishGBIF$Julian<-as.numeric(FishGBIF$date_diff + 2444917)

##### CREATE BLANK MATRIX
temp.matrix<-matrix(0,dim(FishGBIF)[1],4)
colnames(temp.matrix)<-c("Survey.T.Match","Mean.Survey.T.4wkPrior","Site.T.Min","Site.T.Max")

###### MATCHING OCCURRENCE LOCATIONS W/ TEMP DATA BY TIME POINTS
for (i in 1:dim(FishGBIF)[1]){
  
  diff.times<-abs(FishGBIF$Julian[i]-time[,1]) #calculates difference in time difference between RLS days and sst days
  time.point<-which(diff.times==min(diff.times)) #chooses time point where time diff is minimal (closest to survey date)
  time.point4wkprior<-time.point-3 #chooses time point for 4 weeks before survey
  time.point1yrprior<-time.point-52 #choose time point 1 year before survey
  #### EUCLIDEAN DISTANCE 
  values<-sqrt(((FishGBIF[i,]$latitude+90)-(all.lat.long[,1]+90))^2+((FishGBIF[i,]$longitude+180)-(all.lat.long[,2]+180))^2)
  lat.long<-all.lat.long[which(values==min(values)),]
  ### second [1] "first element" is added here because there's two identitical spatial points for lat.long[,1]
  lat1<-which(lat[,1]==lat.long[,1][1])
  long1<-which(long[,1]-180==lat.long[,2][1])
  temp.matrix[i,1]<-sst.array[long1,lat1,time.point]
  temp.matrix[i,2]<-mean(sst.array[long1,lat1,time.point4wkprior:time.point],na.omit=T)
  temp.matrix[i,3]<-min(sst.array[long1,lat1,time.point:time.point1yrprior]) # Min temp in 1 year prior survey period
  temp.matrix[i,4]<-max(sst.array[long1,lat1,time.point:time.point1yrprior]) # max temp in 1 year prior 
  
}

GBIF_TempsFish<-cbind(FishGBIF,temp.matrix)
#### use this snipped to remove no-longer-needed columns (Origin, days since origin)
GBIF_TempsFish<-GBIF_TempsFish[,-14]
write.csv(GBIF_TempsFish,file="GBIF_Temps_Fish.csv",col.names=TRUE)

### MERGE GBIF AND RLS OCCURRENCE DATA
GBIF_TempsFish<-GBIF_TempsFish %>% dplyr::select(SPECIES_NAME,latitude,longitude,Survey.T.Match,Site.T.Min,Site.T.Max) #global occurrences of species based on GBIF
colnames(GBIF_TempsFish)[which(names(GBIF_TempsFish) == "longitude")] <- "SiteLong" #rename to match RLS col name
colnames(GBIF_TempsFish)[which(names(GBIF_TempsFish) == "latitude")] <- "SiteLat" #rename to match RLS col name
GBIF_TempsFish$Record <- rep('GBIF',nrow(GBIF_TempsFish)) 

RLS_Temps<-ALLSPP %>% dplyr::select(SPECIES_NAME,SiteLat,SiteLong,Survey.T.Match,Site.T.Min,Site.T.Max) #global occurrences of species based on RLS
RLS_Temps$Record <- rep('RLS',nrow(RLS_Temps)) 
RLS_GBIF_OCC<-rbind(GBIF_TempsFish,RLS_Temps) # combine RLS and GBIF occurrence data

### THERMAL RANGES OF GBIF+RLS SPECIES OCCURRENCE
# Summarize minimum and maximum temperature experienced across geographic range for each species
RLS_GBIF_T_Ranges<-plyr::ddply(RLS_GBIF_OCC,c("SPECIES_NAME"),summarize,Min_T=quantile(Site.T.Min,probs = 0.05,na.rm=TRUE),Max_T=quantile(Site.T.Max,probs = 0.95,na.rm=TRUE))
# Calculate Species Thermal Index (STI)
RLS_GBIF_T_Ranges$T_MidPoint<-rowMeans(RLS_GBIF_T_Ranges[c("Min_T","Max_T")])
# Calculate Species Thermal Range (STR)
RLS_GBIF_T_Ranges$T_Range<-(RLS_GBIF_T_Ranges$Max_T-RLS_GBIF_T_Ranges$Min_T)
write.xlsx(RLS_GBIF_T_Ranges,file="RLS_GBIF_T_Ranges.xlsx",col.names=TRUE) # this was merged back with survey metadata to create spp-level dfs (FishSiteAbundance & FishAffinities)

####### FISH DATA (Species-level and community level)
### SPECIES DATA
# Subset fish data from survey sites where sea urchins are present (for fish community analysis at urchin sites)
Fish_Communities<-ALLSPP %>% filter(SurveyID %in% UrchinSiteAbundance$SurveyID) # subset fish survey data to relevant urchin sites

# Calculate total site abundance for each fish species:
# raw data has two separate blocks (left and right of transect) for each survey; i.e. sum across blocks for total abundance per 500m^2 fish survey transect
FishSiteAbundance<-aggregate(TotalAbundance~SiteCode+SurveyID+SurveyDate+SiteLat+SiteLong+ECOregion+realm+Country+Location+Depth+SPECIES_NAME+Site.T.Min+Site.T.Max,data=Fish_Communities,FUN=sum)
# Merge with Thermal range data for each spp
FishSiteAbundance<-left_join(FishSiteAbundance, RLS_GBIF_T_Ranges,by="SPECIES_NAME")

# Cleaning up species data frame (remove unidentified spp)
library(stringr)
FishSiteAbundance<- FishSiteAbundance %>% filter(str_detect(SPECIES_NAME, "spp.") == FALSE) #removes any record of unidentified species
FishSiteAbundance<-FishSiteAbundance %>% filter(str_detect(SPECIES_NAME, "sp.") == FALSE) #removes any record of unidentified species
FishSiteAbundance<-FishSiteAbundance %>% filter(str_detect(SPECIES_NAME, "Unidentified") == FALSE) #removes any record of unidentified species

# Limit to temperate latitudes
FishSiteAbundance$absSiteLat<-abs(FishSiteAbundance$SiteLat) # absolute latitude
FishSiteAbundance$r.SiteLat<-round(FishSiteAbundance$SiteLat) # rounded latitude
FishSiteAbundance<-subset(FishSiteAbundance,absSiteLat> 23.5 & absSiteLat< 50)
FishSiteAbundance$SurveyID<-as.factor(FishSiteAbundance$SurveyID)

### COMMUNITY DATA (for each survey site)
# First, calculate total abundance of fishes (all fish species summed) per survey site
FINAL_URCHIN_DF<-aggregate(TotalAbundance~SiteCode+SurveyID+SurveyDate+SiteLat+SiteLong+ECOregion+absSiteLat+r.SiteLat+realm+Country+Location+Depth+Site.T.Min+Site.T.Max,data=FishSiteAbundance,FUN=sum)
FINAL_URCHIN_DF<-left_join(FINAL_URCHIN_DF,UrchinSiteAbundance,by='SurveyID') # add sea urchin site abundance (density)

# Urchin density score (barren if >2 urchins per m^2, kelp if <2 urchins per m^2; i.e. threshold is 200 urchins per 100m^2 transect area)
FINAL_URCHIN_DF$Site_Score200[FINAL_URCHIN_DF$Barren_former_AB>=200]<-'High' #barren
FINAL_URCHIN_DF$Site_Score200[FINAL_URCHIN_DF$Barren_former_AB<200]<-'Low'  #kelp

FINAL_URCHIN_DF$Site_Score200<-
  relevel(FINAL_URCHIN_DF$Site_Score200,ref="Low") # relevel site score so 'low' (kelp) is first level

##### HABITAT CLASSIFICATION (PQs + urchin density) 

### FINAL_URCHIN_DF uses the highest quality data, where barren/kelp habitat status is varified with PQ data
### FULL_URCHIN_DATA is full data set of temperate barren/kelp sites, but habitat status classified from urchin density only (2 per m^2 threshold; see methods)

# PQ scoring completed externally, final score = COMBINED_SCORE which classifies habitat based on:
# Kelp = barren former urchin density <2 per m^2 and kelp cover > 50%
# Barren = barren former urchin density >2 per m^2 and kelp cover < 50%
PQ_scored<-read.csv('PQ_scored.csv')
PQ_scored<-PQ_scored %>% filter(COMBINED_SCORE!='Neither') %>% droplevels() # remove sites that are neither barren nor kelp
PQ_scored<-PQ_scored %>% filter(COMBINED_SCORE!='Uncertain') %>% droplevels() # remove sites with fuzzy classification

COMB_SCORE<-PQ_scored %>% select(SurveyID,COMBINED_SCORE)
COMB_SCORE$COMBINED_SCORE<-as.factor(COMB_SCORE$COMBINED_SCORE)
COMB_SCORE$SurveyID<-as.factor(COMB_SCORE$SurveyID)
COMB_SCORE$COMBINED_SCORE<-
  relevel(COMB_SCORE$COMBINED_SCORE,ref="Kelp")

FINAL_URCHIN_DF<-left_join(FINAL_URCHIN_DF,COMB_SCORE,by='SurveyID')

# remove ecoregions with <5 total surveys & ecoregions that only have kelp sites but no barrens, or vice versa)
FINAL_URCHIN_DF %>% group_by(ECOregion,COMBINED_SCORE) %>% tally()
FINAL_URCHIN_DF<-FINAL_URCHIN_DF %>% filter(ECOregion!="Puget Trough/Georgia Basin" & 
                                              ECOregion!='Oyashio Current' & ECOregion!='Gulf of Maine/Bay of Fundy' &
                                              ECOregion!='Sea of Japan/East Sea'& ECOregion!='Leeuwin'& 
                                              ECOregion!='Houtman'& ECOregion!='Shark Bay') %>% droplevels()


####### COMMUNITY THERMAL COMPOSITION ############
# Calculate CTI, CTDiv and CTR (=AvrgBreadth) of fish communities
CTI_Fish<-plyr::ddply(FishSiteAbundance,c("SurveyID"),summarize,CTI=mean(T_MidPoint),CTDiv=sd(T_MidPoint),AvrgBreadth=mean(T_Range),CTImax=mean(Max_T))
FULL_URCHIN_DATA<-left_join(FULL_URCHIN_DATA,CTI_Fish,by='SurveyID') # add community thermal composition data to full community df
FINAL_URCHIN_DF<-left_join(FINAL_URCHIN_DF,CTI_Fish,by='SurveyID') # add community thermal composition data to highest quality community df


##################################
#### Simulation of CTI change under different scenarios of CTI sensitivity and habitat change
### Calculate thermal diversity of ecoregional species pools
ecoregion_pool<- aggregate(TotalAbundance~ECOregion+SPECIES_NAME+Class+T_MidPoint+T_Range,data=FishSiteAbundance,FUN=sum)
Fish_ecoregionPool<-ecoregion_pool %>% group_by(ECOregion) %>% dplyr::summarize(sdSTI=sd(T_MidPoint),STR=mean(T_Range),sdSTR=sd(T_Range))
FULL_URCHIN_DATA<-left_join(FULL_URCHIN_DATA,Fish_ecoregion_Pool,by='ECOregion')

m1sim=gamm(data=FULL_URCHIN_DATA,CTImax_fish_ECS~s(Depth)+STR+sdSTI+Survey.T.Match*Site_Score200*STR+Survey.T.Match*Site_Score200*sdSTI,random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit)
summary(m1sim$lme)

# Scenario A: Predict for average depth, narrow STR, high response diversity (sdSTI): predict high CTI sensitivity
# This is plot pHigh
newdat=data.frame(Survey.T.Match=seq(8,28,length.out=1000),Site_Score200="High",Depth=10,STR=5,sdSTI=4)
newdat<-cbind(newdat,predict(m1sim$gam, newdat, se.fit=TRUE))
newdat<-transform(newdat,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals
newdat2=data.frame(Survey.T.Match=seq(8,28,length.out=1000),Site_Score200="Low",Depth=10,STR=5,sdSTI=4)
newdat2<-cbind(newdat2,predict(m1sim$gam, newdat2, se.fit=TRUE))
newdat2<-transform(newdat2,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals

# Plot scenario A
## ggplot without legend 
pHigh<-ggplot(NULL,aes(x=Survey.T.Match,y=fit))+xlim(c(8,28))+ylim(c(-3,3))+
  geom_ribbon(data=newdat,aes(ymin=lwr,ymax=upr),fill='mediumslateblue',alpha=0.2)+geom_ribbon(data=newdat2,aes(ymin=lwr,ymax=upr),fill='lawngreen',alpha=0.2)+
  geom_line(data=newdat,colour='#330066',linetype='dotted',size=1.4)+geom_line(data=newdat2,colour='#669933',size=1.4)+xlab('Temperature °C')+ylab(expression(Delta*'CTI'))+
  theme_bw()+theme(panel.background = element_blank(),panel.grid=element_blank(),legend.position='none')

# Scenario B: Predict for average depth, narrow STR, low response diversity (sdSTI): predict med CTI sensitivity
# this is plot pMed1
newdat=data.frame(Survey.T.Match=seq(8,28,length.out=1000),Site_Score200="High",Depth=10,STR=5,sdSTI=1)
newdat<-cbind(newdat,predict(m1sim$gam, newdat, se.fit=TRUE))
newdat<-transform(newdat,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals
newdat2=data.frame(Survey.T.Match=seq(8,28,length.out=1000),Site_Score200="Low",Depth=10,STR=5,sdSTI=1)
newdat2<-cbind(newdat2,predict(m1sim$gam, newdat2, se.fit=TRUE))
newdat2<-transform(newdat2,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals

# Plot scenario B
pMed1<-ggplot(NULL,aes(x=Survey.T.Match,y=fit))+xlim(c(8,28))+ylim(c(-3,3))+
  geom_ribbon(data=newdat,aes(ymin=lwr,ymax=upr),fill='mediumslateblue',alpha=0.2)+geom_ribbon(data=newdat2,aes(ymin=lwr,ymax=upr),fill='lawngreen',alpha=0.2)+
  geom_line(data=newdat,colour='#330066',linetype='dotted',size=1.4)+geom_line(data=newdat2,colour='#669933',size=1.4)+xlab('Temperature °C')+ylab(expression(Delta*'CTI'))+
  theme_bw()+theme(panel.background = element_blank(),panel.grid=element_blank(),legend.position='none')

# Scenario C: Predict for average depth, wide STR, high response diversity (sdSTI): predict med-high CTI sensitivity
# this is plot pMed2
newdat=data.frame(Survey.T.Match=seq(8,28,length.out=1000),Site_Score200="High",Depth=10,STR=19,sdSTI=4)
newdat<-cbind(newdat,predict(m1sim$gam, newdat, se.fit=TRUE))
newdat<-transform(newdat,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals
newdat2=data.frame(Survey.T.Match=seq(8,28,length.out=1000),Site_Score200="Low",Depth=10,STR=19,sdSTI=4)
newdat2<-cbind(newdat2,predict(m1sim$gam, newdat2, se.fit=TRUE))
newdat2<-transform(newdat2,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals

# Plot scenario C
pMed2<-ggplot(NULL,aes(x=Survey.T.Match,y=fit))+xlim(c(8,28))+ylim(c(-3,3))+
  geom_ribbon(data=newdat,aes(ymin=lwr,ymax=upr),fill='mediumslateblue',alpha=0.2)+geom_ribbon(data=newdat2,aes(ymin=lwr,ymax=upr),fill='lawngreen',alpha=0.2)+
  geom_line(data=newdat,colour='#330066',linetype='dotted',size=1.4)+geom_line(data=newdat2,colour='#669933',size=1.4)+xlab('Temperature °C')+ylab(expression(Delta*'CTI'))+
  theme_bw()+theme(panel.background = element_blank(),panel.grid=element_blank(),legend.position='none')

# Scenario D: Predict for average depth, wide STR, low response diversity (sdSTI): predict low CTI sensitivity (Burrows et al fiture 1 g)
# This is plot pLow
newdat=data.frame(Survey.T.Match=seq(8,28,length.out=1000),Site_Score200="High",Depth=10,STR=19,sdSTI=1)
newdat<-cbind(newdat,predict(m1sim$gam, newdat, se.fit=TRUE))
newdat<-transform(newdat,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals
newdat2=data.frame(Survey.T.Match=seq(8,28,length.out=1000),Site_Score200="Low",Depth=10,STR=19,sdSTI=1)
newdat2<-cbind(newdat2,predict(m1sim$gam, newdat2, se.fit=TRUE))
newdat2<-transform(newdat2,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals

# Plot scenario D
pLow<-ggplot(NULL,aes(x=Survey.T.Match,y=fit))+xlim(c(8,28))+ylim(c(-3,3))+
  geom_ribbon(data=newdat,aes(ymin=lwr,ymax=upr),fill='mediumslateblue',alpha=0.2)+geom_ribbon(data=newdat2,aes(ymin=lwr,ymax=upr),fill='lawngreen',alpha=0.2)+
  geom_line(data=newdat,colour='#330066',linetype='dotted',size=1.4)+geom_line(data=newdat2,linetype='solid',colour='#669933',size=1.4)+xlab('Temperature °C')+ylab(expression(Delta*'CTI'))+
  theme_bw()+theme(panel.background = element_blank(),panel.grid=element_blank(),legend.position='none')


##### MANUSCRIPT FIG 2
####### MAP of CTI sensitivity scenarios & simulation outputs
# Data for mapping of spatial polygons
# regions.df was created using OGR vector maps as spatial objects (from MEOW public library containing spatial
# information on Marine Ecoregions Of the World)
regions<-readOGR("/Users/MEOW/", "meow_ecos") # reads OGR vector maps into spatial objects

#Turn SpatialPolygons into data frames
#from https://github.com/hadley/ggplot2/wiki/plotting-polygon-shapefiles
regions@data$id = rownames(regions@data)
regions.points = fortify(regions, ECOREGION="id")
regions.df = plyr::join(regions.points, regions@data, by="id")

# Fill SpatialPolygons
mything <- data.frame(ECOREGION = regions$ECOREGION,
                      score = runif(nrow(regions), -100, 100))

names(Fish_ecoregionPool)[names(Fish_ecoregionPool) == "ECOregion"] <- "ECOREGION" #match spelling of ECOregion across spatial regions df and RLS data

mything <- left_join(mything, Fish_ecoregionPool,by="ECOREGION") # merge with RLS ecoregion data (CTI sensitivity)
mything$ECOREGION<-as.factor(mything$ECOREGION)

regions.df <- merge(regions.df, mything)
#

regions.df<-read.csv('regions.df.csv')
regions.df$CTIsensitivity<-as.factor(regions.df$CTIsensitivity)
regions.df$CTIsensitivity <- ordered(regions.df$CTIsensitivity, levels = c("Low", "Medium", "High"))

# Adds coloured outlines of each scenario plot 
pLow<-pLow+theme(panel.background = element_blank(),panel.border = element_rect(colour = "yellowgreen", fill=NA,size=1.5),panel.grid=element_blank())
pHigh<-pHigh+theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA,size=1.5),panel.grid=element_blank(),axis.title.x=element_blank())
pMed1<-pMed1+theme(panel.background = element_blank(),panel.border = element_rect(colour = "turquoise3", fill=NA,size=1.5),panel.grid=element_blank(),axis.title.x=element_blank())
pMed2<-pMed2+theme(panel.background = element_blank(),panel.border = element_rect(colour = "turquoise3", fill=NA,size=1.5),panel.grid=element_blank())

### Annotates CTI sensitivity scenario plots (adding labels and setting location)
pMed1a<-pMed1+ annotate(geom="text", x=8, y=2, label="Medium\n(narrow STR)",size=3,hjust=0,fontface=1)
pLowa<-pLow+ annotate(geom="text", x=8, y=2.3, label="Low",size=3,hjust=0,fontface=1)
pHigha<-pHigh+ annotate(geom="text", x=8, y=2.3, label="High",size=3,hjust=0,fontface=1)
pMed2a<-pMed2+ annotate(geom="text", x=8, y=2, label="Medium\n(wide STR)",size=3,hjust=0,fontface=1)
text <- ggplot() +theme_void() +geom_text(aes(0,0,label='CTI sensitivity scenarios'),cex=4,angle=90) + xlab(NULL)

# Map plot:
library(viridis); library(ggpubr);library(sf);library("rnaturalearth"); library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")

map_STR<-ggplot(data=regions.df,na.rm=TRUE) + theme_bw()  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),axis.text = element_blank(),plot.margin = unit(c(0,8,0,8), "lines"))+
  geom_sf(data = world, fill = "grey80", colour = "grey80", size = 0.25)+
  geom_polygon(data=regions.df3,mapping=aes(long,lat,group=group,fill=CTIsensitivity2),colour='black',size=0.05) +xlab("Longitude °")+ylab("Latitude °")+ scale_y_continuous(expand = c(0,0), limits = c(-84, 84))+
  scale_x_continuous(expand = c(0,0), limits = c(-160, 180))+scale_fill_manual(values=c("turquoise3",'black'))+
  theme(legend.position = "bottom", legend.box = "vertical")+ labs(fill = "CTI sensitivity")

# Add side panel plots with 4 scenarios of CTI sensitivity
MAP<-map_STR+
  annotation_custom(ggplotGrob(pMed1a), xmin = 180, xmax = 290, ymin = -5, ymax=90)+ 
  annotation_custom(ggplotGrob(pMed2a), xmin = 180, xmax = 290, ymin = -115, ymax=-5)+
  annotation_custom(ggplotGrob(pHigha), xmin = -270, xmax = -160, ymin = -5, ymax=90)+
  annotation_custom(ggplotGrob(pLowa), xmin = -270, xmax = -160, ymin = -115, ymax=-5)#+
annotation_custom(ggplotGrob(text),xmin=-280,xmax=-255,ymin=-100,ymax=100)

# Export final map
pdf(width = 8, useDingbats=TRUE,height = 4, bg="white", file="Fig1_Map")
MAP
dev.off()

###############################################
# Determine CTI sensitivity for each Ecoregion (based on distribution of sdSTI and STR)
# Where might we expect a CTI response?
FINAL_URCHIN_DF<-FINAL_URCHIN_DF %>% mutate(CTIsensitivity = case_when(
  sdSTI >= 2 & STR<13 ~ 'High',
  sdSTI >= 2 & STR>13 ~ 'Medium',
  sdSTI < 2 & STR<13 ~ 'Medium',
  sdSTI<2 & STR>13 ~ 'Low'))


################################################
#### MODELS: CTI, CTDiv, CTR
FINAL_URCHIN_DF$COMBINED_SCORE<-as.factor(FINAL_URCHIN_DF$COMBINED_SCORE)
FINAL_URCHIN_DF$COMBINED_SCORE<-
  relevel(FINAL_URCHIN_DF$COMBINED_SCORE,ref="Kelp") # order Habitat factor

# SCALE Community Thermal Composition Metrics for each ecoregion
FINAL_URCHIN_DF <- FINAL_URCHIN_DF %>% group_by(ECOregion) %>% mutate(CTImax_fish_ECS=scale(CTImax),CTDiv_fish_ECS=scale(CTDiv),AvrgBreadth_fish_ECS=scale(AvrgBreadth))

# CTI
m1High=gamm(data=subset(FINAL_URCHIN_DF,CTIsensitivity=='High'),CTImax_fish_ECS~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) 
m1Medium=gamm(data=subset(FINAL_URCHIN_DF,CTIsensitivity=='Medium'),CTImax_fish_ECS~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) 

# CTDiv
m2High=gamm(data=subset(FINAL_URCHIN_DF,CTIsensitivity=='High'),CTDiv_fish_ECS~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) 
m2Medium=gamm(data=subset(FINAL_URCHIN_DF,CTIsensitivity=='Medium'),CTDiv_fish_ECS~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) 

# CTR
m3High=gamm(data=subset(FINAL_URCHIN_DF,CTIsensitivity=='High'),AvrgBreadth_fish_ECS~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) 
m3Medium=gamm(data=subset(FINAL_URCHIN_DF,CTIsensitivity=='Medium'),AvrgBreadth_fish_ECS~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit)

##### MANUSCRIPT FIG 3
### PLOTTING (CTI change)
pdf(width = 6, useDingbats=TRUE,height = 8, bg="white", file="Fig3_CTI")
par(mfrow=c(3,2),oma=c(4,1,1,1))

#CTImax
# High
par(mai=c(0.2,0.6,0.2,0))
plot(m1High$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"CTI")),cex.main=1.6,rug=FALSE,ylim=c(-1.5,1.1))
axis(side=2,at=c(-1,-0,1),cex.axis=1.5)
par(new=TRUE)
plot(m1High$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=3,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.8,rug=FALSE,ylim=c(-1.5,1.1))
text(9,0.8,"a",cex=1.8)
mtext('High', side=3, line=0.5, font=1,cex=1.5)
text(25,-1.4,"***",col='#669933',cex=2)
text(25,-1.2,"....",col='#330066',cex=2)

# Medium
par(mai=c(0.2,0.5,0.2,0.1))
plot(m1Medium$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, yaxt='n',cex.axis=1.5, cex.lab=1.6, xlab="",ylab='',cex.main=1.6,rug=FALSE,ylim=c(-2.1,1.2))
axis(side=2,at=c(-1,-0,1),cex.axis=1.5)
par(new=TRUE)
plot(m1Medium$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=3,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.8,rug=FALSE,ylim=c(-2.1,1.2))
text(11,0.9,"d",cex=1.8)
mtext('Medium', side=3, line=0.5, font=1,cex=1.5)
text(27,-1.8,"***",col='#669933',cex=2)
text(27,-1.6,"....",col='#330066',cex=2)

#CTDiv
# High
par(mai=c(0.2,0.6,0.2,0))
plot(m2High$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, yaxt='n',cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"CTDiv")),cex.main=1.6,rug=FALSE,ylim=c(-0.7,1))
axis(side=2,at=c(-0.5,0.1,0.6),cex.axis=1.5)
par(new=TRUE)
plot(m2High$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=3,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.7,1))
text(9,0.8,"b",cex=1.8)
text(25,0.7,"***",col='#669933',cex=2)
text(25,0.8,"....",col='#330066',cex=2)

# Medium
par(mai=c(0.2,0.5,0.2,0.1))
plot(m2Medium$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, xlab="",ylab='',cex.main=1.6,rug=FALSE,ylim=c(-1.1,1.2))
axis(side=2,at=c(-0.5,-0,0.5),cex.axis=1.5)
par(new=TRUE)
plot(m2Medium$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=3,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1.1,1.2))
text(11,1,"e",cex=1.8)
text(27,0.9,"***",col='#669933',cex=2)
# text(27,0.7,"....",col='#330066',cex=2) no signif barren effect


# CTR
# High
par(mai=c(0.2,0.6,0.2,0))
plot(m3High$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, yaxt='n',cex.axis=1.5, cex.lab=1.6, xlab="Survey Temperature °C",ylab=expression(paste(~Delta,"CTR")),cex.main=1.6,rug=FALSE,ylim=c(-0.5,1.1))
axis(side=2,at=c(-0.4,0.1,0.6),cex.axis=1.5)
par(new=TRUE)
plot(m3High$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=3,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.5,1.1))
text(9,0.9,"c",cex=1.8)
text(25,0.8,"***",col='#669933',cex=2)
text(25,0.9,"....",col='#330066',cex=2)

# Medium
par(mai=c(0.2,0.5,0.2,0.1))
plot(m3Medium$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, yaxt='n',cex.axis=1.5, cex.lab=1.6, xlab="Survey Temperature °C",ylab='',cex.main=1.6,rug=FALSE,ylim=c(-1,1.5))
axis(side=2,at=c(-1,0,1),cex.axis=1.5)
par(new=TRUE)
plot(m3Medium$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=3,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1,1.5))
legend("bottomright", c("Barren","Kelp"),lty=c(3,1), lwd=3,col=c("#330066", "#669933"), bty="n", cex=1.5)
text(11,1.3,"f",cex=1.8)
mtext("Survey Temperature °C", side=1, outer=T, at=0.5,cex=1.5,line=2)
text(27,1.3,"***",col='#669933',cex=2)


dev.off()

### THERMAL AFFINITY CLASSIFICATION (cold vs warm fishes)
# Upper and Lower Quantiles of each ecoregional species pool thermal midpoint:
FishAffinitiesQuartile<-plyr::ddply(FishAffinities,c("ECOregion"),summarize,coldestQ=quantile(T_MidPoint,probs = 0.25,na.rm=TRUE),warmestQ=quantile(T_MidPoint,probs = 0.75,na.rm=TRUE)) #Quartile (Main Analysis)
FishAffinitiesDecile<-plyr::ddply(FishAffinities,c("ECOregion"),summarize,coldestQ=quantile(T_MidPoint,probs = 0.10,na.rm=TRUE),warmestQ=quantile(T_MidPoint,probs = 0.90,na.rm=TRUE)) # Decile (Sensitivity Analysis)

FishAffinities<-left_join(FishAffinities,FishAffinitiesQuartile,by="ECOregion")

# Classify thermal affinity of fishes relative to ecoregional species pool thermal midpoints (here based on quartiles)
FishAffinities<-FishAffinities %>% mutate(T_Affinity = case_when(
  ECOregion =='Alboran Sea' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Alboran Sea' & T_MidPoint <= coldestQ ~ 'cold', 
  ECOregion =='Azores Canaries Madeira' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Azores Canaries Madeira' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Bassian' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Bassian' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Cape Howe' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Cape Howe' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Central Kuroshio Current' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Central Kuroshio Current' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Chiloense' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Chiloense' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Kermadec Island' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Kermadec Island' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Lord Howe and Norfolk Islands' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Lord Howe and Norfolk Islands' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Manning-Hawkesbury' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Manning-Hawkesbury' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='North Patagonian Gulfs' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='North Patagonian Gulfs' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Northeastern New Zealand' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Northeastern New Zealand' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Northern California' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Northern California' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='South Australian Gulfs' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='South Australian Gulfs' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='South European Atlantic Shelf' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='South European Atlantic Shelf' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Southern California Bight' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Southern California Bight' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Three Kings-North Cape' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Three Kings-North Cape' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Tweed-Moreton' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Tweed-Moreton' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Western Bassian' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Western Bassian' & T_MidPoint <= coldestQ ~ 'cold',
  ECOregion =='Western Mediterranean' & T_MidPoint >= warmestQ ~ 'warm',
  ECOregion =='Western Mediterranean' & T_MidPoint <= coldestQ ~ 'cold'))

# Calculate site-level richness and abundance of warm & cold affinity fishes
FishAffinities<-FishAffinities %>% group_by(SurveyID, T_Affinity) %>% summarise(Richness=length(SPECIES_NAME),Abundance=sum(TotalAbundance)) %>% ungroup() 
coldFishQuart<-FishAffinities %>% filter(T_Affinity=='cold') %>% droplevels() # extract cold affinity fish
warmFishQuart<-FishAffinities %>% filter(T_Affinity=='warm') %>% droplevels() # extract warm affinity fish
# richness and AB of cold and warm affinity fishes then merged back to final df (FINAL_URCHIN_DF)

#### Models: Thermal Affinities
# showing results when modelling with 'percentage of cold/warm Fish' per site here, but same result with raw data
mcold=gamm(data=subset(FINAL_URCHIN_DF,CTIsensitivity=='High'),scale(Richness_coldF_perc)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) 
mwarm=gamm(data=subset(FINAL_URCHIN_DF,CTIsensitivity=='High'),scale(Richness_warmF_perc)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit)

#### MANUSCRIPT FIG 4
### PLOT T-affinity Fish
pdf(width = 4, useDingbats=TRUE,height = 7, bg="white", file="Fig3_ThermalAffinity")
par(mfrow=c(2,1))
par(mar=c(4,5,1,0.5),oma=c(1,1,1,1))

# Warm Affinity Fishes
plot(mwarm$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-2,1.1))
axis(side=2,at=c(-1.5,-0.5,0.5),cex.axis=1.5)
par(new=TRUE)
plot(mwarm$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=3,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-2,1.1))
text(10,0.8,"a",cex=1.8)
mtext('Warm affinity fishes', side=2, line=4.5,font=2,cex=1.5)
text(11,-1.5,'....',cex=2,col='#330066')
text(11,-1.7,'***',cex=2,col='#669933')

# Cold Affinity Fishes
plot(mcold$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, xlab="Survey Temperature °C",ylab=expression(paste(~Delta,"Richness")),rug=FALSE,ylim=c(-0.8,1.5))
axis(side=2,at=c(-0.6,0.1,0.8),cex.axis=1.5)
par(new=TRUE)
plot(mcold$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=3,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.8,1.5))
text(10,1.3,"b",cex=1.8)
mtext('Cold affinity fishes', side=2, line=4.5, font=2,cex=1.5)
legend("topright", c("Barren","Kelp"),lty=c(3,1), lwd=3,col=c("#330066", "#669933"), bty="n", cex=1.5)
text(11,-0.5,'....',cex=2,col='#330066')
text(11,-0.7,'***',cex=2,col='#669933')

dev.off()

######## END OF MAIN RESULTS/ANALYSIS #########


#### SUPPLEMENTARY ANALYSIS
## 1 ## FISH SPECIES ABUNDANCE, RICHNESS AND FUNCTIONAL RICHNESS

### FUNCTIONAL RICHNESS
library(FD)

# create trait matrix (spp = rows, traits =col)
FishSiteAbundance$SPECIES_NAME2<-FishSiteAbundance$SPECIES_NAME
FishSiteAbundance$SPECIES_NAME2 <- make.names(FishSiteAbundance$SPECIES_NAME)
FishSiteAbundance$sqrt_Abundance<-sqrt(FishSiteAbundance$TotalAbundance) # take sqrt of fish abundance

#FishSiteAbundanceID$SurveyID<-as.factor(FishSiteAbundanceID$SurveyID)
#FishSiteAbundance<-FishSiteAbundance %>% mutate_if(is.numeric,scale) #scales all numeric columns

site.fish.matrix <- reshape::cast(FishSiteAbundance, SurveyID ~ SPECIES_NAME2, value='sqrt_Abundance') #makes site x species matrix with number of occurrences per trait group and site
site.fish.matrix[is.na(site.fish.matrix)] <- 0 # add zeros for 'trait absence'
site.fish.matrix <- as.data.frame(site.fish.matrix) # make a data frame

rownames(site.fish.matrix)<-site.fish.matrix$SurveyID # turn SurveyID column into rownames and 
site.fish.matrixSQRT<-site.fish.matrixSQRT[-c(1)] # remove SurveyID col

# Do this to turn empty/blank cells in NAs (otherwise dbFD doesn't work; do here bc if doing with TRAITS df, rownames as spp get lost)
FishSiteAbundance_2<-FishSiteAbundance %>% mutate_all(na_if,"")

# Select trait data only
# note: Water.column = Activity Group (fish position in the water column [benthic, demersal, pelagic])
TRAITS<-FishSiteAbundance_2 %>% dplyr::select(SPECIES_NAME2,Trophic.group,Water.column,Diel.Activity,MaxLength,Gregariousness)  %>% distinct(SPECIES_NAME2,.keep_all=TRUE) %>% droplevels()
TRAITS$SPECIES_NAME2<-as.character(TRAITS$SPECIES_NAME2)

# use this to make sure col and row names are identical in the 2 df
colnames(site.fish.matrix)<-TRAITS$SPECIES_NAME2

# then finish TRAITS df
rownames(TRAITS)<-TRAITS$SPECIES_NAME2
TRAITS<-TRAITS[-c(1)]
TRAITS$MaxLength<-log(TRAITS$MaxLength+1) 
TRAITS$Gregariousness<-log(TRAITS$Gregariousness+1) 
#TRAITS$Gregariousness<-as.factor(TRAITS$Gregariousness)

FD_FishSQRT<-dbFD(TRAITS,site.fish.matrix,w.abun = TRUE,corr='none',stand.x=TRUE) 

FRic<-as.data.frame(FD_FishSQRT$FRic) # Functional Richness
Fnsp<-as.data.frame(FD_FishSQRT_log$nbsp) # Number of Species (Species Richness)

FD_Fish_dfSQRT<-cbind(FRic,Fnsp)

FD_Fish_dfSQRT$SurveyID<-rownames(FD_Fish_dfSQRT) # Add SurveyID back as a column
FD_Fish_dfSQRT$SurveyID<-as.numeric(FD_Fish_dfSQRT$SurveyID)

colnames(FD_Fish_dfSQRT)[which(names(FD_Fish_dfSQRT) == "FD_FishSQRT$FRic")] <- "FRic"
colnames(FD_Fish_dfSQRT)[which(names(FD_Fish_dfSQRT) == "FD_FishSQRT$Fnsp")] <- "Site_Rich"

FD_Fish_dfSQRT_log <- FD_Fish_dfSQRT_log %>% filter(Fnsp>2) %>% droplevels()
FINAL_URCHIN_DF<-left_join(FINAL_URCHIN_DF,FD_Fish_dfSQRT,by="SurveyID") # Merge back with Community df

## SITE LEVEL FISH ABUNDANCE, RICHNESS AND FUNCTIONAL RICHNESS MODELS 
# Standardizing FRic by site species richness (at level of x species richness, level of functional richness is...)
res.FRic_nsp=lme(FRic~Site_Rich,random=~1|ECOregion/r.SiteLat,data=FINAL_URCHIN_DF,na.action=na.exclude)
res.FRic_nsp=data.frame(residuals(res.FRic_nsp));colnames(res.FRic_nsp)=c("res.FRic_nsp")
FINAL_URCHIN_DF=bind_cols(FINAL_URCHIN_DF,res.FRic_nsp)

m0r=gamm(data=FINAL_URCHIN_DF,scale(res.FRic_nsp)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit)
m1r=gamm(data=FINAL_URCHIN_DF,scale(Site_Rich)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) 
m2r=gamm(data=FINAL_URCHIN_DF,scale(Site_AB)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1))

#### SUPP FIG S1 #####
# Use this 'newdat' snipet to predict for each model (m0r to make plot pres.FRic_nsp; m1r to make plot pRic; m2r to make plot pAB) )
newdat=data.frame(Survey.T.Match=seq(10,26,length.out=3),COMBINED_SCORE="Barren",Depth=10)
newdat<-cbind(newdat,predict(m0r$gam, newdat, se.fit=TRUE))
newdat<-transform(newdat,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals
newdat2=data.frame(Survey.T.Match=seq(10,26,length.out=3),COMBINED_SCORE="Kelp",Depth=10)
newdat2<-cbind(newdat2,predict(m0r$gam, newdat2, se.fit=TRUE))
newdat2<-transform(newdat2,upr=fit+(2*se.fit),lwr=fit-(2*se.fit)) # adds upr and lower Confidence Intervals

newdat3<-rbind(newdat,newdat2)
newdat3$Survey.T.Match<-as.factor(newdat3$Survey.T.Match)

#
pres.FRic_nsp<-ggplot(newdat3,aes(x=Survey.T.Match,y=fit,group=COMBINED_SCORE,colour=COMBINED_SCORE))+ 
  geom_pointrange(aes(ymin=fit-se.fit, ymax=fit+se.fit))+ylab(expression(Delta*'Functional Richness'))+xlab('Survey Temperature °C')+theme_bw(base_size = 15)+
  theme(panel.background = element_blank(),panel.grid=element_blank(),legend.position=c(0.7, 0.25),legend.box.background = element_rect(colour = "black"))+
  labs(colour = "Urchin Class")+ 
  scale_color_manual(values=c("#330066", "#669933"),name = "Habitat", labels = c("Barren", "Kelp"))+
  scale_y_continuous(breaks=c(-0.90,-0.4,0.10),limits=c(-1.2,0.4))

pRic<-ggplot(newdat3,aes(x=Survey.T.Match,y=fit,group=COMBINED_SCORE,colour=COMBINED_SCORE,fill=COMBINED_SCORE))+ 
  geom_pointrange(aes(ymin=fit-se.fit, ymax=fit+se.fit),position=position_dodge(width=0.3))+ylab(expression(Delta*'Species Richness'))+xlab('Survey Temperature °C')+theme_bw(base_size = 15)+
  theme(panel.background = element_blank(),panel.grid=element_blank(),legend.position="none")+ scale_color_manual(values=c("#330066", "#669933"), labels = c("Kelp", "Barren"))+scale_fill_manual(values=c('white','white'))+
  scale_y_continuous(breaks=c(-0.90,-0.30,0.30),limits=c(-1.2,0.8))

pAB<-ggplot(newdat3,aes(x=Survey.T.Match,y=fit,group=COMBINED_SCORE,colour=COMBINED_SCORE))+ 
  geom_pointrange(aes(ymin=fit-se.fit, ymax=fit+se.fit),position=position_dodge(width=0.3))+ylab(expression(Delta*'Site Abundance'))+xlab('Survey Temperature °C')+theme_bw(base_size = 15)+
  geom_text(x=15,y=0.5,label="a",fontface=1,cex=3)+
  theme(panel.background = element_blank(),panel.grid=element_blank(),legend.position="none")+ scale_color_manual(values=c("#330066", "#669933"),name = "Habitat", labels = c("Kelp", "Barren"))+
  scale_y_continuous(breaks=c(-0.30,0.2,0.70),limits=c(-0.6,0.9))

dev.new()
FunRic_plot<-ggarrange(pAB, pRic, pres.FRic_nsp,
                       labels = c("a", "b", "c"), hjust=-8,vjust=2.5,
                       font.label = list(size = 16, face = "plain"),
                       ncol = 3, nrow = 1,
                       align = "v")

pdf(width = 9.5, useDingbats=TRUE,height = 3, bg="white", file="Fig4_Funct_Richn")
FunRic_plot
dev.off()


## 2 ## INDIVIDUAL FISH TRAITS: TROPHIC GROUPS AND ACTIVITY GROUPS
benthic.invertivoreF<-read.csv("benthic.invertivoreF.csv") # richness & abundance
omnivoreF<-read.csv("omnivoreF.csv") # richness & abundance
carnivoreF2<-read.csv("carnivoreF2.csv") # richness
carnivoreF3<-read.csv("carnivoreF3.csv") # abundance
planktivoreF<-read.csv("planktivoreF.csv") # richness & abundance
herbivoreF2<-read.csv("herbivoreF2.csv") # richness
herbivoreF3<-read.csv("herbivoreF3.csv") # abundance

# Before running models, make sure Habitat Classification (COMBINED_SCORE) is factor & ordered correctly
benthic.invertivoreF$COMBINED_SCORE<-as.factor(benthic.invertivoreF$COMBINED_SCORE)
benthic.invertivoreF$COMBINED_SCORE<-
  relevel(benthic.invertivoreF$COMBINED_SCORE,ref="Kelp") # order Habitat factor; kelp first

### TROPHIC GROUP MODELS
# Benthic invertivore fihes
mBIF=gamm(data=benthic.invertivoreF,scale(TG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) #gaussian fits here, no difference, scaling, ECS, T>10 and modelling each separate = same pattern
mBIFA=gamm(data=benthic.invertivoreF,scale(log(TG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # scale(log()) is best model (resids), quassi poiss pulls up barrens at cold end because of one Low val; same pattern with ECS
# Planktivore fishes
mPl=gamm(data=subset(planktivoreF,TG_Richness>0),scale(TG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) #same pattern with gaussian, scale, scale(log()) and when subsetting for T>10deg
mPlA=gamm(data=subset(planktivoreF,TG_Abundance>0),scale(log(TG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # pattern consistent for subset, scaled, separate model, residuals best with scale(log), quasipois doesn't converge,negbin shows same pattern but doesn't fit well, ECS pulls cold Low sites above High
# Omnivore fishes
mOm=gamm(data=subset(omnivoreF,TG_Richness>0),scale(TG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) #looks good, pattern consistent with scale, separate models and T>10, same pattern with ECS
mOmA=gamm(data=subset(omnivoreF,TG_Abundance>0),scale(log(TG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) #pattern same with scale(log)), same pattern with ECS, separate model and poiss/negbin don't converge
# Herbivore fishes
mHer=gamm(data=subset(herbivoreF2,TG_Richness>0),scale(TG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) #gaussian is good, same pattern with scale, ECS ,transform and subset
mHerA=gamm(data=subset(herbivoreF3,TG_Abundance>0),scale(log(TG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth), random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1))#pattern is same for untransformed, separate and negbin, resid only look good for scale(log); but reduces significance, ECS shows slightly bigger diff for cold barrens but otherwise same
# Carnivore fishes
mCa=gamm(data=subset(carnivoreF2,TG_Richness>0),scale(TG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # gaussian fits, pattern consistent with transformation/scale/subset, ECS
mCaA=gamm(data=subset(carnivoreF3,TG_Abundance>0),scale(log(TG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) #same pattern with scale(log), ECS, separate and subset models, but only diff driven by cold High reg, plots look like no diff but model says yes?

#### SUPPL Fig S2 #####
### TROPHIC GROUP RICHNESS & ABUNDANCE PLOTS
pdf(width = 8, useDingbats=TRUE,height = 11, bg="white", file="SuppFig_Fish_TrophicGroups")
par(mfrow=c(5,2))
par(mar=c(4,5.5,1.25,0.5),oma=c(0.5,4,1,1))

plot(mBIF$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE)
par(new=TRUE)
plot(mBIF$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE)
text(9,0.4,"a",cex=1.8)
mtext('Benthic \ninvertivores', side=2, line=4.5,cex=1.5)

plot(mBIFA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-1,1))
par(new=TRUE)
plot(mBIFA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1,1))
text(9,0.8,"b",cex=1.8)

plot(mPl$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-1,1))
par(new=TRUE)
plot(mPl$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1,1))
text(11,0.8,"c",cex=1.8)
mtext('Planktivores', side=2, line=6,cex=1.5)

plot(mPlA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.6))
par(new=TRUE)
plot(mPlA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.6))
text(11,0.5,"d",cex=1.8)

plot(mOm$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-1.5,1))
par(new=TRUE)
plot(mOm$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1.5,1))
text(10,0.55,"e",cex=1.8)
mtext('Omnivores', side=2, line=6,cex=1.5)

plot(mOmA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-1.2,0.8))
par(new=TRUE)
plot(mOmA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1.2,0.8))
text(10,0.6,"f",cex=1.8)

plot(mHer$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-1,1.2))
par(new=TRUE)
plot(mHer$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1,1.2))
text(11,0.9,"g",cex=1.8)
mtext('Herbivores', side=2, line=6,cex=1.5)

plot(mHerA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-1,0.5))
par(new=TRUE)
plot(mHerA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1,0.5))
text(11,0.35,"h",cex=1.8)

plot(mCa$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="Survey Temperature °C",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-1.2,1))
par(new=TRUE)
plot(mCa$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1.2,1))
text(11,0.75,"i",cex=1.8)
mtext('Carnivores', side=2, line=6,cex=1.5)

plot(mCaA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="Survey Temperature °C",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.8))
par(new=TRUE)
plot(mCaA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.8))
text(11,0.6,"j",cex=1.8)
legend("bottomright", c("Barren","Kelp"), lwd=3,col=c("#330066", "#669933"), bty="n", cex=1.5)
dev.off()

#### ACTIVITY GROUP MODELS
benthicF<-read.csv("benthicF.csv")
demersalF<-read.csv("demersalF.csv")
pelagicSAF<-read.csv("pelagicSAF.csv") # site-attached pelagic fishes
pelagicF<-read.csv("pelagicF.csv") # all pelagic fishes
pelagicNSAF<-read.csv("pelagicNSAF.csv") # non site-attached pelagic fishes

mBenthic=gamm(data=benthicF,scale(AG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overalp); model check plots look good; 
mBenthicA=gamm(data=benthicF,scale(log(AG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overlap); model check plots look good

mDemersal=gamm(data=demersalF,scale(AG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overalp); model check plots look good; 
mDemersalA=gamm(data=demersalF,scale(log(AG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overlap); model check plots look good

mPelagicSAF=gamm(data=pelagicSAF,scale(AG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overlap)
mPelagicSAFA=gamm(data=pelagicSAF,scale(log(AG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overlap)

mPelagic=gamm(data=pelagicF,scale(AG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # not enough data
mPelagicA=gamm(data=pelagicF,scale(log(AG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overlap)

mPelagicNSA=gamm(data=pelagicNSAF,scale(AG_Richness)~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overlap)
mPelagicNSAA=gamm(data=pelagicNSAF,scale(log(AG_Abundance+1))~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1)) # no difference (CI overlap)

#### ACTIVITY GROUP RICHNESS & ABUNDANCE PLOTS
### SUPPL FIG S3 #####
pdf(width = 8, useDingbats=TRUE,height = 11, bg="white", file="SuppFig_Fish_ActivityGroup")
par(mfrow=c(5,2))
par(mar=c(4,5.5,1.25,0.5),oma=c(0.5,4,1,1))

plot(mBenthic$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-0.5,0.5))
par(new=TRUE)
plot(mBenthic$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.5,0.5))
text(9,0.4,"a",cex=1.8)
mtext('Benthic', side=2, line=6,cex=1.5)

plot(mBenthicA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-0.7,0.7))
par(new=TRUE)
plot(mBenthicA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.7,0.7))
text(9,0.6,"b",cex=1.8)

plot(mDemersal$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-0.9,1.2))
par(new=TRUE)
plot(mDemersal$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.9,1.2))
text(9,0.95,"c",cex=1.8)
mtext('Demersal', side=2, line=6,cex=1.5)

plot(mDemersalA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-1,0.7))
par(new=TRUE)
plot(mDemersalA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1,0.7))
text(9,0.5,"d",cex=1.8)

plot(mPelagic$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.5))
par(new=TRUE)
plot(mPelagic$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.5))
text(10,0.35,"e",cex=1.8)
mtext('Pelagic', side=2, line=6,cex=1.5)

plot(mPelagicA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.4))
par(new=TRUE)
plot(mPelagicA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.4))
text(10,0.25,"f",cex=1.8)

plot(mPelagicSAF$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-1,0.8))
par(new=TRUE)
plot(mPelagicSAF$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1,0.8))
text(10,0.6,"g",cex=1.8)
mtext('Pelagic \nSA', side=2, line=4.5,cex=1.5)

plot(mPelagicSAFA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.6))
par(new=TRUE)
plot(mPelagicSAFA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.6))
text(10,0.4,"h",cex=1.8)

plot(mPelagicNSA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="Survey Temperature °C",ylab=expression(paste(~Delta,"Richness")),cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.5))
par(new=TRUE)
plot(mPelagicNSA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-0.8,0.5))
text(12,0.35,"i",cex=1.8)
mtext('Pelagic \nNSA', side=2, line=4.5,cex=1.5)

plot(mPelagicNSAA$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="Survey Temperature °C",ylab=expression(paste(~Delta,"log Abundance")),cex.main=1.6,rug=FALSE,ylim=c(-1,0.5))
par(new=TRUE)
plot(mPelagicNSAA$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.6,rug=FALSE,ylim=c(-1,0.5))
text(12,0.3,"j",cex=1.8)
legend("bottomright", c("Barren","Kelp"), lwd=3,col=c("#330066", "#669933"), bty="n", cex=1.5)
dev.off()


####### SENSITIVITY ANALYSIS USING FULL DATASET (i.e. also including ecoregions/sites where no photo-quadrats
# are available; 19 ecoregions, 5596 sites included)
# presented in supplementary table S12 and supplementary figure S9
FULL_URCHIN_DATA$Site_Score200<-as.factor(FULL_URCHIN_DATA$Site_Score200)
FULL_URCHIN_DATA$Site_Score200<-
  relevel(FULL_URCHIN_DATA$Site_Score200,ref="Low") # order Habitat factor (low urchin density [= kelp] is first level)

m1Full=gamm(data=subset(FULL_URCHIN_DATA,CTIsensitivity=='High'),
            CTImax_fish_ECS~s(Survey.T.Match,by=Site_Score200,k=3)+scale(Depth),
            random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) 

# Northern Hemisphere only
m1NH=gamm(data=subset(FULL_URCHIN_DATA,CTIsensitivity!='Low' & r.SiteLat >0),
          CTImax_fish_ECS~s(Survey.T.Match,by=Site_Score200,k=3)+scale(Depth),
          random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) # only 2 ecoregions with high CTI sensitivty in NH, 
# so modelling high & medium combined due to limited temp range of the 2 high CTI sensitivity regions.
# pattern qualitatively the same but weaker, given limited number of sensitive ecoregions with data available.

# Southern Hemisphere only
m1SH=gamm(data=subset(FULL_URCHIN_DATA,CTIsensitivity=='High' & r.SiteLat <0),
          CTImax_fish_ECS~s(Survey.T.Match,by=Site_Score200,k=3)+scale(Depth),
          random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit) 

plot(m1Full$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"CTI")),cex.main=1.6,rug=FALSE,ylim=c(-1.5,1.1))
axis(side=2,at=c(-1,-0,1),cex.axis=1.5)
par(new=TRUE)
plot(m1Full$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.8,rug=FALSE,ylim=c(-1.5,1.1))
text(9,0.8,"a",cex=1.8)
mtext('High', side=3, line=0.5, font=1,cex=1.5)

plot(m1NH$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"CTI")),cex.main=1.6,rug=FALSE,ylim=c(-1.5,1.1))
axis(side=2,at=c(-1,-0,1),cex.axis=1.5)
par(new=TRUE)
plot(m1NH$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.8,rug=FALSE,ylim=c(-1.5,1.1))
text(9,0.8,"a",cex=1.8)
mtext('Northern Hemisphere', side=3, line=0.5, font=1,cex=1.5)

plot(m1SH$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, xlab="",ylab=expression(paste(~Delta,"CTI")),cex.main=1.6,rug=FALSE,ylim=c(-1.5,1.1))
axis(side=2,at=c(-1,-0,1),cex.axis=1.5)
par(new=TRUE)
plot(m1SH$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.8,rug=FALSE,ylim=c(-1.5,1.1))
text(9,0.8,"a",cex=1.8)
mtext('Southern Hemisphere', side=3, line=0.5, font=1,cex=1.5)


##### SENSITIVITY ANALYSIS ### 
# COMPARING TO NATURALLY BARE AREAS (with not kelp and < urchins per m^2)
# Supplementary Figure S5 and Table S7

#CTImax
# High
m1Natbare=gamm(data=subset(NATBARE,CTIsensitivity=='High'),CTImax_fish_ECS~s(Survey.T.Match,by=COMBINED_SCORE,k=3)+scale(Depth),random=list(ECOregion=~1,r.SiteLat=~1,SiteCode=~1),na.action=na.omit)

pdf(width = 12, useDingbats=TRUE,height = 6, bg="white", file="Suppl_FigS5")
par(mfrow=c(1,2),oma=c(0,0,0.5,0.5))
par(mai=c(1,1,0.5,0.5))

# SUPP FIG S5 panel a
plot(m1Natbare$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, xlab="Survey Temperature °C",ylab=expression(paste(~Delta,"CTI")),cex.main=1.6,rug=FALSE,ylim=c(-1.5,1.1))
axis(side=2,at=c(-1,-0,1),cex.axis=1.5)
par(new=TRUE)
plot(m1Natbare$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.8,rug=FALSE,ylim=c(-1.5,1.1))
par(new=TRUE)
plot(m1Natbare$gam,shade=TRUE,select=3,shade.col=rgb(10,0,55,40,maxColorValue=655),col="black",lty=3,lwd=3,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.8,rug=FALSE,ylim=c(-1.5,1.1))
legend("bottomright", c('Barren',"Naturally bare","Kelp"),lty=c(1,3,1), lwd=3,col=c("#330066","black", "#669933"), bty="n", cex=1.5)
text(10,0.8,"a",cex=1.8)

# SUPP FIG S5 panel b
par(mai=c(0.5,0.6,0.5,0.1))

col_test<-c('#330066', '#669933')[as.numeric(FINAL_URCHIN_DF$COMBINED_SCORE)]
col_transp<-adjustcolor(col_test, alpha.f = 0.5)

plot(CTImax_fish_ECS~Survey.T.Match,data=subset(FINAL_URCHIN_DF,CTIsensitivity=='High'),yaxt='n',xaxt='n',ylab='',xlab='',col=col_transp,pch=16,ylim=c(-4.5,3.5))
legend("bottomright",legend=levels(FINAL_URCHIN_DF$COMBINED_SCORE),inset=.02,col=c('#330066', '#669933')[as.numeric(FINAL_URCHIN_DF$COMBINED_SCORE)], title="Habitat", horiz=TRUE, cex=0.8,pch=16,c('Kelp','Barren'))
par(new=TRUE)
plot(m1High$gam,select=1,shade=TRUE,shade.col=rgb(100,255,0,40,maxColorValue=255),col="#669933",lty=1,lwd=4,cex=1.8,yaxt='n', cex.axis=1.5, cex.lab=1.6, ,xlab="Survey Temperature °C",ylab=expression(paste(~Delta,"CTI")),cex.main=1.6,rug=FALSE,ylim=c(-4,3.5))
axis(side=2,at=c(-4,-2,-0,2,4),cex.axis=1.5)
par(new=TRUE)
plot(m1High$gam,shade=TRUE,select=2,shade.col=rgb(10,0,255,40,maxColorValue=255),col="#330066",lty=1,lwd=4,cex=1.8, cex.axis=1.5, cex.lab=1.6, xlab="",ylab="",yaxt="n", xaxt="n",cex.main=1.8,rug=FALSE,ylim=c(-4,3.5))
text(9.5,2.5,"b",cex=1.8)

dev.off()


##### END OF SCRIPT