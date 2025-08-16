# Add data and packages

library("ggplot2")  #package for plotting stacked bar charts, dendrogram, nmds
library('RColorBrewer') #package for plotting
library(patchwork) #package for multipanel plots
library(Polychrome) # for making color pallates with createpalette() .. makes new random palettes of distinct colors 
library(reshape2) # for data management stuff like wide-long
library(knitr) # for ??
library(raster) # for the cv function .. why not in base?!?!
library(plyr) # turn turn an array into a list .. 
library(gtools) # for the combinations() function used to get all pairwise combinations of focal taxa
library(wesanderson) #package for plotting
library(tidyverse) # package for data management (nmds, clustering plots)
library(ggdendro) # for the hierachical clustering plot
library(vegan) # for the NMDS
library(dplyr)
library(tidyr)
library('stats4')
library('viridis')
## For pocillo code
library("tools") # for "file_path_sans_ext()"
library('multcompView')
library('lmtest')

# Make a fuction to paste data directly into excel
pastey <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,file = paste0("clipboard-", object.size(x)),sep="\t",row.names=row.names,col.names=col.names,...)
}

# Make a function to copy data directly FROM excel
copey <- function(header=TRUE,...) {
  read.table("clipboard",sep="\t",header=header,...)
}



# additional instructions for multi panel plots - not used but maybe some good stuff in here
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

# Data Preparation
# Data are exported as .csv from TL from each site by transition (e.g. y1-y2, y2-y3 ...), combine into a single data frame 

# Directory path

file_path <- "~/Documents/Coral_Model/data/raw_coral_data"

# Get all CSV file names in the directory
taglab <- dir(path = file_path, full.names = TRUE, recursive = FALSE)

# Extract metadata from file names
file.meta <- do.call("rbind", strsplit(file_path_sans_ext(basename(taglab)), "_"))

# Initialize an empty list to store data
raw.data <- list()

# Loop through the files and read them, adding metadata as columns
for (i in 1:length(taglab)) {
  # Read the CSV file
  temp_data <- read.csv(taglab[i], header = TRUE, sep = ",")
  
  # Add metadata columns as the first columns
  temp_data <- temp_data %>%
    mutate(Island = file.meta[i, 2],  # PAL
           Site = file.meta[i, 3],    # FR132, FR14, etc.
           transition = file.meta[i, 4]) %>%
    dplyr::select(Island, Site, transition, everything())  # Move metadata to the front
  
  # Append to the list
  raw.data[[i]] <- temp_data
}

# Combine all data frames in the list into one data frame
#final_data <- bind_rows(raw.data)

names(raw.data)=file_path_sans_ext(taglab)
raw.combined=data.frame(do.call("rbind",raw.data))


# Due to the structure of the data resulting from FnF, colonies are repeated on rows, which creates some issues. The solution is to separate data into each transition,  pull out the unique entries for each ramet (eg unique blob in a given year) each year, sum genet for each year, then merge genets across years.  

# Repeated rows explanation: a colony that fuses will have two IDs in the first year of a transition, and one ID in the second year of a transition. Colonies that split will have one ID in the first year and two IDs in the second year. Explore one of the transition.split[[1]], filter to just split and fuse in action2 to get an idea of what is going on (I think that transition split is now blobs[[]])

# 1- Data prep
# 2- separate the data into transitions
# 3- collect the unique Genet-blob combos for each year
# 4- sum blobs within genetID
# 5- Merge genet IDs across years


############ 1- Combine IDs (blob IDs are repeated within sites across years, and GenIDs are repeated across sites)
############ x
#RID is the ramet or blob ID

# RID_ini uses the ID of a blob in the first year of a transition, RID_fin uses the ID in the second year  
raw.combined$RID_ini=paste(raw.combined$Site,raw.combined$Genet,raw.combined$Blob1,sep="_")
raw.combined$RID_fin=paste(raw.combined$Site,raw.combined$Genet,raw.combined$Blob2,sep="_")
raw.combined$SiteGenID=paste(raw.combined$Site,raw.combined$Genet,sep="_")
### subset to just Pocillopora
poc.raw=raw.combined[which(raw.combined$Class=="Pocillopora"),]
rownames(poc.raw)=NULL



############ 2- sum genets within years

transitions=c("t1","t2","t3","t4","t5","t6")
trans_years=unique(poc.raw$transition)
yearz=2013:2019

# data are stored as transitions, so pull out the 2013 sizes via initial sizes and IDs of the first year in that transition , then pull out sizes and IDs of subsequent years as the final sizes and IDs of that transition
# Store unique genet-blob combos for each year while keeping Class column
blobs <- list()
blobs[[1]] <- unique(poc.raw[which(poc.raw$transition == trans_years[1]), 
                             c("Genet", "Blob1", "RID_ini", "Area1", "Site", "SiteGenID", "Class")])

for (i in 1:length(trans_years)) {
  blobs[[i + 1]] <- unique(poc.raw[which(poc.raw$transition == trans_years[i]), 
                                   c("Genet", "Blob2", "RID_fin", "Area2", "Site", "SiteGenID", "Class")])
}

names(blobs) <- yearz

# Rename columns dynamically while keeping Class
for (i in 1:length(blobs)) {
  names(blobs[[i]]) <- c(
    paste(yearz[i], "GenID", sep = "_"),
    paste(yearz[i], "RID", sep = "_"),
    paste(yearz[i], "SiteGenRID", sep = "_"),
    paste(yearz[i], "size", sep = "_"),
    paste(yearz[i], "site", sep = "_"),
    paste(yearz[i], "SiteGenID", sep = "_"),
    "Class")  # Retaining Class column
}


############ 3- sum ramets within genets for each year
############

# 3- sum ramets/blobs within genets for each year
# sum size by GenID and SiteGenID (GenID is redundant, but keeping it for)
genetsums=list()
for (i in 1:length(blobs)){
  genetsums[[i]]=aggregate(blobs[[i]][,4],by=list(blobs[[i]][,6],blobs[[i]][,1],blobs[[i]][,5]),sum)
}


# Add names not to make column names after merging a little cleaner
for (i in 1:length(genetsums)){
  names(genetsums[[i]])=c(
    paste(yearz[i],"GenSiteID",sep="_"),
    paste(yearz[i],"GenID",sep="_"),
    paste(yearz[i],"Site",sep="_"),
    paste("size",yearz[i],sep="_"))
}

############ 4- Merge genets across year
############

merge1=merge(genetsums[[1]],genetsums[[2]],by.x="2013_GenSiteID",by.y="2014_GenSiteID",all=T)
merge2=merge(merge1,genetsums[[3]],by.x="2013_GenSiteID",by.y="2015_GenSiteID",all=T)
merge3=merge(merge2,genetsums[[4]],by.x="2013_GenSiteID",by.y="2016_GenSiteID",all=T)
merge4=merge(merge3,genetsums[[5]],by.x="2013_GenSiteID",by.y="2017_GenSiteID",all=T)
merge5=merge(merge4,genetsums[[6]],by.x="2013_GenSiteID",by.y="2018_GenSiteID",all=T)
final_merge_raw=merge(merge5,genetsums[[7]],by.x="2013_GenSiteID",by.y="2019_GenSiteID",all=T)


# Turn NAs into 0s for easier fate tracking. maybe not needed, but also helps to remove some weird rows with all 0's
final_merge_raw[is.na(final_merge_raw)] <- 0

# there are some entries that have a gen and site ID but all the size columns are NA. for these colonies the fate for every transition is alive. not sure what is going on here. filtering these out by hand for the time being as they can mess up sample size sums
size_columns <- grep("size_", names(final_merge_raw), value = TRUE)
# Remove rows where all size columns are 0
final_merge_raw2 <- final_merge_raw[-which(rowSums(final_merge_raw[, size_columns]) == 0), ]



# Extract the site names (should prob be done when reading in the original files .. )
merge_meta=data.frame(do.call("rbind",strsplit(final_merge_raw2$'2013_GenSiteID',"_")))



# Pull out just the columns of interest
final_merge=data.frame(site=merge_meta[,1],final_merge_raw2[,c('2013_GenSiteID','2013_GenID','size_2013','size_2014','size_2015','size_2016','size_2017','size_2018','size_2019')])


# Rename the first three columns
names(final_merge)[1:3]=c("Site","GenSiteID","GenID")




############ 5- Fate tracking
############

#For various reasons we want a matrix that has the size of each genet for each year, as well as a fate status which are as follows (this should prob be split into two parallel dataframes, but for now is stored in one)
#1- Alive: alive adult colony (not a new recruit)
#2- Recruit: new colony for that year
#3- Died: Colony that is alive in year ti, but which is dead by the following year (ti+1) 
#4- AD: First year after death. the associated size is 0.
#5- BC: temp fate for colonies that recruit later in the time series. Associated size is 0

# AD and BC are alter converted to NA for most purposes


####### 1- AD - Alive - Dead 


## Add "AD", "alive" and "dead" fates. (almost certainly a better way to do this)
#AD is for things that had already died earlier in the time series


#2013 is easy
final_merge$Fate_2013<-ifelse(final_merge$size_2013>0 & final_merge$size_2014 == 0,"died","alive")

#2014 is ok
final_merge$Fate_2014<-ifelse(final_merge$size_2014>0 & final_merge$size_2015 == 0,"died","alive")
final_merge$Fate_2014<-ifelse(final_merge$Fate_2013=="died" | final_merge$Fate_2013=="AD","AD",final_merge$Fate_2014)
#final_merge$Fate_2014<-ifelse(final_merge$size_2014>0 & final_merge$size_2015 == 0,"died",final_merge$Fate_2014)


#2015 is harder
final_merge$Fate_2015<-ifelse(final_merge$size_2015>0 & final_merge$size_2016 == 0,"died","alive")
final_merge$Fate_2015<-ifelse(final_merge$Fate_2014=="AD"| final_merge$Fate_2014=="AD","AD",final_merge$Fate_2015)
final_merge$Fate_2015<-ifelse(final_merge$Fate_2014=="died"| final_merge$Fate_2014=="AD","AD",final_merge$Fate_2015)


#2016 is harder
final_merge$Fate_2016<-ifelse(final_merge$size_2016>0 & final_merge$size_2017 == 0,"died","alive")
final_merge$Fate_2016<-ifelse(final_merge$Fate_2015=="AD"| final_merge$Fate_2015=="AD","AD",final_merge$Fate_2016)
final_merge$Fate_2016<-ifelse(final_merge$Fate_2015=="died"| final_merge$Fate_2015=="AD","AD",final_merge$Fate_2016)

#2017 is harder
final_merge$Fate_2017<-ifelse(final_merge$size_2017>0 & final_merge$size_2018 == 0,"died","alive")
final_merge$Fate_2017<-ifelse(final_merge$Fate_2016=="AD"| final_merge$Fate_2016=="AD","AD",final_merge$Fate_2017)
final_merge$Fate_2017<-ifelse(final_merge$Fate_2016=="died"| final_merge$Fate_2016=="AD","AD",final_merge$Fate_2017)

#2018 is too
final_merge$Fate_2018<-ifelse(final_merge$size_2018>0 & final_merge$size_2019 == 0,"died","alive")
final_merge$Fate_2018<-ifelse(final_merge$Fate_2017=="AD" | final_merge$Fate_2017=="AD","AD",final_merge$Fate_2018)
final_merge$Fate_2018<-ifelse(final_merge$Fate_2017=="died"| final_merge$Fate_2017=="AD","AD",final_merge$Fate_2018)



####### 1- Add 'recruit' and 'BC' fates  (I think that the above fates have to be added before these or else need a different approach)

#2018-2019
final_merge$Fate_2018<-ifelse(final_merge$size_2018==0 & final_merge$size_2019 > 0,"recruit",final_merge$Fate_2018)

#2017-2018
final_merge$Fate_2017<-ifelse(final_merge$size_2017==0 & final_merge$size_2018 > 0,"recruit",final_merge$Fate_2017)
final_merge$Fate_2017<-ifelse(final_merge$Fate_2018 == "recruit" | final_merge$Fate_2018 == "BC","BC",final_merge$Fate_2017)

#2016-2017
final_merge$Fate_2016<-ifelse(final_merge$size_2016==0 & final_merge$size_2017 > 0,"recruit",final_merge$Fate_2016)
final_merge$Fate_2016<-ifelse(final_merge$Fate_2017 == "recruit" | final_merge$Fate_2017 == "BC","BC",final_merge$Fate_2016)

#2015-2016
final_merge$Fate_2015<-ifelse(final_merge$size_2015==0 & final_merge$size_2016 > 0,"recruit",final_merge$Fate_2015)
final_merge$Fate_2015<-ifelse(final_merge$Fate_2016 == "recruit" | final_merge$Fate_2016 == "BC","BC",final_merge$Fate_2015)


#2014-2015
final_merge$Fate_2014<-ifelse(final_merge$size_2014==0 & final_merge$size_2015 > 0,"recruit",final_merge$Fate_2014)
final_merge$Fate_2014<-ifelse(final_merge$Fate_2015 == "recruit" | final_merge$Fate_2015 == "BC","BC",final_merge$Fate_2014)

#2013-2014
final_merge$Fate_2013<-ifelse(final_merge$size_2013==0 & final_merge$size_2014 > 0,"recruit",final_merge$Fate_2013)
final_merge$Fate_2013<-ifelse(final_merge$Fate_2014 == "recruit" | final_merge$Fate_2014 == "BC","BC",final_merge$Fate_2013)

### WRITE CSV HERE IF NECESSARY
#write.csv(final_merge, "~/Documents/Coral_Model/data/processed_data/species_data/{FILENAME}, row.names = FALSE)

# I repeated this with process for all coral species (not automated) and saved it into "~/Documents/Coral_Model/data/processed_data/species_data"
# under the name "final_merge_{species name}"
