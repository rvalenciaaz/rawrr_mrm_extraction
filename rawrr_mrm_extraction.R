library(tidyverse)
library(rawrr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

#define the folder path where the RAW files are present 
folder_path<-"14_Sep_Bacillus_Iter_0/"
#remove the last character
f2<-str_sub(folder_path,end=-2)

#create a empty list for the RAW filenames
names<-c()
#append the names to the list
for (i in  list.files(path = folder_path)){
    if (endsWith(i,".raw")){
        names<-append(names,i)
    }
}
#the data extraction of one file is performed by calling the extract function 
#for the file (insert the filename in the function argument)
extract<-function(name){
	#generate a path for the file
    rawfile<-paste(folder_path,name, sep="/")
    #use the rawrr package to give an index for the RAW file, the index is like a summary, where scan information 
    #is displayed for the whole acquisition time
    index<-rawrr::readIndex(rawfile = rawfile)
    
    #get the unique scan names in the file, the unique scan names correspond simply to a precursor mass and the defined fragments
    #for that precursor in the MS method
    scans<-unique(index$scanType)

    #create an empty dataframe to allocate data, define column names
	megatable <- data.frame(matrix(ncol = 3, nrow = 0))
    x <- c("rt", "tic","scan_type")
    colnames(megatable) <- x
    megatable$scan_type<-as.character(megatable$scan_type)
    #iterate over the unique scan names list and get the specific information for each scan associated with that name. 
    #The number of scans in the file associated with a precursor/fragments mass is defined in the MS method, as cycles or dwell time.
    for (i in scans){
    	#retrieve information of all the scans matching the scan name	
    	subi<-rawrr::readIndex(rawfile = rawfile) |>
    	subset(scanType == i)
    	#get a list of the scans ids    
    	listscan<-subi$scan
    	#get the MS2 spectrum associated to that scan using rawrr
    	S<- rawrr::readSpectrum(rawfile = rawfile, scan = listscan)
    	#iterate over the spectrum and retrieve the TIC (sum of counts for the fragments of a specific precursor) and the recorder time for that scan.
    	for (j in 1:length(listscan)) {
    		datos<-S[[j]]
    		ticn<-datos$TIC
    		sec<-datos$rtinseconds
    		#add to table
    		megatable<-megatable %>% add_row(rt = sec, tic = ticn,scan_type = i)
    	}

    }
    #add the filename to table as an additional column
    megatable$file<-name
    #return the table
    return(megatable)
}
#create a empty list to keep the the extracted data as tables for each file
probar<-list()
#iterate over filenames and add the extracted data table to the list
for (j in 1:length(names)){
  probar[[j]] <- extract(names[j])
}

#merge the tables vertically
hypertable<-do.call("rbind", probar)
#write to computer
write_csv(hypertable, paste(f2,"_hypertable.csv",sep=""))