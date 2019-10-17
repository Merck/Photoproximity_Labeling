# Copyright (c) 2019 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
#
# This file is part of the Protein Proximity Labeling program.
#
# Protein Proximity Labeling is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#  Filter_by_PeptideCount.r
#  Cory Haley White

#  Code for filtering data by a peptide count.
#  This code is designed to work with data formatted as those in the data folder.  Other formats will need adjustment of the code.  

#  Set variables
Arg<-commandArgs(trailingOnly=TRUE)
if(length(Arg)==0){
	#  No arguments, using defaults
	Name<-"CD45_1137_p22"
	Target<-"CD45"
	#  TargetName must come from Accession.  
	#  Warning:  Some proteins don't have accession numbers and are replaced by other identifiers.  
	TargetName="P08575"
	Filter=1
	File="CD45_1137_p22.csv"
}else{
	Name=Arg[1]
	Target=Arg[2]
	TargetName=Arg[3]
	Filter=as.numeric(Arg[4])
	File=Arg[5]
}

print(Arg)

#  Open libraries
library(reshape)

#  Set working directory
bindir<-paste0(getwd(),"/")
setwd("..")
WD<-paste0(getwd(),"/")

datdir=paste0(WD,"data/")
outdir=paste0(datdir)

#####  Read data  #####
Data_0<-read.csv(paste0(datdir,File),stringsAsFactors=FALSE,check.names=FALSE)

#  Remove peptides not found in Uniprot
if(any(grepl("\\|",Data_0[1,]))==TRUE){
	ANcol<-which(grepl("\\|",Data_0[1,]))
	checkpep<-as.matrix(colsplit(as.character(Data_0[,ANcol]), split="\\|", names=c("P1","P2","P3")))
	checkpep<-gsub(".*:","",checkpep)
	checkpep<-gsub(" .*","",checkpep)
	checkpep<-checkpep[,2]
}else{
	print("Peptide names don't follow the string|string|string format.  Quitting!")
	stop()
}

Data<-data.frame(Data_0,checkpep,stringsAsFactors=FALSE,check.names=FALSE)

#  Count up peptides
PepCounts<-table(Data$checkpep)
Keep<-PepCounts>Filter

#  Check if target is present
#  Keep target even if it only has a single peptide  
Check<-Data[Data[,"checkpep"]==TargetName,"checkpep"]
if(length(Check)==0){
	cat("Your correlation target protein is not found in your dataset!",sep="\n")
	cat(paste("Your target protein was",TargetName),sep="\n")
	TargetKeep<-rep(FALSE,length(PepCounts))
}else{
	TargetKeep<-names(PepCounts)==unique(Check)
}
Keep<-Keep|TargetKeep

output<-Data[Data[,"checkpep"]%in%names(PepCounts)[Keep],1:(ncol(Data)-1)]

write.csv(output,file=paste0(outdir,Name,"_FilterPep_",Filter,".csv"),row.names=FALSE)
