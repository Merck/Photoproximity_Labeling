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

#  Make_Sup_Table.r
#  Cory Haley White

#  Code to add in membrane information to supplemental table

#  Set variables
Arg<-commandArgs(trailingOnly=TRUE)
if(length(Arg)==0){
	#  No arguments, using defaults
	Metric="Median"
	Norm="Total"
}else{
	Metric=Arg[1]
	Norm=Arg[2]
	Database=Arg[3]
}

print(Arg)
#  Open libraries
library(xlsx)

#  Set directories
setwd("../data")
WD<-getwd()

if(substr(WD,nchar(WD),nchar(WD))!="/"){
	WD<-paste0(WD,"/")
}

#  Check if database was given
if(!exists("Database")){
	print("No Database was given, defaulting to data directory")
	Database=paste0(WD,"/Uniprot_9-21-18_wHMPAS.tab")
}

Dirs<-list.dirs(WD,recursive=FALSE)
Dirs<-Dirs[grep("_p",Dirs)]

#  Set names and targets
Names<-gsub(WD,"",Dirs)
Names<-gsub("1137_","",Names)
Names<-gsub("/","",Names)

Targets<-strsplit(Names,"_")
Targets<-do.call(rbind,Targets)
Targets<-Targets[,1]

#  Read membrane data
Uniprot<-read.delim(paste0(Database),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

#  Fix PTPRCAP calls to be membrane.  Database update.  
Fix<-grep("PTPRCAP",Uniprot[,"Gene names"])
Uniprot[Fix,"Membrane"]<-"Known Membrane Protein"

#  Splitting data into two files and merging manually due to RAM problems writing all to a single file.  
for(i in 1:(length(Dirs)-9)){
	Files<-list.files(paste0(Dirs[i],"/",Targets[i],"/",Metric,"_",Norm),full.names=TRUE)
	File<-Files[grep("AllAnalyses.csv",Files)]
	Data<-read.csv(File,stringsAsFactors=FALSE)
	Utmp<-Uniprot[match(Data$Accession,Uniprot$Entry),]
	Utmp<-Utmp[,"Membrane"]
	#  Restricting membrane call to only those IDed in Uniprot and not other possible sources.  
	Utmp[Utmp=="Membrane_HMPAS"]<-""
	Utmp[Utmp=="Potential_MP"]<-""
	Data<-cbind(Data,Membrane=Utmp)
	Col<-grep("_Cor",colnames(Data))
	Row<-which(Data[,Col]==1)
	if(length(Row)>0){
		Dat2<-Data[Row,]
		Data<-Data[-Row,]
		Data<-Data[order(Data$logFC,decreasing=TRUE),]
		Data<-rbind(Dat2,Data)
	}else{
		Data<-Data[order(Data$logFC,decreasing=TRUE),]
	}
	output<-Data[,c("Gene.Symbol","Accession","logFC","P.Value","adj.P.Val","Membrane")]
	Sheet=paste0(Names[i])
	print(Sheet)
	write.xlsx(output,file=paste0(WD,"Supplemental_Data_1.xlsx"),sheetName=Sheet,append=TRUE,col.names=TRUE,row.names=FALSE)
}

for(i in (length(Dirs)-8):length(Dirs)){
	Files<-list.files(paste0(Dirs[i],"/",Targets[i],"/",Metric,"_",Norm),full.names=TRUE)
	File<-Files[grep("AllAnalyses.csv",Files)]
	Data<-read.csv(File,stringsAsFactors=FALSE)
	Utmp<-Uniprot[match(Data$Accession,Uniprot$Entry),]
	Utmp<-Utmp[,"Membrane"]
	#  Restricting membrane calls to only those IDed in Uniprot and not other possible sources.  
	Utmp[Utmp=="Membrane_HMPAS"]<-""
	Utmp[Utmp=="Potential_MP"]<-""
	Data<-cbind(Data,Membrane=Utmp)
	Col<-grep("_Cor",colnames(Data))
	Row<-which(Data[,Col]==1)
	if(length(Row)>0){
		Dat2<-Data[Row,]
		Data<-Data[-Row,]
		Data<-Data[order(Data$logFC,decreasing=TRUE),]
		Data<-rbind(Dat2,Data)
	}else{
		Data<-Data[order(Data$logFC,decreasing=TRUE),]
	}
	output<-Data[,c("Gene.Symbol","Accession","logFC","P.Value","adj.P.Val","Membrane")]
	Sheet=paste0(Names[i])
	print(Sheet)
	write.xlsx(output,file=paste0(WD,"Supplemental_Data_2.xlsx"),sheetName=Sheet,append=TRUE,col.names=TRUE,row.names=FALSE)
}
