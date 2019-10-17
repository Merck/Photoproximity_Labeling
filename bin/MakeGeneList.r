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

#  MakeGeneList.r
#  Cory Haley White

#  Create a gene list based upon filtering criteria for input into MetaCore and IPA
#  Requires that you generated data with the correlation and differential protein abundance scripts.  

#  Set variables
Arg<-commandArgs(trailingOnly=TRUE)
if(length(Arg)==0){
	#  No arguments, using defaults
	Name<-"CD45_1137_p21"
	Target<-"CD45"
	Metric="Median"
	Normalization="Total"
	LogCutoff=1
}else{
	Name=Arg[1]
	Target=Arg[2]
	Metric=Arg[3]
	Normalization=Arg[4]
	LogCutoff=as.numeric(Arg[5])
}

print(Arg)

#  Set working directory
bindir<-paste0(getwd(),"/")
setwd("..")
WD<-paste0(getwd(),"/")
#  Get directories
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=datadir

#  Get file from correlation script
FilesList<-list.files(paste0(datadir))
File<-FilesList[grep("AllAnalyses.csv",FilesList)]

#  Read data.  
Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE)
#  Remove contaminants if present
Data_0<-Data_0[!grepl("contaminant", Data_0[,"Protein.ID"]),]

#  Filter dataset by criteria
#  Log2FC Cutoff
Data<-Data_0[Data_0$logFC>=LogCutoff,]

#  Get genes.
Genes<-data.frame(Genes=Data$Gene.Symbol,stringsAsFactors=FALSE)

write.table(Genes,file=paste0(outdir,Name,"_lg2=",LogCutoff,".txt"),sep="/t",row.names=FALSE,quote=FALSE)
