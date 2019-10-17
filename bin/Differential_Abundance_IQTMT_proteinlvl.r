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

#  Differential_Abundance_IQTMT_proteinlvl.r
#  Cory Haley White

#  Performs differential abundance measurements given an experiment with biological replicates  
#  This script is designed to work with data formatted as that in the data folder.  For other formats, this script needs to be altered.  
#  This script must be run after correlations as it requires that information to build the output files.  
#  This script also requires an experimental design file to indicate which samples belong to which group (See Experimental_Design_Files folder for examples).  

#  Set variables
Arg<-commandArgs(trailingOnly=TRUE)
if(length(Arg)==0){
	#  No arguments, using defaults
	Name<-"CD45_1137_p22"
	Target<-"CD45"	
	#  TargetName must come from Accession.  
	#  Warning:  Some proteins don't have accession numbers and are replaced by other identifiers.  
	TargetName="P08575"
	Group2="Target"
	Group1="Iso"
	Metric="Median"
	Normalization="Total"
	Paired="NotPaired"
}else{
	Name=Arg[1]
	Target=Arg[2]
	TargetName=Arg[3]
	Group2=Arg[4]
	Group1=Arg[5]
	Metric=Arg[6]
	Normalization=Arg[7]
	Paired=Arg[8]
}

print(Arg)
#  Open libraries
library(limma)
library(stringr)

#  Set working directory
bindir<-paste0(getwd(),"/")
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=datadir
EDdir<-paste0(WD,"data/Experimental_Design_Files/")
datadir2=paste0(WD,"data/")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(imgdir,recursive=TRUE,showWarnings = FALSE)

#  Get file from correlation script
FilesList<-list.files(paste0(datadir))
File<-FilesList[grep("AllAnalyses.csv",FilesList)]

#  Read data.  
Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE)
ED<-read.csv(paste0(EDdir,"ED_",Name,".csv"))

#  Select data columns
Start<-grep("Warning",colnames(Data_0))+1
End<-grep("Original_Order",colnames(Data_0))-1
Data<-Data_0[,c(Start:End)]

#  Replace 0 with min
Data[Data==0]<-min(Data[Data!=0])
rownames(Data)<-Data_0[,1]

Datalg2<-log(Data,2)

#  Data was previously normalized in correlation script.  

Group<-factor(ED$Group,levels=c(Group1,Group2))

if(Paired=="Paired"){
	Donor=factor(ED$Donor)
	Design=model.matrix(~Group+Donor)
}else{
	Design=model.matrix(~Group)
}

print("Factors made, fitting next.")

fit<-lmFit(Datalg2,Design)
p.fit<-eBayes(fit)
LMout<-topTable(p.fit,n=nrow(p.fit),coef=2,sort.by="none")

print("Fitting done")

#  Merge with correlation data
output<-cbind(Data_0,LMout)

#  write output
outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
write.csv(output,paste0(outdir,outName),row.names=FALSE)
