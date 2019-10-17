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

#  Correlations_IQTMT_proteinlvl.R
#  Cory Haley White

#  Code to run normalization and correlation for peptide level data 
#  This script is designed to work with data formatted as that in the data folder, in particular peptide level data that is to be compressed to protein level data.

#  Set variables
Arg<-commandArgs(trailingOnly=TRUE)
if(length(Arg)==0){
	#  No arguments, using defaults
	Name<-"CD45_1137_p22"
	Target<-"CD45"
	#  TargetName must come from Accession.  
	#  Warning:  Some proteins don't have accession numbers and are replaced by other identifiers.  
	TargetName="P08575"
	Metric="Median"
	Normalization="Total"
	File="CD45_1137_p22_FilterPep_1.csv"
	RemoveAB=FALSE
}else{
	Name=Arg[1]
	Target=Arg[2]
	TargetName=Arg[3]
	Metric=Arg[4]
	Normalization=Arg[5]
	File=Arg[6]
	RemoveAB=Arg[7]
}

print(Arg)
#  Intensity Filter
#  For intensity data, filtering happens before switching to log2.
IntFilter=0
#  Open libraries
library(reshape)

#  Set working directory
bindir<-paste0(getwd(),"/")
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
datdir=paste0(WD,"data/")
outdir=paste0(datdir,Name,"/",Target,"/",Metric,"_",Normalization,"/")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(imgdir,recursive=TRUE,showWarnings = FALSE)

#####  Read data  #####
Data_0<-read.csv(paste0(datdir,File),stringsAsFactors=FALSE)

#  Checking peptides for proper structure
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

#  Separate into data and header info
Datapep<-as.matrix(Data_0[,8:ncol(Data_0)])
Headerpep<-as.matrix(Data_0[,4:7])

#  Normalize data, only total is implemented for now.  
if(Normalization=="Total"){
	totals<-apply(Datapep,2,sum)
	Avg_totals<-mean(totals)
	Scale_factors=Avg_totals/totals
	Datapep<-sweep(Datapep,2,Scale_factors,"*")
}else if(Normalization=="None"){
	print("Using no normalization")
}else{
	err_message<-paste0(Normalization, " normalization method not implemented, quitting")
	print(err_message)
	stop("Methods allowed are Total and None.")
}

#  Compress data from peptide level to protein level.  

Split_Data<-split(data.frame(Datapep,stringsAsFactors=FALSE),Headerpep[,1])
Split_Headerpep<-split(data.frame(Headerpep,stringsAsFactors=FALSE),Headerpep[,1])
Data<-data.frame(matrix(NA,nrow=length(Split_Data),ncol=ncol(Datapep)))
Header<-data.frame(matrix(NA,nrow=length(Split_Data),ncol=ncol(Headerpep)-1))
colnames(Data)<-colnames(Data_0)[8:ncol(Data_0)]
colnames(Header)<-colnames(Data_0)[4:6]

#  Default method is merging with median.  
if(Metric=="Geom_Avg"){
	#  For geometric average, we'll add 1 to all 0s.  
	#  The reason for this is that we want it to behave in a similar fashion to the average.
	#  We want the 0s to shift the trend towards a smaller value, but not collapse the entire average to neg inf or 0 (product and root method).  
	#  However, if there is only one entry, and it is 0, We want the old 0 back for the average and not 1.  
	for (i in 1:length(Split_Data)){
		Header[i,]<-Split_Headerpep[[i]][1,1:3]
		df1<-Split_Data[[i]]
		Fix<-df1==0			
		df1[df1==0]<-1
		df1<-log(df1,2)
		df2<-2^apply(df1,2,mean,na.rm=TRUE)
		if(nrow(df1)==1){
			df2[Fix]<-0
		}
		Data[i,]<-df2
	}
}else{
	for (i in 1:length(Split_Data)){
		Header[i,]<-Split_Headerpep[[i]][1,1:3]
		df1<-Split_Data[[i]]
		Data[i,]<-apply(df1,2,median,na.rm=TRUE)
	}
}

#  Create column of only accession number for later use in Uniprot matching.  
#  Check columns of header for string|string|string format
if(any(grepl("\\|",Header[1,]))==TRUE){
	ANcol<-which(grepl("\\|",Header[1,]))
	AC<-as.matrix(colsplit(as.character(Header[,ANcol]), split="\\|", names=c("P1","P2","P3")))
	AC<-gsub(".*:","",AC)
	AC<-gsub(" .*","",AC)
	AC<-AC[,2]
	Header<-cbind(Header,AC)
	colnames(Header)[ncol(Header)]<-"Accession"
}else{
	print("Protein names don't follow the string|string|string format.  Quitting!")
	stop()
}

#  Check if target is present.  
if(length(Header[Header[,"Accession"]==TargetName,"Accession"])==0){
	cat("Your correlation target protein is not found in your dataset!",sep="\n")
	cat(paste("Your target protein was",TargetName),sep="\n")
	cat("Proceeding without correlation measures...",sep="\n")
}

#  Filter proteins that don't reach at least the filter setting in half of the samples.  
SampFilter=ncol(Data)/2

keep<-rowSums(Data>IntFilter)>=SampFilter
Warning<-rowSums(Data<2*IntFilter)>=SampFilter
Warning<-ifelse(Warning==TRUE,"Low_Intensities","")

Header<-cbind(Header,Warning)
#  Filter low count data
Data<-Data[keep,]
Header<-Header[keep,]
output<-data.frame(Header,Data)

#  Check if target is actually present in the filtered data.  
if(length(Header[Header[,"Accession"]==TargetName,"Accession"])==0){
	cat("Your target protein did not survive filtering...",sep="\n")
	cat(paste("Your target protein was",TargetName),sep="\n")
	cat("Proceeding without correlation measures...",sep="\n")
}

#  Remove contaminants before any analyses
output<-output[!grepl("_contaminant", output[,"Protein.ID"]),]

#  Remove known antibodies if setting is given
if(RemoveAB==TRUE){
	Remove1<-grepl("IGKV",output[,"Gene.Symbol"])&grepl("Immunoglobulin",output[,"Description"])
	Remove2<-grepl("IGLV",output[,"Gene.Symbol"])&grepl("Immunoglobulin",output[,"Description"])
	Remove<-Remove1|Remove2
	output<-output[!Remove,]
}

GeneNamesCol<-which(grepl("Gene",colnames(output))&grepl("Symbol",colnames(output)))
DataCols<-grep("SUM",colnames(output))
output2<-output[,c(colnames(Header)[GeneNamesCol],colnames(Data))]
write.table(output2,file=paste(outdir,Name,"_",Target,"Target_Filtered_Short.txt",sep=""),sep="\t",row.names=FALSE)

#  Correlations
c_out<-data.frame(matrix(NA,ncol=1,nrow=(nrow(output))))
p_out<-data.frame(matrix(NA,ncol=1,nrow=(nrow(output))))

if(length(Header[Header[,"Accession"]==TargetName,"Accession"])==0){
	print("Target not found, using NAs as placeholders")
}else{
	for (i in 1:nrow(output)){
		outlist<-cor.test(as.numeric(output[i,6:ncol(output)]),as.numeric(output[output[,"Accession"]==TargetName,6:ncol(output)]))[4:3]
		c_out[i,1]<-outlist[[1]][[1]]
		p_out[i,1]<-outlist[[2]]
	}
}
colnames(c_out)[1]<-paste0(Target,"_Cor")
colnames(p_out)[1]<-paste0(Target,"_PVal")

# Merge together
s <- rep(1:ncol(c_out), each = 2) + (0:1) * ncol(c_out)
cor_out<-cbind(c_out,p_out)[s]
cor_out<-cbind(c(1:nrow(output)),cor_out)
output<-cbind(output,cor_out)
colnames(output)[ncol(output)-2]<-"Original_Order"

#  Order by P-value
output<-output[order(output[,ncol(output)]),]
p_adjust<-p.adjust(output[,ncol(output)])
output<-cbind(output,p_adjust)
colnames(output)[ncol(output)]<-"P.adjust"

#  Order by correlations
output<-output[order(output[,ncol(output)-2],decreasing=TRUE),]

#  write output
outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
write.csv(output,paste0(outdir,outName),row.names=FALSE)
