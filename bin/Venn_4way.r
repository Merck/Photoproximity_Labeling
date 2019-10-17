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

#  Venn_4way.r
#  Cory Haley White

#  Script to compare four experiments and get overlapping results.   
#  This script is designed to work with data formatted by the correlation and differential abundance scripts.  
#  The design of this script is to take each original file for all results and filter based around a particular criteria such as Log2FC. 
#  If you have already filtered your results, just set the filter to none and pick your pre-filtered lists.  
#  Filtering by P-value is not implemented

#  Set variables

Arg<-commandArgs(trailingOnly=TRUE)

File1=Arg[1]
File2=Arg[2]
File3=Arg[3]
File4=Arg[4]
Out1<-Arg[5]
Out2<-Arg[6]
Out3<-Arg[7]
Out4<-Arg[8]
Lg2Filter=as.numeric(Arg[9])
Exclude_List<-Arg[10]

print(Arg)

#  Open libraries
library(VennDiagram)
#  Stop generating logger messages...
flog.threshold(ERROR)

#  Set working directory
bindir<-paste0(getwd(),"/")
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/Results_Comparison/Venn_4way/")
dir.create(imgdir,recursive=TRUE,showWarnings=FALSE)

#  Read data
Files<-c(File1,File2,File3,File4)
Out<-c(Out1,Out2,Out3,Out4)

if(file.exists(paste0(WD,"data/",Exclude_List))){
	Exclude_List<-read.table(paste0(WD,"data/",Exclude_List),stringsAsFactors=FALSE,sep="\t")
}else{
	Exclude_List<-NULL
}

Dlist<-list()

for(i in 1:length(Files)){
	Dat<-read.csv(Files[i],stringsAsFactors=FALSE)
	if(Lg2Filter!="none"){
		Dat<-Dat[Dat$logFC>=Lg2Filter,]
	}	
	Dat<-Dat[!Dat$Gene.Symbol%in%Exclude_List[,1],]
	Dlist[[Out[i]]]<-Dat$Gene.Symbol
}

#  Fix names for venns later
Fixme1<-which(names(Dlist)==Out1)
Fixme2<-which(names(Dlist)!=Out1)
names(Dlist)[Fixme1]<-"rep1"
Fixed<-c("rep4","rep2","rep3")
names(Dlist)[Fixme2]<-Fixed

#  Make venn
Alpha=rep(0.5,n=length(Dlist))

#  Get overlap for center
Overlap<-calculate.overlap(Dlist)$a6
Overlap<-paste(Overlap,collapse="\n")
Overlap<-paste0(Overlap,"\n")

Cex=c(rep(3.5,5),1.7,rep(3.5,9))

venn.plot<-venn.diagram(Dlist,filename=NULL,alpha=Alpha,fill=c("#FFF2DE","#EEFBEE","#E7EEFF","#FFFFE3"),ext.text=FALSE,cex=Cex,main.cex=3.5,cat.cex=3.5)
#  Fix label of center, should be position 14 for 4-way venn.  
venn.plot[[14]]$label<-Overlap
venn.plot[[14]]$y<-unit(0.422,"npc")

pdf(paste0(imgdir,"4_way_Venn.pdf"),width=12,height=12)
grid.draw(venn.plot)
dev.off()

#  Add in creation of tiff file if in windows.  
if(Sys.info()[[1]]=="Windows"){
	tiff(paste0(imgdir,"4_way_Venn.tiff"),width=12,height=12,units="in",compression="lzw",res=600)
	grid.draw(venn.plot)
	dev.off()
}
