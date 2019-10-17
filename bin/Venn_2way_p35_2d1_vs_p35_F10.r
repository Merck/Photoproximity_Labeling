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

#  Venn_2way.r
#  Cory Haley White

#  Script to compare two experiments and get overlapping results.   
#  This script is designed to work with data formatted by the correlation and differential abundance scripts.  
#  The design of this script is to take each original file for all results and filter based around a particular criteria such as Log2FC. 
#  If you have already filtered your results, just set the filter to none and pick your pre-filtered lists.  
#  Filtering by P-value is not implemented

#  Set variables

Arg<-commandArgs(trailingOnly=TRUE)

File1=Arg[1]
File2=Arg[2]
Out1<-Arg[3]
Out2<-Arg[4]
Lg2Filter=as.numeric(Arg[5])
Exclude_List<-Arg[6]
if(length(Arg)>6){
	Sub1<-Arg[7]
	Sub2<-Arg[8]
}else{
	Sub1<-NULL
	Sub2<-NULL
}

print(Arg)

#  Open libraries
library(VennDiagram)
#  Stop generating logger messages...
flog.threshold(ERROR)

#  Set working directory
bindir<-paste0(getwd(),"/")
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/Results_Comparison/")
dir.create(imgdir,recursive=TRUE,showWarnings=FALSE)

#  Read data
Files<-c(File1,File2)

if(is.null(Sub1)){
	Name1<-Out1
}else{
	Name1<-paste0(Out1,"\n",Sub1)
}

if(is.null(Sub2)){
	Name2<-Out2
}else{
	Name2<-paste0(Out2,"\n",Sub2)
}

Out<-c(Name1,Name2)

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

#  Make venn
Alpha=rep(0.5,n=length(Dlist))

#  Get overlap for center
All<-calculate.overlap(Dlist)
Overlap<-All$a3
Overlap<-paste(Overlap,collapse="\n")
Overlap<-paste0(Overlap,"\n")

#  Get differences
Diff1<-setdiff(Dlist[[1]],Dlist[[2]])
Diff2<-setdiff(Dlist[[2]],Dlist[[1]])

Cex=c(2,2,2)

#  Set rotation so that the first entry is always first in the venn.  
if(length(Dlist[[2]])>length(Dlist[[1]])){
	Rotation=180
}else{
	Rotation=0
}
#  Set position of labels
Position=c(330,30)

venn.plot<-venn.diagram(Dlist,filename=NULL,alpha=Alpha,fill=c("#E2E6FC","#D3E4D6"),ext.text=TRUE,main.cex=3.5,cat.cex=3,euler=FALSE,euler.d=FALSE,scaled=FALSE,cex=Cex,cat.pos=Position,rotation.degree=Rotation,)

#  Fix label of center, should be position 7 for 2-way venn.  
venn.plot[[7]]$label<-Overlap
#  Fix labels of left and right, should be positions 6 and 5 respectively for 2-way venn.  
if(length(Diff1)>0){
	venn.plot[[6]]$label<-Diff1
}
if(length(Diff2)>0){
	venn.plot[[5]]$label<-Diff2
}

#  Move labels further away from venn.  Should be positions 8 and 9 for 2-way venn.  
if(is.null(Sub1)){
	venn.plot[[8]]$y<-unit(0.83,"npc") 
	venn.plot[[9]]$y<-unit(0.83,"npc")
}else{
	venn.plot[[8]]$y<-unit(0.87,"npc") 
	venn.plot[[9]]$y<-unit(0.87,"npc")
}

pdf(paste0(imgdir,Out1,"_vs_",Out2,"_Venn.pdf"),width=12,height=12)
grid.draw(venn.plot)
dev.off()

#  Add in creation of tiff file if in windows.  
if(Sys.info()[[1]]=="Windows"){
	tiff(paste0(imgdir,Out1,"_vs_",Out2,"_Venn.tiff"),width=12,height=12,units="in",compression="lzw",res=600)
	grid.draw(venn.plot)
	dev.off()
}
