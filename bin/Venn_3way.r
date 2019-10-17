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

#  Venn_3way.r
#  Cory Haley White

#  Script to compare three experiments and get overlapping results.   
#  This script is designed to work with data formatted by the correlation and differential abundance scripts.  
#  The design of this script is to take each original file for all results and filter based around a particular criteria such as Log2FC. 
#  If you have already filtered your results, just set the filter to none and pick your pre-filtered lists.  
#  Filtering by P-value is not implemented

#  Set variables

Arg<-commandArgs(trailingOnly=TRUE)

File=Arg[1]
Sheet1<-Arg[2]
Sheet2<-Arg[3]
Sheet3<-Arg[4]
Out1<-Arg[5]
Out2<-Arg[6]
Out3<-Arg[7]
Lg2_1<-as.numeric(Arg[8])
Lg2_2<-as.numeric(Arg[9])
Lg2_3<-as.numeric(Arg[10])
Exclude_List<-Arg[11]
Include_List<-Arg[12]

print(Arg)

#  Open libraries
library(VennDiagram)
library(xlsx)
#  Stop generating logger messages...
flog.threshold(ERROR)

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd("..")
WD<-paste0(getwd(),"/")
datadir<-paste0(WD,"data/")

imgdir<-paste0(WD,"figures/Results_Comparison/Venn_3way/")
dir.create(imgdir,recursive=TRUE,showWarnings=FALSE)

#  Read data
Sheets<-c(Sheet1,Sheet2,Sheet3)
Out<-c(Out1,Out2,Out3)
Lgs<-c(Lg2_1,Lg2_2,Lg2_3)

if(file.exists(paste0(WD,"data/",Exclude_List))){
	Exclude_List<-read.table(paste0(WD,"data/",Exclude_List),stringsAsFactors=FALSE,sep="\t")
}else{
	Exclude_List<-NULL
}

if(file.exists(paste0(WD,"data/",Include_List))){
	Include_List<-read.table(paste0(WD,"data/",Include_List),stringsAsFactors=FALSE,sep="\t")
}else{
	Include_List<-NULL
}

Dlist<-list()

for(i in 1:length(Sheets)){
	Dat<-read.xlsx(paste0(datadir,File),sheetName=Sheets[i],stringsAsFactors=FALSE)
	Lg2Filter<-Lgs[i]
	#  If filter given, use it.  
	if(Lg2Filter!="none"){
		Dat<-Dat[Dat$logFC>=Lg2Filter,]
	}	
	Dat<-Dat[!Dat$Gene.Symbol%in%Exclude_List[,1],]
	Keep1<-Dat[Dat$Membrane=="Known Membrane Protein","Gene.Symbol"]
	Keep2<-Dat[Dat$Gene.Symbol%in%Include_List[,1],"Gene.Symbol"]
	Genes<-unique(c(Keep1,Keep2))
	Dlist[[i]]<-Genes
	names(Dlist)[i]<-Out[i]
}

#  Make venn
Alpha=rep(0.5,n=length(Dlist))

#  Get overlap for center (a5 in 3-way venn)
All<-calculate.overlap(Dlist)
Overlap<-All$a5
Overlap<-paste(Overlap,collapse="\n")
Overlap<-paste0(Overlap,"\n")

Cex=c(3,3,3,3,1.7,3,3)

venn.plot<-venn.diagram(Dlist,filename=NULL,alpha=Alpha,fill=c("#E2E6FC","#D3E4D6","#FFFFE3"),ext.text=TRUE,main.cex=3.5,cat.cex=3,euler=FALSE,euler.d=FALSE,scaled=FALSE,cex=Cex,)

#  Fix label of center, should be position 11 for 3-way venn.  
venn.plot[[11]]$label<-Overlap
venn.plot[[11]]$y<-unit(0.53,"npc")

pdf(paste0(imgdir,Out1,"_",Out2,"_",Out3,"_Venn.pdf"),width=12,height=12)
grid.draw(venn.plot)
dev.off()

#  Add in creation of tiff file if in windows.  
if(Sys.info()[[1]]=="Windows"){
	tiff(paste0(imgdir,Out1,"_",Out2,"_",Out3,"_Venn.tiff"),width=12,height=12,units="in",compression="lzw",res=600)
	grid.draw(venn.plot)
	dev.off()
}
