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

#  VolcanoPlots_IQTMT_proteinlvl.r
#  Cory Haley White

#  Create volcano plots with different custom settings.  
#  See examples below for required arguments.  
#  This script is designed to work with data formatted as that in the data folder.  For other formats, this script needs to be altered.  
#  This script must be run after correlations and differential abundance measurements as it requires that information to build the output files.  

#  Set variables
set.seed(1945)
Arg<-commandArgs(trailingOnly=TRUE)
if(length(Arg)==0){
	#  No arguments, using defaults
	Name<-"CD45_1137_p22"
	Target<-"CD45"	#  TargetName must come from Accession.  
	#  Warning:  Some proteins don't have accession numbers and are replaced by other identifiers.  
	TargetName="P08575"
	Group2="Target"
	Group1="Iso"
	Metric="Median"
	Normalization="Total"
	LogCutoff=1
	Assoc_List<-"CD45_Associated_Proteins.txt"
	Extra_List<-"CD45_Extra_Labels.txt"
	Exclude_List<-"CD45_Exclude.txt"
}else{
	Name=Arg[1]
	Target=Arg[2]
	TargetName=Arg[3]
	Group2=Arg[4]
	Group1=Arg[5]
	Metric=Arg[6]
	Normalization=Arg[7]
	LogCutoff=as.numeric(Arg[8])
	if(length(Arg)>8){
		Assoc_List<-Arg[9]
	}
	if(length(Arg)>9){
		Extra_List<-Arg[10]
	}
	if(length(Arg)>10){
		Exclude_List<-Arg[11]
	}
}

print(Arg)
#  Open libraries
library(ggplot2)
library(ggrepel)
library(extrafont)
#  This command below only needs to be run once on the system after installation of the extrafonts package.  
#loadfonts()

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/Vplot/")
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=datadir

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(imgdir,recursive=TRUE,showWarnings = FALSE)

#  Get file from correlation script
FilesList<-list.files(paste0(datadir))
File<-FilesList[grep("AllAnalyses.csv",FilesList)]

#  Read data.  
Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE)
if(file.exists(paste0(WD,"data/",Assoc_List))){
	AssocP<-read.table(paste0(WD,"data/",Assoc_List),stringsAsFactors=FALSE,sep="\t")
}else{
	print("Association file not found")
	AssocP<-NULL
}

if(file.exists(paste0(WD,"data/",Extra_List))){
	Extra_List<-read.table(paste0(WD,"data/",Extra_List),stringsAsFactors=FALSE,sep="\t")
}else{
	print("Extra label file not found")
	Extra_List<-NULL
}

if(file.exists(paste0(WD,"data/",Exclude_List))){
	Exclude_List<-read.table(paste0(WD,"data/",Exclude_List),stringsAsFactors=FALSE,sep="\t")
}else{
	print("Exclusion file not found")
	Exclude_List<-NULL
}

Labels_List<-rbind(AssocP,Extra_List)

#  Filter dataset by suspected contaminants in removal list
Data_0<-Data_0[!Data_0$Gene.Symbol%in%Exclude_List[,1],]

#  Filter dataset by criteria
#  Log2FC Cutoff
Data<-Data_0[Data_0$logFC>=LogCutoff,]

if(nrow(Data)==0){
	Message<-paste0("WARNING!  Nothing survived your logFC cutoff of ", LogCutoff,".  Quitting!")
}

#  Write output data.  
Dataout<-Data[,1:5]
write.csv(Dataout,file=paste0(datadir,Name,"_",Target,"_lg2=",signif(LogCutoff,3),"_Vplot.csv"),row.names=FALSE)

#  Make volcano plots of log2FC values.  
print("Making volcano plots")
#  Generate volcano plot for all data with selected logFC cutoffs.  
Title=paste0(Group2," vs. ", Group1)
logFC_All<-Data_0$logFC
Log10_All<--log(Data_0$P.Value,10)

#  Set fonts
System<-Sys.info()[[1]]
if(System=="Linux"){
	FONT<-"Liberation Sans"
}else{
	FONT<-"Arial"
}

#  Get antibody information if present
AB1<-grepl("IGKV",Data_0[,"Gene.Symbol"])&grepl("Immunoglobulin",Data_0[,"Description"])
AB2<-grepl("IGLV",Data_0[,"Gene.Symbol"])&grepl("Immunoglobulin",Data_0[,"Description"])
ABs<-AB1|AB2

#  Colors
Colors<-c("#800000","#68BC44","#BE29EC","black")
#  Select p-value cutoff by FDR<=0.05
PVal<-Data_0[which(abs(Data_0$adj.P.Val-0.05)==min(abs(Data_0$adj.P.Val-0.05))),"P.Value"]
#  Select max P-value if there's multiple adjusted p-value matches as can happen with FDR method.  
PVal<-max(PVal)
Color<-ifelse(Data_0$logFC>LogCutoff&Data_0$P.Value<=PVal,"Enriched","NS")
Color<-ifelse(Color=="Enriched"&Data_0[,"Gene.Symbol"]%in%AssocP[,1],"Associated",Color)
Color<-ifelse(Color=="Enriched"&ABs,"Antibody",Color)

#  Set up plot for inclusion of antibodies
#  Get labels
Labels_All=Data_0$Gene.Symbol

if(sum(ABs)>0){
	#  Add in antibodies to labels if present  
	Labels_List<-rbind(Labels_List,data.frame(V1=Data_0[ABs,"Gene.Symbol"],stringsAsFactors=FALSE))
	ABNames<-Data_0[ABs,"Gene.Symbol"]
}

Labels_All<-ifelse(Labels_All%in%Labels_List[,1],as.character(Labels_All),"")
#  Remove label if hit isn't significant
Labels_All<-ifelse(Color=="NS","",Labels_All)

#  Make data frame
DFAll<-data.frame(logFC=logFC_All,Log10=Log10_All,Color=Color,GeneLabels=Labels_All,stringsAsFactors=FALSE,check.names=FALSE)

#  Get range for logFC and -log10p with antibody labels
if(sum(ABs)>0){
	Max<-max(DFAll$logFC)
	Min<-min(DFAll$logFC)

	if(abs(Max)>abs(Min)){
		xlim=Max
	}else{
		
	}

	xlim<-max(abs(Min),abs(Max))
	xlim<-round(xlim)

	if(xlim<Max|xlim>Min){
		xlim<-xlim+0.5
	}
	xlimf<-floor(xlim)
	ylim<-ceiling(max(DFAll$Log10))
	ylim<-ylim+1

	#  Gene labeling for all above filter
	ggvol<-ggplot(data=DFAll,aes(x=logFC,y=Log10,color=Color))+
	geom_point(alpha=0.8, size=8) +
	labs(title=NULL,family="Arial")+
	coord_cartesian(ylim = c(0, ylim), xlim=c(Min,Max+0.1*Max),expand = TRUE,)+
	xlab(bquote("log"[2]~'fold change')) + ylab(bquote("-log"[10]~ "p-value"))+
	scale_color_manual(values=c("Antibody"="#800000","Associated"="#68BC44","Enriched"="#BE29EC","NS"="black"))+
	scale_x_continuous(breaks=c(-xlimf:xlimf))+
	scale_y_continuous(breaks=c(0:ylim))+
	geom_label_repel(size=12,aes(label=GeneLabels),color="black",box.padding=0.3,point.padding=.8,segment.color="grey50",fill=NA,label.size=NA,family=FONT,fontface="bold")+
	theme_classic(base_size=50)+
	theme(legend.position = "none",axis.title.y=element_text(family=FONT),axis.title.x=element_text(family=FONT))+
	guides(fill=guide_legend(title=""))

	ggsave(paste0(imgdir,Name,"_",Group2,"_vs_",Group1,"_Lg2FC_",signif(LogCutoff,3),"_VPlot_wAB.pdf"),width=16,height=10,units="in",device="pdf",dpi=600)
}

#  Generate plot without antibody data

DFAll<-DFAll[DFAll$Color!="Antibody",]

Max<-max(DFAll$logFC)
Min<-min(DFAll$logFC)

if(abs(Max)>abs(Min)){
	xlim=Max
}else{
	
}

xlim<-max(abs(Min),abs(Max))
xlim<-round(xlim)

if(xlim<Max|xlim>Min){
	xlim<-xlim+0.5
}
xlimf<-floor(xlim)
ylim<-ceiling(max(DFAll$Log10))
ylim<-ylim+1

#  Gene labeling for all above filter
ggvol<-ggplot(data=DFAll,aes(x=logFC,y=Log10,color=Color))+
geom_point(alpha=0.8, size=8) +
labs(title=NULL,family="Arial")+
coord_cartesian(ylim = c(0, ylim), xlim=c(Min,Max+0.1*Max),expand = TRUE,)+
xlab(bquote("log"[2]~'fold change')) + ylab(bquote("-log"[10]~ "p-value"))+
scale_color_manual(values=c("Associated"="#68BC44","Enriched"="#BE29EC","NS"="black"))+
scale_x_continuous(breaks=c(-xlimf:xlimf))+
scale_y_continuous(breaks=c(0:ylim))+
geom_label_repel(size=12,aes(label=GeneLabels),color="black",box.padding=0.3,point.padding=.8,segment.color="grey50",fill=NA,label.size=NA,family=FONT,fontface="bold")+
theme_classic(base_size=50)+
theme(legend.position = "none",axis.title.y=element_text(family=FONT),axis.title.x=element_text(family=FONT))+
guides(fill=guide_legend(title=""))

pdf(NULL)
ggsave(paste0(imgdir,Name,"_",Group2,"_vs_",Group1,"_Lg2FC_",signif(LogCutoff,3),"_VPlot.pdf"),width=16,height=10,units="in",device="pdf",dpi=600)
