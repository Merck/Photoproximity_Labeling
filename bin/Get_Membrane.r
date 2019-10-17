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

#  Get_Membrane.r
#  Cory Haley White

#  Code for generating membrane information given a database downloaded from Uniprot. 
#  This code assumes you put the Uniprot database into your data directory.  

#  Set variables
Arg<-commandArgs(trailingOnly=TRUE)

if(length(Arg)==0){
	#  No arguments, using defaults
	Database<-"Uniprot_9-21-18.tab"
	OutName<-"Uniprot_9-21-18.tab"
}else{
	Database=Arg[1]
	OutName=Arg[2]
}

print(c(Database,OutName))

#  Set directories
bindir<-paste0(getwd(),"/")
setwd("../")
WD<-getwd()
setwd(paste0(WD,"/data"))

datdir=paste0(WD,"data/")

#####  Read data  #####
Data<-read.delim(paste0(Database),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

#  Extract membrane from Uniprot list
keep<-grepl("Extracellular",Data$"Topological domain",ignore.case=TRUE)
Membrane_Uniprot<-ifelse(keep,"Known Membrane Protein","")

output<-cbind(Data,Membrane=Membrane_Uniprot)
#  Write fasta file
write.table(output,file=paste0(datdir,OutName),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
