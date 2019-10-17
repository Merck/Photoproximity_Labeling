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

#  Cleanup_IQTMT.r
#  Cory Haley White
#  Clean up script to compress output for sending.  
#  This script must be run at the end as it requires prior information to build the output files.  

#  Set variables
Arg<-commandArgs(trailingOnly=TRUE)

Name=Arg[1]
Target=Arg[2]
TargetName=Arg[3]
Metric=Arg[4]
Normalization=Arg[5]

print(Arg)

#  Set working directory
bindir<-paste0(getwd(),"/")

if(length(Arg[6])==0){
	print("You didn't provide a variable for your shared directory, defaulting to git data directory")
	setwd("../data")
	Tempdir<-getwd()
	Sharedir=paste0(Tempdir,"/Compressed_Output/")
}
Sharedir=Arg[6]

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/Vplot/")
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=paste0(Sharedir)

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

#  Compress final output for sending.  
setwd(datadir)
files1<-dir(datadir)
files1<-files1[grepl(".csv",files1)|grepl(".txt",files1)]
zipfile1<-paste(Name,"_",Target,"Target_",Metric,"_",Normalization,"_Datafiles",sep="")
zip(zipfile=zipfile1,files=files1)
file.rename(from=paste(datadir,zipfile1,".zip",sep=""),to=paste(outdir,zipfile1,".zip",sep=""))

setwd(imgdir)
files2<-dir(imgdir)
files2<-files2[grepl(".pdf",files2)|grepl(".tiff",files2)]
zipfile2<-paste0(Name,"_",Metric,"_",Normalization,"_Images")
zip(zipfile=zipfile2,files=files2)
file.rename(from=paste0(imgdir,zipfile2,".zip"),to=paste0(outdir,zipfile2,".zip"))
