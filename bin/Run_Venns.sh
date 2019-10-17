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

#  Run_Venns.sh
#  Cory Haley White

#  Shell script to run venn diagram code
#  Please note that the tiff device on Linux can be problematic depending on your setup, so if on Linux, it defaults to producing a pdf file instead of a both a pdf and tiff file.  Tiff creation is done in this script using the convert command instead of within R.  

#  Need to source modules.sh if running as a shell script.  
source /etc/profile.d/modules.sh

module load R/3.5.2

export bindir=$(pwd)
cd ../data/
export WD=$(pwd)
cd ../bin

Rscript Venn_2way_p29_1_vs_p29_2.r $WD/CD45_1137_p29_1/CD45/Median_Total/CD45_1137_p29_1_CD45Target_AllAnalyses.csv $WD/CD45_1137_p29_2/CD45/Median_Total/CD45_1137_p29_2_CD45Target_AllAnalyses.csv p29_1 p29_2 1 CD45_Exclude_Venn.txt "clone HI30" "clone 2D1"

Rscript Venn_2way_p33_vs_p22.r $WD/CD45_1137_p33/CD45/Median_Total/CD45_1137_p33_CD45Target_AllAnalyses.csv $WD/CD45_1137_p22/CD45/Median_Total/CD45_1137_p22_CD45Target_AllAnalyses.csv HRP Diazirine 1.32 CD45_Exclude_Venn.txt

Rscript Venn_2way_p35_2d1_vs_p35_F10.r $WD/CD45_1137_p35_2d1/CD45/Median_Total/CD45_1137_p35_2d1_CD45Target_AllAnalyses.csv $WD/CD45_1137_p35_F10/CD45/Median_Total/CD45_1137_p35_F10_CD45Target_AllAnalyses.csv p35_2d1 p35_F10 1 CD45_Exclude_Venn.txt

Rscript Venn_3way.r Supplemental_Data_9-13-19_3wayVenn.xlsx CD30_p69 CD300A_p69 PDL1_p36_10m CD30 CD300A PDL1 1.32 1 2.6 PDL1_Exclude.txt Venn_3way_Include.txt

Rscript Venn_4way.r $WD/CD45_1137_p22/CD45/Median_Total/CD45_1137_p22_CD45Target_AllAnalyses.csv $WD/CD45_1137_p27/CD45/Median_Total/CD45_1137_p27_CD45Target_AllAnalyses.csv $WD/CD45_1137_p29_1/CD45/Median_Total/CD45_1137_p29_1_CD45Target_AllAnalyses.csv $WD/CD45_1137_p32/CD45/Median_Total/CD45_1137_p32_CD45Target_AllAnalyses.csv p22 p27 p29 p32 1 CD45_Exclude_Venn.txt

#  Make high def tiff images

cd ../figures/Results_Comparison
for PDF in $(ls *.pdf); do
	echo $PDF
	TMPNAME=$(echo $PDF|awk '{print substr($0, 1, length()-4)}')
	echo $TMPNAME
	convert -density 600 -compress lzw $PDF $TMPNAME.tiff
done

cd Venn_4way
for PDF in $(ls *.pdf); do
	echo $PDF
	TMPNAME=$(echo $PDF|awk '{print substr($0, 1, length()-4)}')
	echo $TMPNAME
	convert -density 600 -compress lzw $PDF $TMPNAME.tiff
done
cd ../Venn_3way
for PDF in $(ls *.pdf); do
	echo $PDF
	TMPNAME=$(echo $PDF|awk '{print substr($0, 1, length()-4)}')
	echo $TMPNAME
	convert -density 600 -compress lzw $PDF $TMPNAME.tiff
done
