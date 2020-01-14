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

#  Run_Analysis.sh
#  Cory Haley White

#  Shell script to run R functions with variable inputs
#  The folder containing the Rscript program (Rscript.exe in windows) must be in your path for this wrapper to run, otherwise give it the full path of Rscript in your wrapper.  

#  Need to source modules.sh if running as a shell script.  
source /etc/profile.d/modules.sh

module load R/3.5.2

#  Check for variable for removing antibodies
if [ -z "$REMOVEAB" ];then 
	echo "Antibody removal was not assigned, setting to FALSE.";
	echo "If you want this changed, add in a 'export REMOVEAB=TRUE' line in your wrapper."
	REMOVEAB=FALSE
fi
#  Set special coloring list for when proteins are below a cutoff but are still colored
SPECIAL=${SPECIAL:-none}

#  Run pipeline
Rscript Filter_by_PeptideCount.r $NAME $TARGET $TARGETNAME $FILTER $NAME.csv
Rscript Correlations_IQTMT_proteinlvl.r $NAME $TARGET $TARGETNAME $METRIC $NORM $NAME"_"FilterPep_$FILTER.csv $REMOVEAB
Rscript Differential_Abundance_IQTMT_proteinlvl.r $NAME $TARGET $TARGETNAME $GROUP2 $GROUP1 $METRIC $NORM $PAIRED
Rscript MakeGeneList.r $NAME $TARGET $METRIC $NORM $LG2
Rscript VolcanoPlots_IQTMT_proteinlvl.r $NAME $TARGET $TARGETNAME $GROUP2 $GROUP1 $METRIC $NORM $LG2 $ASSOC $EXTRA $EXCLUDE $SPECIAL

export BINDIR=$(pwd)

cd ..
cd figures/$NAME/$TARGET/$METRIC"_"$NORM/Vplot/

#  Make high def tiff images
for PDF in $(ls *.pdf); do
	echo $PDF
	TMPNAME=$(echo $PDF|awk '{print substr($0, 1, length()-4)}')
	echo $TMPNAME
	convert -density 600 -compress lzw $PDF $TMPNAME.tiff
done

cd $BINDIR

Rscript Cleanup_IQTMT.r $NAME $TARGET $TARGETNAME $METRIC $NORM $SHARE
