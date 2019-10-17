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

#  Run_Setup.sh
#  Cory Haley White

#  Shell script to extract data and experimental design files 

#  Need to source modules.sh if running as a shell script.  
source /etc/profile.d/modules.sh

cd ../data/

tar -xzvf Data_Files.tar.gz
cd Peptide_Files
mv *.csv ..
cd ..

cd Vplot_Label_Files
mv *.txt ..
cd ..

rm -r Peptide_Files
rm -r Vplot_Label_Files 

cd ../bin
