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

#  Run_All.sh
#  Cory Haley White

#  Large wrapper to run all experiments.  
#  Run from the /bin directory of the git repo.  

#  Need to source modules.sh if running as a shell script.  
source /etc/profile.d/modules.sh

module load R/3.5.2

#  Set data directory and directory for storing compressed output
cd ../data
export DATDIR=$(pwd)
export SHARE=$DATDIR/TempSend/
mkdir -p $SHARE
cd ../bin

#  Run setup script to decompress data files
bash Run_Setup.sh

#  CD45 Runs
bash Run_All_CD45_1137_p22.sh
bash Run_All_CD45_1137_p27.sh
bash Run_All_CD45_1137_p29_1.sh
bash Run_All_CD45_1137_p29_2.sh
bash Run_All_CD45_1137_p32.sh
bash Run_All_CD45_1137_p33.sh
bash Run_All_CD45_1137_p35_2d1.sh
bash Run_All_CD45_1137_p35_F10.sh
bash Run_All_CD45_1137_p66.sh
bash Run_All_CD45_1137_p70.sh

#  PDL1 Runs
bash Run_All_PDL1_1137_p30.sh
bash Run_All_PDL1_1137_p31_LG2_2.sh
bash Run_All_PDL1_1137_p36_2m.sh
bash Run_All_PDL1_1137_p36_10m.sh

#  CD29 Runs
bash Run_All_CD29_1137_p32.sh

#  CD44 Runs
bash Run_All_CD44_1137_p56.sh

#  CD47 Runs
bash Run_All_CD47_1137_p27.sh

#  EGFR Runs
bash Run_All_EGFR_1137_p56.sh

#  CD30 Runs
bash Run_All_CD30_1137_p69.sh

#  CD300A Runs
bash Run_All_CD300A_1137_p69.sh

#  Make Venn diagrams
bash Run_Venns.sh

#  Make supplemental table
Rscript Make_Sup_Table.r Median Total $DATDIR/Uniprot_9-21-18_wHMPAS.tab
