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

#  Cory Haley White

#  Shell script to run R functions with variable inputs
#  The folder containing Rscript.exe must be in your path for this file to run, otherwise give it the full path of Rscript.  

export NAME=CD45_1137_p22
export METRIC=Median
export NORM=Total
export GROUP2=Target
export GROUP1=Iso
export PAIRED=NotPaired
export FILTER=1
export ASSOC="CD45_Associated_Proteins.txt"
export EXTRA="CD45_Extra_Labels.txt"
export EXCLUDE="CD45_Exclude.txt"
export LG2=1.32

#  Correlation Targets
export TARGET=CD45
export TARGETNAME=P08575

bash Run_Analysis.sh

#  Adding in more log2FCs
export LG2=2
bash Run_Analysis.sh

export LG2=1
bash Run_Analysis.sh