# New Version Available:
hello a new version of this code base has been released as a new project available here: https://github.com/Merck/MicroMap_Pipeline, this was done due to the fact the code base has vast changes from this repo. This repo will be archived and new work will be done on the new repo.

#  Introduction

###  Written by Cory Haley White

This code is provided to replicate the analysis performed for the publication "Exploring Cellular Microenvironments via Photocatalytic Carbene Generation".  A link will replace this line upon acceptance of the publication.  In order to run the code as is, you must have R (tested with version 3.5.2), and all dependences installed.  


#  R Libraries Used

```xlsx```
```ggrepel```
```Venndiagram```
```ggplot2```
```extrafont```
```limma```
```stringr```
```reshape```


#  Creation of Membrane Protein Database

A membrane protein database was used for this work and is available upon request.  This database was downloaded from Uniprot and curated for membrane proteins using the following command.  

```
#  This script assumes the database is stored in the data directory of the GitHub repository.  
Rscript Get_Membrane.r Uniprot_9-21-18.tab Uniprot_9-21-18.tab
```


#  Running Code

The code may be run with the following command in a Linux environment.  The creation of the supplemental table requires a membrane protein database to be downloaded to the data directory.  

```
bash Run_All.sh
```

If this code crashes upon reaching the use of the supplemental table with membrane information, the most likely cause is that you have not provided a proper membrane database file.  Please see the "Get_Membrane.r" script and the "Protein_Proximity_Labeling_Code_File_Description.md" document in the doc directory 


##  NOTES

Data files are provided for each experiment in a tar.gz file which is unzipped in the wrapper script Run_All.sh.  Please see the files "Protein_Proximity_Labeling_Code_and_Directory_Structure.md" and "Protein_Proximity_Labeling_Code_File_Description.md" for a detailed description of the directory structure and code used, and the file "Protein_Proximity_Labeling_Background_Methods.md" for a detailed description and rational of the methods.  These markdown files are located in the doc directory.  
