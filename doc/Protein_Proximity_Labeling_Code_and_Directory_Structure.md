#  Proximity Labeling Directory Structure

###  Written by Cory Haley White

##  Running Code

All code is set to run from the bin directory in the bitbucket repository.  You will need Rscript.exe installed and in your PATH statement to run this on windows.  R code will work in Linux or windows, but you may need to build a wrapper in either batch or shell depending on your environment instead of using the provided wrappers.  You can also run the code in the proper order manually, but this is not recommended as it will be set to debugging defaults unless you specify the variables directly.  For converting pdf files to tiff, you will need to install Imagemagick in windows to use the convert command.  

##  File Locations and Directory Structure

All data files but the downloaded Uniprot database are present in the data folder in a tar.gz file.  A description of the directory and structure is given below.  

bin - Code is located here.  R code serves as the primary workhorse scripts while .sh files serves to run R code with inputs representing the type of analysis performed.  All code contains a short description of the purpose.  

data - Data directory containing small data files for input and subdirectories named according to the experiment which contain output files.  Beneath each experiment directory is the target (e.g. CD45), and then the merging method (median) and normalization method (total) e.g. "CD45_1137_p22\CD45\Median_Total".  Within this directory is the output in the form of .csv or .txt files.  Output files are named based upon the experiment, the target, a name representing the type of analysis, and the filter settings e.g. ( "Name"_"Target"_"AnalysisType"_"Filters"  CD45_1137_p22_CD45Target_AllAnalyses.csv).  

data/Experimental_Design_Files - Directory within data directory that contains structure of experiment required for fold change measures

data/TempSend - Directory for compressed output for ease of sending to others.  

doc - RMarkdown directory which contains this document and other documents explaining the methods and/or code.  

figures - Contains subdirectories for each experiment where graphical output is stored.  Graphical output is named according to the experiment name, Target, conditions compared, and analysis (e.g. CD45_1137_p22_Target_vs_Iso_Lg2FC_1_VPlot.tiff).  

figures/Results_Comparison - Contains any comparisons between experiments present in the repository.  Two-way comparisons are stored in this folder (e.g. p35_2d1_vs_p35_F10_Venn.pdf) while three-way and four-way venns are stored in sub folders.  Venn diagrams are generated for the intersection at the gene level only for log2FC cutoffs indicated in the name of the file.  
