#  Proximity Labeling Code File Description

###  Written by Cory Haley White

#  Introduction and Background

The goal of this project is to apply protein proximity labeling techniques where a protein is targeted by an antibody or small molecular ligand that contains an attached catalytic moiety which can generate reactive small molecules.  These reactive molecules may then covalently bind to proteins in close proximity to the original target protein.  The labeled proteins are then analyzed through LC-MS/MS based proteomics to produce a list of identified proteins and their TMT-labeled intensity measures.  This data is then subjected to correlation profiling and fold change analysis to search for proteins within the experiment with highly similar profiles under the assumption that as the reactive molecules spread away from the original target protein, they will label proteins in close proximity in a similar fashion as to how the original target protein was labeled whereas ones not in close proximity will have a different profile or lower fold changes due to aspects such as time of spread of the reactive molecules.  For a further description of the chemical/biological methodology behind this work, please see the paper entitled "Exploring Cellular Microenvironments via Photocatalytic Carbene Generation".  

#  Profiling Methodology

Under the assumption that a protein in close proximity to the target protein will have a similar protein intensity profile, we chose to attempt two separate analysis methods.  The first method generates a fold change for the labeling when using an antibody target vs. a generic "isotype" which does not target any particular protein.  The second method examines the target protein profile, and correlates the profile to that of all proteins surviving filtering.  Both methods are applied in the code, but for the correlation methodology, having more than two conditions is highly recommended.  

#  Data Filtering and Normalization

Peptide level is data obtained from a vendor in an excel file and is converted to a .csv file for input into R.  Data is filtered to remove all proteins with only a single (adjustable setting) peptide detected unless that peptide is the target.  Peptide data is normalized according to the summed total intensity value of each sample and then multiplied by the average summed total across samples to rescale values.  Peptide data is converted to protein data by taking the median (default) or geometric average of all peptides belonging to a Uniprot accession value.  Data is then subset to only proteins surviving a filter cutoff of 0 in half the samples.  

#  Correlation Methodology

Protein level data is correlated against the target protein abundance values using Pearson correlation.  A p-value is generated for this and corrected for multiple comparisons.  

#  Fold Change Methodology

Normalized protein data is converted to log2 values with replacement of 0 by the minimum non-zero value.  A design file (located in the Experimental_Design folder) is used to generate a group variable.  This group variable is used by Limma to assess log2FC values and a level of significance to determine if a protein is significantly different between groups.  Log2FC, p-values, and FDR adjusted p-values are generated and stored in the output file with name tail "AllAnalyses.csv".  All data but contaminants indicated as such by the vendor are present in this file including suspected contaminants not indicated by the vendor and antibodies.  

#  Volcano Plot Generation

Proteins are filtered to remove known contaminants such as antibody contaminants or anything listed as a contaminant by the vendor.  An exclusion list of proteins can be designated in the code wrappers and any proteins present in this list are removed prior to plot generation.  An associated proteins list can also be specified in the wrapper and those are given a label in the volcano plot.  An extra label text file can be given in the wrapper and proteins present in that file are labeled if they are above the log2FC cutoff limit.  An additional "Special_List" may be given to color and label proteins regardless of if they are above or below the cutoff.  Colors for each designation are set in the .r files themselves and may be adjusted there.  Finally, for some experiments, a volcano plot with antibodies present was generated.  These files contains "wAB" in the tail of the name and antibodies are colored in maroon.  

#  Membrane Protein Determination

Membrane information generated from the Uniprot database (downloaded in Sept. 2018) is used to assess membrane protein count in the total proteins surviving filtering.  Membrane proteins are determined by having an Extracellular domain in the topological domain category of Uniprot.  As databases are continually updated, some proteins were not properly identified as membrane, e.g. PTPRCAP.  This protein was manually corrected as a membrane protein in the downloaded database.  

#  Venn Diagram Generation

When comparing two experiments or two different conditions vs. isotype in one experiment in a Venn diagram, proteins are filtered to remove known contaminants such as antibody contaminants.  This is done by providing an exclusion list to the wrapper.  Occasionally this wrapper was used to remove non-membrane proteins or suspected contaminants in visualizations.  

#  Data Availability

Data for this project is stored in a tar.gz file in the data directory.  The membrane database was downloaded from Uniprot and is not included in this repository.  To take advantage of the membrane calling code, you will need to download such a database from Uniprot and include the subcellular location and topological domain as columns in the database.  
