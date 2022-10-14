## Data Analysis and Optimization for CCLE and TCGA databases platform.

It facilitates data analysis for the well known databases CCLE and TCGA with two main avenues of modelling schemes implemented:

* GLM with ReNoiR (Generalized Linear Regression models)
* Optimization with CARNIVAL (Mixed-Integer Linear Programming optimization for gene regulatory networks)

The user can start the platform by running the "app.R" function in R (version 4.0 and above recommended) using a command like "shiny::runApp()".

In order to read a new version of CCLE the user needs to create a folder (such as the ones existing in \inputs) with the version of the CCLE (e.g 20Q1) and place there all the necessary files for mutation and expression profiles (the exact ones as in the example folders). Then , when launching the platform the user can select the option "read new version" from the bottom and input the name of the version in the box, which will be used to wrap all matrices and data in the .R  object later on. 

For the main basic analysis the platform will use only CCLE. TCGA data sets can be used to perform GLM and CARNIVAL. Each run of the platform creates a folder named after the main database used (CCLE, TCGA) and the type of cancer selected by the user. All the generated files from the platform are stored therein, with the "opt" folder including all the optimized networks per condition in the .DOT format, and the "GLM" folder containing the .html repor for the analysis of the regression results.

In the current version, the emphasis was given in analysing gene TP53. If the user selectes a different gene, some extra input files might be required, such as mutation profiles and CNV from cBioPortal specifically for the TCGA pancancer altas. 



### References
Reconstructing the functional effect of TP53 somatic mutations on its regulon using causal signalling network modelling
Charalampos P. Triantafyllidis, Alessandro Barberis, Ana Miar Cuervo, Philip Charlton, Fiona Hartley, Linda Van Bijsterveldt, Enio Gjerga, Julio Saez Rodriguez, Francesca M. Buffa
bioRxiv 2022.06.23.497293; doi: https://doi.org/10.1101/2022.06.23.497293
