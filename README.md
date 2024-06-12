## Data Analysis and Optimization for CCLE and TCGA databases platform.

It facilitates data analysis for the well known databases CCLE and TCGA with two main avenues of modelling schemes implemented:

* GLM with ReNoiR (Generalized Linear Regression models)
* Optimization with CARNIVAL (Mixed-Integer Linear Programming optimization for gene regulatory networks)

The user can start the platform by running the "app.R" function in R (version 4.0 and above recommended) using a command like "shiny::runApp()".

In order to read a new version of CCLE the user needs to create a folder (such as the ones existing in \inputs) with the version of the CCLE (e.g 20Q1) and place there all the necessary files for mutation and expression profiles (the exact ones as in the example folders). Then , when launching the platform the user can select the option "read new version" from the bottom and input the name of the version in the box, which will be used to wrap all matrices and data in the .R  object later on. 

For the main basic analysis the platform will use only CCLE. TCGA data sets can be used to perform GLM and CARNIVAL. Each run of the platform creates a folder named after the main database used (CCLE, TCGA) and the type of cancer selected by the user. All the generated files from the platform are stored therein, with the "opt" folder including all the optimized networks per condition in the .DOT format, and the "GLM" folder containing the .html repor for the analysis of the regression results.

In the current version, the emphasis was given in analysing gene TP53. If the user selectes a different gene, some extra input files might be required, such as mutation profiles and CNV from cBioPortal specifically for the TCGA pancancer altas. 



### References
Charalampos P. Triantafyllidis, Alessandro Barberis, Fiona Hartley, Ana Miar Cuervo, Enio Gjerga, Philip Charlton, Linda van Bijsterveldt, Julio Saez Rodriguez, Francesca M. Buffa,
A machine learning and directed network optimization approach to uncover TP53 regulatory patterns,
iScience,
Volume 26, Issue 12,
2023,
108291,
ISSN 2589-0042,
https://doi.org/10.1016/j.isci.2023.108291.
(https://www.sciencedirect.com/science/article/pii/S2589004223023684)

Abstract: Summary
TP53, the Guardian of the Genome, is the most frequently mutated gene in human cancers and the functional characterization of its regulation is fundamental. To address this we employ two strategies: machine learning to predict the mutation status of TP53 from transcriptomic data, and directed regulatory networks to reconstruct the effect of mutations on the transcipt levels of TP53 targets. Using data from established databases (Cancer Cell Line Encyclopedia, The Cancer Genome Atlas), machine learning could predict the mutation status, but not resolve different mutations. On the contrary, directed network optimization allowed to infer the TP53 regulatory profile across: (1) mutations, (2) irradiation in lung cancer, and (3) hypoxia in breast cancer, and we could observe differential regulatory profiles dictated by (1) mutation type, (2) deleterious consequences of the mutation, (3) known hotspots, (4) protein changes, (5) stress condition (irradiation/hypoxia). This is an important first step toward using regulatory networks for the characterization of the functional consequences of mutations, and could be extended to other perturbations, with implications for drug design and precision medicine.
Keywords: Regulatory networks; directed networks; causal inference; mutations; cancer systems biology; machine learning; TP53; trascriptomics; regulon
