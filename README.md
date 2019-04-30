# linc_PCPG
codes and metadata related to lincRNA profiles in PCPG

########################################################################################################################################

R code files:

comparemodels_PCPG_molsubtypes_classification_by_lincRNAs.r
description: R script for comparing four machine learning models (CART, Elastic net, Ridge and LASSO) for classification of five molecular subtypes of PCPG: SDHx-mutated pseudohypoxia, non SDHx-mutated pseudohypoxia, Wnt-altered, kinase signaling and cortical admixture

Comparemodels_rf_nn_svm_xgb_for_PCPG_metastatic_linc.r
description: R script for comparing four machine learning models (random forest, neural network, support vector machine and xtreme gradient boosting) for classification of five metastatic and non-metastatic PCPG

CoxRegression_Prognostic_Index_CV.r
description: R script for building the prognostic index combining the expression of 18 lincRNAs, SDHB germline mutation, ATRX somatic mutation, TERT gene expression, tumor location, and hormone secretary profiles

########################################################################################################################################

Metadata files:

TCGA_PCPG_clinical_molSubtype.txt
clinical data containing driver mutation status and molecular subtype information for TCGA PCPG cohort.

gencode22_lincRNA_GeneID.txt
lincRNA gene symbols and gene Ensembl ids from Gencode22 included in the analysis

PCPG_GDC_CountData_normal_nr.txt
RNA-Seq v2 count data for 7649 lincRNAs for TCGA PCPG cohort and 3 normal adrenal medulla samples

SelectedlincRNAs_by_ElasticNet_for_PCPG_MolSubtypes.txt
list of 106 lincRNAs chosen for classification of PCPG molecular subtypes by elastic net model

Clinical_Aggressive_lincRNA.txt
18 lincRNAs having elevated expression signature in aggressive/metastatic PCPG compared to non aggressive/metastatic PCPG from TCGA cohort

SDHX-mut_Pseudohypoxia_PCPG_lincRNA_Survival.txt
lincRNAs whose over-expression/under-expression has significant difference in metastasis-free survival in PCPG patients with germline mutation in SDHB/SDHD (Kaplan-Meier analysis)

Pseudohypoxia_PCPG_lincRNA_Survival.txt
lincRNAs whose over-expression/under-expression has significant difference in metastasis-free survival in PCPG patients in non-SDHx-mutated pseudohypoxia subtype (Kaplan-Meier analysis)

Kinase_signaling_PCPG_lincRNA_Survival.txt
lincRNAs whose over-expression/under-expression has significant difference in metastasis-free survival in PCPG patients in kinase-signaling subtype (Kaplan-Meier analysis)
