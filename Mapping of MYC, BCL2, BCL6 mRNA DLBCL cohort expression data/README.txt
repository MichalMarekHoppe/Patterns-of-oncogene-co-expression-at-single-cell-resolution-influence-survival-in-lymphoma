https://github.com/MichalMarekHoppe/Mapping-of-MYC-BCL2-BCL6-mRNA-DLBCL-cohort-expression-data-into-extent
author: Choi Hyung Won
maintained by: Michal Marek Hoppe (mmlhoppe@gmail.com)
source: https://www.medrxiv.org/content/10.1101/2020.10.20.20216101

Mapping of MYC, BCL2, BCL6 mRNA DLBCL cohort expression data into % extent based on empirical IHC protein quantification. 

IHC.csv - empirical ground truth measurement of MYC, BCL2, BCL6 protein percentage extent by IHC in five DLBCL cohorts (n=712).
mRNA_example.csv - MYC, BCL2, BCL6 mRNA expression values to be pasted by the user for a cohort of patients. Here values are from Reddy et al. 
required_functions.R - source for 'gauss.kernel.smooth', 'eCDF' and 'predict.eCDF' functions.
code_smooth_eCDF.R - source code to perform % extent mapping.
 