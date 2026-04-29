# PWAS-EXCEED
External validation of protein-wide association study of antidepressant exposure

Shared files required:
1) GS_lasso_weights_MDD_newProts.rds: Proteins and their LASSO weights trained on MDD status in GS cohort
2) GS_lasso_weights_likert_D_newProts.rds: Proteins and their LASSO weights trained on Depressive symptoms subscale of GHQ in GS cohort
3) GS_lasso_weights_likert_total_newProts.rds: Proteins and their LASSO weights trained on total GHQ in GS cohort
4) protein_list_newProts.rds: List of proteins with non-zero coefficient (LASSO) used in calculating the proteins scores

## Overall summary:
PWAS_pipeline.r carries out validation of protein scores in one script by:
1) Performing z-score standardisation to scale the Proteins
2) Applying the weights from our training model to your cohort and produce a protein score for each participant in your cohort
3) Run predictive models to use the protein scores created in the previous step to predict MDD and/or depressive symptoms in your cohort

As well as the weights files for made uisng binary and continuous input variables provided by us, you may supply the script with one binary outcome variable and/or one continuous outcome variable using the appropriate flags. You may also supply the script with a file containing appropriate covariates. If weights and outcome phenotypes are provided for both continuous and binary phenotypers the script will perform cross validation for each weights file and outcome phenotype. e.g.
  - binary_phenotype ~ scale(binary_weights) + covariates
  - binary_phenotype ~ scale(continuous_weights) + covariates
  - continuous_phenotype ~ scale(binary_weights) + covariates
  - continuous_phenotype ~ scale(continuous_weights) + covariates


Please refer to the following information for further details and please let us know if you have any questions, thank you so much for your help!

## pwas_pipeline.r
The R script pwas_pipeline.r will read :
1) The dataframe (as .rds format) containing your cohort's protein levels, it should have rows as participant ID and columns as protein names. 
2) The protein_list_newProts.rds file which contains a list of probe Proteins provided by us. The R script will filter your cohort's Proteins to the Proteins from our training model and scale them

Arguments:
--cohort : Cohort name, e.g 'GS' or 'EXCEED' \
--proteins : The file path to the protein file (rds format) \
--probe : The file path for the list of proteins from LASSO(non-zero coefficients) provided by us (protein_list_newProts.rds)  \
--id_column : The column name of the identifier column (default == ID) \
--binary_weights : \
--binary_weights_name : \
--continuous_weights : \
--continuous_weights_name : \
--binary_pheno : \
--binary_pheno_name : \
--continuous_pheno : \
--continuous_pheno_name : \ 
--covs : \
--outdir : The directory where the outputs will be saved

Example:
```bash
Rscript PWAS_pipeline.r 
    --cohort "GS" 
    --proteins "test_prot_data_newProts.rds" 
    --probe "protein_list_newProts.rds" 
    --id_column "id" 
    --binary_weights "GS_lasso_weights_MDD_newProts.rds" 
    --binary_weights_name "MDD" 
    --continuous_weights "GS_lasso_weights_likert_D_newProts.rds" 
    --continuous_weights_name "GHQ_D" 
    --binary_pheno "mdd_pheno.rds" 
    --binary_pheno_name "MDD" 
    --continuous_pheno "likert_D_pheno.rds" 
    --continuous_pheno_name "GHQ_D" 
    --covs "basic_covs.rds" 
    --outdir "PWAS_pipeline_out/"
```