# PWAS-EXCEED
External validation of protein-wide association study of major depressive disorder

Shared files required:
1) GS_lasso_weights_MDD.rds: proteins and their LASSO weights trained using GS cohort
2) protein_list.rds: List of proteins with non-zero coefficient (LASSO) used in calculating the protein scores

## Overall summary:
1) Perform rank transformation using preprocessing_proteins.R to standardise the proteins
2) Run protS_calc.R to apply the weights from our training model to your cohort and produce a protein score for each participant in your cohort
3) Run predictive models to use the protein scores created in the previous step to predict major depressive disorder status in your cohort
  - MDD ~ scale(MDD_protS) + age + sex + techinical covariates in ALL participants
*Technical covariates can be variables such as assessment centre, spectrometer (if applicable)

Please refer to the following information for further details and please let us know if you have any questions, thank you so much for your help!

## protein preprocessing: Rank standardisation
The protein scores were trained using standardised protein levels using rank standardisation. Therefore the protein scores should also be calculated using standardised protein levels. 
This script will carry out rank transformation on the protein data.

The R script preprocessing_proteins.R will read :
1) The dataframe (as .rds format) containing your cohort's protein levels, it should have rows as participant ID and columns as protein names.
2) The protein_list.rds file which contains a list of proteins provided by us
Then the R script will filter your cohort's proteins to the proteins from our training model and scale them

Arguments:
--cohort : Cohort name, e.g 'GS' or 'EXCEED' \
--proteins : The file path to the protein file (rds format) \
--probe : The file path for the list of proteins from LASSO (non-zero coefficients) provided by us \
--id_column : The column name of the identifier column (default == ID) \
--outdir : The directory where the outputs will be saved

Example:
```bash
Rscript preprocessing_proteins.R \
  --cohort "GS" \
  --proteins "/Users/Desktop/test_prot_data.rds" \
  --probe "/Users/Desktop/protein_list.rds" \
  --id_column "ID" \
  --outdir "/Users/Desktop/"
```

## Calculate protein scores
The R script protS_calc.R will read: 
1) The _mrs_preproc_proteins.rds file from preprocessing_proteins.R,
2) The LASSO weights provided by us
3) The major depressive disorder phenotype in your cohort (0 = controls, 1 = cases), this file should be .rds format with two columns, ID and MDD
Then, it will calculate protein score for each participant.

Arguments:
--cohort : Cohort name, e.g 'GS' or 'EXCEED' \
--std_prot : The file path for the standardised protein file from preprocessing_proteins.R \
--id_column : The column name of the identifier column (default == ID) \
--GS_weights : The file path for the protein weights file provided by us \
--pheno : The file path to the major depressive disorder phenotype file for your cohort (rds format). The file should have two columns: ID and MDD with MDD coded as 0 (controls) and 1 (cases) \
--outdir : The directory where the results and graphs will be saved 

Example:
```bash
Rscript protS_calc.R \
--cohort "GS" \
--std_prot "/Users/Desktop/GS_mrs_preproc_proteins.rds" \
--id_column "ID" \
--GS_weights "/Users/Desktop/lasso_coef.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--outdir "/Users/Desktop/"
```

## Predictive model 1 : basic model in all participants
The R script basic_model_all.R will read:
1) The _MDD_protS.rds file from protS_calc.R output
2) The major depressive disorder phenotype in your cohort (0 = controls, 1 = cases), this file should be .rds format with two columns, ID and MDD
3) A .rds file with ID and the **basic covariates** including sex, age, techinical covariates such as assessment_centre and/or spectrometer (if applicable) as columns.

Then the protein scores created in the previous step will be used to predict major depressive disorder status in your cohort. 

Key elements of the model: Phenotype (major depressive disorder phenotype), predictor (protein scores output from protS_calc.R), basic covariates (should have rows as participant ID and columns as covariates)

General Output: We would like the coefficients of the model, alongside the standard errors, t values, P values. We will assess the model's performance using AUC, including ROC and PR curves.

Arguments:
--cohort : Cohort name, e.g 'GS' or 'EXCEED' \
--id_column : The column name of the identifier column (default == ID) \
--ps : file path to the _MDD_protS.rds file from protS_calc.R output
--pheno : The file path to the major depressive disorder phenotype file for your cohort (rds format). The file should have two columns: ID and MDD with MDD coded as 0 (controls) and 1 (cases) \
--basic_covs : The file path to the basic covariate file for your cohort (rds format). It should include ID and each covariate as column \
--outdir : The directory where the results and graphs will be saved 


```bash
Rscript basic_model_all.R \
--cohort "GS" \
--id_column "ID" \
--ps "/Users/Desktop/GS_MDD_protS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--basic_covs "/Users/Desktop/basic_covs.rds" \
--outdir "/Users/Desktop/"
```