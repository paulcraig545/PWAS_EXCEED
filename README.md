# MetWAS-EXCEED
External validation of protein-wide association study of major depressive disorder

Shared files required:
1) GS_lasso_weights_MDD.rds: proteins and their LASSO weights trained using GS cohort
2) protein_list.rds: List of proteins with non-zero coefficient (LASSO) used in calculating the protein scores

## Overall summary:
1) Perform z-score standardisation using preprocessing_proteins.R to scale the proteins
2) Run protS_calc.R to apply the weights from our training model to your cohort and produce a protein score for each participant in your cohort
3) Run four predictive models to use the protein scores created in the previous step to predict major depressive disorder status in your cohort
  - MDD ~ scale(MDD_protS) + age + sex + techinical covariates in ALL participants
  - MDD ~ scale(MDD_protS) + age + sex + techinical covariates in MDD cases only (please add mdd as a column to the covariate file)
  - MDD ~ scale(MDD_protS) + age + sex + techinical covariates + lifestyle covariates in ALL participants
  - MDD ~ scale(MDD_protS) + age + sex + techinical covariates + lifestyle covariates in MDD cases only \
*Technical covariates can be variables such as assessment centre, spectrometer (if applicable) \
*Lifestyle covariates can be variables such as BMI, smoking status, socioeconomic status, educational level, MDD diagnosis, ethnicity, depressive symptom scores, alcohol drinking (only when applicable) 
4) Run cohort_demograph_info.R to generate a demographic table about your cohort 

Please refer to the following information for further details and please let us know if you have any questions, thank you so much for your help!

## protein preprocessing: Z-score standardisation
The protein scores were trained using standardised protein levels using rank standardisation. Therefore we would like the protein scores to be calculated also using standardised protein levels.

The R script preprocessing_proteins.R will read :
1) The dataframe (as .rds format) containing your cohort's protein levels, it should have rows as participant ID and columns as protein names.
2) The protein_list.rds file which contains a list of probe proteins provided by us
Then the R script will filter your cohort's proteins to the probe proteins from our training model and scale them

Arguments:
--cohort : Cohort name, e.g 'GS' or 'EXCEED' \
--proteins : The file path to the protein file (rds format) \
--probe : The file path for the list of probe proteins from LASSO (non-zero coefficients) provided by us \
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

## Predictive model 2 : basic model in MDD cases only
The basic_model_mdd.R script is similar to model 1, but with additional mdd column (MDD diagnosis) and it will filter out participants with no MDD diagnosis
1) The _MDD_protS.rds file from protS_calc.R output
2) The major depressive disorder phenotype in your cohort (0 = controls, 1 = cases), this file should be .rds format with two columns, ID and MDD
3) A .rds file with ID and the **basic covariates** including sex, age, mdd, techinical covariates such as assessment_centre and/or spectrometer (if applicable) as columns. The mdd column should be coded into 0/1 (1 = mdd, 0 = controls), the basic_model_mdd.R script will filter out the controls.
   
```bash
Rscript basic_model_mdd.R \
--cohort "GS" \
--id_column "ID" \
--ps "/Users/Desktop/GS_MDD_protS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--basic_covs "/Users/Desktop/basic_covs.rds" \
--outdir "/Users/Desktop/"
```

## Predictive model 3 : complex model in all participants
This is also similar to model 1 but with additional lifestyle-related covariates such as BMI, smoking status, socioeconomic status, educational level, MDD diagnosis, ethnicity, depressive symptom scores, alcohol drinking (only when applicable)

The R script complex_model_all.R will read:
1) The _MDD_protS.rds file from protS_calc.R output
2) The major depressive disorder phenotype in your cohort (0 = controls, 1 = cases), this file should be .rds format with two columns, ID and MDD
3) A .rds file with ID and the **complex covariates** including sex, age, techinical covariates such as assessment_centre and/or spectrometer, and lifestyle covariates (if applicable) as columns.

```bash
Rscript complex_model_all.R \
--cohort "GS" \
--id_column "ID" \
--ps "/Users/Desktop/GS_MDD_protS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--complex_covs "/Users/Desktop/complex_covs.rds" \
--outdir "/Users/Desktop/"
```

## Predictive model 4 : complex model in MDD cases only
This is also similar to model 3 (with lifestyle-related covariates) but in MDD cases only

The R script complex_model_mdd.R will read:
1) The _MDD_protS.rds file from protS_calc.R output
2) The major depressive disorder phenotype in your cohort (0 = controls, 1 = cases), this file should be .rds format with two columns, ID and MDD
3) A .rds file with ID and the **complex covariates** including sex, age, techinical covariates such as assessment_centre and/or spectrometer, and lifestyle covariates (if applicable) as columns. The mdd column should be coded into 0/1 (1 = mdd, 0 = controls), the complex_model_mdd.R script will filter out the controls.

```bash
Rscript complex_model_mdd.R \
--cohort "GS" \
--id_column "ID" \
--ps "/Users/Desktop/GS_MDD_protS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--complex_covs "/Users/Desktop/complex_covs.rds" \
--outdir "/Users/Desktop/"
```

## Demographic information about your cohort
The file cohort_demograph_info.R will generate a table of demographic information of your cohort (useful for the manuscript and interpretating our results). It formats information on age, sex, bmi, MDD diagnosis, and AD MetS (generated from protS_calc.R).

--cohort: Cohort name, e.g 'GS' or 'EXCEED' \
--id_column: The column name of the identifier column (default == ID) \
--ps: The filepath to the _MDD_protS.rds file from protS_calc.R output (colnames: ID, MDD_protS) \
--pheno: The filepath to the AD phenotype file (colnames: ID, MDD(0 = no exposure/1 = major depressive disorder)) \
--demo: Filepath to the file containing demographic variables: age(numeric), sex(factor: Male/Female), bmi(numeric), mdd(factor: 0 = controls/1 = cases), MDD_protS)
--outdir : The directory where the results will be saved 


```bash
Rscript cohort_demograph_info.R \
--cohort "GS" \
--id_column "ID" \
--ps "/Users/Desktop/GS_MDD_protS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--demo "/Users/Desktop/covs.rds" \
--outdir "/Users/Desktop/"
```