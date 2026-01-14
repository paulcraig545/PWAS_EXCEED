###############################################################################

# Set up libraries and options/files

###############################################################################
# install.packages("PRROC")
library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tidyr)
library(ggplot2)
library(tools)
library(lme4)
library(tibble)
library(pROC)
library(caret)
library(PRROC)             


parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort, ideally no spaces (for graphs and documentation)", action='store'),
  make_option('--id_column', type = 'character', default="ID", help = "Column names of identifier column", action = 'store'),
  make_option('--ps', type = 'character', help = 'File path to protein score file made using ProtS_calc.R'),
  make_option('--pheno', type = 'character', help = 'File path to MDD phenotype file'),
  make_option('--basic_covs', type = 'character', help = 'File path to covariate file for basic model'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
id_col <- opt$id_column # Vector of identifier columns 
pheno_fp=opt$pheno # MDD exposure (phenotype of cohort)
ProtS_fp=opt$ps # MDD ProtS (predictor) 
basic_covs_fp=opt$basic_covs # Just age and sex + cohort's technical covs
outdir <- opt$outdir # File path of output directory

args <- commandArgs(trailingOnly = FALSE) # get script name
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_name <- tools::file_path_sans_ext(basename(script_path))
outdir_fp <- file.path(outdir, script_name)

# sinking all output to a log file 

sink(paste0(cohort, script_name, "_ProtS_AD_assoc.log"))

###############################################################################

# Read in covariates, phenotype and ProtS files 

###############################################################################

ProtS <- readRDS(ProtS_fp) # File with ID and MDD_protS
print(colnames(ProtS))

# check that there is a ProtS column in the file 

if('MDD_protS' %in% colnames(ProtS) == FALSE){
  stop('No MDD_protS column in the ProtS file')
} else {
  print('MDD_protS column in the ProtS file')
}

# support phenotype .rds file

if (endsWith(pheno_fp, '.rds')){
  mdd_pheno <- readRDS(pheno_fp)
  mdd_pheno <- mdd_pheno |> rename(!!id_col:= ID)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .rds file')
}

# check that there is an antidep column in the file 

if('MDD' %in% colnames(mdd_pheno) == FALSE){
  stop('No MDD column in the phenotype file')
} else {
  print('MDD column in the phenotype file')
}
 
# remove missing values (if any?)
mdd_pheno <- mdd_pheno %>% filter(!is.na(MDD)) 

# logging phenotype characteristics 
print(paste0('Read in the MDD phenotype for ', cohort, script_name, ' : Number of cases: ',
             nrow(mdd_pheno %>% 
                    filter(MDD==1)), 
             'Number of controls: ',
             nrow(mdd_pheno%>% 
                    filter(MDD==0))))

all_covs <- readRDS(basic_covs_fp)

print(paste0("Covariates read in ", paste(colnames(all_covs %>% dplyr::select(-all_of(id_col))), collapse = ", ")))
print(paste0("Covariates data type:\n", paste(capture.output(str(all_covs)), collapse = "\n")))
 
#merge the phenotype and ProtS file together 

ProtS_pheno <- merge(ProtS, mdd_pheno, by = id_col)
ProtS_pheno_covs <- merge(ProtS_pheno, all_covs, by = id_col)


# logging phenotype characteristics after merging 

print(paste0('Read in the MDD phenotype for ', cohort, script_name, ' after merging with ProtS and pheno: Number of cases: ',
             nrow(ProtS_pheno_covs %>% 
                    filter(MDD==1)), 
             ' \n Number of controls: ',
             nrow(ProtS_pheno_covs %>% 
                    filter(MDD==0))))

###############################################################################
  
  # Generalised linear model (GLM)
  
###############################################################################

# Create the formula
# Fit a logistic model 
# logit link function and binomial family 
# Outcome - MDD phenotype 
# Predictor - Antidepressant ProtS (from ProtS_calc.R)

covs_ls <- colnames(all_covs)[colnames(all_covs) != id_col]
covs_formu = paste0(covs_ls, collapse = " + ")

formula <- as.formula(
    paste0("as.factor(MDD) ~ scale(MDD_protS) + ", paste(covs_formu, collapse = " + "))
)

assoc_mod <- glm(formula, 
                 family=binomial (link=logit), 
                 data = ProtS_pheno_covs)

# Extract the effect estimates, standard errors and p-value 
warnings()
print(assoc_mod)
print(summary(assoc_mod)$coefficients %>% as.data.frame())

assoc_coefs <- summary(assoc_mod)$coefficients %>% as.data.frame()
assoc_coefs <- rownames_to_column(assoc_coefs, var = "Coefficient")

# save the coefficients 
outfile <- file.path(outdir, paste0(cohort, "_", script_name, "_ProtS_AD_coefficients.rds"))
print(paste0('Saving the model coefficients to ', outfile))
saveRDS(assoc_coefs, outfile)

###############################################################################

# Calculating the AUC and ROC graph

###############################################################################

predicted_probs <- predict(assoc_mod, type = 'response')

# Take the true outcomes (ProtS_pheno_covs$MDD)
# and the predicted probabilities for the 1 ('case') class
# returns false positive and true positive rates for different classification thresholds

roc_curve <- roc(ProtS_pheno_covs[complete.cases(ProtS_pheno_covs),]$MDD, predicted_probs)
auc_value <- auc(roc_curve)

# save ROC curve object for plotting all cohorts together
outfile_roc <- file.path(outdir, paste0(cohort, "_", script_name, "_roc_curve.rds"))
print(paste0('Saving the ROC curve object for plotting all cohorts together to rds object: ', outfile_roc))
saveRDS(roc_curve, outfile_roc)

# ROC Graph 
outfile_roc_graph <- file.path(outdir, paste0(cohort, "_", script_name, "_assoc_ROC_curve.pdf"))
print(paste0('Saving the ROC curve for the cohort alone to ', outfile_roc_graph))
cairo_pdf(file = outfile_roc_graph, width = 8, height = 6)
plot.roc(roc_curve, col = "blue", lwd =2, main = paste0('ROC Curve: ', cohort, "_", script_name))
dev.off()

###############################################################################

# Precision Recall Curve 

###############################################################################

pr_curve <- pr.curve(ProtS_pheno_covs$MDD, predicted_probs, curve = T)

outfile_pr <- file.path(outdir, paste0(cohort, "_", script_name, "_assoc_precision_recall.pdf"))
print(paste0('Saving the PR curve for the cohort alone to ', outfile_pr))
cairo_pdf(file = outfile_pr, width = 8, height = 6)
plot(pr_curve, col = "red", main= paste0('Precision Recall Curve: ', cohort, "_", script_name))
dev.off()

sink()