###### Set up libraries and options/files

library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tools)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(pROC)
library(caret)
library(PRROC)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort", action='store'),
  make_option('--proteins', type='character', help="The filepath for proteins file", action='store'), 
  make_option('--probe', type = 'character', help= "The filepath for the list of proteins", action = 'store'),
  make_option('--id_column', type = 'character', default="ID", help = "Column name for identifier column", action = 'store'),
  make_option('--binary_weights', type = 'character', default = "NULL", help = "Weights file from binary phenotype for calculating scores"),
  make_option('--binary_weights_name', type = 'character', default = "MDD", help = "Name of phenotype used for making binary weights"),
  make_option('--continuous_weights', type = 'character', default = "NULL", help = "Weights file from continuous phenotype for calculating scores"),
  make_option('--continuous_weights_name', type = 'character', default = "GHQ", help = "Name of phenotype used for making continuous weights"),
  make_option('--binary_pheno', type = 'character', default = "NULL", help = 'File path to binary phenotype file'),
  make_option('--binary_pheno_name', type = 'character', default = "MDD", help = 'Name of binary phenotype'),
  make_option('--continuous_pheno', type = 'character', default = "NULL", help = 'File path to continuous phenotype file'),
  make_option('--continuous_pheno_name', type = 'character', default = "GHQ", help = 'Name of continuous phenotype'),
  make_option('--covs', type = 'character', default = "NULL", help = 'File path to covariates file'),  
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store'),
  make_option('--missingness_cov', type = 'character', default = FALSE, help = 'If TRUE use N missing proteins as a covariate. Default False')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

# setting up arguments from the options 
print('Setting up the options')
cohort <- opt$cohort
prot_filepath=opt$proteins # protein file
probe_filepath=opt$probe # List of protein names
id_col <- opt$id_column # Vector of identifier column
covs_fp <- opt$covs
out_dir <- opt$outdir
missingness_cov <- opt$missingness_cov

dir.create(out_dir, showWarnings = FALSE)
sink(file.path(out_dir, "protein_score.log"))


print(paste0('protein file from : ', prot_filepath))
print(paste0('List of proteins from : ', probe_filepath))
print(paste0('ID column : ', id_col))
print(paste0('Covariates from : ', opt$covs))
print(paste0('Output to be saved in : ', out_dir))


proteins <- readRDS(prot_filepath)
prot_list <- readRDS(probe_filepath)

ms_prot <- proteins %>%
    rename(id = all_of(id_col)) |>
    select(id, intersect(names(proteins), prot_list))

missingness <- data.frame(id = ms_prot$id, missing_prots = rowSums(is.na(ms_prot)))
saveRDS(missingness, file.path(out_dir, "missing_proteins.rds"))

rm(proteins) # remove large protein object 

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

prot_std <- mapply(
              transform,
              ms_prot[, prot_list]
            ) |>
    as.data.frame() |>
    mutate(id = ms_prot$id) |>
    select(c(id, everything())) |>
    filter(!is.na(id))

long_prot_std <- prot_std %>% 
  pivot_longer(-"id", 
               names_to = "proteins", 
               values_to= "Protval")

# Covariates formula
all_covs <- readRDS(covs_fp)

if(missingness_cov==TRUE){
  all_covs <- all_covs |> rename(id = all_of(id_col)) |> merge(missingness, by = "id")
}
covs_ls <- colnames(all_covs)[colnames(all_covs) != id_col]
covs_formu = paste0(covs_ls, collapse = " + ")

print(opt$binary_weights)
if(opt$binary_weights!="NULL"){

  print("This shouldn't print")
  
  ##### MDD weights #####
  binary_weights_filepath <- opt$binary_weights
  binary_weights <- readRDS(binary_weights_filepath)
  binary_weights_name <- opt$binary_weights_name
  dir.create(file.path(out_dir, binary_weights_name))
  colnames(binary_weights) = c('proteins','weights')

  long_prot_std_binary_weights <- merge(long_prot_std, binary_weights, by = 'proteins')

  ProtS_binary <- long_prot_std_binary_weights %>% group_by(id) %>%
    summarise(binary_protS = sum(weights*Protval, na.rm = T)) %>% as.data.frame()

  binary_protS_outfile <- file.path(out_dir, binary_weights_name, paste0(cohort, "_binary_PS.rds"))
  print(paste0('Saving the protein score to ', binary_protS_outfile))
  saveRDS(ProtS_binary, binary_protS_outfile)

  PS_dist_binary <- ggplot(ProtS_binary, aes(x = binary_protS)) + 
    geom_histogram() + 
    theme_minimal() + 
    labs(x = 'Protein Profile Score', y = 'Count')+
    ggtitle(cohort)

  PS_density_binary <- ggplot(ProtS_binary, aes(x=binary_protS)) +
    geom_density() +
    theme_minimal() + 
    labs(x = 'Protein Profile Score', y = 'Count')+
    ggtitle(cohort)

  ggsave(file.path(out_dir, binary_weights_name, paste0(cohort, "_", binary_weights_name, "_protein_score_overall_distribution.png")), PS_dist_binary, width = 8, height = 6, device='png', dpi=300)
  ggsave(file.path(out_dir, binary_weights_name, paste0(cohort, "_", binary_weights_name, "_protein_score_overall_density.png")), PS_density_binary, width = 8, height = 6, device='png', dpi=300)
}

if(opt$continuous_weights!="NULL"){
  ##### GHQ weights #####
  continuous_weights_filepath <- opt$continuous_weights
  continuous_weights <- readRDS(continuous_weights_filepath)
  continuous_weights_name <- opt$continuous_weights_name
  dir.create(file.path(out_dir, continuous_weights_name))
  colnames(continuous_weights) = c('proteins','weights')

  long_prot_std_continuous_weights <- merge(long_prot_std, continuous_weights, by = 'proteins')

  ProtS_continuous <- long_prot_std_continuous_weights %>% group_by(id) %>%
    summarise(continuous_protS = sum(weights*Protval, na.rm = T)) %>% as.data.frame()

  continuous_protS_outfile <- file.path(out_dir, continuous_weights_name, paste0(cohort, "_continuous_PS.rds"))
  print(paste0('Saving the protein score to ', continuous_protS_outfile))
  saveRDS(ProtS_continuous, continuous_protS_outfile)

  PS_dist_continuous <- ggplot(ProtS_continuous, aes(x = continuous_protS)) + 
    geom_histogram() + 
    theme_minimal() + 
    labs(x = 'Protein Profile Score', y = 'Count')+
    ggtitle(cohort)

  PS_density_continuous <- ggplot(ProtS_continuous, aes(x=continuous_protS)) +
    geom_density() +
    theme_minimal() + 
    labs(x = 'Protein Profile Score', y = 'Count')+
    ggtitle(cohort)

  ggsave(file.path(out_dir, continuous_weights_name, paste0(cohort, "_", continuous_weights_name, "_protein_score_overall_distribution.png")), PS_dist_continuous, width = 8, height = 6, device='png', dpi=300)
  ggsave(file.path(out_dir, continuous_weights_name, paste0(cohort, "_", continuous_weights_name, "_protein_score_overall_density.png")), PS_density_continuous, width = 8, height = 6, device='png', dpi=300)
}


if(file.exists(file.path(out_dir, "descriptive_statistics.csv"))){
  desc_stats <- read.csv(file.path(out_dir, "descriptive_statistics.csv"))[,-1]
}else{
desc_stats <- data.frame(
  weights = character(), 
  Phenotype = character(), 
  beta = numeric(), 
  std = numeric(), 
  p = numeric(), 
  R2 = numeric(), 
  Incremental_R2 = numeric(), 
  AUC = numeric(), 
  Incremental_AUC = numeric()
  )
}

##### Scores #####
if(opt$binary_pheno!="NULL"){

  ##### Read in binary pheno ###############################################################################################################

  binary_pheno_filepath <- opt$binary_pheno
  binary_pheno <- readRDS(binary_pheno_filepath)
  binary_pheno_name <- opt$binary_pheno_name
  
  colnames(binary_pheno) <- c("id", "pheno")
  binary_pheno <- binary_pheno %>% filter(!is.na(pheno)) # remove missing values if there are any 

  binary_pheno_covs <- merge(binary_pheno, all_covs, by = id_col)

  binary_covs_formula <- as.formula(
      paste0("as.factor(pheno) ~ ", paste(covs_formu, collapse = " + "))
      )

  binary_pheno_assoc_mod <- glm(
      binary_covs_formula, 
      family=binomial (link=logit), 
      data = binary_pheno_covs
      )

  aucR <- auc(binary_pheno_covs[complete.cases(binary_pheno_covs),]$pheno, binary_pheno_assoc_mod$linear.predictors)  ### AUC for reduced module

  if(opt$binary_weights!="NULL"){

    dir.create(file.path(out_dir, binary_weights_name, binary_pheno_name))
    
    ##### Binary Protein Score #####
    
    binary_pheno_binary_PS <- merge(binary_pheno, ProtS_binary, by = "id")
    
    binary_pheno_binary_PS_dists <- ggplot(binary_pheno_binary_PS, aes(x = binary_protS, fill = as.factor(pheno))) + 
      geom_histogram(alpha = 0.8) + 
      theme_minimal() + 
      labs(x = 'protein Profile Score', y = 'Count', fill = 'pheno') +
      ggtitle(cohort)
    
    binary_pheno_binary_PS_densities <- ggplot(binary_pheno_binary_PS, aes(x=binary_protS, group=as.factor(pheno), color=as.factor(pheno))) +
        geom_density() +
        labs(title="Protein score distribution per pheno",x="Protein Score", y = "Density", color="pheno") +
        theme_classic()
    
    ggsave(file.path(out_dir, binary_weights_name, binary_pheno_name, paste0(cohort, "_", binary_weights_name, "_protein_score_case_control_distribution.png")), binary_pheno_binary_PS_dists, width = 8, height = 6, device='png', dpi=300)
    ggsave(file.path(out_dir, binary_weights_name, binary_pheno_name, paste0(cohort, "_", binary_weights_name, "_protein_score_case_control_density.png")), binary_pheno_binary_PS_densities, width = 8, height = 6, device='png', dpi=300)

    binary_pheno_binary_PS_covs <- merge(binary_pheno_binary_PS, all_covs, by = id_col)
    
    binary_covs_PS_formula <- as.formula(
        paste0("as.factor(pheno) ~ ", paste(covs_formu, collapse = " + "), " + scale(binary_protS)")
        )    
    
    binary_pheno_binary_PS_assoc_mod <- glm(
        binary_covs_PS_formula, 
        family=binomial (link=logit), 
        data = binary_pheno_binary_PS_covs
        )
    
    binary_pheno_binary_PS_assoc_coefs <- summary(binary_pheno_binary_PS_assoc_mod)$coefficients %>% as.data.frame() |> rownames_to_column(var = "Coefficient")
    
    # save the coefficients 
    binary_pheno_binary_PS_assoc_coefs_outfile <- file.path(out_dir, binary_weights_name, binary_pheno_name, paste0(cohort, "_", binary_pheno_name, "_", binary_weights_name, "_PS_assoc_coefs.rds"))
    saveRDS(binary_pheno_binary_PS_assoc_coefs, binary_pheno_binary_PS_assoc_coefs_outfile)
    
    aucF <- auc(binary_pheno_binary_PS_covs[complete.cases(binary_pheno_binary_PS_covs),]$pheno, binary_pheno_binary_PS_assoc_mod$linear.predictors)  ### AUC for full module

    print(paste0("AUC value for binary phenotype binary weights: ", aucF))
    print(paste0("Incremental AUC value for binary phenotype binary weights: ", aucF - aucR))

    binary_pheno_binary_PS_predicted_probs <- predict(binary_pheno_binary_PS_assoc_mod, type = 'response')
    roc_curve <- roc(binary_pheno_binary_PS_covs[complete.cases(binary_pheno_binary_PS_covs),]$pheno, binary_pheno_binary_PS_predicted_probs)    
    
    # save ROC curve object for plotting all cohorts together
    saveRDS(roc_curve, file.path(out_dir, binary_weights_name, binary_pheno_name, paste0(cohort, "_", binary_pheno_name, "_", binary_weights_name, "_roc_curve.rds")))
    
    # ROC Graph 
    cairo_pdf(file = file.path(out_dir, binary_weights_name, binary_pheno_name, paste0(cohort, "_", binary_pheno_name, "_", binary_weights_name, "_binary_pheno_binary_PS_assoc_ROC_curve.pdf")), width = 8, height = 6)
    plot.roc(roc_curve, col = "blue", lwd =2, main = paste0('ROC Curve: ', cohort, "_", binary_pheno_name, "_", binary_weights_name))
    dev.off()
    
    pr_curve <- pr.curve(binary_pheno_binary_PS_covs$pheno, binary_pheno_binary_PS_predicted_probs, curve = T)
    
    cairo_pdf(file = file.path(out_dir, binary_weights_name, binary_pheno_name, paste0(cohort, "_", binary_pheno_name, "_", binary_weights_name, "_PS_assoc_precision_recall.pdf")), width = 8, height = 6)
    plot(pr_curve, col = "red", main= paste0('Precision Recall Curve: ', cohort, "_", binary_pheno_name, "_", binary_weights_name))
    dev.off()

    binary_pheno_binary_weights_desc_stats <- data.frame(
      weights = binary_weights_name,
      Phenotype = binary_pheno_name,
      beta = summary(binary_pheno_binary_PS_assoc_mod)$coefficients["scale(binary_protS)", "Estimate"],
      std = summary(binary_pheno_binary_PS_assoc_mod)$coefficients["scale(binary_protS)", "Std. Error"],
      p = summary(binary_pheno_binary_PS_assoc_mod)$coefficients["scale(binary_protS)", "Pr(>|z|)"],
      R2 = NA,
      Incremental_R2 = NA,
      AUC = aucF,
      Incremental_AUC = aucF - aucR
      )

    desc_stats <- desc_stats |> rbind(binary_pheno_binary_weights_desc_stats)
  }
  
  if(opt$continuous_weights!="NULL"){

    dir.create(file.path(out_dir, continuous_weights_name, binary_pheno_name))

    ##### Binary Protein Score #####
    
    binary_pheno_continuous_PS <- merge(binary_pheno, ProtS_continuous, by = "id")
    
    binary_pheno_continuous_PS_dists <- ggplot(binary_pheno_continuous_PS, aes(x = continuous_protS, fill = as.factor(pheno))) + 
      geom_histogram(alpha = 0.8) + 
      theme_minimal() + 
      labs(x = 'protein Profile Score', y = 'Count', fill = 'pheno') +
      ggtitle(cohort)
    
    binary_pheno_continuous_PS_densities <- ggplot(binary_pheno_continuous_PS, aes(x=continuous_protS, group=as.factor(pheno), color=as.factor(pheno))) +
        geom_density() +
        labs(title="Protein score distribution per pheno",x="Protein Score", y = "Density", color="pheno") +
        theme_classic()
    
    ggsave(file.path(out_dir, continuous_weights_name, binary_pheno_name, paste0(cohort, "_", continuous_weights_name, "_protein_score_case_control_distribution.png")), binary_pheno_continuous_PS_dists, width = 8, height = 6, device='png', dpi=300)
    ggsave(file.path(out_dir, continuous_weights_name, binary_pheno_name, paste0(cohort, "_", continuous_weights_name, "_protein_score_case_control_density.png")), binary_pheno_continuous_PS_densities, width = 8, height = 6, device='png', dpi=300)

    binary_pheno_continuous_PS_covs <- merge(binary_pheno_continuous_PS, all_covs, by = id_col)
    
    binary_covs_PS_formula <- as.formula(
        paste0("as.factor(pheno) ~ ", paste(covs_formu, collapse = " + "), " + scale(continuous_protS)")
        )
    
    binary_pheno_continuous_PS_assoc_mod <- glm(
        binary_covs_PS_formula, 
        family=binomial (link=logit), 
        data = binary_pheno_continuous_PS_covs
        )
    
    binary_pheno_continuous_PS_assoc_coefs <- summary(binary_pheno_continuous_PS_assoc_mod)$coefficients %>% as.data.frame() |> rownames_to_column(var = "Coefficient")
    
    # save the coefficients 
    binary_pheno_continuous_PS_assoc_coefs_outfile <- file.path(out_dir, continuous_weights_name, binary_pheno_name, paste0(cohort, "_", binary_pheno_name, "_", continuous_weights_name, "_PS_assoc_coefs.rds"))
    saveRDS(binary_pheno_continuous_PS_assoc_coefs, binary_pheno_continuous_PS_assoc_coefs_outfile)
    
    aucF <- auc(binary_pheno_continuous_PS_covs[complete.cases(binary_pheno_continuous_PS_covs),]$pheno, binary_pheno_continuous_PS_assoc_mod$linear.predictors)  ### AUC for full module

    print(paste0("AUC value for binary phenotype continuous weights: ", aucF))
    print(paste0("Incremental AUC value for binary phenotype continuous weights: ", aucF - aucR))

    binary_pheno_continuous_PS_predicted_probs <- predict(binary_pheno_continuous_PS_assoc_mod, type = 'response')
    roc_curve <- roc(binary_pheno_continuous_PS_covs[complete.cases(binary_pheno_continuous_PS_covs),]$pheno, binary_pheno_continuous_PS_predicted_probs)    
    
    # save ROC curve object for plotting all cohorts together
    saveRDS(roc_curve, file.path(out_dir, continuous_weights_name, binary_pheno_name, paste0(cohort, "_", binary_pheno_name, "_", continuous_weights_name, "_roc_curve.rds")))
    
    # ROC Graph 
    cairo_pdf(file = file.path(out_dir, continuous_weights_name, binary_pheno_name, paste0(cohort, "_", binary_pheno_name, "_", continuous_weights_name, "_binary_pheno_continuous_PS_assoc_ROC_curve.pdf")), width = 8, height = 6)
    plot.roc(roc_curve, col = "blue", lwd =2, main = paste0('ROC Curve: ', cohort, "_", binary_pheno_name, "_", continuous_weights_name))
    dev.off()
    
    pr_curve <- pr.curve(binary_pheno_continuous_PS_covs$pheno, binary_pheno_continuous_PS_predicted_probs, curve = T)
    
    cairo_pdf(file = file.path(out_dir, continuous_weights_name, binary_pheno_name, paste0(cohort, "_", binary_pheno_name, "_", continuous_weights_name, "_PS_assoc_precision_recall.pdf")), width = 8, height = 6)
    plot(pr_curve, col = "red", main= paste0('Precision Recall Curve: ', cohort, "_", binary_pheno_name, "_", continuous_weights_name))
    dev.off()

    binary_pheno_continuous_weights_desc_stats <- data.frame(
      weights = continuous_weights_name,
      Phenotype = binary_pheno_name,
      beta = summary(binary_pheno_continuous_PS_assoc_mod)$coefficients["scale(continuous_protS)", "Estimate"],
      std = summary(binary_pheno_continuous_PS_assoc_mod)$coefficients["scale(continuous_protS)", "Std. Error"],
      p = summary(binary_pheno_continuous_PS_assoc_mod)$coefficients["scale(continuous_protS)", "Pr(>|z|)"],
      R2 = NA,
      Incremental_R2 = NA,
      AUC = aucF,
      Incremental_AUC = aucF - aucR
      )

    desc_stats <- desc_stats |> rbind(binary_pheno_continuous_weights_desc_stats)
  }
}

if(opt$continuous_pheno!="NULL"){

  ##### Read in continuous pheno ###############################################################################################################

  continuous_pheno_filepath <- opt$continuous_pheno
  continuous_pheno <- readRDS(continuous_pheno_filepath)
  continuous_pheno_name <- opt$continuous_pheno_name
  
  colnames(continuous_pheno) <- c("id", "pheno")
  continuous_pheno <- continuous_pheno %>% filter(!is.na(pheno)) # remove missing values if there are any 

  if(opt$binary_weights!="NULL"){

    dir.create(file.path(out_dir, binary_weights_name, continuous_pheno_name))

    ##### Binary Protein Score #####
    
    continuous_pheno_binary_PS <- merge(continuous_pheno, ProtS_binary, by = "id")
    
    continuous_pheno_binary_PS_scatter <- ggplot(continuous_pheno_binary_PS, aes(pheno, binary_protS))+
        geom_jitter() +
        theme_classic()+
        labs(title="Protein score per continuous phenotype", x="Continuous phenotype", y="Protein Score") +
        geom_smooth(method = "lm", se = TRUE)
    
    ggsave(file.path(out_dir, binary_weights_name, continuous_pheno_name, paste0(cohort, "_", continuous_pheno_name, "_", binary_weights_name, "_scatter.png")), continuous_pheno_binary_PS_scatter, width = 8, height = 6, device='png', dpi=300)
  
    continuous_pheno_binary_PS_covs <- merge(continuous_pheno_binary_PS, all_covs, by = id_col)
    
    binary_PS_formula <- as.formula(
        paste0("scale(pheno) ~ scale(binary_protS) + ", paste(covs_formu, collapse = " + "))
    )

    continuous_pheno_binary_PS_assoc_mod <- glm(binary_PS_formula, 
                  #    family=binomial (link=logit), 
                    data = continuous_pheno_binary_PS_covs)
    
    continuous_pheno_binary_PS_assoc_coefs <- summary(continuous_pheno_binary_PS_assoc_mod)$coefficients %>% as.data.frame() |> rownames_to_column(var = "Coefficient")
    
    # save the coefficients 
    continuous_pheno_binary_PS_assoc_coefs_outfile <- file.path(out_dir, binary_weights_name, continuous_pheno_name, paste0(cohort, "_", continuous_pheno_name, "_", binary_weights_name, "_PS_assoc_coefs.rds"))
    saveRDS(continuous_pheno_binary_PS_assoc_coefs, continuous_pheno_binary_PS_assoc_coefs_outfile)

    rsq_value <- summary(lm(binary_PS_formula, data = continuous_pheno_binary_PS_covs))$r.squared
    incr_rsq_value <- summary(lm(binary_PS_formula, data = continuous_pheno_binary_PS_covs))$r.squared - summary(lm(as.formula(paste0("scale(pheno) ~ ", paste(covs_formu, collapse = " + "))), data = continuous_pheno_binary_PS_covs))$r.squared

    print(paste0("R squared value for continuous phenotype binary weights: ", rsq_value))
    print(paste0("Incremental R squared value for continuous phenotype binary weights: ", incr_rsq_value))

    png(filename=file.path(out_dir, binary_weights_name, continuous_pheno_name, paste0(cohort, "_", continuous_pheno_name, "_", binary_weights_name, "_PS_scatter_.rds")))
      ggplot(continuous_pheno_binary_PS, aes(pheno, binary_protS))+
        geom_jitter() +
        theme_classic()+
        labs(title="Protein score per continuous phenotype", x="Continuous phenotype", y="Binary Protein Score") +
        geom_smooth(method = "lm", se = TRUE)
    dev.off()

    continuous_pheno_binary_weights_desc_stats <- data.frame(
      weights = binary_weights_name,
      Phenotype = continuous_pheno_name,
      beta = summary(continuous_pheno_binary_PS_assoc_mod)$coefficients["scale(binary_protS)", "Estimate"],
      std = summary(continuous_pheno_binary_PS_assoc_mod)$coefficients["scale(binary_protS)", "Std. Error"],
      p = summary(continuous_pheno_binary_PS_assoc_mod)$coefficients["scale(binary_protS)", "Pr(>|t|)"],
      R2 = rsq_value,
      Incremental_R2 = incr_rsq_value,
      AUC = NA,
      Incremental_AUC = NA
      )
    
    desc_stats <- desc_stats |> rbind(continuous_pheno_binary_weights_desc_stats)
  }
  
  if(opt$continuous_weights!="NULL"){

    dir.create(file.path(out_dir, continuous_weights_name, continuous_pheno_name))

    ##### Continuous Protein Score #####
    
    continuous_pheno_continuous_PS <- merge(continuous_pheno, ProtS_continuous, by = "id")
    
    continuous_pheno_continuous_PS_scatter <- ggplot(continuous_pheno_continuous_PS, aes(pheno, continuous_protS))+
        geom_jitter() +
        theme_classic()+
        labs(title="Protein score per continuous phenotype", x="Continuous phenotype", y="Protein Score") +
        geom_smooth(method = "lm", se = TRUE)
    
    ggsave(file.path(out_dir, continuous_weights_name, continuous_pheno_name, paste0(cohort, "_", continuous_pheno_name, "_", continuous_weights_name, "_scatter.png")), continuous_pheno_continuous_PS_scatter, width = 8, height = 6, device='png', dpi=300)
  
    continuous_pheno_continuous_PS_covs <- merge(continuous_pheno_continuous_PS, all_covs, by = id_col)
    
    continuous_PS_formula <- as.formula(
        paste0("scale(pheno) ~ scale(continuous_protS) + ", paste(covs_formu, collapse = " + "))
    )

    continuous_pheno_continuous_PS_assoc_mod <- glm(continuous_PS_formula, 
                  #    family=binomial (link=logit), 
                    data = continuous_pheno_continuous_PS_covs)
    
    continuous_pheno_continuous_PS_assoc_coefs <- summary(continuous_pheno_continuous_PS_assoc_mod)$coefficients %>% as.data.frame() |> rownames_to_column(var = "Coefficient")
    
    # save the coefficients 
    continuous_pheno_continuous_PS_assoc_coefs_outfile <- file.path(out_dir, continuous_weights_name, continuous_pheno_name, paste0(cohort, "_", continuous_pheno_name, "_", continuous_weights_name, "_PS_assoc_coefs.rds"))
    saveRDS(continuous_pheno_continuous_PS_assoc_coefs, continuous_pheno_continuous_PS_assoc_coefs_outfile)

    rsq_value <- summary(lm(continuous_PS_formula, data = continuous_pheno_continuous_PS_covs))$r.squared
    incr_rsq_value <- summary(lm(continuous_PS_formula, data = continuous_pheno_continuous_PS_covs))$r.squared - summary(lm(as.formula(paste0("scale(pheno) ~ ", paste(covs_formu, collapse = " + "))), data = continuous_pheno_continuous_PS_covs))$r.squared

    print(paste0("R squared value for continuous phenotype continuous weights: ", rsq_value))
    print(paste0("Incremental R squared value for continuous phenotype continuous weights: ", incr_rsq_value))

    png(filename=file.path(out_dir, continuous_weights_name, continuous_pheno_name, paste0(cohort, "_", continuous_pheno_name, "_", continuous_weights_name, "_PS_scatter_.rds")))
      ggplot(continuous_pheno_continuous_PS, aes(pheno, continuous_protS))+
        geom_jitter() +
        theme_classic()+
        labs(title="Protein score per continuous phenotype", x="Continuous phenotype", y="Continuous Protein Score") +
        geom_smooth(method = "lm", se = TRUE)
    dev.off()

    continuous_pheno_continuous_weights_desc_stats <- data.frame(
      weights = continuous_weights_name,
      Phenotype = continuous_pheno_name,
      beta = summary(continuous_pheno_continuous_PS_assoc_mod)$coefficients["scale(continuous_protS)", "Estimate"],
      std = summary(continuous_pheno_continuous_PS_assoc_mod)$coefficients["scale(continuous_protS)", "Std. Error"],
      p = summary(continuous_pheno_continuous_PS_assoc_mod)$coefficients["scale(continuous_protS)", "Pr(>|t|)"],
      R2 = rsq_value,
      Incremental_R2 = incr_rsq_value,
      AUC = NA,
      Incremental_AUC = NA
      )

    desc_stats <- desc_stats |> rbind(continuous_pheno_continuous_weights_desc_stats)
  }
}

write.csv(desc_stats, file.path(out_dir, "descriptive_statistics.csv"))

sink()
