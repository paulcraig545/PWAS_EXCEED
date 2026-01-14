###############################################################################

# Set up libraries and options/files

###############################################################################
library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tidyr)
library(ggplot2)
library(tools)
library(ggpubr)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
# GS_weights contain two columns: proteins and weights
# pheno should have two columns: ID and MDD

option_list <- list(
  make_option('--cohort', type='character', help="Cohort", action='store'),
  make_option('--std_prot', type='character', help="The filepath for standardised protein file", action='store'),
  make_option('--id_column', type = 'character', default="ID", help = "Column names of participant identifier column", action = 'store'),
  make_option('--GS_weights', type = 'character', default = 'GS_MDD_prot_weights.rds', help = "Weights file provided for calculating PS"),
  make_option('--pheno', type = 'character', help = 'File path to MDD phenotype file'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

# cohort <- "GS" 
# std_prot_filepath <- "/Users/pcraig4/Desktop/PWAS_EXCEED/GS__std_proteins.rds" 
# id_col <- "id" 
# GS_weights_filepath <- "/Users/pcraig4/Desktop/data_for_lasso/GS_lasso_weights_MDD.rds" 
# pheno_filepath <- "/Users/pcraig4/Desktop/data_for_lasso/mdd_pheno.rds"
# outdir <-  "/Users/pcraig4/Desktop/PWAS_EXCEED/"

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
std_prot_filepath=opt$std_prot
GS_weights_filepath=opt$GS_weights # Weights file from GS
id_col <- opt$id_column # Vector of identifier columns 
pheno_filepath=opt$pheno
outdir <- opt$outdir

# Create a log file
sink(paste0(outdir, cohort, "_MDD_ps.log"))
print(paste0('Calculating protein scores for ', cohort))
print(paste0('Read in the processed protein file from: ', std_prot_filepath))
print(paste0('Read in the weights file from: ', GS_weights_filepath))
print(paste0('Read in the pheno file from: ', pheno_filepath))
print(paste0('Output to be saved in: ', outdir))

###############################################################################

# Read the files

###############################################################################

std_prot <- readRDS(std_prot_filepath) 
print(paste0('The protein file has data for ', nrow(std_prot), ' participants and ', std_prot %>% select(-all_of(id_col)) %>% ncol(), " proteins"))  # Remove ID col
GS_weights <- readRDS(GS_weights_filepath) 
colnames(GS_weights) = c('proteins','weights')
print(paste0('The weights file has weights for ', nrow(GS_weights), ' proteins'))

###############################################################################

# Assess missingness of proteins

###############################################################################

if(std_prot %>% select(-all_of(id_col)) %>% ncol() != nrow(GS_weights)){
  print(paste0('Number of non-zero coef proteins read in for protein file (', std_prot %>% select(-all_of(id_col)) %>% ncol(),
               ') does not match the number of weights provided (', nrow(GS_weights), ')'))
  # If not all weights included, save the proteins included in the PS for the cohort 
  readRDS(std_prot %>% select(-all_of(id_col)) %>% colnames(), paste0(outdir,cohort, '_proteinPS.txt'))
  print(paste0('proteins used in the PS are saved in: ', outdir, cohort, '_proteinPS.txt'))
} else {
  print(paste0('Number of probes read in for std_prot matches the number of weights provided: n = ', nrow(GS_weights)))
}

missing_percentage <- std_prot %>% select(-all_of(id_col)) %>%
  summarise_all(~ mean(is.na(.)) * 100) %>%
  gather(proteins, MissingPercentage) %>% arrange(desc(MissingPercentage))

print(paste0('The protein with the highest level of missingness is: ', missing_percentage$proteins[1]))
print(paste0('There are ', nrow(missing_percentage %>% filter(MissingPercentage > 5)), ' proteins with more than 5% missingness and ',
             nrow(missing_percentage %>% filter(MissingPercentage > 50)), ' with more than 50% missingness'))

print(missing_percentage)

missing_hist <- ggplot(missing_percentage, aes(x = MissingPercentage, fill = MissingPercentage < 50)) + 
  geom_histogram() + 
  theme_minimal() + 
  labs(x = '% of Missingness', y = 'Number of protein')+
  geom_vline(xintercept = 50, color = "red", linetype = 'dashed')+
  ggtitle(cohort)

missing_plot <- ggplot(missing_percentage %>% filter(MissingPercentage > 50),
                           aes(x = reorder(proteins, -MissingPercentage), y = MissingPercentage)) +
    geom_bar(stat='identity', fill = 'skyblue', color = 'black') + 
    labs(title = paste0(cohort, ': protein missingness in PS'), x = "protein", y = "% Missing") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# Plot the % of missingness against the absolute value of the PS coefficient 

missingness_weights <- merge(missing_percentage, GS_weights, by = 'proteins')

missing_weights_plt <- ggplot(missingness_weights, aes(x = MissingPercentage, y = abs(weights))) + 
  geom_point() + 
  theme_minimal() + 
  labs(x = '% of Missingness', y ='PS weight', 
       title = paste0(cohort, ': % missing vs PS weights for proteins > 50% missingness'))

all_missing_plots <- ggarrange(missing_hist, missing_weights_plt, missing_plot, nrow = 3, ncol = 1)

ggsave(paste0(outdir, cohort, "_ps_prot_missingness.png"), missing_plot, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_ps_prot_missingness_hist.png"), missing_hist, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_ps_prot_missingness_weights.png"), missing_weights_plt, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_ps_prot_missingness_allplots.png"), all_missing_plots, width = 8, height = 6, device='png', dpi=300)

write.table(missing_percentage, paste0(outdir, cohort, '_prot_missingness.txt'), row.names = F, quote = F)

###############################################################################

# Calculate the protein scores 

###############################################################################
print('Converting std_prot to long format')
long_std_prot <- std_prot %>% 
  pivot_longer(-c(all_of(id_col)), 
               names_to = "proteins", 
               values_to= "Protval")
print('Merging with the weights')
long_std_prot <- merge(long_std_prot, GS_weights, by = 'proteins')

print('Calculating the protein scores')

# Group by ID
# and then calculate a weighted sum of all the proteins per ID

ProtS <- long_std_prot %>% group_by(!!sym(id_col)) %>%
  summarise(weighted_sum = sum(weights*Protval, na.rm = T)) %>% as.data.frame()

print(paste0('protein scores calculated for ', nrow(ProtS), ' participants in the ', cohort, ' cohort'))

###############################################################################

# Distribution of the number of proteins within the ProtS (non-missing)

###############################################################################

# Showing as similar metric to the above missing plots, but another way of looking at it 

num_prot <- long_std_prot %>%
  group_by(!!sym(id_col)) %>%
  summarise(prot_ps = sum(!is.na(Protval)))

prot_ps_hist <- ggplot(num_prot, aes(x = prot_ps)) + 
  geom_histogram() +
  theme_minimal() + 
  labs(x = 'Number of proteins included in the ProtS', y = 'Frequency') + 
  ggtitle(paste0(cohort, ': Histogram of number of proteins included in the ProtS'))

ggsave(paste0(outdir, cohort, "_indiv_numprot_ps.png"), prot_ps_hist, width = 8, height = 6, device='png', dpi=300)

###############################################################################

# protein Score Distributions

###############################################################################

print(paste0('Plotting distributions across the whole ', cohort, ' sample'))

# Distribution of the protein score across the cohort 

PS_dist <- ggplot(ProtS, aes(x = weighted_sum)) + 
  geom_histogram() + 
  theme_minimal() + 
  labs(x = 'protein Profile Score', y = 'Count')+
  ggtitle(cohort)

ggsave(paste0(outdir, cohort, "_MDD_ps_overalldist.png"), PS_dist, width = 8, height = 6, device='png', dpi=300)

# look at Distribution in MDD cases and controls (violin plots)

if (endsWith(pheno_filepath, '.rds')){
mdd_pheno <- readRDS(pheno_filepath)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .rds file')
}

if('MDD' %in% colnames(mdd_pheno) == FALSE){
  stop('No MDD column in the phenotype file')
} else {
  print('MDD column in the phenotype file')
}

mdd_pheno <- mdd_pheno %>% filter(!is.na(MDD)) # remove missing values if there are any 

print(paste0('Read in the MDD phenotype for ', cohort, ': Number of cases: ',
             nrow(mdd_pheno %>% 
                    filter(MDD==1)), 
             ' Number of controls: ',
             nrow(mdd_pheno%>% 
                    filter(MDD==0))))

mdd_pheno <- mdd_pheno |> rename(!!id_col:= ID )
mdd_pheno_PS <- merge(mdd_pheno, ProtS, by = id_col)

print('Plotting ProtS distributions for MDD exposure cases and controls ')
PS_pheno_dists <- ggplot(mdd_pheno_PS, aes(x = weighted_sum, fill = as.factor(MDD))) + 
  geom_histogram(alpha = 0.8) + 
  theme_minimal() + 
  labs(x = 'protein Profile Score', y = 'Count', fill = 'MDD') +
  ggtitle(cohort)

ggsave(paste0(outdir, cohort, "_MDD_ps_phenodist.png"), PS_pheno_dists, width = 8, height = 6, device='png', dpi=300)

###############################################################################

# Saving the protein score 

###############################################################################
outfile <- file.path(outdir, paste0(cohort, "_MDD_protS.rds"))
print(paste0('Saving the protein score to ', outfile))

colnames(ProtS)[2] <- 'MDD_protS'
saveRDS(ProtS, outfile)
sink()
