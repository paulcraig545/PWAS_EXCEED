
###### Set up libraries and options/files

library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tools)
library(ggplot2)
library(tidyr)
library(tidyverse)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort", action='store'),
  make_option('--proteins', type='character', help="The filepath for proteins file", action='store'), 
  make_option('--probe', type = 'character', help= "The filepath for the list of proteins", action = 'store'),
  make_option('--id_column', type = 'character', default="ID", help = "Column name for identifier column", action = 'store'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)


# setting up arguments from the options 
print('Setting up the options')
cohort <- opt$cohort
prot_filepath=opt$proteins # protein file
probe_filepath=opt$probe # List of protein names
id_col <- opt$id_column # Vector of identifier column
out_dir <- opt$outdir

sink(paste0(out_dir, cohort, "_", "_std_proteins.log"))  #log file

print(paste0('protein file from : ', prot_filepath))
print(paste0('List of proteins from : ', probe_filepath))
print(paste0('ID column : ', id_col))
print(paste0('Output to be saved in : ', out_dir))

###############################################################################

# Read the files

###############################################################################

if (endsWith(prot_filepath, ".rds")){
  proteins <- readRDS(prot_filepath)
} else {
  stop("Unsupported file format. Please provide a file with .rds format")
}

prot_list <- readRDS(probe_filepath)
print('Read in files')

###############################################################################

# Filter the proteins to just the proteins used in the LASSO training 

###############################################################################
print('Filter to just the LASSO proteins')
ms_prot <- proteins %>% select(c(all_of(id_col), intersect(names(proteins), prot_list)))

print(paste0('Filtered to ', ms_prot %>% select(-id_col) %>% ncol())) 
rm(proteins) # remove large protein object 

###############################################################################

# Apply quantile normalisation the protein levels 

###############################################################################
print('Scaling the protein columns')

quantile_normalisation <- function(df){
  
  # Find rank of values in each column
  df_rank <- map_df(df,rank,ties.method="average")
  # Sort observations in each column from lowest to highest 
  df_sorted <- map_df(df,sort)
  # Find row mean on sorted columns
  df_mean <- rowMeans(df_sorted)
  
  # Function for substiting mean values according to rank 
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  # Replace value in each column with mean according to rank 
  df_final <- map_df(df_rank,index_to_mean, my_mean=df_mean)
  
  return(df_final)
}

# Quantile normalisation of new protein data
prot_std <- ms_prot |> 
  select(-c(id)) |> 
  quantile_normalisation() |>
  mutate(id = ms_prot$id) |>
  select(c(id, everything())) |>
  as.data.frame()

# Plotting the distribution of unstandardised and standardised values
# for 3 randomly selected proteins

prot_std_compare <- sample(setdiff(names(prot_std), id_col), 3)
prot_both <- rbind(
  ms_prot %>% 
  select(c(all_of(id_col), all_of(prot_std_compare))) %>%
    pivot_longer(cols = -c(all_of(id_col)), 
    names_to = "proteins", 
  values_to = "protval") %>% 
  as.data.frame() %>% 
  mutate(Values = 'Unstandardised'),
prot_std %>% 
  select(c(all_of(id_col), all_of(prot_std_compare))) %>%
  pivot_longer(cols = -c(all_of(id_col)), 
               names_to = "proteins", 
               values_to = "protval") %>%
  mutate(Values = 'Standardised')
)

prot_dists <- ggplot(prot_both, aes(x = protval, fill = proteins)) +
  geom_histogram() + 
  facet_grid(Values~proteins) +
  ggtitle(paste0(cohort, ': Random sample of proteins - standardisation'))

ggsave(filename=paste0(out_dir, cohort, "_", "_std_unstd_proteins.png"),prot_dists, 
       width = 8, height = 6, device='png', dpi=300)

###############################################################################

# Write out the filtered and standardised protein data 

###############################################################################

outfile <- paste0(out_dir, cohort, "_", "_std_proteins.rds")
saveRDS(prot_std, outfile)
sink()

