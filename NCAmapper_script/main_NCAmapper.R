# Non-Contributing Area mapping model
# NCAmapper (v 1.10)
# © Mohamed Ismaiel Ahmed 2024 @ UCalgary
# This software is provided under the MIT license
# A copy of this script can be found on GitHub at https://github.com/MIsmlAhmed/NCAmapper

# NCAmapper main script
# clear envirnment variables
cat("\014") 
rm(list = ls())
# load NCAmapper functiond
source("functions_NCAmapper.R")
# main code #
##########################################
# Specify inputs (location where NCAmapper_config.ini and other inputs are located)
inp_dir <- '/Users/mohamed/Documents/scratch/SCRB'
# read config file and initialize NCAmapper
config_file <- init_NCAmapper_config(inp_dir)
##########################################
#preprocess the data using Whiteboxtools
preprocessing_wbt(config_file)
####################
#get summary of depressions
dep_smr <- get_depressions_summary(config_file)

####################
# Start the depressional storage & NCA component
depression_state <- simulate_depStorage_NCA(config_file, dep_smr)

