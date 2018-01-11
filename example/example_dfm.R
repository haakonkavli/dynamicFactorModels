# Example use of dfm on South African Macro data


rm(list = ls())
source('../R/dfm.R')

# Load packages:
libs <- c('data.table','MARSS')
lapply(libs, require, character.only = TRUE)

# Load example data:
load('dat_input_example.Rdata')

# dat_input is an (m, 1+Nx) data.table, where:
#   m is number of examples
#   Nx is number of variables
#   The added column holds datestamps


# Must transpose dat_input to (Nx, m) and remove datestamp column:

dat_input <- t(dat_input[,-c('datestamp'),with = F])


# Specify hyper parameters:
n_factors <- 1
model_structure <- 'diagonal and equal'


# Estimate dfm using dfm function: (Will take a while to run)
pea_output <- dfm(dat_input = dat_input,
                  n_factors = 1,
                  model_structure = model_structure)



plot(ts(pea_output$trends, pea_output$fitted_data), type = 'l')





