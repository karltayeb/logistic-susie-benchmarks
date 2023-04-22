# A collection of benchmarks for logistic SuSiE

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint

library(future)

library(future.callr)
plan(callr)

# library(future.batchtools)
# plan(batchtools_slurm, resources=list(account='pi-mstephens', time='1:0:0'))

# Set target options:
tar_option_set(
  packages = c("tibble", "logisticsusie", "dplyr", "tidyr", "purrr", "ggplot2", "susieR", "stringr"), # packages that your targets need to run
  format = "rds", # default storage format
  workspace_on_error = TRUE # Save a workspace file for a target that errors out.
  # Set other options as needed.
)

tar_source()

# All targets -------
list(
  make_X_targets(),
  half_normal_target,
  constant_sim_target,
  # realX_constant_sim_target,
  website,
  yusha_example_targets
)
