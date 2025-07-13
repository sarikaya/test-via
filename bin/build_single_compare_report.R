#!/usr/bin/env Rscript
#build_single_compare_report.R
suppressMessages(library(optparse))
suppressMessages(library(aws.s3))
suppressMessages(library(htmltools))
suppressMessages(library(plotrix))
suppressMessages(library(plotly))
suppressMessages(library(ggseqlogo))
suppressMessages(library(tidyverse))
suppressMessages(library(knitr))
suppressMessages(library(kableExtra))
suppressMessages(library(cowplot))
theme_set(theme_cowplot())
knitr::opts_chunk$set(echo = TRUE)
options(knitr.table.format = "markdown")
options(warn=-1)

# ARGUMENTS ---------------------------------------------------------------

option_list = list(
  make_option(c("--input_table"), type="character", default='comparisons_table.csv', 
              help="Input groups table (output from prepare_for_plotting.R)"),
  make_option(c("--c_num"), type="integer",
              help="Comparison number for this report ")
  )
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_table = opt$input_table
c_num = opt$c_num

# FUNCTIONS ---------------------------------------------------------------

source(system("which end_seq_functions.R",  intern = TRUE))

# LOAD DATA ---------------------------------------------------------------

ll = load('data_prepared_for_plotting.rda')
full_compare_df = read_csv(input_table,
                           show_col_types = FALSE) #this table is build by prepare_for_plotting.R

# APPLY PARAMETERS --------------------------------------------------

#INPUTS BASED ON DN RUN
source('run_parameters.R')

#S3 VARS
source(system("which s3_vars.R",  intern = TRUE))

#GET COMPARISON VARIABLES FROM SINGLE ROW OF full_compare_df BASED ON c_num
comparison_row = full_compare_df %>% 
  dplyr::slice(c_num)
comparison_num = comparison_row$comparison_num
x_group = comparison_row$x_group
y_group = comparison_row$y_group
x_library_type = comparison_row$x_library_type
y_library_type = comparison_row$y_library_type

# GLOBAL VARIABLES --------------------------------------------------------

#BUILD LOCAL OUTPUT DIRS
#dir.create("report_data", showWarnings = FALSE)
#dir.create(file.path(comparison_output_dir), showWarnings = FALSE)


# BUILD COMPARISON REPORTS -------------------------------------

#WRITE REPORT OBJECT
report_object = build_comparison_report_obj(comparison_num,
                            x_group,
                            y_group,
                            x_library_type,
                            y_library_type)

#RENDER REPORT
report_written = render_individual_report(report_object)


# SEND TO S3 --------------------------------------------------------------

#create the folder
# aws.s3::put_folder(folder = paste0(s3_comparisons_prefix),
#                   bucket = s3_bucket)

# #write
# report_to_s3(report_written, s3_bucket, s3_comparisons_prefix)
# print(str_c(basename(report_written), ' written to ', str_c(s3_bucket, '/', s3_comparisons_prefix), basename(report_written)))
