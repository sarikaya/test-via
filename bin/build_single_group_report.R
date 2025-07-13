#!/usr/bin/env Rscript
#build_single_group_report.R
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
  make_option(c("--s_transcript"), type = "character", help = "Select transcript (SQ-number)")
  )
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
s_transcript = opt$s_transcript

# FUNCTIONS ---------------------------------------------------------------
source(system("which end_seq_functions.R",  intern = TRUE))

# LOAD DATA ---------------------------------------------------------------
ll = load("data_prepared_for_plotting.rda")

# APPLY PARAMETERS --------------------------------------------------

#PARAMETERS FOR THIS RUN
source('run_parameters.R')

#S3 VARS
source(system("which s3_vars.R",  intern = TRUE))

# GLOBAL VARIABLES --------------------------------------------------------

#FOR GROUP-LEVEL TRANSCRIPT REPORTS
seq_column = 'seq_context' #column to use for sequence logos
# dir.create("report_data", showWarnings = FALSE)
# dir.create(file.path(groups_output_dir), showWarnings = FALSE)
# BUILD GROUP-LEVEL REPORT FOR S_TRANSCRIPT ------------------------------------------

#BUILD REPORT OBJECT
if (s_transcript == 'All'){
  report_title = 'All transcripts - group report'
  report_written = build_all_group_report(s_transcript)
} else {
  report_object = build_group_report_obj(s_transcript)
  report_written = render_individual_report(report_object)
}
# SEND TO S3 --------------------------------------------------------------

#create the folder
#aws.s3::put_folder(folder = paste0(s3_groups_prefix),
#                   bucket = s3_bucket)
#write

#report_to_s3(report_written, s3_bucket, s3_groups_prefix)
#print(str_c(basename(report_written), ' written to ', str_c(s3_bucket, '/', s3_groups_prefix), basename(report_written)))
