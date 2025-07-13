#!/usr/bin/env Rscript
#build_combined_group_report.R
suppressMessages(library(optparse))
suppressMessages(library(aws.s3))
suppressMessages(library(htmltools))
suppressMessages(library(tidyverse))
suppressMessages(library(knitr))
suppressMessages(library(kableExtra))
knitr::opts_chunk$set(echo = TRUE)
options(knitr.table.format = "markdown")
print('Running build_combined_group_report.R')

# ARGUMENTS ---------------------------------------------------------------

option_list = list(
  make_option(c("--input_table"), type="character", default='groups_table.csv', 
              help="Input groups table (output from prepare_for_plotting.R)")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_table = opt$input_table

# FUNCTIONS ---------------------------------------------------------------

source(system("which end_seq_functions.R",  intern = TRUE))


# LOAD DATA ---------------------------------------------------------------

ll = load('data_prepared_for_plotting.rda')

# APPLY PARAMETERS ---------------------------------------------------------------

#parameters for this run
source('run_parameters.R')

#s3 vars
source(system("which s3_vars.R",  intern = TRUE))

#input groups table with transcripts to write into combined report
dat = read_csv(input_table, show_col_types = FALSE)

#dolphinnext run information for report
dn_run_url = "{{DNEXT_RUN_URL}}"
dn_tag_list = list(link = a(href = dn_run_url,
                            str_c('DolphinNext Run'),
                            target = "_blank"))

#format significance parameters
sig_param_rep = format_sig_params(sig_param_df)

#vectors for running report loop below
transcripts_written = dat$transcript
group_report_urls = str_c(s3_url_prefix, dat$s3)
names(group_report_urls) = transcripts_written


# BUILD COMBINED REPORT ---------------------------------------------------

#set report variables
report_title = 'Group reports'
report_description = "This report contains links to results for each individual group organized by transcript."
#set up tag_list with entry for All sqs first
all_id = 'All'
all_url = str_c(s3_url_prefix, s3_bucket, '/', s3_groups_prefix, all_id, '.html')
tag_list = list('All - header' = h3(all_id),
                'All - link' = a(href = all_url,
                                 all_id,
                                 target = "_blank"))

#set up report description
#building urls into tag_list:
for (tw in transcripts_written){
  tag_list[[paste(tw, 'header', sep='-')]] = h3(tw)
  tag_list[[paste(tw, 'link', sep='-')]] = a(href = group_report_urls[tw],
                                             tw,
                                             target = "_blank")
}

#render the report
output_file = 'combined_groups.html'
rmarkdown::render(input = combined_markdown_script, 
                  output_format = "html_document",
                  output_file = output_file,
                  output_dir = groups_output_dir)


# SEND TO S3 --------------------------------------------------------------

# local_file_path = paste(groups_output_dir, output_file, sep='/')
# report_to_s3(local_file_path, s3_bucket, s3_groups_prefix)
# combined_groups_url = str_c(s3_url_prefix, s3_bucket, '/', s3_groups_prefix, basename(local_file_path))
# print('link to combined groups report:')
# print(combined_groups_url)

