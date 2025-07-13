#!/usr/bin/env Rscript
#build_combined_comparison_report.R
suppressMessages(library(optparse))
suppressMessages(library(htmltools))
suppressMessages(library(tidyverse))
suppressMessages(library(knitr))
suppressMessages(library(kableExtra))
knitr::opts_chunk$set(echo = TRUE)
options(knitr.table.format = "markdown")
print('Running build_combined_comparison_report.R...')

# ARGUMENTS ---------------------------------------------------------------

option_list = list(
  make_option(c("--input_table"), type="character", default='comparisons_table.csv', 
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


#compile dataframe and revise group names to include library types for clarity

get_written_path = function(comparison_num){
  str_c(comparison_output_dir, '/', comparison_num, '.html')
}
get_compare_url = function(f){
  str_c(s3_url, s3_comparisons_folder, '/', basename(f))
}

written_compare_df = full_compare_df %>% 
  mutate(x_group_libtype = str_c(x_group, x_library_type, sep='-'),
         y_group_libtype = str_c(y_group, y_library_type, sep='-'),
         written = get_written_path(comparison_num),
         url = get_compare_url(written))


# BUILD COMBINED COMOPARISON REPORT ---------------------------------------------------
#use the dataframe compare_df to build named links to each comparison in a single rmarkdown report

#set up empty taglist
tag_list = list()

#set report variables
report_title = 'Comparison reports'

#report description
report_description = "Comparisons between treatment groups are organized based on the library types compared.
                             Subheaders indicate the group used for the X axis in the comparison.
                             Link names give the group used for the Y axis in the comparison."

#building urls into tag_list:
for (libtype_pair in unique(written_compare_df$library_type_pair)){
  lt_compare_df = written_compare_df %>% 
    filter(library_type_pair == libtype_pair)
  libtype_pair_header = str_replace(libtype_pair, '_', ' vs ')
  libtype_pair_id = paste(libtype_pair, 'header', sep='-')
  tag_list[[libtype_pair_id]] = h2(libtype_pair_header)
  for (xg in unique(lt_compare_df$x_group_libtype)){
    xg_id = paste0(libtype_pair, 'X=', xg)
    xg_header = paste0(xg)
    tag_list[[xg_id]] = h3(xg_header)
    xg_compare_df = lt_compare_df %>% 
      filter(x_group_libtype == xg)
    for (i in 1:nrow(xg_compare_df)){
      yg = xg_compare_df$y_group_libtype[i]
      url = xg_compare_df$url[i]
      link_id = paste(libtype_pair, xg, yg, 'link')
      link_text = paste(yg)
      tag_list[[link_id]] = a(href = url,
                              link_text,
                              target = "_blank")
      tag_list[[paste(link_id, '-')]] = p('')
    }
  }
}

#render the report
output_file = 'Combined_comparisons.html'
rmarkdown::render(input = combined_markdown_script, 
                  output_format = "html_document",
                  output_file = output_file,
                  output_dir = comparison_output_dir)

# SEND TO S3 --------------------------------------------------------------
# local_file_path = paste(comparison_output_dir, output_file, sep='/')
# report_to_s3(local_file_path, s3_bucket, s3_comparisons_prefix)
# combined_comparisons_url = str_c(s3_url_prefix, s3_bucket, '/', s3_comparisons_prefix, basename(local_file_path))
# print('link to combined comparisons report:')
# print(combined_comparisons_url)

