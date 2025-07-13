#set_S3_vars.R
#set varaibles for copying to S3
#requires s3_prefix set in run parameters

#S3 STRINGS USED AT MULTIPLE STEPS
DNEXT_PUBLISH_DIR = "{{DNEXT_PUBLISH_DIR}}" 
s3_bucket = 's3-informatics'
s3_prefix =str_c(sub(paste("s3://",s3_bucket,"/", sep=""), "", DNEXT_PUBLISH_DIR),"/")

#handling network fails
s3_path = str_c('s3://', s3_bucket, '/', s3_prefix)
s3_attempts = 8 #number of times to try copying to S3 before giving up
s3_sleep_seconds = 3
#url strings
s3_url_prefix = 'https://s3.modernatx.net/'
s3_non_url_prefix = 's3://'
s3_url = str_replace(s3_path, s3_non_url_prefix, s3_url_prefix)
#groups strings
s3_groups_folder = 'downstream_groups'
s3_groups_prefix = str_c(s3_prefix, s3_groups_folder, '/')
#comparisons strings
s3_comparisons_folder = 'downstream_comparisons'
s3_comparisons_prefix = str_c(s3_prefix, s3_comparisons_folder, '/')