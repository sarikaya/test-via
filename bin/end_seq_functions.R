

# GLOBAL VARS -------------------------------------------------------------

#R SCRIPT FOR S3 VARIABLES
s3_vars = 's3_vars.R'

#PATH STRINGS USED AT MULTIPLE STEPS
groups_output_dir = '.'
comparison_output_dir = '.'

#MARKDOWN SCRIPT NAMES
individual_markdown_script = system("which individual_report.Rmd",  intern = TRUE) #the markdown script to build individual group- and comparison reports
combined_markdown_script = system("which combined_report.Rmd",  intern = TRUE)    #markdown script to build combined group- and comparison reports

#DN RUN VARIABLES
dn_run_url_prefix = 'https://dnext.modernatx.net/index.php?np=3&id='

# PROCESSING BIGWIGS ------------------------------------------------------
#these functions will replace those based on bedgraphs below
#the functions were tested and produce the same output as the bedgraph functions
#reading in from bigwigs is better


bigwig_to_df = function(bw_file, score_column){
  rtracklayer::import(bw_file) %>% 
    data.frame() %>% 
    tibble() %>% 
    dplyr::select(-width, -strand) %>% 
    set_names(c('chr', 'start', 'stop', score_column))
}

read_lib_bigwig = function(lib) {
  # print(lib)
  pos_file = paste0( lib, '.pos.raw.bw')
  neg_file = paste0( lib, '.neg.raw.bw')
  p = bigwig_to_df(pos_file, 'pos.depth')
  n = bigwig_to_df(neg_file, 'neg.depth')
  d = p %>%
    left_join(n, by = c('chr', 'start', 'stop')) %>%
    replace_na(list(pos.depth = 0,
                    neg.depth = 0)) %>%
    mutate(library = lib,
           pos = stop) %>%
    dplyr::select(-start, -stop) %>%
    arrange(chr, pos)
  return(d)
}

#function to get df for set of libs from bedgraphs
get_bw_trace_df = function(libs, metadata_df) {
  depth_df = map(libs, read_lib_bigwig) %>%
    purrr::reduce(rbind) %>%
    left_join(metadata_df, by = 'library') %>%
    mutate(
      library = factor(library, levels = libs),
      group = factor(group, levels = unique(metadata_df$group)),
      pos = as.numeric(pos)
    )
  return(depth_df)
}

# PROCESSING BEDGRAPHS ----------------------------------------------------

#function to read in the positive and negative strand bedGraphs for a library
read_lib_trace = function(lib) {
  pos_bed = paste0(indir, lib, '.pos.sorted.bedGraph')
  neg_bed = paste0(indir, lib, '.neg.sorted.bedGraph')
  p = read_tsv(
    pos_bed,
    col_names = c('chr', 'start', 'stop', 'pos.depth'),
    show_col_types = FALSE
  )
  n = read_tsv(
    neg_bed,
    col_names = c('chr', 'start', 'stop', 'neg.depth'),
    show_col_types = FALSE
  )
  d = p %>%
    left_join(n, by = c('chr', 'start', 'stop')) %>%
    replace_na(list(pos.depth = 0,
                    neg.depth = 0)) %>%
    mutate(library = lib,
           pos = stop) %>%
    dplyr::select(-start, -stop) %>%
    arrange(chr, pos)
  return(d)
}

#function to get df for set of libs from bedgraphs
get_trace_df = function(libs, metadata_df) {
  depth_df = map(libs, read_lib_trace) %>%
    purrr::reduce(rbind) %>%
    left_join(metadata_df, by = 'library') %>%
    mutate(
      library = factor(library, levels = libs),
      group = factor(group, levels = unique(metadata_df$group))
    )
  return(depth_df)
}

#strand counting
#for now just have combining, but could make this more sophisticated based
#on the gene's strand if needed
count_strand_selection = function(dat, count_strand = 'both') {
  if (count_strand == 'both') {
    print('Using depth summed accross both strands')
    dat = dat %>%
      mutate(depth = pos.depth + neg.depth)
  }
  if (count_strand == 'pos') {
    print('Using depth from only positive strand')
    dat = dat %>%
      mutate(depth = pos.depth)
  }
  dat %>%
    dplyr::select(-pos.depth, -neg.depth)
}

##function that will report % endogenous ends from output of read_depth()
get_percent_endogenous = function(dat) {
  dat %>%
    mutate(RNA_source = if_else(grepl("^SQ", chr),
                                'SQ',
                                'endogenous')) %>%
    group_by(library, library_type, group, RNA_source) %>%
    summarize(summed_reads = sum(depth)) %>%
    pivot_wider(names_from = RNA_source,
                values_from = summed_reads) %>%
    ungroup() %>%
    mutate(
      pct.endogenous = endogenous / (endogenous + SQ) * 100,
      pct.endogenous = round(pct.endogenous, digits = 3)
    )
}


# PROCESS METADATA --------------------------------------------------------

get_long_sq_df = function(metadata_df){
  long_sq_df = metadata_df %>% 
    select(library, group, library_type, sq) %>% 
    separate_rows(sq, sep=';')
  return(long_sq_df)
}

get_dosed_lib_list = function(metadata_df){
  long_sq_df = get_long_sq_df(metadata_df)
  unique_sqs = unique(long_sq_df$sq)
  libs_dosed = list()
  for (my_sq in unique_sqs){
    libs_dosed[[my_sq]] = long_sq_df %>% 
      filter(sq==my_sq) %>% 
      pull(library)
  }
  return(libs_dosed)
}

# SET UP DATA FOR TRANSCRIPT ----------------------------------------------

#HANDLING FASTA INFORMAITON
#convert a fasta object to a dataframe
fasta_to_tibble = function(fa_object){
  to_df = function(seq_id){
    seq_rec = fa_object[[seq_id]]
    seq_vec = getSequence(seq_rec)
    df = tibble('seq_id' = seq_id,
                'pos' = 1:length(seq_vec),
                'seq' = seq_vec)
  }
  map(names(fa_object), function(seq_id) to_df(seq_id)) %>% 
    bind_rows()
}

#get list of sequence vectors from a fasta object
fasta_to_vector_list = function(fa_object){
  res = getSequence(fa_object)
  names(res) = names(fa_object)
  return(res)
}

#FUNCTIONS FOR WRANGING BIGWIG DATA BASED ON TRANSCRIPT/EXON COORDINATES
#subset for gene
subset_trace_for_gene = function(dat, gchr, gstart, gstop, trim_left = 0, trim_right = 0) {
  left_cut = gstart + trim_left
  right_cut = gstop - trim_right
  dat %>%
    dplyr::filter(pos >= left_cut,
                  pos <= right_cut) %>%
    dplyr::filter(if_any(starts_with("chr"), ~ .x == gchr))
}


build_zero_list = function(start, stop) {
  pos = start:stop
  res = data.frame('pos' = pos,
                   'zero_depth' = 0)
}

build_zero_tibble = function(select_exon_df, keep = FALSE) {
  if (!keep) {
    # this was the original code
    starts = select_exon_df$start
    stops = select_exon_df$stop
    map2(starts, stops, function(start, stop) build_zero_list(start, stop)) %>% 
    bind_rows() %>% 
    tibble() %>% 
    distinct(pos, .keep_all = TRUE)
  } else {
    # this is to maintain chromosome info etc.
    group_split(select_exon_df, gene, chrom, start, stop) %>% map_dfr(~ tibble(
      gene = .x$gene,
      chrom = .x$chrom,
      pos = .x$start:.x$stop,
      zero_depth = 0
    )) %>% distinct()
  }
}

merge_in_zeros = function(l, e, g, zdat) {
  zdat$lib = l
  zdat$end = e
  zdat$group = g
}


add_zero_positions = function(s_depth_df, select_exon_df, keep = FALSE) {
  #-----------------------------------------------------------------------------
  # commented out to work with the normalization strategy -- dont want
  #   to lose the metadata information
  #subset for minimal variables
  # s_depth_df = s_depth_df %>%
  # dplyr::select(c('library', 'pos', 'depth'))
  
  #build set of zeros for each exon position (originally based on gtf)
  zdat = build_zero_tibble(select_exon_df, keep)
  #build dataframe with each unique library paired with all zero positions
  libs = as.character(unique(s_depth_df$library))
  lib_zdat_list = list()
  for (l in libs) {
    lib_zdat_list[[l]] = zdat %>%
      mutate(library = l) %>%
      data.frame()
  }
  lib_zdat = lib_zdat_list %>%
    bind_rows() %>%
    tibble()
  #merge zero positions back with original data, keeping the true counts and zeros for positions without counts
  if(!keep) {
    # Groves original code
    z_depth_df = s_depth_df %>%
    full_join(lib_zdat, by = c('pos', 'library')) %>%
    arrange(library, pos) %>%
      replace_na(list(depth = 0)) %>%
      dplyr::select(-zero_depth)  
  } else {
    # adapting to case when we want to retain a more information
    z_depth_df = s_depth_df %>%
      full_join(lib_zdat, by = c('chrom', 'pos', 'library')) %>%
      arrange(library, pos) %>%
      replace_na(
        list(
          depth = 0,
          group = .$group[!is.na(.$group)][1],
          library_type = .$library_type[!is.na(.$library_type)][1],
          sq = .$sq[!is.na(.$sq)][1]
        )
      ) %>%
      dplyr::select(-zero_depth)
  }
  
  return(z_depth_df)
}

add_pseudocount = function(depth_df, pseudocount) {
  depth_df %>%
    mutate(pseudo_depth = depth + pseudocount)
}

# Prepare transcript depth by library for a vector of genes
# (this is an extension of the example written by Groves in subset_transcript.R)
prep_transcript_depth_by_library <-
  function(gvec,
           gene_df,
           exon_df,
           depth_df,
           position_shift_options = c('+' = -1, '-' = 1),
           PSEUDOCOUNT = 1) {
    map_dfr(gvec, ~ {
      # pull transcript variables
      gchr <- gene_df[.x, 'chrom']
      gstart <- gene_df[.x, 'start'] - 1
      gstop <- gene_df[.x, 'stop']
      glength <- gene_df[.x, 'summed_exon_length']
      gstrand <- gene_df[.x, 'strand']
      
      select_exon_df <- exon_df %>%
        dplyr::filter(gene == .x)
      
      # check if we have any reads in the gene
      idx <- rownames(gene_df) == .x
      
      # only continue if there are any reads, otherwise return NULL for this gene
      if (any(idx)) {
        transcript_vector <- gene_df[idx,]
        
        # subset for transcript
        tmp <- depth_df %>%
          subset_trace_for_gene(gchr, gstart, gstop)
        
        if (nrow(tmp) > 0) {
          tmp %>%
            group_split(library, group, library_type, sq) %>%
            # add zeros
            map_dfr(~ add_zero_positions(.x, select_exon_df, keep = TRUE)) %>%
            # add psuedocount column
            add_pseudocount(PSEUDOCOUNT) %>%
            shift_position_by_lib_type(gstrand, position_shift_options) %>%
            group_by(library, group, library_type, sq, gene) %>%
            nest() %>%
            ungroup()
        } else {
          # ignore exons right now where there aren't any reads
          NULL
        }
        
      } else {
        NULL
      }
    })
  }


#function to shift position of 5p libraries by 1 to match 3p libraries
shift_position_by_lib_type = function(depth_df,
                                      gstrand,
                                      position_shift_options = c('+' = -1, '-' = 1)) {
  position_shift = position_shift_options[gstrand]
  depth_df %>%
    mutate(pos = if_else(library_type == '5p',
                         pos + position_shift,
                         pos))
}


#function to combine paired libraries in transcript_depth_df based on pair information in pairs_df
combine_pair = function(pair_id, group, p3_library, p5_library){
  transcript_depth_df %>% 
    dplyr::select(library, library_type, pos, depth) %>% 
    mutate(library = factor(library, levels = c(p3_library, p5_library))) %>% 
    filter(library %in% c(p3_library, p5_library)) %>% 
    tidyr::complete(library_type, pos, fill = list(depth = 0)) %>% 
    pivot_wider(id_cols = c('pos'),
                names_from = library_type,
                values_from = depth,
                values_fill = NA,
                names_prefix = "depth_") %>% 
    arrange(pos) %>% 
    #filter out the zero and final position created by subtracting from 5' position to make them align
    filter(pos > 0,
           ! ( is.na(depth_5p) & (pos == max(pos)) ) ) %>% 
    mutate(pair_id = pair_id,
           group = group,
           p3_library = p3_library,
           p5_library = p5_library)
}


# VISUALIZE DATA ----------------------------------------------------------

plot_group_scatter_matrix = function(lib_scale_depth_df,
                                     grp,
                                     stat_to_plot = 'depth',
                                     log_transform = TRUE) {
  if (log_transform) {
    plot_stat = paste0('log_', stat_to_plot)
    plot_df = lib_scale_depth_df %>%
      mutate(plot_stat = log(lib_scale_depth_df[[stat_to_plot]], 2)) %>%
      filter(lib_scale_depth_df[[stat_to_plot]] > 0)
  } else {
    plot_df = lib_scale_depth_df
    plot_stat = stat_to_plot
  }
  plot_df %>%
    filter(group == grp) %>%
    dplyr::select(c('pos', 'library', plot_stat)) %>%
    pivot_wider(names_from = 'library',
                values_from = plot_stat) %>%
    select(-pos) %>%
    ggpairs(title = paste(grp, plot_stat, sep = ' - '))
}

plot_all_group_scatter_matrix = function(lib_scale_depth_df,
                                         stat_to_plot = 'depth',
                                         log_transform = TRUE) {
  groups_to_plot = lib_scale_depth_df %>%
    group_by(group) %>%
    summarize(N = length(unique(library))) %>%
    filter(N > 1) %>%
    pull(group)
  mat_list = list()
  for (grp in groups_to_plot) {
    mat_list[[grp]] = plot_group_scatter_matrix(lib_scale_depth_df,
                                                grp,
                                                stat_to_plot,
                                                log_transform = log_transform) %>%
      ggmatrix_gtable()
  }
  plot_grid(plotlist = mat_list)
}


# COMPUTE STATISTICS ------------------------------------------------------

get_scaled_stats = function(dat) {
  dat %>%
    mutate(
      per_total = depth / sum(depth),
      zscore = (depth - mean(depth)) / sd(depth),
      medians = pseudo_depth / median(pseudo_depth)
    )
}

get_percentile = function(dat) {
  dat %>%
    mutate(rank = rank(depth),
           percentile = rank / max(rank) * 100) %>%
    dplyr::select(-rank)
}

get_ppois_pvalue = function(dat) {
  dat %>%
    mutate(pois_pvalue = ppois(depth, lambda = mean(depth), lower = FALSE))
}

get_exp_pvalue = function(dat) {
  dat %>%
    mutate(exp_pvalue = pexp(depth, rate = 1 / mean(depth), lower = FALSE))
}

run_grouped_stats = function(stat_df) {
  stat_df %>%
    get_percentile() %>%
    get_ppois_pvalue() %>%
    get_exp_pvalue() %>%
    ungroup()
}

# SEQUENCE CONTEXT --------------------------------------------------------


old_pull_seq_context = function(sq_vector, pos, to_left, to_right) {
  tot_len = length(sq_vector)
  l = pos - to_left
  l[l < 1] <- 1
  r = pos + to_right
  r[r > tot_len] <- tot_len
  seq_slice = paste(sq_vector[l:r], collapse = '')
  return(seq_slice)
}

pull_seq_context = function(left, right) {
  paste(sq_vector[left:right], collapse = '')
}

add_seq_context_to_df = function(stat_df,
                                 sq_vector,
                                 to_left = 0,
                                 to_right = 1) {
  min_pos = 1
  max_pos = as.numeric(length(sq_vector))
  slice_coords = stat_df %>%
    mutate(
      left = pos - to_left,
      right = pos + to_right,
      left = if_else(left < min_pos,
                     min_pos,
                     left),
      right = if_else(right > max_pos,
                      max_pos,
                      right)
    ) %>%
    dplyr::select(left, right)
  #add sequence context windows
  seq_context = pmap(slice_coords, function(left, right)
    pull_seq_context(left, right)) %>%
    unlist()
  stat_df %>%
    mutate(
      seq_context = seq_context,
      dinucleotide = str_sub(seq_context, start = to_left + 1, end = to_left + 2)
    )
}

add_bp_prop = function(stat_df, pair_prob_df) {
  stat_df %>%
    left_join(pair_prob_df, by = 'pos')
}


# DOUBLE-CHECKING FUNCTIONS -----------------------------------------------
#quick re-coding to check that computations are correct
check_single_end_lib_data = function(my_lib, lib_depth_df) {
  check_results = c()
  for (my_lib in libs) {
    sub_lib_depth_df = lib_depth_df %>%
      filter(library == my_lib)
    sub_lib_depth_df
    #manually compute stats
    mn = mean(sub_lib_depth_df$depth)
    sd = sd(sub_lib_depth_df$depth)
    med = median(sub_lib_depth_df$pseudo_depth)
    tot = sum(sub_lib_depth_df$depth)
    #build reference for checking by subsetting output dataset
    checkref0 = sub_lib_depth_df %>%
      slice_sample(n = 5)
    checkref = sub_lib_depth_df %>%
      filter(depth > 0) %>%
      slice_sample(n = 5) %>%
      rbind(checkref0) %>%
      dplyr::relocate(library_type) %>%
      select(depth:last_col())
    #build manual check df
    check = checkref %>%
      select(depth, pseudo_depth) %>%
      mutate(
        per_total = depth / tot,
        zscore = (depth - mn) / sd,
        medians = pseudo_depth / med
      )
    #confirm they match
    res = all.equal(checkref, check)
    check_results = append(check_results, res)
  }
  lib_check = data.frame('library' = libs,
                         'lib_data_check' = check_results)
  print(lib_check)
  return(lib_check)
}

check_library_group_stats = function(lib_depth_df, lib_group_stat_df) {
  check_results = c()
  dist_check_results = c()
  ugroups = unique(lib_group_stat_df$group)
  for (my_group in ugroups) {
    sub_lib_depth_df = lib_depth_df %>%
      filter(group == my_group)
    
    #build reference for checking by subsetting output dataset
    checkref0 = lib_group_stat_df %>%
      filter(group == my_group) %>%
      slice_sample(n = 5)
    checkref = lib_group_stat_df %>%
      filter(group == my_group,
             depth > 0) %>%
      slice_sample(n = 5) %>%
      rbind(checkref0) %>%
      arrange(pos) %>%
      distinct(pos, .keep_all = TRUE)
    
    #build manual check df for the summary columns
    check = sub_lib_depth_df %>%
      semi_join(checkref, by = c('group', 'pos')) %>%
      arrange(pos)
    
    #check single position
    single_pos = sample(checkref$pos, 1)
    check %>%
      filter(pos == single_pos)
    single_pos_check_in = check %>%
      filter(pos == single_pos)
    single_pos_check = tibble('pos' = single_pos,
                              'depth' = sum(single_pos_check_in$depth),
                              'per_total' = mean(single_pos_check_in$per_total),
                              'zscore' = mean(single_pos_check_in$zscore),
                              'medians' = mean(single_pos_check_in$medians))
    single_pos_checkref = checkref %>% 
      filter(pos == single_pos) %>% 
      dplyr::select(colnames(single_pos_check))
    res = all.equal(single_pos_check, single_pos_checkref)
    check_results = append(check_results, res)
    
    #get input for checking the distributions
    dist_check = sub_lib_depth_df %>%
      group_by(group, pos) %>%
      summarize(depth = sum(depth),
                .groups = 'drop') %>%
      mutate(
        rank = rank(depth),
        percentile = rank / max(rank) * 100,
        pois_pvalue = ppois(depth, lambda = mean(depth), lower = FALSE),
        exp_pvalue = pexp(depth, rate = 1 / mean(depth), lower = FALSE)
      ) %>%
      dplyr::select(-rank) %>%
      semi_join(checkref, by = c('group', 'pos')) %>%
      arrange(pos) %>%
      ungroup()
    dist_check_ref = checkref %>%
      dplyr::select(colnames(dist_check)) %>%
      arrange(pos) %>%
      ungroup()
    dist_res = all.equal(dist_check, dist_check_ref)
    dist_check_results = append(dist_check_results, dist_res)
  }
  check_result_df = data.frame(
    'group' = ugroups,
    'check_summary_cols' = check_results,
    'check_dist_cols' = dist_check_results
  )
  return(check_result_df)
}


check_pair_data = function(transcript_depth_df, paired_depth_df, pairs_df){
  res_list = c()
  stat_res_list = c()
  for (i in 1:nrow(pairs_df)){
    p3_lib = pairs_df$p3_library[i]
    p5_lib = pairs_df$p5_library[i]
    p3_df = transcript_depth_df %>% 
      filter(library == p3_lib) %>% 
      dplyr::rename(p3_library = library)
    p5_df =  transcript_depth_df %>% 
      filter(library == p5_lib) %>% 
      dplyr::rename(p5_library = library)
    check = p3_df %>% 
      inner_join(p5_df, by = c('pos')) %>% 
      mutate(depth = depth.x * depth.y) %>% 
      mutate(pseudo_depth = depth + PSEUDOCOUNT,
             per_total = depth / sum(depth),
             zscore = (depth - mean(depth)) / sd(depth),
             medians = pseudo_depth / median(pseudo_depth)) %>% 
      dplyr::select(pos, depth, pseudo_depth, per_total, zscore, medians)
    check_ref = paired_depth_df %>% 
      filter(p3_library == p3_lib,
             p5_library == p5_lib) %>% 
      dplyr::select(colnames(check))
    res = all.equal(check, check_ref)
    res_list = append(res_list, res)
  }
  res_df = pairs_df
  res_df$match = res_list
  return(res_df)
}

check_pair_stats = function(paired_depth_df, paired_group_stat_df){
  pair_groups = unique(paired_depth_df$group)
  res_list = c()
  for (my_group in pair_groups){
    check_paired_depth_df = paired_depth_df %>% 
      filter(group == my_group) %>% 
      group_by(pos) %>% 
      summarize(depth = sum(depth),
                per_total = mean(per_total),
                zscore = mean(zscore),
                medians = mean(medians)) %>% 
      mutate(group = my_group) %>% 
      dplyr::relocate(group) %>% 
      mutate(rank = rank(depth),
             percentile = rank/max(rank)*100) %>% 
      dplyr::select(-rank) %>% 
      mutate(pois_pvalue = ppois(depth, lambda=mean(depth), lower=FALSE)) %>% 
      mutate(exp_pvalue = pexp(depth, rate=1/mean(depth), lower=FALSE))
    ref_check_paired_depth_df = paired_group_stat_df %>% 
      filter(group == my_group) %>% 
      dplyr::select(colnames(check_paired_depth_df))
    res = all.equal(check_paired_depth_df, ref_check_paired_depth_df)
    res_list = append(res_list, res)
  }
  res_df = data.frame('group' = pair_groups,
                      'match' = res_list)
  return(res_df)
}



# REPORT FUNCTIONS --------------------------------------------------------



# DATA MANIPULATION FOR PLOTTING ------------------------------------------

trim_for_plot = function(stat_df, trim_left, trim_right){
  stat_df %>%
    group_by(transcript) %>%
    filter(Position > trim_left,
           Position < (max(Position) - trim_right)) %>%
    ungroup()
}


correct_pvalues = function(stat_df, false_discovery_method){
  correct_pvalue = function(x){
    p.adjust(x, method = false_discovery_method)
  }
  if (false_discovery_method != FALSE){
    res = stat_df %>% 
      mutate_at(str_subset(colnames(stat_df), '_pvalue'), correct_pvalue)
  }
  else{
    res = stat_df
  }
  return(res)
}

add_significance_calls = function(stat_df, sig_param_df){
  #build dataframe of significance calls for each position in stat_df
  sig_df = tibble(pos = stat_df$pos)
  for (stat in rownames(sig_param_df)){
    stat_min = sig_param_df[stat, 'min']
    stat_max =  sig_param_df[stat, 'max']
    sig_df[[stat]] = stat_df[[stat]] > stat_min & stat_df[[stat]] < stat_max
  }
  sig_df = sig_df %>% 
    dplyr::select(-pos)
  #apply accross each to get significance for all included stats in min_vals and max_vals
  sig_call = apply(sig_df, 1, function(x) all(x, na.rm = TRUE))
  #output with significance calls
  stat_df %>% 
    mutate(sig = sig_call)
}

factorize_significance = function(stat_df){
  stat_df %>% 
    mutate(Significant = factor(Significant, levels = c(TRUE, FALSE)))
}

#make column adjustments for plotting
prep_for_plotting = function(dat, var_label_df){
  #swap out for nice column names for plotting and filtering
  sub_var_label_df = var_label_df %>% 
    filter(var %in% colnames(dat))
  #add hover variables
  dat %>% 
    rename_at(vars(sub_var_label_df$var), ~ sub_var_label_df$label) %>% 
    mutate(context = paste0('Position: ', Position, '\n',
                            'Dinucleotide: ', dinucleotide, '\n',
                            'Motif: ', seq_context),
           stat_labs = paste0('Raw count: ', `Count`, '\n',
                              'Proportion: ', formatC(`Proportion`, format = "e", digits = 0), '\n',
                              'Percentile: ', round(Percentile, digits = 1), '\n',
                              'Z-score: ', round(`Zscore`, digits=1), '\n',
                              'Times median: ', round(`Medians`, digits=1), '\n',
                              'Poisson p-value: ', formatC(`pois_pvalue`, format = "e", digits = 0), '\n',
                              'Exponential p-value: ', formatC(`Pvalue`, format = "e", digits = 0)))
}


# PLOT SINGLE GROUP TRACE -------------------------------------------------


#function to plot ggplot scatterplot
gg_scatterplot = function(stat_df, ycol, facet_by_group = TRUE, yaxis_scale = 'log2'){
  xcol = 'Position'
  plt = stat_df %>% 
    mutate(Significant = factor(Significant, levels = c(TRUE, FALSE))) %>% 
    ggplot(aes(x = .data[[xcol]], y = .data[[ycol]],
               color = Significant, fill=Significant, alpha=Significant)) +
    geom_point(pch = 21, color = 'black') +
    scale_fill_manual(values = c('red', 'grey')) + 
    scale_color_manual(values = c('black', 'grey')) +
    scale_alpha_manual(values = c(1, 0.1)) +
    scale_y_continuous(trans = yaxis_scale)
  if (facet_by_group){
    plt = plt + 
      facet_wrap(~group, nrow = length(unique(stat_df$group)))
  }
  return(plt)
}

add_axis_buffer = function(lim, left_buffer_prop, right_buffer_prop){
  res = lim
  left_buffer_add = round((lim[2] - lim[1]) * left_buffer_prop, digits=0)
  right_buffer_add = round((lim[2] - lim[1]) * right_buffer_prop, digits=0)
  res[1] = lim[1]-left_buffer_add
  res[2] = lim[2]+right_buffer_add
  return(res)
}

#pull out library strings from either metadata_df or pairs_df based on library_type and group
lib_strs_by_group_libtype = function(lib_type, grp, pairs_df, metadata_df){
  if (lib_type == 'paired'){
    libs = pairs_df %>% 
      filter(group == grp) %>% 
      mutate(pair_string = str_c(p3_library, p5_library, sep= ' - ')) %>% 
      pull(pair_string) %>% 
      str_c(collapse = "\n")
  } else{
    libs = metadata_df %>%
      filter(group == grp,
             library_type == lib_type) %>%
      pull(library) %>%
      str_c(collapse = "\n")
  }
  return(libs)
}

plotly_single_scatter = function(stat_df, ycol, xlim = c(FALSE, FALSE),log2Y = FALSE, depth_min_rm = -1){
  #handle Y axis scale
  ytitle = ycol
  if (log2Y){
    stat_df[[ycol]] = log(stat_df[[ycol]], 2)
    ytitle = paste('log2', ycol)
  }
  #handle minimum depth
  stat_df = stat_df %>% 
    filter(Count > depth_min_rm)
  
  #prepare the hover labels and factorize Significance
  plt_df = stat_df %>% 
    factorize_significance() %>% 
    mutate(hover_text = paste(stat_labs, context, sep='\n'))
  
  #handle setting x axis
  if (xlim[1] == FALSE){
    xlim = c(min(plt_df$Position, na.rm = TRUE),
             max(plt_df$Position, na.rm = TRUE))
  }
  xlim_list = as.list(xlim)
  
  #build the plot
  plt_df %>% 
    plot_ly(x = ~Position,
            y = as.formula(paste0('~', ycol)),
            color = ~Significant,
            colors = c("red","grey50"),
            text = ~I(hover_text),
            type = 'scatter',
            mode = 'markers',
            marker = list(size = 10,
                          line = list(color = 'black',
                                      width = 1)),
            showlegend=FALSE,
            hoverlabel = list(align = 'left'),
            hovertemplate = "%{text}<extra></extra>") %>% 
    layout(yaxis = list(title = ytitle),
           xaxis = list(range = xlim_list))
}


# MEAN BY DINUCLEOTIDE FUNCTIONS --------------------------------------------

get_dnt_err_bar_df = function(stat_df, ycol){
  plt_df = stat_df %>% 
    group_by(dinucleotide) %>% 
    summarize(mn = mean(.data[[ycol]]),
              stdv = sd(.data[[ycol]]),
              se = std.error(.data[[ycol]])) %>% 
    filter(str_length(dinucleotide) == 2) %>% 
    arrange(desc(mn)) 
  order = plt_df$dinucleotide
  plt_df$dinucleotide = factor(plt_df$dinucleotide, levels = order)
  levels(plt_df$dinucleotide)
  return(plt_df)
}

plot_dnt_errorbars = function(dnt_plt_df, title = ''){
  dnt_plt_df %>% 
    ggplot(aes(x = dinucleotide,
               y = mn)) +
    geom_point() +
    geom_errorbar(aes(ymin = mn-se,
                      ymax = mn+se), width = 0.2) +
    labs(title = title)
}

build_dnt_errorbar_plot = function(stat_df, ycol){
  n_non_zero = stat_df %>% 
    filter(Count > 0) %>% 
    nrow()
  title = paste0('N = ', n_non_zero, ' non-zero positions')
  get_dnt_err_bar_df(stat_df, ycol) %>% 
    plot_dnt_errorbars(title = title) +
    labs(y = paste('mean', ycol))
}

# SEQ LOGO FUNCTIONS ------------------------------------------------------

filter_for_same_seq_len = function(stat_df, seq_column = 'seq_context'){
  seq_lengths = str_length(stat_df[[seq_column]])
  stat_df %>% 
    mutate(seq_length = seq_lengths) %>% 
    filter(seq_length == max(seq_length)) 
}

seqs_by_group_libtype = function(in_df, s_group, s_library_type, seq_column = 'seq_context'){
  in_df %>% 
    filter(group == s_group,
           library_type == s_library_type) %>% 
    filter_for_same_seq_len() %>% 
    pull(seq_column)
}

plot_seq_logo = function(seqs, title_color = 'black', title_prefix = ''){
  title = paste0(title_prefix, 'N = ', length(seqs), ' peaks')
  ggplot() + geom_logo( seqs,  method='p' ) + 
    theme_logo() + 
    labs(title = title) +
    theme(plot.title = element_text(colour = title_color))
}


plot_group_seq_logos = function(stat_df, seq_column = 'seq_context'){
  seq_lengths = str_length(stat_df[[seq_column]])
  stat_df = stat_df %>% 
    mutate(seq_length = seq_lengths) %>% 
    filter(seq_length == max(seq_length))
  plotlist = list()
  # all_seqs = stat_df[[seq_column]]
  # all_lens = str_length(all_seqs)
  # all_seqs = all_seqs[all_lens == max(all_lens)]
  # all_pos_lab = 'all positions'
  # plotlist[[all_pos_lab]] = plot_seq_logo(all_seqs, all_pos_lab)
  peak_df = stat_df %>% 
    filter(Significant == TRUE)
  for (grp in unique(peak_df$group)){
    seqs = peak_df %>% 
      filter(group == grp) %>% 
      pull(seq_column)
    plotlist[[grp]] = plot_seq_logo(seqs, grp)
  }
  plot_grid(plotlist = plotlist, nrow = length(plotlist))
}

#ENRICHMENT TEST
get_p = function(x){
  x$p.value
}
get_odds = function(x){
  x$estimate
}
sequence_fisher_test = function(stat_df, seq_column){
  seqs = stat_df[[seq_column]]
  useqs = unique(seqs)
  fisher_y = stat_df$sig
  test_list = list()
  for (s in useqs){
    fisher_x = factor(seqs == s, levels = c(TRUE, FALSE))
    test_list[[s]] = fisher.test(x = fisher_y, y = fisher_x)
  }
  pvals = unlist(map(test_list, get_p))
  odds = unlist(map(test_list, get_odds))
  res_df0 = tibble(seq_column = useqs,
                   'odds ratio' = odds,
                   'p-value' = pvals) %>% 
    set_names(c(seq_column, 'odds ratio', 'p-value'))
  #get counts
  res_df = stat_df %>% 
    group_by_at(c(seq_column, 'sig'), .drop = FALSE) %>% 
    summarize(N = n(), .groups = 'drop') %>% 
    mutate(sig = as.logical(sig),
           sig = if_else(sig, 'significant', 'not significant')) %>% 
    pivot_wider(id_cols = c('dinucleotide'),
                names_from = sig,
                values_from = N) %>% 
    left_join(res_df0, by = seq_column) %>% 
    arrange(desc(`odds ratio`), seq_column)
  return(res_df)
}



# NEGATIVE BINOMIAL FUNCTIONS ---------------------------------------------

#pull model coefficients from glm.nb() into a tibble
mod_to_tibble = function(mod){
  mod_res = coef(summary(mod)) %>% 
    data.frame() %>% 
    set_names('estimate', 'stderr', 'zvalue', 'pvalue') %>% 
    rownames_to_column('covariate') %>% 
    mutate(padj = p.adjust(pvalue, method = 'fdr')) %>% 
    tibble()
  return(mod_res)
}


get_intercept = function(mod_df){
  mod_df %>% 
    filter(covariate=='(Intercept)') %>% 
    pull(estimate) %>% 
    round(digits = 2)
}

get_mod_plot_df = function(mod_df, seq_column){
  mod_df %>% 
    filter(str_detect(covariate, seq_column)) %>% 
    arrange(estimate) %>% 
    mutate(covariate = str_remove(covariate, seq_column),
           covariate = factor(covariate, levels = unique(covariate)))
}
plot_mod_estimates = function(mod_plt_df, intercept, seq_column){
  mod_plt_colors = inferno(length(mod_plt_df$covariate))
  mod_plt_df %>% 
    ggplot(aes(x = covariate, y = estimate, fill = covariate)) +
    geom_bar(stat = 'identity', color = 'black') +
    scale_fill_manual(values = mod_plt_colors) +
    coord_flip() +
    labs(title = paste('intercept =', intercept),
         x = seq_column) +
    theme(legend.position = 'none')
}

run_and_plot_nb_mod = function(stat_df, ycol, seq_column){
  mod_formula = paste(ycol, seq_column, sep = ' ~ ')
  m_nb <- glm.nb(mod_formula, data = stat_df)
  mod_df = mod_to_tibble(m_nb)
  intercept = mod_df %>% 
    get_intercept()
  mod_df %>% 
    get_mod_plot_df(seq_column) %>% 
    plot_mod_estimates(intercept, seq_column)
}


# GROUP-LEVEL REPORT --------------------------------------------

format_sig_params = function(sig_param_df){
  sig_param_df %>% 
    rownames_to_column('Statistic') %>% 
    mutate(min = na_if(abs(min), Inf),
           max = na_if(abs(max), Inf),
           min = as.character(min),
           max = as.character(max)) %>% 
    replace_na(list(min = '--',
                    max = '--')) %>% 
    dplyr::rename(Minimum = min,
                  Maximum = max)
}

build_group_report_obj = function(s_transcript){
  # print(paste0('adding tags for transcript=', s_transcript, '...'))
  #set up report title, the html taglist, data, and x limits for this transcript
  report_title = paste(s_transcript, '- group report')
  tag_list = list()
  st_df = stat_df %>%
    filter(transcript == s_transcript)
  st_long_sq_df = long_sq_df %>% 
    filter(sq == s_transcript)
  st_ugroups = st_long_sq_df %>% 
    pull(group) %>% 
    unique()
  st_ulib_types = st_df %>% 
    pull(library_type) %>% 
    unique()
  xlim = c(min(st_df$Position, na.rm = TRUE), max(st_df$Position, na.rm = TRUE)) #set x limits for full set to match axis scales across plots
  xlim = add_axis_buffer(xlim, 0.01, 0.01)
  #for each group
  for (grp in st_ugroups){
    # print(paste0('adding tags for transcript=', s_transcript, ', group=', grp,'...'))
    #subset for this group
    st_sg_df = st_df %>% 
      filter(group == grp)
    seq_logo_list = list() #new seq logo list to be populated with 3p and 5p data for this group
    
    #set up the group header, subset the data, and set up list for logo plots
    group_header_id = paste(s_transcript, grp, 'header')
    tag_list[[group_header_id]] = h3(grp) #level 3 header
    
    #for each library type within group
    for (lib_type in st_ulib_types){
      # print(paste0('adding tags for transcript=', s_transcript, ', group=', grp,', lib=', lib_type, '...'))
      #set up plot title
      plot_title = paste(s_transcript, grp, lib_type, 'peaks:', sep='; ')
      plt_id = paste(s_transcript, grp, lib_type, 'plot')
      tag_list[[plot_title]] = h4(plot_title) #level 4 title for scatterplot
      #subset data for this trascript/group/lib_type
      st_sg_sl_df = st_sg_df %>% 
        filter(library_type == lib_type)
      
      #add line for the libraries included in this group/lib_type
      lib_strings = lib_strs_by_group_libtype(lib_type, grp, pairs_df, metadata_df)
      libs_id = paste(s_transcript, grp, lib_type, 'libs')
      tag_list[[libs_id]] = p(lib_strings)
      
      #check this subset has sufficient data for plotting
      n_nonzero = st_sg_sl_df %>% 
        filter(Count > depth_min_rm) %>% 
        nrow()
      has_data = n_nonzero > 1
      
      #build trace plot or save as a no data string
      if (has_data){
        #build plot
        plt = st_sg_sl_df %>%
          plotly_single_scatter(ycol = ycol, xlim = xlim, log2Y = log2Y, depth_min_rm = depth_min_rm)
        tag_list[[plt_id]] = plt
        #build the seq logo figure if there are significant peaks
        sig_seqs = st_sg_sl_df %>%
          filter(Significant == TRUE) %>% 
          filter_for_same_seq_len() %>%
          pull(seq_column)
        if (length(sig_seqs)>0){
          logo_plt = st_sg_sl_df %>%
            filter(Significant == TRUE) %>% 
            filter_for_same_seq_len() %>%
            pull(seq_column) %>%
            plot_seq_logo() %>%
            ggplotly()
        }
        #if there were no significant peaks for logo
        else{
          logo_plt = p('no significant peaks')
        }
        seq_logo_list[[lib_type]] = logo_plt
      }
      #if there were no data for
      else {
        no_obs_message = p('not enough observations')
        tag_list[[plt_id]] = no_obs_message
        seq_logo_list[[lib_type]] = no_obs_message
      }
    }
    #add the two logo figures to the tag_list
    for (lib_type in st_ulib_types){
      #add title
      logo_title = paste(s_transcript, grp, lib_type, 'logo:', sep='; ')
      tag_list[[logo_title]] = h4(logo_title)
      #add logo figure
      logo_id = paste(s_transcript, grp, lib_type, 'logo')
      tag_list[[logo_id]] = seq_logo_list[[lib_type]]
    }
  }
  report_outfile = paste0(s_transcript, '.html')
  report_title = paste0(s_transcript, '- group report')
  report_list = list(input = individual_markdown_script,
                     output_file = report_outfile,
                     output_dir = groups_output_dir,
                     report_title = report_title,
                     tag_list = tag_list)
  object_path = paste0(groups_output_dir, '/', s_transcript, '_rep.Rdata')
  object_path
  report_list
  save(report_list, file = object_path)
  return(object_path)
}


render_individual_report = function(r){
  #get required packages for rendering
  require(rmarkdown)
  require(kableExtra)
  require(htmltools)
  #load report_list and re-assign variables for rendering
  load(r)
  tag_list = report_list$tag_list
  report_title = report_list$report_title
  #set up temp dir for intermediates
  tf <- tempfile()
  dir.create(tf)
  #render
  res = rmarkdown::render(input = report_list$input,
                          output_file = report_list$output_file,
                          output_dir = report_list$output_dir,
                          intermediates_dir=tf,
                          quiet = TRUE)
  unlink(tf)
  return(res)
}



build_all_group_report = function(s_transcript = 'All'){
  tag_list = list()
  ugroups = unique(stat_df$group)
  ulib_types = unique(stat_df$library_type)
  for (grp in ugroups){
    # print(paste0('adding tags for transcript=', s_transcript, ', group=', grp,'...'))
    #set up the group header, subset the data, and set up list for logo plots
    group_header_id = paste(s_transcript, grp, 'header')
    tag_list[[group_header_id]] = h3(grp) #level 3 header
    st_sg_df = stat_df %>% 
      filter(group == grp)
    seq_logo_list = list() #new seq logo list to be populated with 3p and 5p data for this group
    #for each library type within group
    for (lib_type in ulib_types){
      # print(paste0('adding tags for transcript=', s_transcript, ', group=', grp,', lib=', lib_type, '...'))
      #set up plot title
      plot_title = paste(s_transcript, grp, lib_type, 'Dinucleotide means:', sep='; ')
      plt_id = paste(s_transcript, grp, lib_type, 'plot')
      tag_list[[plot_title]] = h4(plot_title) #level 4 title for error bars
      #subset data for this trascript/group/lib_type
      st_sg_sl_df = st_sg_df %>% 
        filter(library_type == lib_type)
      n_nonzero = st_sg_sl_df %>% 
        filter(Count > depth_min_rm) %>% 
        nrow()
      has_data = n_nonzero > 0
      #build trace plot or save as a no data string
      if (has_data){
        # #build plot
        plt = st_sg_sl_df %>%
          build_dnt_errorbar_plot(error_bar_ycol) %>% 
          ggplotly()
        tag_list[[plt_id]] = plt
        #build the seq logo figure if there are significant peaks
        sig_seqs = st_sg_sl_df %>%
          filter(Significant == TRUE) %>% 
          filter_for_same_seq_len() %>%
          pull(seq_column)
        if (length(sig_seqs)>0){
          logo_plt = st_sg_sl_df %>%
            filter(Significant == TRUE) %>% 
            filter_for_same_seq_len() %>%
            pull(seq_column) %>%
            plot_seq_logo() %>%
            ggplotly()
        }
        #if there were no significant peaks for logo
        else{
          logo_plt = p('no significant peaks')
        }
        seq_logo_list[[lib_type]] = logo_plt
      }
      #if there were no data for
      else {
        no_obs_message = p('not enough non-zero positions')
        tag_list[[plt_id]] = no_obs_message
        seq_logo_list[[lib_type]] = no_obs_message
      }
    }
    #add the two logo figures to the tag_list
    for (lib_type in ulib_types){
      #add title
      logo_title = paste(s_transcript, grp, lib_type, 'logo:', sep='; ')
      tag_list[[logo_title]] = h4(logo_title)
      #add logo figure
      logo_id = paste(s_transcript, grp, lib_type, 'logo')
      tag_list[[logo_id]] = seq_logo_list[[lib_type]]
    }
  }
  #render the html
  report_outfile = paste0(s_transcript, '.html')
  res = rmarkdown::render(input = individual_markdown_script, 
                    output_format = "html_document",
                    output_file = report_outfile,
                    output_dir = groups_output_dir)
  return(res)
}

#function to copy group report to s3
report_to_s3 = function(local_file_path, s3_bucket, s3_prefix){
  attempts = 0
  check = FALSE
  while (check != TRUE & attempts < s3_attempts){
    attempts = attempts + 1
    # print(str_c('attempt ', attempts, '...'))
    check = try(aws.s3::put_object(local_file_path,
                                   object = paste0(s3_prefix, basename(local_file_path)),
                                   bucket = s3_bucket,
                                   multipart = FALSE))
    Sys.sleep(s3_sleep_seconds)
  }
  return(check)
}

# FUNCTIONS FOR COMPARING GROUPS ------------------------------------------

#
add_sig_agree_col = function(widened_df){
  widened_df %>% 
    mutate(Overlap = 'none',
           Overlap = if_else(Significant.x & !Significant.y,
                             'X only',
                             Overlap),
           Overlap = if_else(!Significant.x & Significant.y,
                             'Y only',
                             Overlap),
           Overlap = if_else(Significant.x & Significant.y,
                             'Both',
                             Overlap),
           Overlap = factor(Overlap, levels = c('Both', 'X only', 'Y only', 'none')))
}

#
widen_for_pair = function(stat_df,
                          stat_to_plot,
                          x_group,
                          y_group,
                          x_library_type,
                          y_library_type){
  xdf = stat_df %>% 
    filter(group == x_group,
           library_type == x_library_type)
  ydf = stat_df %>% 
    filter(group == y_group,
           library_type == y_library_type)
  plt_df = xdf %>% 
    left_join(ydf, by = c('transcript', 'Position')) %>% 
    add_sig_agree_col()
  return(plt_df)
}


#
gg_compare_group_stat = function(stat_df,
                                 stat_to_plot,
                                 x_group,
                                 y_group,
                                 x_library_type,
                                 y_library_type,
                                 log2_axes = FALSE,
                                 facet_by_transcript = TRUE){
  plt_df = stat_df %>% 
    filter(Position < max(Position)) %>% #remove final position which doesn't appear in 5' because of subtracting 1
    widen_for_pair(stat_to_plot,
                   x_group,
                   y_group,
                   x_library_type,
                   y_library_type)
  xcol = paste(stat_to_plot, 'x', sep='.')
  ycol = paste(stat_to_plot, 'y', sep='.')
  corrs = plt_df %>% 
    group_by(transcript) %>% 
    dplyr::select(all_of(c('transcript', 'Position', xcol, ycol))) %>% 
    # na.omit() %>% 
    summarize(corr = cor(.data[[xcol]],
                         .data[[ycol]]))
  mn_corr = round(mean(corrs$corr), digits = 2)
  # std_err = round(std.error(corrs$corr), digits = 2)
  subtitle = str_c('mean cor = ', mn_corr)
  plt = plt_df %>% 
    ggplot(aes(x = .data[[xcol]],
               y = .data[[ycol]],
               color = Overlap)) +
    geom_point() +
    scale_color_manual(values = c('firebrick', 'dodgerblue', 'orchid', 'grey20')) +
    labs(x = paste0(x_group, ' (', x_library_type, ')'),
         y = paste0(y_group, ' (', y_library_type, ')'),
         title = stat_to_plot,
         subtitle = subtitle)
  if (log2_axes){
    plt = plt + 
      scale_x_continuous(trans='log2') +
      scale_y_continuous(trans='log2')
  }
  if (facet_by_transcript){
    plt = plt + facet_wrap(~transcript, scales = 'free')
  }
  return(plt)
}


overlap_colors = c('firebrick', 'dodgerblue', 'orchid', 'grey20')
names(overlap_colors) = c('Both', 'X only', 'Y only', 'none')
plotly_compare_group_stat = function(wide_pair_df,
                                     stat_to_plot,
                                     x_group,
                                     y_group,
                                     x_library_type,
                                     y_library_type,
                                     log2_axes = FALSE,
                                     overlap_colors = overlap_colors){
  xcol = paste(stat_to_plot, 'x', sep='.')
  ycol = paste(stat_to_plot, 'y', sep='.')
  
  title_text = stat_to_plot
  if(log2_axes){
    wide_pair_df[[xcol]] = log(wide_pair_df[[xcol]], 2)
    wide_pair_df[[ycol]] = log(wide_pair_df[[ycol]], 2)
    title_text = paste('log2', stat_to_plot)
  }
  x_left = min(wide_pair_df[[xcol]], na.rm = TRUE)
  y_top = max(wide_pair_df[[ycol]], na.rm = TRUE) * 0.98
  corrs = wide_pair_df %>% 
    filter(.data[[xcol]] != -Inf,
           .data[[ycol]] != -Inf) %>% 
    group_by(transcript) %>% 
    dplyr::select(all_of(c('transcript', 'Position', xcol, ycol))) %>% 
    na.omit() %>%
    summarize(corr = cor(.data[[xcol]],
                         .data[[ycol]]))
  mn_corr = round(mean(corrs$corr), digits = 2)
  std_err = round(std.error(corrs$corr), digits = 2)
  subtitle = str_c('cor = ', mn_corr)
  
  xtitle = paste0(x_group, ' (', x_library_type, ')')
  ytitle = paste0(y_group, ' (', y_library_type, ')')
  wide_pair_df %>% 
    plot_ly(x = as.formula(paste0('~', xcol)),
            y = as.formula(paste0('~', ycol)),
            color = ~Overlap,
            colors = overlap_colors,
            text = ~I(context.x),
            type = 'scatter',
            mode = 'markers',
            marker = list(size = 10,
                          line = list(color = 'black',
                                      width = 1)),
            showlegend = TRUE,
            hoverlabel = list(align = 'left'),
            hovertemplate = "X: %{x}<br>Y: %{y}<br>%{text}<extra></extra>") %>% 
    layout(xaxis = list(title = xtitle),
           yaxis = list(title = ytitle),
           title = list(text = title_text),
           annotations = list(
             x = x_left,
             xanchor = "left",
             y = y_top,
             text = subtitle,
             showarrow = FALSE ,
             font = list(size = 14)
           ))
}

#function to get list of sequence logos based on overlapping significance between two groups
#assumed overlap_colors is named vector with names matching the levels in 'Overlap' column generated by add_sig_agree_col()
#see plotly_compare_group_stat()


build_overlap_logos = function(wide_pair_df, overlap_colors){
  logo_list = list()
  wide_pair_df = wide_pair_df %>% 
    rename(seq_context = seq_context.x) %>% 
    filter_for_same_seq_len()
  for (n in names(overlap_colors)){
    seqs =  wide_pair_df %>% 
      filter(Overlap == n) %>% 
      pull(seq_context)
    #handle cases of empty sets
    if (length(seqs)>0){
      res = seqs %>% 
        plot_seq_logo(title_color = overlap_colors[n], title_prefix = paste0(n, '; ')) %>% 
        ggplotly()
      
    } else{
      res = htmltools::p(paste(n, '(no data)'))
    }
    logo_list[[n]] = res
  }
  return(logo_list)
}



# LARGE PANEL VISUALIZATIONS ----------------------------------------------


compare_two_groups = function(lib_stat_df,
                              stat_to_plot,
                              x_group,
                              y_group,
                              x_library_type,
                              y_library_type,
                              seq_column = 'seq_context'){
  scat_plt = lib_stat_df %>% 
    plot_group_stat_by_sq(stat_to_plot,
                          x_group,
                          y_group,
                          x_library_type,
                          y_library_type)
  seq_df = lib_stat_df %>% 
    filter_for_same_seq_len()
  null_logo = seq_df %>% 
    distinct(transcript, Position, .keep_all = TRUE) %>% 
    pull(seq_column) %>% 
    plot_seq_logo('all positions')
  x_seqs = seq_df %>% 
    filter(Significant) %>% 
    seqs_by_group_libtype(x_group, x_library_type, seq_column)
  x_logo = x_seqs %>% 
    plot_seq_logo(title = paste(x_group, x_library_type, sep='; '))
  y_seqs = seq_df %>% 
    filter(Significant) %>% 
    seqs_by_group_libtype(y_group, y_library_type, seq_column)
  y_logo = y_seqs %>% 
    plot_seq_logo(title = paste(y_group, y_library_type, sep='; '))
  both_seqs = lib_stat_df %>% 
    filter_for_same_seq_len() %>% 
    widen_for_pair(stat_to_plot,
                   x_group,
                   y_group,
                   x_library_type,
                   y_library_type) %>% 
    filter(Significant.x & Significant.y) %>% 
    pull(paste(seq_column, 'x', sep='.'))
  both_logo = both_seqs %>% 
    plot_seq_logo(title = 'Both')
  logo_list = list(x_logo,
                   y_logo,
                   both_logo)
  logo_plt = plot_grid(plotlist = logo_list, nrow = length(logo_list))
  plot_grid(scat_plt, logo_plt, nrow=1,
            rel_widths = c(3,2))
  
}


# PREPARE COMPARISONS -----------------------------------------------------

pull_libtype_groups = function(metadata_df, libtype){
  metadata_df %>% 
    arrange(library) %>% 
    filter(library_type == libtype) %>% 
    pull(group)
}

get_compare_df_by_libtype = function(metadata_df,
                                     xlibtype,
                                     ylibtype){
  xgroups = pull_libtype_groups(metadata_df, xlibtype)
  ygroups = pull_libtype_groups(metadata_df, ylibtype)
  
  #get sq assignments for each group (to only compare groups dosed with same sqs)
  long_sq_group_df = metadata_df %>% 
    select(group, sq) %>% 
    separate_rows(sq, sep=';')
  
  sub_df_list = list()
  for (xg in xgroups){
    #select sqs dosed for this X group
    xg_sqs = long_sq_group_df %>% 
      filter(group == xg) %>% 
      pull(sq)
    #select the other groups dosed with sqs dosed for this X group
    sq_groups = long_sq_group_df %>% 
      filter(sq %in% xg_sqs) %>% 
      pull(group)
    #set up options for Y groups (exclude if it's the same as X, only keep those dosed with same SQ)
    if (xlibtype==ylibtype){
      pair_options = ygroups[ygroups != xg]
    } else{
      pair_options = ygroups
    }
    pair_options = pair_options[pair_options %in% sq_groups]
    
    sub_df = tibble('x_group' = xg,
                    'x_library_type' = xlibtype,
                    'y_group' = pair_options,
                    'y_library_type' = ylibtype)
    sub_df_list[[xg]] = sub_df
  }
  compare_df = bind_rows(sub_df_list) %>% 
    mutate(library_type_pair = paste(xlibtype, ylibtype, sep='_')) %>% 
    data.frame()
  return(compare_df)
}

get_paired_compare_df = function(pairs_df){
  sub_df_list = list()
  xgroups = pairs_df %>% 
    pull(group) %>% 
    unique()
  for (xg in xgroups){
    pair_options = xgroups[xgroups != xg]
    sub_df = tibble('x_group' = xg,
                    'x_library_type' = 'paired',
                    'y_group' = pair_options,
                    'y_library_type' = 'paired')
    sub_df_list[[xg]] = sub_df
  }
  compare_df = bind_rows(sub_df_list) %>% 
    mutate(library_type_pair = 'pair_pair') %>% 
    data.frame()
  return(compare_df)
}


# BUILD COMPARISON REPORT -------------------------------------------------


build_comparison_report_obj = function(comparison_num, x_group, y_group, x_library_type, y_library_type){
  report_title = paste0('Comparison report<br>',
                        'X: ', x_group, ' (', x_library_type, ')<br>',
                        'Y: ', y_group, ' (', y_library_type, ')')
  
  #set up tag_list for this comparison
  tag_list = list()
  
  #sqs for these groups
  x_sqs = long_sq_df %>% 
    filter(group == x_group) %>% 
    pull(sq) %>% 
    unique()
  y_sqs = long_sq_df %>% 
    filter(group == x_group) %>% 
    pull(sq) %>% 
    unique()
  grp_utranscripts = x_sqs[x_sqs %in% y_sqs] 
  
  #BUILD FIGURES FOR ALL TRANSCRIPTS
  #set up scatterplot
  # print('plotting all transcripts...')
  s_transcript = 'All transcripts'
  tag_list[[s_transcript]] = h3(s_transcript)
  plt_id = paste(s_transcript, 'scatterplot')
  #get widened df for this comparison
  wide_pair_df = stat_df %>% 
    widen_for_pair(compare_ycol,
                   x_group,
                   y_group,
                   x_library_type,
                   y_library_type)
  #plot scatterplot
  plt = wide_pair_df %>% 
    plotly_compare_group_stat(compare_ycol,
                              x_group,
                              y_group,
                              x_library_type,
                              y_library_type,
                              log2_axes = compare_log2_axes,
                              overlap_colors = overlap_colors)
  tag_list[[plt_id]] = plt
  
  #SEQUENCE LOGOS
  #build logos
  logo_list = wide_pair_df %>% 
    build_overlap_logos(overlap_colors)
  #add the logos to tag_list
  for (n in names(overlap_colors)){
    logo_id = paste(s_transcript, n, 'logo')
    tag_list[[logo_id]] = logo_list[[n]]
  }
  
  # #REPEAT FOR EACH TRANSCRIPT
  for (s_transcript in grp_utranscripts){
    # print(paste0('plotting ', s_transcript, '...'))
    #set up title and subset data
    tag_list[[s_transcript]] = h3(s_transcript)
    plt_id = paste(s_transcript, 'scatterplot')
    st_wide_pair_df = wide_pair_df %>% 
      filter(transcript == s_transcript)
    
    #plot scatterplot
    plt = st_wide_pair_df %>% 
      plotly_compare_group_stat(compare_ycol,
                                x_group,
                                y_group,
                                x_library_type,
                                y_library_type,
                                log2_axes = compare_log2_axes,
                                overlap_colors = overlap_colors)
    tag_list[[plt_id]] = plt
    #SEQUENCE LOGOS
    #build logos
    logo_list = st_wide_pair_df %>% 
      build_overlap_logos(overlap_colors)
    #add the logos to tag_list
    for (n in names(overlap_colors)){
      logo_id = paste(s_transcript, n, 'logo')
      tag_list[[logo_id]] = logo_list[[n]]
    }
  }
  # #render the report for this comparison
  report_outfile = paste0(comparison_num, '.html')
  report_list = list(input = individual_markdown_script,
                     output_file = report_outfile,
                     output_dir = comparison_output_dir,
                     report_title = report_title,
                     tag_list = tag_list)
  
  object_path = paste0(comparison_output_dir, '/', comparison_num, '_rep.Rdata')
  save(report_list, file = object_path)
  return(object_path)
}
