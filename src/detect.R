#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-a", "--in_signals"), type="character", help="Signals input file (.tsv(.gz/.bz2/.xz))", metavar="FILE"),
  make_option(c("-b", "--in_tracks"), type="character", help="Tracks input file (.tsv(.gz/.bz2/.xz))", metavar="FILE"),
  make_option(c("-c", "--in_reads"), type="character", help="Reads input file (.rds)", metavar="FILE"),
  make_option(c("-y", "--out_signals"), type="character", help="Signals output file (.tsv(.gz/.bz2/.xz))", metavar="FILE"),
  make_option(c("-z", "--out_tracks"), type="character", help="Tracks output file (.tsv(.gz/.bz2/.xz)); contains detected‐segment parameters", metavar="FILE"),
  make_option("--legacy", action = "store_true", default = FALSE, help = "Run in legacy mode"),
  make_option("--duration_pulse", type="double", default=4, help="BrdU pulse duration (min) [default: %default]"),
  make_option("--frequency_pulse", type="double", default=30, help="BrdU pulse frequency (min) [default: %default]"),
  make_option("--subsampling_step", type="integer", default=100, help="Signal subsampling step [default: %default]"),
  make_option("--size_margin", type="integer", default=3000, help="Margin size around segments (bp) [default: %default]")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print parsed arguments for verification
print(opt)

# suppressWarnings(suppressMessages(library(tidyverse)))
library(tibble)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(reticulate)
library(readr)



################### Helper functions ###################
# Paste operator
`%+%` = paste0

check_empty = function(tib, exec_step_message="Table is empty.") {
  if (nrow(tib)==0) {
    message(exec_step_message %+% " Stopping execution.")
    quit(save="no")
  } else {
    return(tib)
  }
}

################### Arguments control ###################
arg_output = function(file) {
  # Control that the output directory exists and is writable
  # input: character
  # output: character
  dir = dirname(file)
  # directory existence
  if ( !(dir.exists(dir)) ) stop("Output directory does not exist")
  # directory write and read permissions
  if ( !(file.access(dir, mode=6)==0) ) stop("Output directory is not writeable/readable")
  # finally
  return(file)
}
arg_duration_pulse = function(x) {
  # Control the range of the argument
  # input: string
  # output: double
  x_num = as.numeric(x)
  if( !(x_num>0) ) stop('Option --duration_pulse needs to be in the range ]0;+inf[, not "' %+% x_num %+% '"')

  return(x_num)
}
arg_size_margin = function(x) {
  # Control the range of the argument
  # input: string
  # output: integer
  x_num = as.numeric(x)
  if( !(x_num>=1) ) stop('Option --size_margin needs to be in the range [1;+inf[, not "' %+% x_num %+% '"')

  return(x_num)
}
get_params = function(args) {
  # Extract and control command line arguments
  # input: list (from docopt)
  # output: list
  p = list()
  p$in_signals = args$in_signals
  p$in_tracks = args$in_tracks
  p$in_reads = args$in_reads
  p$out_signals = arg_output(args$out_signals)
  p$out_tracks = arg_output(args$out_tracks)
  p$is_legacy = args$legacy
  p$t_pulse = arg_duration_pulse(args$duration_pulse)
  p$f_pulse = args$frequency_pulse
  p$subsampling_step = arg_duration_pulse(args$subsampling_step)
  p$size_neutral = arg_size_margin(args$size_margin)

  return(p)
}

################### Input data control ###################
# Input tibble
sig_input_class = function(signals) {
  # Check that the input object is a tibble.
  if( !is_tibble(signals) ) stop("Input data needs to be a tibble, not a " %+% class(signals))
}
sig_input_names = function(signals) {
  # Check that all needed columns are present in the input data
  # input: tibble
  # output: NULL
  # check columns
  cols_needed = c("read_id", "read_chrm", "read_strand", "signal")
  not_included = !(cols_needed %in% colnames(signals))
  # if at least one is not in
  if(any(not_included)) {
    # paste column names
    cols_needed[not_included] %>%
      map_chr(~ paste0("\"", .x, "\"")) %>%
      paste(collapse=", ") -> cols_not_included
    # print error message
    err_msg = 'Columns ' %+% cols_not_included %+% ' do not exist in the input data table.'
    stop(err_msg)
  }
}
sig_input_types = function(signals) {
  # Check that all needed columns have the right type
  # input: tibble
  # output: NULL
  if ( !is.character(signals$read_id) ) {
    stop('Column "read_id" in the input data needs to be of type "character", not "' %+% class(signals$read_id) %+% '".')
  }
  if ( !is.character(signals$read_chrm) ) {
    stop('Column "read_chrm" in the input data needs to be of type "character", not "' %+% class(signals$read_chrm) %+% '".')
  }
  if ( !is.character(signals$read_strand) ) {
    stop('Column "read_strand" in the input data needs to be of type "character", not "' %+% class(signals$read_strand) %+% '".')
  }
  if ( !is.list(signals$signal) ) {
    stop('Column "signal" in the input data needs to be of type "list", not "' %+% class(signals$signal) %+% '".')
  }
}
sig_input_empty = function(signals) {
  # Check that the input data table has at least 1 row
  # input: tibble
  # output: NULL
  if (nrow(signals)==0) stop("Input data table is empty.")
}

# Column "signal" of nested signals
sig_input_signal_class = function(sig) {
  # Check that the signal column is a list of tibbles
  if( !is_tibble(sig) ) stop('Input "signal" column needs to be a list of tibbles, not a list of ' %+% class(sig))
}
sig_input_signal_names = function(sig) {
  # Check that the nested signals have the right column names
  # input: tibble
  # output: NULL
  cols_needed = c("x", "y")
  not_included = !(cols_needed %in% colnames(sig))
  # if at least one is not in
  if(any(not_included)) {
    # paste column names
    cols_needed[not_included] %>%
      map_chr(~ paste0("\"", .x, "\"")) %>%
      paste(collapse=", ") -> cols_not_included
    # print error message
    err_msg = 'Columns ' %+% cols_not_included %+% ' do not exist in the nested signal tibbles.'
    stop(err_msg)
  }
}
sig_input_signal_types = function(sig) {
  # Check the types of the nested signals
  # input: tibble
  # output: NULL
  if( !is.numeric(sig$x) ) {
    stop('Column "x" in the nested signals must be of type "numeric", not "' %+% class(sig$x) %+% '".')
  }
  if( !is.numeric(sig$y) ) {
    stop('Column "y" in the nested signals must be of type "numeric", not "' %+% class(sig$y) %+% '".')
  }
}
sig_input_signal_contains_na = function(list_sig) {
  # Check if some signals have NAs/NaNs
  # input: tibble
  # output: NULL
  contains_na = list_sig %>% 
    map_lgl(function(sig) { any(is.na(sig$x), is.na(sig$y)) }) %>%
    any()
  if (contains_na) warning("Signals contain missing values (NAs/NaNs)")
}
get_signals = function(file) {
  # Load input data table and perform controls
  # input: <character> file path
  # output: <tibble> signals
  signals = read_tsv(
    file,
    col_types = cols(
      read_id = col_character(),
      read_chrm = col_character(),
      read_strand = col_character(),
      read_start = col_double(),
      read_end = col_double(),
      read_inf = col_double(),
      read_sup = col_double(),
      sample = col_character(),
      x = col_double(),
      y = col_double()
    )
  ) %>% nest(signal = c(x, y))
  sig_input_class(signals)
  sig_input_names(signals)
  sig_input_types(signals)
  sig_input_empty(signals)
  walk(signals$signal, sig_input_signal_class)
  walk(signals$signal, sig_input_signal_names)
  walk(signals$signal, sig_input_signal_types)
  sig_input_signal_contains_na(signals$signal)
  signals
}

# Input tibble
trk_input_class = function(tracks) {
  # Check that the input object is a tibble.
  if( !is_tibble(tracks) ) stop("Input data needs to be a tibble, not a " %+% class(tracks))
}
trk_input_names = function(tracks) {
  # Check that all needed columns are present in the input data
  # input: tibble
  # output: NULL
  # check columns
  cols_needed = c("read_id", "read_chrm", "read_strand", "direction", "pred_inf", "pred_sup")
  not_included = !(cols_needed %in% colnames(tracks))
  # if at least one is not in
  if(any(not_included)) {
    # paste column names
    cols_needed[not_included] %>%
      map_chr(~ paste0("\"", .x, "\"")) %>%
      paste(collapse=", ") -> cols_not_included
    # print error message
    err_msg = 'Columns ' %+% cols_not_included %+% ' do not exist in the input data table.'
    stop(err_msg)
  }
}
trk_input_types = function(tracks) {
  # Check that all needed columns have the right type
  # input: tibble
  # output: NULL
  if ( !is.character(tracks$read_id) ) {
    stop('Column "read_id" in the input data needs to be of type "character", not "' %+% class(tracks$read_id) %+% '".')
  }
  if ( !is.character(tracks$read_chrm) ) {
    stop('Column "read_chrm" in the input data needs to be of type "character", not "' %+% class(tracks$read_chrm) %+% '".')
  }
  if ( !is.character(tracks$read_strand) ) {
    stop('Column "read_strand" in the input data needs to be of type "character", not "' %+% class(tracks$read_strand) %+% '".')
  }
  if ( !is.numeric(tracks$direction) ) {
    stop('Column "direction" in the input data needs to be of type "numeric", not "' %+% class(tracks$direction) %+% '".')
  }
  if ( !is.numeric(tracks$pred_inf) ) {
    stop('Column "pred_inf" in the input data needs to be of type "numeric", not "' %+% class(tracks$pred_inf) %+% '".')
  }
  if ( !is.numeric(tracks$pred_sup) ) {
    stop('Column "pred_sup" in the input data needs to be of type "numeric", not "' %+% class(signals$pred_sup) %+% '".')
  }   
}
trk_input_empty = function(tracks) {
  # Check that the input data table has at least 1 row
  # input: tibble
  # output: NULL
  if (nrow(tracks)==0) stop("Input data table is empty.")
}
get_tracks = function(file) {
  # Load input data table and perform controls
  # input: <character> file path
  # output: <tibble> tracks
  tracks = read_tsv(
    file,
    col_types = cols(
      read_id = col_character(),
      read_chrm = col_character(),
      read_strand = col_character(),
      read_start = col_double(),
      read_end = col_double(),
      read_inf = col_double(),
      read_sup = col_double(),
      sample = col_character(),
      track_id = col_double(),
      direction = col_double(),
      pred_inf = col_double(),
      pred_sup = col_double()
    )
  )
  trk_input_class(tracks)
  trk_input_names(tracks)
  trk_input_types(tracks)
  trk_input_empty(tracks)
  tracks
}


################### Processing ###################
# Helper
pre_filter_empty = function(signals) {
  # Remove NAs/NaNs from the signals, then remove the empty signals
  # input: <n-elts list> signals
  # output: <m-elts list> signals (m<=n)
  signals %>%
    # remove NAs/NaNs in signals
    mutate(signal = map(
      signal,
      function(signal) {
        signal %>%
          filter(!is.na(x) & !is.na(y))
      }
    )) %>%
    # remove reads with empty signals
    filter(map_lgl(signal, ~ nrow(.x)>0)) %>%
    # return
    return()
}
pre_extend_bounds = function(tracks, size_neutral) {
  # Extend boundaries of tracks when possible
  # input: <tibble>
  # output: <tibble>
  tracks %>%
    group_by(read_id) %>%
    arrange(read_id, pred_inf) %>%
    mutate(
      pred_mid_inf = case_when(
        row_number()==1 ~ read_inf, 
        T               ~ (pred_inf + lag(pred_sup)) / 2
      ),
      pred_mid_sup = case_when(
        row_number()==n() ~ read_sup,
        T                 ~ (pred_sup + lead(pred_inf)) / 2
      )
    ) %>%
    ungroup() %>%
    mutate(
      pred_ext_inf = pmax(pred_inf-size_neutral, pred_mid_inf),
      pred_ext_sup = pmin(pred_sup+size_neutral, pred_mid_sup)
    )
}
pre_join_signals = function(tracks, signals) {
  # Join signals to the tracks tibble
  # input: <tibble>, <tibble>
  # output: <tibble>
  tracks %>%
    left_join(
      signals %>% select(read_id, signal), 
      by="read_id"
    ) %>%
    mutate(signal = pmap(
      .,
      function(signal, pred_ext_inf, pred_ext_sup, ...) {
        signal %>%
          filter(x>=pred_ext_inf & x<=pred_ext_sup)
      }
    ))
}

preprocessing = function(tracks, signals, size_neutral) {
  # Preprocess the input data
  signals %>% 
    pre_filter_empty() %>% 
    check_empty("No more signals after NAs/NaNs filtering.") -> signals_nonempty
  tracks %>%
    pre_extend_bounds(size_neutral=size_neutral) %>%
    pre_join_signals(signals_nonempty) %>%
    filter(map_lgl(signal, ~ nrow(.x)>0)) %>%
    check_empty("No more tracks after joining signals & filtering out empty signals.") %>%
    return()
}

post_bind_models_with_signals = function(signals, tracks) {
  signals %>%
    semi_join(tracks, by="read_id") %>%
    unnest(signal) %>%
    mutate(type="raw", group=0) %>%
    group_by(read_id) %>%
    nest(signal = c(x, y, type, group)) %>%
    ungroup()
}

postprocessing_signals = function(signals, tracks, p) {
  # Postprocess the signals by adding the models of the detection and the tracks.

  # bind to signals
  signals %>%
    select(read_id, signal) %>%
    post_bind_models_with_signals(tracks) %>%
    return()
}

convert.mapping.coordinate <- function(tracks, reads, subsampling_step, t_pulse, f_pulse, is_legacy=FALSE){
  # subsampling_step <- p$subsampling_step
  # t_pulse <- p$t_pulse
  # f_pulse <- p$f_pulse
  # is_legacy <- TRUE
  
  # Define read <-> genome mapping table
  if(is_legacy){
    col_bin <- "ref_pos"
    col_match <- "read_pos"
  }else{    
    col_bin <- "read_pos"
    col_match <- "ref_pos"
  }
  
  mapping_df <- reads %>%
    mutate(mapping = pmap(., function(mapping, subsampling_step, col_bin, col_match, ...){
        # res1 <- mapping %>% # For full bins
        #   # dplyr::filter(read_pos > max(read_pos) - subsampling_step*2) %>% # half_subsampling_step is intermediate but not exact. and missing last bin sometimes
        #   mutate(read_bin = round(read_pos/subsampling_step) * subsampling_step) %>%
        #   group_by(read_bin) %>%
        #   summarize(approx_ref_pos=round(mean(ref_pos), 0)) %>%
        #   arrange(read_bin)
        # res2 <- mapping %>% # For half bins, full bin computed twice
        #   mutate(read_bin = round(read_pos/(subsampling_step/2)) * (subsampling_step/2)) %>%
        #   group_by(read_bin) %>%
        #   summarize(approx_ref_pos=round(mean(ref_pos), 0)) %>%
        #   arrange(read_bin)
        # res <- rbind(res1, res2 %>% dplyr::filter(!read_bin %in% res1$read_bin)) %>%
        #   arrange(read_bin)

        res1 <- mapping %>% # For full bins
          mutate(bin_pos = !!sym(col_bin), matching_pos = !!sym(col_match)) %>%
          mutate(bin = round(bin_pos/subsampling_step) * subsampling_step) %>%
          group_by(bin) %>%
          summarize(approx_matching_pos = round(mean(matching_pos), 0)) %>%
          arrange(bin)
        res2 <- mapping %>% # For half bins, full bin computed twice
          mutate(bin_pos = !!sym(col_bin), matching_pos = !!sym(col_match)) %>%
          mutate(bin = round(bin_pos/(subsampling_step/2)) * (subsampling_step/2)) %>%
          group_by(bin) %>%
          summarize(approx_matching_pos = round(mean(matching_pos), 0)) %>%
          arrange(bin)
        res <- rbind(res1, res2 %>% dplyr::filter(!bin %in% res1$bin)) %>%
          arrange(bin)

        return(res)
      }, subsampling_step=subsampling_step)) %>%
    dplyr::select(read_id, mapping) %>%
    unnest(mapping) %>%
    # dplyr::select(read_id, read_bin, approx_ref_pos)
    dplyr::select(read_id, bin, approx_matching_pos)

  # Convert coordinate
  pred_cols <- c("pred_inf","pred_sup","pred_mid_inf","pred_mid_sup","pred_ext_inf","pred_ext_sup")

  # tracks <- tracks %>%
  #   pivot_longer(cols = all_of(pred_cols), names_to = "pred_field", values_to = "read_bin") %>%
  #   left_join(mapping_df, by = c("read_id", "read_bin")) %>%
  #   pivot_wider(names_from = pred_field, values_from = c(read_bin, approx_ref_pos), names_glue = "{pred_field}_{.value}") %>%
  #   rename_with(.fn = ~ sub("pred_(.*)_approx_ref_pos", "pred_\\1", .x), .cols = ends_with("approx_ref_pos")) %>%
  #   mutate(
  #     pred_inf_flipped = if_else(read_strand == "-", pred_sup, pred_inf),
  #     pred_sup = if_else(read_strand == "-", pred_inf, pred_sup),
  #     pred_inf = pred_inf_flipped,

  #     pred_mid_inf_flipped = if_else(read_strand == "-", pred_mid_sup, pred_mid_inf),
  #     pred_mid_sup = if_else(read_strand == "-", pred_mid_inf, pred_mid_sup),
  #     pred_mid_inf = pred_mid_inf_flipped,

  #     pred_ext_inf_flipped = if_else(read_strand == "-", pred_ext_sup, pred_ext_inf),
  #     pred_ext_sup = if_else(read_strand == "-", pred_ext_inf, pred_ext_sup),
  #     pred_ext_inf = pred_ext_inf_flipped
  #   ) %>%
  #   dplyr::select(-pred_inf_flipped, -pred_mid_inf_flipped, -pred_ext_inf_flipped)

  tracks <- tracks %>%
    pivot_longer(cols = all_of(pred_cols), names_to = "pred_field", values_to = "bin") %>%
    left_join(mapping_df, by = c("read_id", "bin")) %>%
    pivot_wider(names_from = pred_field, values_from = c(bin, approx_matching_pos), names_glue = "{pred_field}_{.value}") %>%
    rename_with(.fn = ~ sub("pred_(.*)_approx_matching_pos", "pred_\\1", .x), .cols = ends_with("approx_matching_pos")) %>%
    mutate(
      pred_inf_flipped = if_else(read_strand == "-", pred_sup, pred_inf),
      pred_sup = if_else(read_strand == "-", pred_inf, pred_sup),
      pred_inf = pred_inf_flipped,

      pred_mid_inf_flipped = if_else(read_strand == "-", pred_mid_sup, pred_mid_inf),
      pred_mid_sup = if_else(read_strand == "-", pred_mid_inf, pred_mid_sup),
      pred_mid_inf = pred_mid_inf_flipped,

      pred_ext_inf_flipped = if_else(read_strand == "-", pred_ext_sup, pred_ext_inf),
      pred_ext_sup = if_else(read_strand == "-", pred_ext_inf, pred_ext_sup),
      pred_ext_inf = pred_ext_inf_flipped
    ) %>%
    dplyr::select(-pred_inf_flipped, -pred_mid_inf_flipped, -pred_ext_inf_flipped)

  tracks$subsampling_step <- subsampling_step
  tracks$t_pulse <- t_pulse
  tracks$f_pulse <- f_pulse

  return(tracks)
}

################### Main ###################
if (sys.nframe()==0) {
  #### Arguments & Input data
  p <- get_params(opt)
  in_signals <- get_signals(p$in_signals)
  in_tracks <- get_tracks(p$in_tracks)
  
  #### Processing
  tracks <- preprocessing(in_tracks, in_signals, p$size_neutral)
  signals <- postprocessing_signals(in_signals, tracks, p)
  reads <- readRDS(p$in_reads)

  #### Convert coordinate: read to genome
  tracks <- convert.mapping.coordinate(tracks, reads, p$subsampling_step, p$t_pulse, p$f_pulse, p$is_legacy)

  #### Save data
  write_tsv(signals %>% unnest(signal), p$out_signals)
  write_tsv(tracks, p$out_tracks)
}
