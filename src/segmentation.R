#!/usr/bin/env Rscript

# global variable
NETWORK_SUBSAMPLING = 100

# Paste operator
`%+%` = paste0

'Segment Fork-seq reads by direction of DNA replication.

Usage:
    segmentation -i=FILE -t=FILE -s=FILE -m=NAME [options]

Mandatory arguments:
    -i FILE, --in_reads FILE    Reads input file, in .tsv(.gz)(.bz2)(.xz) format.
    -t FILE, --tracks FILE      Tracks output file, in .tsv(.gz)(.bz2)(.xz) format.
                                Contains info about the detected replication tracks.
    -s FILE, --signals FILE     Signals output file in .tsv(.gz)(.bz2)(.xz) format.
                                Contains BrdU signals of the reads where at least 
                                one replication track was detected.
    -m NAME, --model_name NAME  Name of ForkML model to use for detection.
                                Models can be stored in <ForkML_dir>/models.

Optional arguments:
    --subsampling X           Subsampling of the input data (divisor of '%+%NETWORK_SUBSAMPLING%+%') [default: '%+%NETWORK_SUBSAMPLING%+%']

' -> doc

# load packages
library(docopt)
args = docopt(doc)

library(tibble)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(reticulate)
library(readr)

################### Helper functions ###################
# Input type check
is_dbl = function(x) {
    # Check if a string is coercible as a double (ex: "1.5")
    # input: string
    # output: logical
    !is.na(suppressWarnings(as.numeric(x)))
}

is_int = function(x) {
    # Check if a string is coercible as an integer (ex: "100")
    # input: string
    # output: logical
    out=is_dbl(x)
    if(out) {
        x_num = as.numeric(x)
        out = x_num==round(x_num)
    }
    out
}

# Find script path automatically
# Will not work through a symbolic link referring directly to the script
getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}

# Exit if tibble is empty
check_empty = function(tib, exec_step_message="Table is empty.") {
    if (nrow(tib)==0) {
        message(exec_step_message %+% " Stopping execution.")
        quit(save="no")
    } else {
        return(tib)
    }
}

################### Arguments check ###################
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

arg_subsampling = function(x) {
    # Control the type and range of the argument
    # input: string
    # output: integer
    # type
    if( !is_int(x) ) stop('Option --subsampling needs to be an integer, not "' %+% x %+% '"')
    # range
    x_num = as.numeric(x)
    if( !(x_num>=1 & x_num<=NETWORK_SUBSAMPLING & NETWORK_SUBSAMPLING%%x_num==0) ) {
        stop('Option --subsampling needs to be a divisor of '%+%NETWORK_SUBSAMPLING%+%' in the range [1;'%+%NETWORK_SUBSAMPLING%+%'], not "'%+%x_num%+%'"')
    } 
    # finally
    return(x_num)
}

get_params = function(args) {
    # Extract and control command line arguments
    # input: list (from docopt)
    # output: list
    p = list()
    p$in_reads = args$in_reads
    p$out_tracks = arg_output(args$tracks)
    p$out_signals = arg_output(args$signals)
    p$subsampling = arg_subsampling(args$subsampling)
    p$model_name = getScriptPath() %+% "/../src/" %+% "../models/" %+% args$model_name

    return(p)
}

################### Input check ###################
# Input tibble
input_class = function(signals) {
    # Check that the input object is a tibble.
    if( !is_tibble(signals) ) stop("Input data needs to be a tibble, not a " %+% class(signals))
}

input_names = function(signals) {
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

input_types = function(signals) {
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

input_empty = function(signals) {
    # Check that the input data table has at least 1 row
    # input: tibble
    # output: NULL
    if (nrow(signals)==0) {
        message("Input data table is empty.")
        quit(save="no")
    }
}

# Column "signal" of nested signals
input_signal_class = function(sig) {
    # Check that the signal column is a list of tibbles
    if( !is_tibble(sig) ) stop('Input "signal" column needs to be a list of tibbles, not a list of ' %+% class(sig))
}

input_signal_names = function(sig) {
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

input_signal_types = function(sig) {
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

input_signal_contains_na = function(list_sig) {
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
    # signals = read_tsv(
    #     file,
    #     col_types = cols(
    #         read_id = col_character(),
    #         read_chrm = col_character(),
    #         read_strand = col_character(),
    #         read_start = col_double(),
    #         read_end = col_double(),
    #         read_inf = col_double(),
    #         read_sup = col_double(),
    #         x = col_double(),
    #         y = col_double(),
    #         sample = col_character()
    #     )
    # ) %>% nest(signal = c(x, y))
    signals <- readRDS(file)

    input_class(signals)
    input_names(signals)
    input_types(signals)
    input_empty(signals)
    walk(signals$signal, input_signal_class)
    walk(signals$signal, input_signal_names)
    walk(signals$signal, input_signal_types)
    input_signal_contains_na(signals$signal)
    
    signals
}

################### Preprocessing ###################
pre_subsample = function(reads, subsampling) {
    # Subsample the signals, if needed, so that it matches NETWORK_SUBSAMPLING
    # input: <tibble>, <numeric>
    # output: <tibble>
    if (subsampling == NETWORK_SUBSAMPLING) {
        # if the signals subsampling already matches NETWORK_SUBSAMPLING, do nothing.
        out = reads
    } else {
        # compute the subsampling step
        subsampling_step = NETWORK_SUBSAMPLING/subsampling
        # subsample
        out = reads %>%
        mutate(signal = map(
            signal, 
            ~ .x %>% slice(seq(1, nrow(.), by=subsampling_step))
        ))
    }
    return(out)
}

pre_reads_lengths = function(reads) {
    # Extract lengths to remove padding after prediction
    # input: <tibble>
    # output: <numeric>
    reads$signal %>%
    map_dbl(nrow)
}

pre_reshape = function(reads) {
    # Reshape from tibble with nested signals, to list of matrices
    # input: <tibble> reads
    # output: <list>
    reads %>% 
        select(read_id, signal) %>%
        unnest(signal) %>%
        # split by read_id
        select(read_id, y) %>%
        nest(data=y) %>%
        pull(data) %>%
        # cast tibbles to matrices
        map(~ .x$y) %>%
        map(as.matrix)
}

################### Segmentation ###################
segmentation = function(reads_lengths, reads_list, model_name) {
    # Perform the segmentation
    # input: <numeric>, <list>
    # output: <array>
    # source python functions

    path_to_R = getScriptPath() %+% "/../src/"
    source_python(path_to_R %+% "segmentation_helper.py")
    
    # pad data
    max_pad_length = compute_padding_length(reads_lengths)
    reads_array = pad_reads(unname(reads_list), as.integer(max_pad_length))
    
    # predict
    preds_array = predict(reads_array, model_name)
    
    return(preds_array)
}

################### Postprocessing ###################
postprocessing = function(preds_array, reads_lengths, reads) {
    # Transform reads info and predictions array to a single tibble,
    # and summarise to 1 row per track
    # input: <array>, <numeric>, <tibble>
    # output: <tibble>

    # transform to tibble
    preds_to_bind <- preds_array %>% 
        apply(c(1,2), which.max) %>% 
        t() %>% 
        as_tibble(.name_repair="unique") %>%
        gather(value="label") %>% 
        # remove padding
        nest(data=-key) %>%
        mutate(length = unlist(reads_lengths)) %>% 
        unnest(data) %>%
        group_by(key) %>% 
        mutate(x = row_number()) %>% 
        filter(x <= length) %>% 
        ungroup() %>%
        select(label)

    # unnest reads
    reads_to_bind <- reads %>% 
        unnest(signal)
    
    # bind reads and predictions
    tracks <- bind_cols(reads_to_bind, preds_to_bind) %>%
        # find segments boundaries
        group_by(read_id) %>%
        mutate(boundary = case_when(
            lag(label, default=3)!=label & lead(label, default=3)!=label ~ as.character(NA), # don't consider segments of length 1
            label!=3 & lag(label, default=3)!=label ~ "pred_inf",
            label!=3 & lead(label, default=3)!=label ~ "pred_sup"
        )) %>%
        # remove other rows
        filter(!is.na(boundary)) %>%
        check_empty("Segmentation : No tracks found.") %>%
        # make track id for spreading, then spread to wide format
        mutate(track_id = ceiling(row_number()/2)) %>% 
        ungroup() %>% 
        mutate(direction=(-1)**(label+1)) %>% 
        select(-y, -label) %>% 
        spread(key=boundary, value=x)
    
    return(tracks)
}

post_filter_reads = function(reads, tracks) {
    # Keep only reads that contain tracks
    # input: <tibble>
    # output: <tibble>
    semi_join(reads, tracks, by="read_id")
}

################### Main ###################
if (sys.nframe()==0) {
    #### Arguments & Input data
    p = get_params(args)    
    reads = get_signals(p$in_reads)
    
    #### Preprocessing
    reads_subsampled = pre_subsample(reads, p$subsampling)
    reads_lengths = pre_reads_lengths(reads_subsampled)
    reads_list = pre_reshape(reads_subsampled)
    
    #### Segmentation
    preds_array = segmentation(reads_lengths, reads_list, p$model_name)
    
    #### Postprocessing
    tracks = postprocessing(preds_array, reads_lengths, reads_subsampled)
    signals = post_filter_reads(reads, tracks)
    
    #### Save
    write_tsv(tracks, p$out_tracks)
    write_tsv(signals %>% unnest(signal), p$out_signals)
}
