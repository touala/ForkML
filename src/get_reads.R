#!/usr/bin/env Rscript

# 'Transform Fork-seq data from .fa and .fa_ratio files to a R table as a .rds file.

library(optparse)

option_list <- list(
  make_option(c("-f", "--fasta"), type = "character", help = "Fasta input file (.fa)", metavar = "FILE"),
  make_option(c("-r", "--ratio"), type = "character", help = "Ratio input file (.fa_ratio)", metavar = "FILE"),
  make_option(c("-c", "--cigar"), type = "character", help = "CIGAR input file (.fa_cigar)", metavar = "FILE"),
  make_option(c("-o", "--out"), type = "character", help = "Reads signals output file (.rds)", metavar = "FILE"),
  make_option("--legacy", action = "store_true", default = FALSE, help = "Run in legacy mode"),
  make_option("--min_length", type = "integer", default = 10000, help = "Minimum length of the read [default: %default]"),
  make_option("--subsampling_step", type = "integer", default = 100, help = "Step for subsampling the signal (1 point every X) [default: %default]"),
  make_option("--sample", type = "character", default = NULL, help = "Sample name [default: %default]")
)

parser <- OptionParser(
  usage = "Usage: %prog -f=FILE -r=FILE -c=FILE -o=FILE [options]",
  option_list = option_list,
  description = "Transform preprocessed Fork-seq data from .fa, .fa_ratio, .fa_cigar files to a R table as a .rds file."
)

opt <- parse_args(parser)


library(tibble)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

extract_headers = function(fasta_lines) {
    fasta_lines[seq(1, length(fasta_lines), by=2)]
}

extract_info = function(fa_info_lines) {
    fa_info_lines[seq(2, length(fa_info_lines), by=2)]
}

parse_header = function(r1) {
    id <- str_match(r1, "^.+\\/([^\\/\\s]+)")[2]
    chrm <- str_match(r1, "'mapped_chrom': '(\\S+)'")[2]
    strand <- str_match(r1, "'mapped_strand': '(.)'")[2]
    start <- str_match(r1, "'mapped_start': (\\d+)")[2] %>% as.numeric() # This is last base without mapping
    end <- str_match(r1, "'mapped_end': (\\d+)")[2] %>% as.numeric()     # This is second to last base with mapping

    tibble(
        "read_id"=id, 
        "read_chrm"=chrm, 
        "read_strand"=strand, 
        "read_start"=start,  # This is last base without mapping
        "read_end"=end       # This is second to last base with mapping
    )
}

parse_ratio = function(ratio_str) {
    ratio_str %>% 
        str_split(' ') %>%  
        unlist %>% 
        as.numeric()
}

parse_cigar <- function(cigar, start_pos){
    # Initialize variables
    ref_position <- start_pos
    query_index <- 1
    mappings <- data.frame(ref_pos = integer(), read_pos = integer())
    
    # Regular expression to capture operations and their lengths from CIGAR
    ops <- gregexpr("\\d+[MIDNSHP=X]", cigar)[[1]]
    lengths <- as.integer(regmatches(cigar, gregexpr("\\d+", cigar))[[1]])
    types <- gsub("\\d+", "", regmatches(cigar, gregexpr("[MIDNSHP=X]", cigar))[[1]])

    # Iterate over operations
    for (i in seq_along(ops)) {
        op_length <- lengths[i]
        op_type <- types[i]
        
        switch(op_type,
            "M" = {
                # Match or mismatch: advance both ref and query indices
                current_mappings <- data.frame(ref_pos = ref_position:(ref_position + op_length - 1), read_pos = query_index:(query_index + op_length - 1))
                mappings <- rbind(mappings, current_mappings)
                ref_position <- ref_position + op_length
                query_index <- query_index + op_length
            },
            "I" = {
                # Insertion to the reference: advance only query index
                query_index <- query_index + op_length
            },
            "D" = {
                # Deletion from the reference: advance only reference position
                ref_position <- ref_position + op_length
            },
            "N" = {
                # Skipped region from the reference: advance only reference position
                ref_position <- ref_position + op_length
            },
            "S" = {
                # Soft clipping: advance only query index
                query_index <- query_index + op_length
            },
            {
                # No operation for H, P, =, X as they are not common or not relevant for query indexing
            }
        )
    }
    
    return(mappings %>% as_tibble())
}

parse_cigar_new <- function(cigar, ref_start_pos, map_strand){
    # cigar is left to right
    decompose_cigar <- data.frame(type=gsub("\\d+", "", regmatches(cigar, gregexpr("[MIDNSHP=X]", cigar))[[1]]), len=as.integer(regmatches(cigar, gregexpr("\\d+", cigar))[[1]])) %>%
        mutate(len_read=ifelse(type=="D", 0, len)) %>%
        mutate(len_ref=ifelse(type %in% c("I","S","H"), 0, len), cs_ref=cumsum(len_ref)) %>%
        mutate(ref_pos_end=(ref_start_pos + 1) + cs_ref - 1) %>% # Correction for offset
        mutate(ref_pos_start=ifelse(type %in% c("I","S","H"), ref_pos_end, ref_pos_end - len_ref + 1))

    if(map_strand=="-"){ # Convert to read orientation
        decompose_cigar <- decompose_cigar[nrow(decompose_cigar):1,]
    }

    # decompose_cigar is read orientation
    mapping <- decompose_cigar %>%
        mutate(cs_read=cumsum(len_read)) %>% # cs=cumsum(len), 
        mutate(read_pos_end=1 + cs_read - 1) %>%
        mutate(read_pos_start=ifelse(type %in% c("D"), read_pos_end, read_pos_end - len_read + 1)) %>%
        mutate(ref_pos=map2(ref_pos_start, ref_pos_end, seq), read_pos=map2(read_pos_start, read_pos_end, function(x, y, s){
            if(s=="+"){
                return(seq(x, y))
            }else{
                return(rev(seq(x, y)))
            }}, s=map_strand))

    mapping <- mapping %>% # map2(read_pos_start, read_pos_end, seq)
        unnest(cols=c(ref_pos, read_pos)) %>% # unnest into one row per aligned base
        select(type, ref_pos, read_pos)

    return(mapping %>% as_tibble())
}

get_xy = function(ratio, strand, start, end) {
    if (strand == '+') x = start + 0:(length(ratio)-1)
    if (strand == '-') x = end - 0:(length(ratio)-1)
    tibble(x, y=ratio) %>% arrange(x)
}

get_xy_new <- function(ratio, mapping) {
    mapping %>%
        full_join(tibble(read_pos=seq(1,length(ratio)), y=ratio), by = join_by(read_pos)) %>%
        dplyr::select(x=read_pos, y) %>%
        arrange(x)
}

subsample_xy = function(signal, subsampling_step) {
    signal <- signal %>%
        mutate(x = round(x/subsampling_step) * subsampling_step) %>% # Binning, TODO Del are duplicated
        group_by(x) %>%
        summarise(y = round(mean(y, na.rm = TRUE), digits=2), .groups='drop')

    return(signal)
}

save_data = function(data_table, path_out) {
    # if there are not reads left, do not save anything
    if(nrow(data_table)==0){
        message("No reads left after filtering.")
    }else{
        saveRDS(data_table, path_out)
    }
}

### MAIN ###
if (sys.nframe()==0) {
    # arguments
    fa_in = opt$fasta
    fa_ratio_in = opt$ratio
    fa_cigar_in = opt$cigar
    out = opt$out

    min_length = as.numeric(opt$min_length)
    subsampling_step = as.numeric(opt$subsampling_step) %>% round()
    sample_name = opt$sample
    
    # load files & extract vectors
    headers_lines = readLines(gzfile(fa_in)) %>% extract_headers()
    ratio_lines = readLines(gzfile(fa_ratio_in)) %>% extract_info()
    cigar_lines = readLines(gzfile(fa_cigar_in)) %>% extract_info()

    # parse headers & make a dataframe
    reads <- map_dfr(headers_lines, parse_header) %>%
        add_column(ratio = map(ratio_lines, parse_ratio)) %>% # add y signal value
        mutate(cigar = cigar_lines) %>% # add cigar value
        mutate(mapping = pmap(., function(cigar, read_strand, read_start, ...){ # add x value & subsample
                parse_cigar_new(cigar, read_start, read_strand)
            }
        ))
    
    if(opt$legacy){
        reads <- reads %>%
            mutate(signal = pmap(., function(ratio, read_strand, read_start, read_end, mapping, ...){ # add x value & subsample
                    get_xy(ratio, read_strand, read_start, read_end) %>%
                        subsample_xy(subsampling_step) %>%
                        filter(!is.na(y))  # remove NaNs
                }
            ))
    }else{
        reads <- reads %>%
            mutate(signal = pmap(., function(ratio, read_strand, read_start, read_end, mapping, ...){ # add x value & subsample
                    get_xy_new(ratio, mapping) %>%
                        subsample_xy(subsampling_step) %>%
                        filter(!is.na(y))  # remove NaNs
                }
            ))
    }

    reads <- reads %>%
        filter(map_lgl(signal, ~ nrow(.x)>0)) %>% # remove empty signals
        mutate(read_inf = map_dbl(signal, ~ min(.x$x)), read_sup = map_dbl(signal, ~ max(.x$x))) %>% # compute signal boundaries
        filter(read_sup - read_inf >= min_length) %>% # filter signal length
        select(read_id, read_chrm, read_strand, read_start, read_end, read_inf, read_sup, signal, mapping) # reorder columns
 
    # add sample 
    if (sample_name!="NULL") {
        reads = reads %>% mutate(sample=sample_name)
    }
    
    # save data
    save_data(reads, out)
}
