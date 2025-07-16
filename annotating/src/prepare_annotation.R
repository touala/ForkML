#!/usr/bin/env Rscript

.libPaths(.libPaths()[1])

library("optparse")

# Option parser setup
option_list <- list(
  make_option(c("-b", "--path_bam"), type = "character", help = "Path to BAM file", metavar = "file"),
  make_option(c("-n", "--name_sample"), type = "character", help = "Sample name", metavar = "name"),
  make_option(c("-o", "--path_output"), type = "character", help = "Path output directory", metavar = "dir"),
  make_option(c("-r", "--name_reviewer"), type = "character", action = "store", help = "Reviewer names (comma-separated)", metavar = "names"),
  make_option(c("-i", "--path_roi"), type = "character", help = "Path to read of interest file", metavar = "file"),
  make_option(c("--path_discard"), type = "character", default = NULL, help = "Path to read to discard file", metavar = "file"),
  make_option(c("--nb_random"), type = "integer", default = 0, help = "Number of random samples (default: 0)"),
  make_option(c("--nb_with_signal"), type = "integer", default = 1000, help = "Number of samples with signal (default: 1000)"),
  make_option(c("--min_read_length"), type = "integer", default = 10000, help = "Minimum read length (default: 10000)"),
  make_option(c("--batch_size"), type = "integer", default = 50, help = "Batch size for processing (default: 50)"),
  make_option(c("--nb_threads"), type = "integer", default = 10, help = "Number of threads (default: 10)"),
  make_option(c("--with_independent_review"), action="store_true", default=FALSE, help="Enable independent review (default: FALSE)"),
  make_option(c("--target_mod_type"), type = "character", default = "N+b?", help = "Target modification type (default: 'N+b?')"),
  make_option(c("--data_type"), type = "character", default = "Megalodon", help = "Data type (default: 'Megalodon', or 'DNAscent')")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print parsed arguments for verification
print(opt)

# Assign values from options
path_annotation <- opt$path_annotation
path_output <- gsub("/$", "", opt$path_output)
name_sample <- opt$name_sample
path_bam <- opt$path_bam
name_reviewer <- unlist(strsplit(opt$name_reviewer, ",")) # Convert comma-separated string to vector
path_list_id_with_signal <- opt$path_roi
path_discard <- opt$path_discard
nb_random <- opt$nb_random
nb_with_signal <- opt$nb_with_signal
min_read_length <- opt$min_read_length
batch_size <- opt$batch_size
nb_threads <- opt$nb_threads
with_independent_review <- opt$with_independent_review
target_mod_type <- opt$target_mod_type
data_type <- opt$data_type

if(path_discard=="null"){path_discard <- NULL}
if(!data_type %in% c("Megalodon", "DNAscent")){quit("data_type unknown.", 2)}
if(data_type=="Megalodon"){target_mod_type <- NULL}

library("future")
library("Rsamtools")
library("tibble")
library("doFuture")
library("dplyr")
library("stringr")
library("tidyr")
library("progressr")
library("furrr")
library("purrr")
library("RcppRoll")
library("data.table")

options(future.rng.onMisuse="ignore") # Warning UNRELIABLE VALUE for RNG

# Load R10 specific signal reading function
list_objects_with_sizes <- function(sort_by="size"){
  # Retrieve a list of objects from the global environment
  object_names <- ls(envir = .GlobalEnv)
  
  # Calculate the size of each object and convert to human-readable format
  object_sizes <- sapply(object_names, function(name) {
    object.size(get(name, envir = .GlobalEnv))
  })
  
  # Create a data frame with object names and sizes
  df <- data.frame(
    Object = object_names,
    Size = object_sizes,
    stringsAsFactors = FALSE
  )
  
  # Convert sizes to a more readable format
  df$Size <- format(df$Size, units = "auto")
  
  # Sort the data frame by size or object name
  if (tolower(sort_by) == "size") {
    df <- df[order(-object_sizes), ]
  } else if (tolower(sort_by) == "name") {
    df <- df[order(df$Object), ]
  } else {
    stop("Invalid sort_by argument; use 'size' or 'name'.")
  }
  
  # Return the data frame
  return(df)
}

find.closest.number <- function(x, n) {
    if (x == 0) {
        return(0)
    }else{
        closest_multiple <- ceiling(x / n) * n
        return(closest_multiple)
    }
}

generate.bam.subset <- function(dataset_to_split, nb_better_target_read, nb_better_target_read_with_signal, min_read_length, path_list_id_to_discard, path_list_id_with_signal, intermediate_bam, intermediate_ori, intermediate_txt){
    if (nb_better_target_read > 0) {
        cat("Select random read_id (long).\n")
        bam <- BamFile(dataset_to_split)

        if (!file.exists(gsub(".bam$", ".bam.bai", dataset_to_split))) {
            stop(paste0(dataset_to_split, " is not indexed yet."))
        }
        param <- ScanBamParam(what = c("qname", "qwidth", "flag", "cigar"))
        bam_data <- scanBam(bam, param = param)

        reads <- as_tibble(bam_data[[1]])
        mapped_reads <- reads[bitwAnd(reads$flag, 0x4) == 0 & reads$qwidth >= min_read_length,]

        # Filtering available read_id
        list_read_id <- mapped_reads$qname # Full list if read_id
        if (!is.null(path_list_id_to_discard)) {
            # Remove blacklisted read_id
            cat("Remove blacklisted read_id.\n")
            blacklisted_ids <- read.table(path_list_id_to_discard, header = FALSE, stringsAsFactors=FALSE)[[1]]
            list_read_id <- list_read_id[!list_read_id %in% blacklisted_ids]
        }
        set.seed(101)
        list_subsampled_read_id <- sample(list_read_id, size = min(length(list_read_id), nb_better_target_read), replace = FALSE) # Random sampling                
    } else {
        cat("Only doped read_id.\n")
        list_subsampled_read_id <- NULL
    }

    list_id_with_signal <- unlist(read.table(path_list_id_with_signal, header = FALSE, stringsAsFactors = FALSE), use.names = FALSE) # Read_id with signal
    list_id_with_signal <- list_id_with_signal[!list_id_with_signal %in% list_subsampled_read_id] # Usable read_id with signal
    if (!is.null(path_list_id_to_discard)) {
        # Remove blacklisted read_id
        blacklisted_ids <- read.table(path_list_id_to_discard, header = FALSE)[[1]]
        list_id_with_signal <- list_id_with_signal[!list_id_with_signal %in% blacklisted_ids]
    }
    set.seed(101)
    list_subsampled_id_with_signal <- sample(list_id_with_signal, size = min(length(list_id_with_signal), nb_better_target_read_with_signal), replace = FALSE) # Random sampling, no length filtering

    # Final list of read_id
    labelled_read_id <- data.frame(
        read_id = c(list_subsampled_read_id, list_subsampled_id_with_signal), 
        source = c(rep("random", length(list_subsampled_read_id)), rep("doped", length(list_subsampled_id_with_signal)))
    )

    write.table(labelled_read_id, file = intermediate_ori, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    list_subsampled_read_id <- c(list_subsampled_read_id, list_subsampled_id_with_signal)

    cat(paste0("Writing temporary read_id list at ", intermediate_txt, " (n=", length(list_subsampled_read_id), ").\n"))
    writeLines(list_subsampled_read_id, intermediate_txt)

    cat("Writing subset of bam file.\n")
    system(sprintf("samtools view -b -N %s %s > %s", intermediate_txt, dataset_to_split, intermediate_bam))
    cat("Indexing subset of bam file.\n")
    system(sprintf("samtools index %s", intermediate_bam))
}

read.bam.batch <- function(bam_file, list_bam_queries, data_type){
    bam_reading <- scanBam(bam_file, param=list_bam_queries)

    if(length(bam_reading[[1]]$qname)==0){
        bam_batch <- FALSE

        return(bam_batch)
    }

    bam_batch <- tibble(
        read_id=bam_reading[[1]]$qname,
        flag=bam_reading[[1]]$flag,
        chrom=bam_reading[[1]]$rname,
        strand=bam_reading[[1]]$strand,
        start=bam_reading[[1]]$pos,
        qwidth=bam_reading[[1]]$qwidth,
        seq=as.character(bam_reading[[1]]$seq)
    )
    
    if(data_type=="DNAscent"){
        bam_batch <- bam_batch %>%
            add_column(cigar=as.character(bam_reading[[1]]$cigar))
    }

    if(is.null(bam_reading[[1]]$tag$Mm)){
        bam_batch <- bam_batch %>%
            add_column(
                pos=bam_reading[[1]]$tag$MM,
                prob=bam_reading[[1]]$tag$ML
            )
    }else{
        bam_batch <- bam_batch %>%
            add_column(
                pos=bam_reading[[1]]$tag$Mm,
                prob=bam_reading[[1]]$tag$Ml
            )
    }
    rm(bam_reading)

    # Filter unmapped reads
    bam_batch <- bam_batch %>%
        filter(flag!=4)

    if(nrow(bam_batch)==0){
        bam_batch <- FALSE

        return(bam_batch)
    }

    return(bam_batch)
}

process.mod.tag <- function(bam_MM_item, bam_ML_item, bam_seq, data_type, target_mod_type = NULL){
    if(data_type=="DNAscent"){
        bam_MM_item <- bam_MM_item %>%
            str_remove(";$") %>%
            str_split(";", simplify=TRUE)

        list_probabilities <- unlist(bam_ML_item)
        list_type_mod <- NULL
        length_ML_processed <- 0
        bam_mod_info <- foreach(idx=seq(1,length(bam_MM_item))) %do% {
            bam_MM_clean <- bam_MM_item[,idx] %>% str_split(",", simplify=TRUE)
            type_mod <- bam_MM_clean[,1]

            list_type_mod <- c(list_type_mod, type_mod)
            if(type_mod == target_mod_type){
                type_base <- substr(bam_MM_clean[,1],1,1)

                mod_rel_pos <- bam_MM_clean[,-1] %>% as.numeric() # Can be longer if supplementary mapping truncated

                if(type_base=="N"){
                    pos_base_in_sequence <- cumsum(mod_rel_pos + 1)
                    mod_idx <- seq(along=pos_base_in_sequence)
                }else{
                    pos_base_in_sequence <- str_locate_all(bam_seq, type_base)[[1]][,1] # Supplementary SEQ are truncated without -Y
                    mod_idx <- cumsum(mod_rel_pos + 1)
                }

                skipped_base_handling <- substr(bam_MM_clean[,1], nchar(bam_MM_clean[,1]), nchar(bam_MM_clean[,1]))
                if(skipped_base_handling == "."){  # Attribute p-value to all non-zero position, [.?]?, . means others are at 0.
                    unreported_value <- 0
                }else if(skipped_base_handling == "?"){
                    unreported_value <- NA
                }else{
                    print("Issue")
                }

                res <- data.frame(
                    mod_type=type_mod %>% as.factor(),
                    mod_pos=pos_base_in_sequence,
                    mod_prob=unreported_value
                )

                mod_prob_idx <- seq(length_ML_processed + 1, length_ML_processed + length(mod_rel_pos))
                res$mod_prob[mod_idx] <- as.numeric(list_probabilities[mod_prob_idx]/255)

                length_ML_processed <- length_ML_processed + length(mod_rel_pos)

                return(res %>% dplyr::filter(!is.na(mod_prob))) 
            }else{
                return(NULL)
            }
        }

        if(all(unlist(lapply(bam_mod_info, is.null)))){
            stop(paste0("Incompatible --target_mod_type. Available mod_type(s): ", paste(list_type_mod, collapse=", ")))
        }

        return(do.call(rbind, bam_mod_info) %>% as_tibble())
    }else{
        bam_MM_item <- bam_MM_item %>%
            str_remove(";$") %>%
            str_split(";", simplify=TRUE)

        bam_mod_info <- foreach(idx=seq(1,length(bam_MM_item)), .combine=rbind) %do% {
            bam_MM_clean <- bam_MM_item[,idx] %>% str_split(",", simplify=TRUE)
            
            res <- data.frame(
                mod_type=bam_MM_clean[,1] %>% as.factor(),
                mod_rel_pos=bam_MM_clean[,-1] %>% as.numeric()
            )
            type_base <- substr(bam_MM_clean[,1],1,1)
            res$mod_pos <- str_locate_all(bam_seq, type_base)[[1]][,1][cumsum(res$mod_rel_pos + 1)]

            return(res) 
        }
        bam_mod_info$mod_prob <- as.numeric(unlist(bam_ML_item)/255)

        return(as_tibble(bam_mod_info))
    }
}

parseCigar <- function(cigar, start_pos){
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
    
    return(mappings)
}

process.signal <- function(start, end, strand, cigar, signal_raw, data_type, p, with_smoothing, min_gap, w_pg, ...){
    p(message = "process.signal")

    # Define list of modification type processed
    list_mod_type <- levels(signal_raw$mod_type)

    if(data_type=="DNAscent"){
        # Convert coordinates
        mapping_table <- parseCigar(cigar, start)

        if(strand=="-"){ # Reverse mod_pos from sequenced direction to bam direction
            read_length <- sum(as.integer(regmatches(cigar, gregexpr("\\d+", cigar))[[1]])) # Approx. TODO refine
            signal_raw <- signal_raw %>%
                mutate(mod_pos=read_length - (mod_pos-1))
        }

        signal_gen <- merge(signal_raw, mapping_table, by.x=c("mod_pos"), by.y=c("read_pos"), all.x=TRUE) %>%
            mutate(positions=ref_pos) %>%
            dplyr::select(-c(mod_pos, ref_pos))
    }else{
        signal_gen <- signal_raw %>%
            mutate(positions=case_when(strand=="+" ~ (start - 1 + mod_pos), strand=="-" ~ end + 1 - mod_pos)) %>%
            dplyr::select(-c(mod_pos, mod_rel_pos))
    }
    rm(signal_raw)
    gc()

    # Find gaps
    size_gaps_tmp <- signal_gen %>%
        group_by(mod_type) %>%
        arrange(positions) %>%
        drop_na()

    if(nrow(size_gaps_tmp)==0){
        summary_gaps=tibble(mod_type=list_mod_type, gap_start=0, gap_end=0, gap_width=0)
    }else{
        size_gaps <- size_gaps_tmp %>%
            mutate(distance=lead(positions, default=as.integer(max(positions)+1)) - positions) %>%
            mutate(is_gap_start=ifelse(distance>min_gap, TRUE, FALSE)) %>%
            mutate(is_gap_end=ifelse(lag(is_gap_start, default=FALSE), TRUE, FALSE)) %>%
            filter(is_gap_start | is_gap_end)

        if(nrow(size_gaps)>0){
            summary_gaps=tibble(
                    mod_type=size_gaps %>% filter(is_gap_start) %>% pull(mod_type),
                    gap_start=size_gaps %>% filter(is_gap_start) %>% pull(distance),
                    gap_end=size_gaps %>% filter(is_gap_end) %>% pull(distance)
                ) %>%
                mutate(gap_width=gap_end - gap_start)
        }else{
            summary_gaps=tibble(mod_type=list_mod_type, gap_start=0, gap_end=0, gap_width=0)
        }
    }

    signal_smoothed <- tibble(positions=rep(start:end, times=length(list_mod_type)), mod_type=rep(as.factor(list_mod_type), each=(end + 1) - start)) %>%
        left_join(signal_gen, by=c("positions", "mod_type")) %>%
        ungroup()
    rm(signal_gen)
    gc()

    if(with_smoothing){
        # Fill in gaps with NA
        # & Smooth signal
        signal_smoothed <- signal_smoothed %>%
            group_by(mod_type) %>%
            mutate(smooth_signal_mean=roll_mean(mod_prob, by=1, align="center", n=100, na.rm=TRUE, fill=NA)) %>%
            mutate(smooth_signal_norm=roll_mean(smooth_signal_mean, by=1, align="center", weights=w_pg, normalize=TRUE, na.rm=TRUE, fill=NA)) %>%
            dplyr::select(-smooth_signal_mean) %>%
            # drop_na() %>%
            ungroup()
    }

    out=list(
        # signal_raw=signal_raw,
        signal_smoothed=signal_smoothed,
        summary_gaps=summary_gaps
    )
}

process.bam.batch <- function(bam_batch, min_read_length, min_gap, w_pg, nb_threads, data_type, with_smoothing=TRUE, target_mod_type = NULL){
    print(paste0("  Extract raw signal (n=",nrow(bam_batch),")"))
    current_processed_bam_batch <- bam_batch %>%
        filter(!is.na(pos)) %>% # Implicitely removing secondary and supplementary mappings 
        mutate(start=as.numeric(start)) %>%
        mutate(end=start + qwidth - 1) %>%  # Approx. for DNAscent TODO refine
        filter((end-start) > min_read_length) # Approx. for DNAscent TODO refine
    rm(bam_batch)

    if(nrow(current_processed_bam_batch)==0){
        print(paste0("No read left after min. length filtering."))
        return(NULL)
    }
    
    grouping_vars <- c("read_id", "flag", "chrom", "strand", "start")
    if ("cigar" %in% colnames(current_processed_bam_batch)) {
        grouping_vars <- c(grouping_vars, "cigar")
    }
    current_processed_bam_batch <- current_processed_bam_batch %>%
        mutate(seq = ifelse(strand == "-", stringi::stri_reverse(chartr("ATGC", "TACG", seq)), seq)) %>%
        group_by(across(all_of(grouping_vars))) %>%
        mutate(signal_raw = list(process.mod.tag(pos, prob, seq, data_type, target_mod_type))) %>%
        select(-c(pos, prob)) %>%
        ungroup()

    if(with_smoothing){
        print("  Smooth raw signal")
    }else{
        print("  Process raw signal")       
    }

    with_progress({
        # handlers("txtprogressbar") # Use a simple text progress bar

        p <- progressor(steps=nrow(current_processed_bam_batch))
        current_processed_bam_batch <- current_processed_bam_batch %>%
            mutate(analysis=future_pmap(., process.signal, data_type=data_type, p=p, with_smoothing=with_smoothing, min_gap=min_gap, w_pg=w_pg)) # %>% dplyr::select(-signal_raw)
        # print(gc())
    }, enable = TRUE)

    print("  Order mappings")
    current_processed_bam_batch <- current_processed_bam_batch %>%
        arrange(chrom, start)

    return(current_processed_bam_batch)
}

extract.processed.signal <- function(analysis, read_id, flag, chrom, strand, p, ...){
    data <- analysis$signal_smoothed
    data$read_id <- read_id
    data$flag <- flag
    data$chrom <- chrom
    data$strand <- strand

    p(message = "extract.processed.signal")

    return(data)
}

extract.local.signal <- function(bam_path, contig_name, contig_start, contig_end, model_name, nb_alignments, nb_threads, data_type, min_read_length=5000, with_smoothing=TRUE, target_mod_type = NULL){
    # Setup data processing
    min_gap <- 100
    w_pg <- dnorm(1:2500,1251,300)

    print("Prepare bam query")
    list_what <- c("qname", "rname", "pos", "qwidth", "strand", "seq", "flag")
    if(data_type=="DNAscent"){
        list_what <- c(list_what, "cigar")
    }
    list_bam_queries <- ScanBamParam(what=list_what, tag=c("Mm", "Ml", "MM", "ML"), which=GRanges(contig_name, IRanges(contig_start, contig_end)))

    # Prepare input file
    print("Open bam input")
    bam_file <- BamFile(bam_path)
    yieldSize(bam_file) <- nb_alignments # Cap maximum number of alignment to read
    open(bam_file)

    print("Read bam batch")
    bam_batch <- read.bam.batch(bam_file, list_bam_queries, data_type) # TODO data_type
    nb_alignments_total <- nrow(bam_batch)
    
    # Close bam file handle
    close(bam_file)
    yieldSize(bam_file) <- NA

    print(paste0("  Nb mappings: ",nb_alignments_total))
    if(is.null(nb_alignments_total)){
        return(NULL)
    }
    if(nb_alignments < nb_alignments_total){
        print(paste0("  Random sampling ",nb_alignments," reads"))
        bam_batch <- bam_batch[sample(seq(1,nb_alignments_total), size=nb_alignments, replace=FALSE),]
    }else{
        print(paste0("  No random sampling"))
    }

    print("Process bam batch")
    current_processed_bam_batch <- process.bam.batch(bam_batch, min_read_length, min_gap, w_pg, nb_threads, data_type, with_smoothing, target_mod_type)
    rm(bam_batch)

    if(is.null(current_processed_bam_batch)){
        return(NULL)
    }

    nb_mappings <- nrow(current_processed_bam_batch)
    print(paste0("Extract smoothed signal (n=",nb_mappings,")"))
    with_progress({
        p <- progressor(steps=nb_mappings)
        chunks_data <- future_pmap(current_processed_bam_batch, extract.processed.signal, p=p) %>% rbindlist()
        rm(current_processed_bam_batch)
        gc()
    }, enable = TRUE)

    chunks_data$model <- model_name
    
    return(chunks_data)
}

extract.local.signal.all <- function(path_bam, min_read_length, nb_threads, data_type, with_smoothing = TRUE, target_mod_type = NULL) {
    max_depth <- 1e6 # 10

    # Get contig names and lengths from BAM header
    header_info <- system(paste0("samtools view -H ", path_bam), intern = TRUE)
    list_contigs <- trimws(sub("^@SQ.*SN:(.*)\t.*", "\\1", grep("^@SQ", header_info, value = TRUE)))
    list_lengths <- as.numeric(sub("^@SQ.*LN:(.*)", "\\1", grep("^@SQ", header_info, value = TRUE)))

    # Process each contig in parallel
    processed_chunks_data_list <- foreach(contig_idx = seq_along(list_contigs), .packages = "dplyr") %do% { # dopar  %>% head(2)
        contig <- list_contigs[contig_idx]
        message("Processing: ", contig)  # Provide feedback during processing

        path_details <- str_split(path_bam, "/", simplify = TRUE)
        tmp_name <- paste0("training_eval_", path_details[ncol(path_details)], "_", contig)

        # Extract local signals for the current contig
        processed_chunks_data_tmp <- extract.local.signal(path_bam, contig, 0, list_lengths[contig_idx], tmp_name, max_depth, nb_threads, data_type, min_read_length, with_smoothing, target_mod_type)

        if(!is.null(processed_chunks_data_tmp)){
            processed_chunks_data_tmp <- processed_chunks_data_tmp %>%
                ungroup() %>%
                dplyr::select(-model) %>%
                dplyr::filter(flag %in% c(0, 16)) %>%
                dplyr::select(-flag) %>%
                dplyr::rename(Bprob=mod_prob, signal=smooth_signal_norm) %>%
                mutate(gap_pos=NA) %>%
                group_by(read_id) %>%
                mutate(start = min(positions), end = max(positions)) %>%
                nest(signalr = c(mod_type, positions, Bprob, signal)) %>%
                dplyr::select(read_id, chrom, start, end, strand, gap_pos, signalr)

            return(processed_chunks_data_tmp)
        }

        return(NULL)
    }

    processed_chunks_data_chromosome <- bind_rows(processed_chunks_data_list)
    rm(processed_chunks_data_list)

    if (nrow(processed_chunks_data_chromosome) == 0) {
        return(NULL)  # Return NULL if there's no data
    }

    return(processed_chunks_data_chromosome)
}

split.reads.for.reviewers <- function(list_subsampled_read_id, nb_reviewer, with_independent_review = TRUE) {
    set.seed(101)

    nb_read_id <- length(list_subsampled_read_id)
    nb_read_id_per_reviewer <- floor(nb_read_id / nb_reviewer)
    list_reviewer <- seq(1, nb_reviewer)

    initial_split <- vector("list", nb_reviewer)
    remaining_read_id <- list_subsampled_read_id

    # Initial split of read IDs for each reviewer
    for (i in list_reviewer) {
        initial_split[[i]] <- sample(remaining_read_id, size = nb_read_id_per_reviewer, replace = FALSE)
        remaining_read_id <- setdiff(remaining_read_id, unlist(initial_split[[i]]))
    }

    # Generate a second set of annotation for independent reviewer if needed
    if(with_independent_review){
        second_split <- list()
        for(i in list_reviewer){
            list_non_reviewer <- list_reviewer[!list_reviewer %in% i]
            remaining_read_id <- sample(unlist(initial_split[i]))
            nb_remaining_read_id <- length(remaining_read_id)
            split_size <- ceiling(nb_remaining_read_id / (nb_reviewer-1))

            tmp_split <- list()
            for (j in 1:(nb_reviewer-1)) {
                start_index <- (j - 1) * split_size + 1
                end_index <- min(j * split_size, nb_remaining_read_id)

                tmp_split[list_non_reviewer[j]] <- list(remaining_read_id[start_index:end_index])
            }
            second_split[[i]] <- tmp_split
        }

        for(i in list_reviewer){
            initial_split[i] <- list(c(unlist(initial_split[i]), unlist(lapply(second_split, function(x) x[i]))))
        }
    }

    return(initial_split)
}

local.smoothing <- function(signal){
    signal <- tibble(x=seq(min(signal$x, na.rm=TRUE), max(signal$x, na.rm=TRUE))) %>%
        left_join(signal, by="x") %>%
        group_by(mod_type) %>%
        arrange(x) %>%
        mutate(avg=roll_mean(raw, by=1, align="center", n=100, na.rm=TRUE, fill=NA)) %>%
        drop_na() # Reads shorted than large smoothing are droped
    
    return(signal)
}

write.forkfit.to.annotate <- function(nfs_data, sample_id, output_file, target_mod_type){
    nfs_data <- nfs_data %>%
        dplyr::rename("read_chrm"="chrom", "read_strand"="strand", "read_start"="start", "read_end"="end", "signal"="signalr") %>%
        dplyr::filter(!is.na(read_id))

    if(!is.null(target_mod_type)){ # Make sure it's NULL for Megalodon
        nfs_data <- nfs_data %>%
            mutate(signal=map(signal, ~ dplyr::filter(.x, mod_type==target_mod_type)))        
    }

    nfs_data <- nfs_data %>%
        mutate(signal=map(signal, ~ dplyr::rename(.x, x=positions, y=signal, raw=Bprob))) %>%
        mutate(signal=map(signal, local.smoothing)) %>%
        rowwise() %>%
        mutate(read_inf=min(unlist(signal$x), na.rm=TRUE), read_sup=max(unlist(signal$x), na.rm=TRUE)) %>%
        mutate(sample=sample_id) %>%
        ungroup() %>%
        dplyr::select("read_id","read_chrm","read_strand","read_start","read_end","read_inf","read_sup","signal","sample")

    saveRDS(nfs_data, output_file)
}

write.batches.for.reviewers <- function(df, initial_split, reviewers, output_dir, batch_size, sample_id, target_mod_type) {
    stifle <- foreach(reviewer_id = seq(1, length(reviewers))) %do% {
        subset_reviewer <- df %>% dplyr::filter(read_id %in% initial_split[[reviewer_id]])

        output_reviewer_dir <- paste0(output_dir, "/data_", reviewers[reviewer_id], "/reads_toannotate")
        latest_batch <- 0
        
        if (!dir.exists(output_reviewer_dir)) {
            dir.create(output_reviewer_dir, recursive = TRUE)
        } else {
            list_existing_batches <- list.files(output_reviewer_dir, pattern = ".batch_[0-9]+.rds")
            if (length(list_existing_batches) > 0) {
                latest_batch <- max(as.integer(str_extract(str_extract(list_existing_batches, "_\\d+.rds"), "\\d+")))
            }
        }

        nb_reads <- nrow(subset_reviewer)
        nb_batches <- ceiling(nb_reads / batch_size)
        cat(paste0("Nb batches: ", nb_batches, "\n"))

        stifle <- foreach(batch_id = seq(1, nb_batches)) %do% {
            offseted_batch_id <- batch_id + latest_batch
            row_idx <- seq((batch_id - 1) * batch_size + 1, batch_id * batch_size)
            batch_file_name <- paste0(gsub("/$", "", output_reviewer_dir), "/", sample_id, ".batch_", offseted_batch_id, ".rds")
            write.forkfit.to.annotate(subset_reviewer[row_idx, ], sample_id, batch_file_name, target_mod_type)
            cat(paste0(output_reviewer_dir, " ", batch_id, " ", offseted_batch_id, "\n"))

            return(NA)
        }

        return(NA)
    }
}

wrapper.prepare.doped.annotation <- function(name_sample, path_bam, path_output, name_reviewer, path_list_id_with_signal, nb_random, nb_with_signal, min_read_length, batch_size, nb_threads, data_type, with_independent_review=TRUE, focus="no", path_discard=NULL, target_mod_type=NULL){

    cat(paste0("Define batches sizes.\n"))
    nb_reviewer <- length(name_reviewer)
    if(nb_reviewer==1){
        nb_better_target_read <- nb_random
        nb_better_target_read_with_signal <- nb_with_signal
    }else{
        nb_better_target_read <- find.closest.number(nb_random, nb_reviewer)
        nb_better_target_read_with_signal <- find.closest.number(nb_with_signal, nb_reviewer)
    }

    output_dir <- paste0(path_output,"/",name_sample,"_all_reviewers")
    cat(paste0("Create/define output directory (",output_dir,").\n"))
    if(!dir.exists(output_dir)){
        dir.create(output_dir, recursive=TRUE)
    }

    intermediate_rds <- paste0(output_dir, "/dataset.subset.rds")
    intermediate_txt <- paste0(output_dir, "/dataset.subset.read_id.txt")
    intermediate_ori <- paste0(output_dir, "/dataset.subset.read_id.source.txt")
    intermediate_bam <- paste0(output_dir, "/dataset.subset.sorted.bam")

    if(!file.exists(intermediate_rds)){
        cat(paste0("Generate bam subset.\n"))
        if(!file.exists(intermediate_bam)){
            generate.bam.subset(path_bam, nb_better_target_read, nb_better_target_read_with_signal, min_read_length, path_discard, path_list_id_with_signal, intermediate_bam, intermediate_ori, intermediate_txt)
        }

        df <- extract.local.signal.all(intermediate_bam, min_read_length, nb_threads, data_type, TRUE, target_mod_type)
        saveRDS(df, intermediate_rds)
    }else{
        cat(paste0("Read bam subset.\n"))
        df <- readRDS(intermediate_rds)
    }

    cat(paste0("Split read set (n=",length(df$read_id),").\n"))
    initial_split <- split.reads.for.reviewers(df$read_id, nb_reviewer, with_independent_review)

    cat(paste0("Writing batches\n"))
    write.batches.for.reviewers(df, initial_split, name_reviewer, output_dir, batch_size, name_sample, target_mod_type)
}

options(future.globals.maxSize = 1.0 * 1e9)
plan(multisession, workers=nb_threads)
wrapper.prepare.doped.annotation(name_sample, path_bam, path_output, name_reviewer, path_list_id_with_signal, nb_random, nb_with_signal, min_read_length, batch_size, nb_threads, data_type, with_independent_review, "no", path_discard, target_mod_type)



