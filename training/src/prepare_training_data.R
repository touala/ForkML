#!/usr/bin/env Rscript

.libPaths(.libPaths()[1])

library("optparse")

# Option parser setup
option_list <- list(
  make_option(c("-b", "--path_bam"), type="character", help="Path mod_mapping.bam"),
  make_option(c("-a", "--path_annotation"), type="character", help="Path annotation"),
  make_option(c("-o", "--path_output"), type="character", help="Path output file"),
  make_option(c("-t", "--nb_threads"), type="integer", default=1, help="Number of threads"),
  make_option(c("-l", "--min_read_length"), type="integer", default=10000, help="Minimum read length"),
  make_option(c("-m", "--nb_alignments"), type="integer", default=2000, help="Maximum alignmentsnumber by contig"),
  # make_option(c("--noise_control"), type = "character", default = NULL, help = "Path to list read_id with noise"),
  make_option(c("--target_mod_type"), type = "character", default = "N+b?", help = "Target modification type (default: 'N+b?')"),
  make_option(c("--data_type"), type = "character", default = "Megalodon", help = "Data type (default: 'Megalodon', or 'DNAscent')")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

path_bam <- opt$path_bam
path_annotation <- opt$path_annotation
path_output <- opt$path_output
nb_threads <- opt$nb_threads
min_read_length <- opt$min_read_length
nb_alignments <- opt$nb_alignments
# noise_control <- opt$noise_control
target_mod_type <- opt$target_mod_type
data_type <- opt$data_type

if(!data_type %in% c("Megalodon", "DNAscent")){quit("data_type unknown.", 2)}
if(data_type=="Megalodon"){target_mod_type <- NULL}

library("Rsamtools")
library("tibble")
library("dplyr")
library("tidyr")
library("stringr")
library("foreach")
library("purrr")
library("progressr")
library("furrr")
library("data.table")

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
                res$mod_prob[mod_idx] <- as.numeric(list_probabilities[mod_prob_idx]) # was /255

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
        bam_mod_info$mod_prob <- as.numeric(unlist(bam_ML_item)) # was /255

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

process.signal <- function(start, end, strand, cigar, signal_raw, data_type, p, ...){
    p(message = "process.signal")

    # Define list of modification type processed
    list_mod_type <- levels(signal_raw$mod_type)

    if(data_type=="DNAscent"){
        # Convert coordinates
        mapping_table <- parseCigar(cigar, start)

        if(strand=="-"){ # Reverse mod_pos from sequenced direction to bam direction
            read_length <- sum(as.integer(regmatches(cigar, gregexpr("\\d+", cigar))[[1]])) # TODO refine
            signal_raw <- signal_raw %>%
                mutate(mod_pos=read_length - (mod_pos-1))
        }

        signal_gen <- merge(signal_raw, mapping_table, by.x=c("mod_pos"), by.y=c("read_pos"), all.x=TRUE) %>%
            mutate(positions=ref_pos) %>%
            dplyr::select(-c(mod_pos, ref_pos))
    }else{
        signal_gen <- signal_raw %>%
            mutate(positions=case_when(strand=="+" ~ (start - 1 + mod_pos), strand=="-" ~ end + 1 - mod_pos)) %>% # TODO refine
            dplyr::select(-c(mod_pos, mod_rel_pos))
    }

    return(as_tibble(signal_gen))
}

process.bam.batch <- function(bam_batch, min_read_length, nb_threads, data_type, target_mod_type = NULL){
    print(paste0("  Extract raw signal (n=",nrow(bam_batch),")"))
    current_processed_bam_batch <- bam_batch %>%
        filter(!is.na(pos)) %>% # Implicitely removing secondary and supplementary mappings 
        dplyr::filter(flag %in% c(0, 16)) %>% # Non-primary mapping are not annotated
        mutate(start=as.numeric(start)) %>%
        mutate(end=start + qwidth - 1) %>% # TODO refine
        filter((end-start) > min_read_length)
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
        select(-c(pos, prob, seq)) %>%
        ungroup()

    print("  Process raw signal")
    with_progress({
        p <- progressor(steps=nrow(current_processed_bam_batch))
        current_processed_bam_batch <- current_processed_bam_batch %>%
            mutate(signal_raw=future_pmap(., process.signal, data_type=data_type, p=p))
    }, enable = TRUE)

    print("  Order mappings")
    current_processed_bam_batch <- current_processed_bam_batch %>%
        arrange(chrom, start)

    return(current_processed_bam_batch)
}

extract.local.signal <- function(bam_path, nb_alignments, nb_threads, data_type, min_read_length=5000, target_mod_type = NULL){
    print("Prepare bam query")
    list_what <- c("qname", "rname", "pos", "qwidth", "strand", "seq", "flag")
    if(data_type=="DNAscent"){
        list_what <- c(list_what, "cigar")
    }
    # list_bam_queries <- ScanBamParam(what=list_what, tag=c("Mm", "Ml", "MM", "ML"), which=GRanges(contig_name, IRanges(contig_start, contig_end)))
    list_bam_queries <- ScanBamParam(what=list_what, tag=c("Mm", "Ml", "MM", "ML"))

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
    current_processed_bam_batch <- process.bam.batch(bam_batch, min_read_length, nb_threads, data_type, target_mod_type)
    rm(bam_batch)

    if(is.null(current_processed_bam_batch)){
        return(NULL)
    }
    
    return(current_processed_bam_batch)
}

extract.local.signal.all <- function(path_bam, min_read_length, nb_alignments, nb_threads, data_type, target_mod_type = NULL) {
    # Get contig names and lengths from BAM header
    processed_chunks_data_chromosome <- extract.local.signal(path_bam, nb_alignments, nb_threads, data_type, min_read_length, target_mod_type)

    if(is.null(processed_chunks_data_chromosome)){
        return(NULL)  # Return NULL if there's no data
    }

    if(nrow(processed_chunks_data_chromosome) == 0){
        return(NULL)  # Return NULL if there's no data
    }

    return(processed_chunks_data_chromosome %>% dplyr::select(read_id, flag, chrom, strand, start, qwidth, end, signal_raw)) # cigar, 
}

processed_chunks_data <- extract.local.signal.all(path_bam, min_read_length, nb_alignments, nb_threads, data_type, target_mod_type)

if(!is.null(processed_chunks_data)){
    manual_annotation_all <- readRDS(path_annotation)

    # if(!is.null(noise_control)){
    #     # Read the list of read_ids from a text file
    #     list_read_id_noise <- readLines(noise_control)

    #     # Get the column names and types from an existing manual_annotation_all$res
    #     template <- manual_annotation_all$res %>% slice(0)

    #     # Create a new tibble with the same columns, filling read_id and leaving others as NA
    #     new_annotations <- template %>%
    #         add_row(!!!set_names(rep(list(NA), length(template)), names(template))) %>%
    #         slice(rep(1, length(list_read_id_noise))) %>%
    #         mutate(read_id=list_read_id_noise)

    #     # Append to manual_annotation_all$res
    #     manual_annotation_all$res <- bind_rows(manual_annotation_all$res, new_annotations)
    # }

    processed_chunks_data <- as_tibble(inner_join(processed_chunks_data, manual_annotation_all$res %>% mutate(read_id=as.character(read_id)) %>% group_by(read_id) %>% nest(.key="annotation"), by="read_id") %>% arrange(read_id))

    saveRDS(processed_chunks_data, file=path_output)    
}
