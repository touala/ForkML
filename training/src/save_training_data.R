#!/usr/bin/env Rscript

.libPaths(.libPaths()[1])

library("optparse")

# Option parser setup
option_list <- list(
  make_option(c("-i", "--input_files"), type="character", default=NULL, help="Files of one-hot encoded annotated data", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Directory for the output ragged tensors", metavar="DIR"),
  make_option("--prop_train", type="numeric", default=0.6, help="Proportion of data for training [default: %default]"),
  make_option("--prop_valid", type="numeric", default=0.2, help="Proportion of data for validation [default: %default]"),
  make_option("--prop_test", type="numeric", default=0.2, help="Proportion of data for test [default: %default]"),
  make_option("--with_augmentation", action="store_true", default=FALSE, help="Perform augmentation of the training set (mirroring only)"),
  make_option("--with_fork_distinction", action="store_true", default=FALSE, help="Keeps distinction between leading/lagging forks"),
  make_option("--nb_pooling", type="integer", default=7, help="Number of downsampling layers in the network [default: %default]"),
  make_option("--max_read_length", type="integer", default=200000, help="Maximum signal/read length, splitted above"),
  make_option("--padding_method", type="character", default="default", help="Choose padding method (rough/clear)"),
  make_option("--for_debug", action="store_true", default=FALSE, help="Run as debugging"),
  make_option(c("-s", "--seed"), type="integer", default=101, help="Random seed")
)
# Parse options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library("tibble")
library("dplyr")
library("tidyr")
library("purrr")
library("reticulate")
library("ggplot2")

if (is.null(opt$input_files) || is.null(opt$output)) {
  cat("Both input files and output directory are required.\n")
  print_help(opt_parser)
  quit(status = 1)
}

# use_python("/usr/bin/python3", required=TRUE)  # Adjust the path as necessary
keras_seq <- import("tensorflow.keras.preprocessing.sequence")
np <- import("numpy")
`%+%` <- paste0

#################### Arguments control ####################
create_dir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir)
    }

    return(dir)
}

is_int <- function(x){
    return(!is.na(suppressWarnings(as.integer(x))))
}

is_dbl <- function(x){
    return(!is.na(suppressWarnings(as.numeric(x))))
}

arg_prop <- function(x){
    if( !is_dbl(x) ){
        stop('Options --prop_train, --prop_valid, and --prop_test need to be doubles, not "' %+% x %+% '"')
    }

    x_num <- as.numeric(x)
    if( !(x_num>=0 & x_num<=1) ){
        stop('Options --prop_train, --prop_valid, and --prop_test need to be in the range [0;1], not "' %+% x_num %+% '"')
    }

    return(x_num)
}

control_sum_prop <- function(prop_train, prop_valid, prop_test){
    # Control that the sum of the proportions is 1
    err_msg = "--prop_train, --prop_test, and --prop_valid do not sum to 1."
    if( (prop_train + prop_valid + prop_test != 1) ){
        stop(err_msg)
    }
}

arg_nb_pooling <- function(x){
    if( !is_int(x) ){
        stop('Option --nb_pooling needs to be an integer, not "' %+% x %+% '"')
    }

    x_num <- as.numeric(x)
    if( !(x_num>=0) ){
        stop('Option --nb_pooling needs to be in the range [0;+inf], not "' %+% x_num %+% '"')
    }
    
    return(x_num)
}

get_params <- function(args){
    p = list()
    p$input_files = args$input_files
    p$output = create_dir(args$output)
    p$prop_train = arg_prop(args$prop_train)
    p$prop_valid = arg_prop(args$prop_valid)
    p$prop_test = arg_prop(args$prop_test)
    control_sum_prop(p$prop_train, p$prop_valid, p$prop_test)

    p$with_augmentation = args$with_augmentation
    p$with_fork_distinction = args$with_fork_distinction
    if(p$with_augmentation & p$with_fork_distinction){
        cat("Both mirroring and fork distinction cannot be set.\n")
        print_help(opt_parser)
        quit(status=1)
    }

    p$nb_pooling = arg_nb_pooling(args$nb_pooling)
    p$max_read_length = args$max_read_length
    p$padding_method = args$padding_method
    p$for_debug = args$for_debug

    return(p)
}

#################### Functions ####################
shuffle_data <- function(df){
    # Shuffle rows of a data frame
    df <- df %>%
        ungroup() %>%  # ungroup just in case
        sample_n(n())  # shuffle

    return(df)
}

compute_padding_length <- function(reads_lengths, nb_pooling, pooling_factor=2){
    # Compute padding length depending on the number of pooling layers
    div_length = (pooling_factor**nb_pooling)
    max_read_length = max(reads_lengths)
    max_pad_length = floor(max_read_length/div_length+1)*div_length
    # print("Number of pooling layers = " %+% nb_pooling)
    # print("Padding length needs to be a multiple of 2**" %+% nb_pooling %+% " = " %+% div_length)
    # print("Max read length = " %+% max_read_length)
    # print("Padding length  = " %+% max_pad_length)
    # print("")
    return(max_pad_length)
}

empty_signals_balancing <- function(data, balancing_factor=0){
    # Flag empty signals
    data <- data %>% mutate(empty=map_lgl(signal_annotated, ~ (sum(.x$fork_right)==0 & sum(.x$fork_left)==0)))

    # Separate empty signals from non empty signals
    data_empty <- data %>% filter(empty)
    data_nonempty <- data %>% filter(!empty)
    
    # Sample the empty reads to have 1/(1+balancing_factor) reads with tracks, and balancing_factor/(1+balancing_factor) empty reads
    N <- nrow(data_nonempty)
    # set.seed(101)
    data_empty <- data_empty %>% sample_n(min(n(), floor(N*balancing_factor)))
    
    # Bind the two back and shuffle
    out <- bind_rows(data_empty, data_nonempty) %>%
        select(-empty) %>%
        shuffle_data()

    return(out)
}

augment_mirroring <- function(data_original){
    # Unnest + mark order of the original reads
    data_original <- data_original %>%
        mutate(signal_annotated=map(signal_annotated, ~ mutate(.x, order=row_number()))) # Define original value order

    # Reverse order + modify read_ids + permute labels of the mirrored reads
    data_mirrored <- data_original %>%
        mutate(read_id = read_id %+% "_mirror") %>%
        mutate(signal_annotated=map(signal_annotated, ~ mutate(.x, tmp_fork_left=fork_left, fork_left=fork_right, fork_right=tmp_fork_left) %>% arrange(desc(order)) %>% dplyr::select(-tmp_fork_left)))
    
    # Bind the two back together + shuffle
    data_augmented <- bind_rows(data_original, data_mirrored) %>%
        shuffle_data()

    return(data_augmented)
}

keep_fork_type_distinction <- function(data_original){
    data_original <- data_original %>%
        mutate(signal_annotated=map(signal_annotated, function(df){
            if(any(strand=="-")){
                df <- df %>%
                    mutate(order=row_number(), tmp_fork_left=fork_left, fork_left=fork_right, fork_right=tmp_fork_left) %>%
                    arrange(desc(order)) %>%
                    select(-tmp_fork_left, -order)
                return(df)
            }else{
                return(df)
            }
        })
    )

    return(data_original)
}

split_by <- function(df, by){
    # Split data frame by column
    df <- df %>%
        nest(data = -by) %>%
        pull(data)

    return(df)
}

make.placeholders <- function(output_name, p){        
    # Create empty placeholders
    empty_array <- array(0, dim = c(1, 1))  # Placeholder array
    empty_list <- list()  # Empty list to mimic the structure
    
    # Save empty signals and labels
    np$save(paste0(p$output, "/", output_name, ".npy"), empty_array)
    np$save(paste0(p$output, "/", gsub("signals", "labels", output_name), ".npy"), empty_array)
    saveRDS(empty_list, paste0(p$output, "/", gsub("signals", "lengths", output_name), ".rds"))
    write.table(empty_list, paste0(p$output, "/", gsub("signals", "read_ids", output_name), ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    if (p$for_debug) {
        return(list(sigs_dataset = empty_array, labs_dataset = empty_array, lengths = empty_list))
    }

    return(NULL)
}

prepare.save.datasets <- function(input_dataset, output_name, p){
    # handle false positive datasets
    if(nrow(input_dataset) == 0){
        message(paste0("Input dataset ",output_name," is empty, maybe for a false positive dataset. Creating empty outputs as placeholders."))

        return(make.placeholders(output_name, p))
    }

    list_read_id <- input_dataset$read_id
    # split into list
    input_dataset = input_dataset %>% pull(signal_annotated)
    # extract signals
    sigs_dataset = input_dataset %>% map(~ select(.x, y)) %>% map(as.matrix) %>% map(unname)
    # extract labels
    labs_dataset = input_dataset %>% map(~ select(.x, fork_right, fork_left, fork_none)) %>% map(as.matrix) %>% map(unname)

    # memorize lengths
    lengths = map_dbl(sigs_dataset, length)

    # compute padding length
    pad_length = compute_padding_length(lengths, nb_pooling=p$nb_pooling, pooling_factor=2) %>% as.integer()

    # combine into array with padding
    if(p$padding_method %in% c("default", "rough")){
        sigs_dataset = keras_seq$pad_sequences(sigs_dataset, dtype="float", padding="post", maxlen=pad_length, value=0.0)
    }else if(p$padding_method == "clean"){
        sigs_dataset = keras_seq$pad_sequences(sigs_dataset, dtype="float", padding="post", maxlen=pad_length, value=-1)
    }
    labs_dataset = keras_seq$pad_sequences(labs_dataset, dtype="int",   padding="post", maxlen=pad_length, value=c(0,0,1))

    # save data
    np$save(paste0(p$output, "/",output_name,".npy"), sigs_dataset)
    np$save(paste0(p$output, "/",gsub("signals","labels",output_name),".npy"), labs_dataset)
    saveRDS(lengths, paste0(p$output,"/",gsub("signals","lengths",output_name),".rds"))
    write.table(list_read_id, paste0(p$output,"/",gsub("signals","read_ids",output_name),".txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)

    if(p$for_debug){
        return(list(sigs_dataset=sigs_dataset, labs_dataset=labs_dataset, lengths=lengths))
    }
}

analyze.threshold.effects <- function(data, threshold_range, split_method="default") {
  results <- data.frame(Threshold = numeric(), Mean = numeric())
  
  for (max_length in threshold_range) {
    adjusted_data <- unlist(sapply(data, function(x) {
      if (x > max_length) {
        if(split_method=="default"){
            div_factor <- ceiling(x/max_length)
            return(rep(x / div_factor, div_factor))  # Split values greater than max_length
        }else if(split_method=="two_stage"){
          div_factor <- ceiling(x/max_length)
          if (x >= max_length*(div_factor-0.5)) {
            return(rep(x / div_factor, div_factor))  # Split values greater than max_length
          }else{
            return(rep(max_length, floor(x/max_length)))
          }
        }else if(split_method=="filled"){
          div_factor <- floor(x/max_length)
          return(rep(max_length, div_factor))  # Split values greater than max_length
        }
      }else{
        return(x)
      }
    }))
    
    # Calculate the score for this threshold
    # score_value <- mean(adjusted_data)
    # score_value <- median(adjusted_data)
    # score_value <- sum(max(adjusted_data)-adjusted_data)
    # score_value <- sum(adjusted_data)/(sum(max(adjusted_data)-adjusted_data)+sum(adjusted_data))
    
    # Append results
    results <- rbind(results, data.frame(Threshold=max_length, q_signal=sum(adjusted_data), q_padding=sum(max(adjusted_data)-adjusted_data), perc_meaningful=sum(adjusted_data)/(sum(max(adjusted_data)-adjusted_data)+sum(adjusted_data))))
  }
  results$split_method <- split_method
  
  return(results)
}

plot.threshold.effect <- function(df){
    threshold_range <- seq(min(df) + 1, max(df), by=10)
    res_default <- analyze.threshold.effects(df, threshold_range, "default")
    res_twostage <- analyze.threshold.effects(df, threshold_range, "two_stage")
    res_filled <- analyze.threshold.effects(df, threshold_range, "filled")

    res_all <- rbind(res_default, res_twostage, res_filled)
    pdf("signal_lengths.processed.pdf")
    hist(df, breaks=100)
    gp <- ggplot(res_all, aes(x=Threshold, y=perc_meaningful, col=split_method)) +
      geom_line() +
      theme_bw() +
      coord_cartesian(ylim=c(0,1)) +
      labs(title="Effect of Max Length Threshold on Score", x="Threshold (max_length)", y="Score of Adjusted Data (% of non-padded data)", col="Split method")
    print(gp)
    gp <- ggplot(res_all %>% group_by(split_method) %>% mutate(pq_signal=q_signal/max(q_signal)), aes(x=Threshold, y=pq_signal, col=split_method)) +
      geom_line() +
      theme_bw() +
      coord_cartesian(ylim=c(0,1)) +
      labs(title="Effect of Max Length Threshold on Score", x="Threshold (max_length)", y="Score of Adjusted Data (% of conserved data)", col="Split method")
    print(gp)
    gp <- ggplot(res_all %>% group_by(split_method) %>% mutate(pq_padding=q_padding/max(q_padding)), aes(x=Threshold, y=pq_padding, col=split_method)) +
      geom_line() +
      theme_bw() +
      coord_cartesian(ylim=c(0,1)) +
      labs(title="Effect of Max Length Threshold on Score", x="Threshold (max_length)", y="Score of Adjusted Data (% of conserved padding)", col="Split method")
    print(gp)
    gp <- ggplot(res_all, aes(x=Threshold, y=perc_meaningful, lty=split_method)) +
      geom_line() +
      geom_line(data=res_all %>% group_by(split_method) %>% mutate(pscore=q_signal/max(q_signal)), aes(x=Threshold, y=pscore, lty=split_method), col=2) +
      geom_line(data=res_all %>% group_by(split_method) %>% mutate(pscore=q_padding/max(q_padding)), aes(x=Threshold, y=pscore, lty=split_method), col=3) +
      coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      labs(title="All metrics", x="Threshold (max_length)", y="Mixed score", lty="Split method")
    print(gp)
    dev.off()
}

split_tibble <- function(data, max_read_length) {
    data_range <- range(data$x) # Independent to smoothing choice
    # data_range <- nrow(data) Assumed no gaps and depends on smoothing
    data_length <- data_range[2] - data_range[1] + 1

    if (data_length > max_read_length) {
        quotient <- data_length %/% max_read_length
        group_size <- data_length %/% (quotient + 1)
        # extra_groups <- data_length %% (quotient + 1)

        start_indices <- seq(1, by = group_size, length.out = quotient + 1)
        end_indices <- c(start_indices[-1] - 1, data_length)

        splitted_data <- lapply(seq_along(start_indices), function(i) {
            range_subset <- start_indices[i]:end_indices[i]
            
            data %>%
                mutate(rel_x=x - (min(x) - 1)) %>%
                dplyr::filter(rel_x >= start_indices[i] & rel_x <= end_indices[i]) %>%
                dplyr::select(-rel_x)
        })

        return(splitted_data)
    } else {
        return(list(data))
    }
}

#### Arguments
p <- get_params(opt)

#### Load data
data <- bind_rows(lapply(strsplit(p$input_files, ",", fixed=TRUE)[[1]], readRDS)) %>%
    arrange(read_id)

#### For padding distinction if not with out of original range
if(p$padding_method=="rough"){
    data <- data %>% mutate(signal_annotated=map(signal_annotated, ~ rowwise(.x) %>% mutate(y=ifelse(y==0, 0.001, y))))
}

#### Shuffle rows
set.seed(opt$seed)
data <- shuffle_data(data) # Reproducible if input in same order?

#### Plot read length properties
plot.threshold.effect(data %>% mutate(nb_points=map_int(signal_annotated, nrow)) %>% pull(nb_points))

#### Split long signal to limit padding
data <- data %>%
    mutate(split_signal_annotated=map(signal_annotated, split_tibble, max_read_length=p$max_read_length)) %>%
    dplyr::select(-signal_annotated) %>%
    unnest(split_signal_annotated) %>%
    group_by(read_id, chrom, strand, start, end) %>%
    mutate(part_index=row_number(), signal_annotated=split_signal_annotated) %>%
    dplyr::select(-split_signal_annotated) %>%
    ungroup()

#### Split into train, validation and test sets
N = nrow(data)
N_train = floor(N * p$prop_train)
N_valid = floor(N * p$prop_valid)
N_test  = floor(N * p$prop_test)

data_train = data %>% slice(1:N_train)
data_valid = data %>% slice(N_train + 1:N_valid)
data_test  = data %>% slice(N_train+N_valid + 1:N_test)
rm(data)

#### Empty signals balancing
data_train_balanced0 = empty_signals_balancing(data_train, 0)
data_train_balanced3 = empty_signals_balancing(data_train, 2.33)

#### Data augmentation
if(p$with_augmentation){
    data_train_balanced0 = augment_mirroring(data_train_balanced0)
    data_train_balanced3 = augment_mirroring(data_train_balanced3)
    data_train = augment_mirroring(data_train)
}

if(p$with_fork_distinction){
    data_train_balanced0 = keep_fork_type_distinction(data_train_balanced0)
    data_train_balanced3 = keep_fork_type_distinction(data_train_balanced3)
    data_train = keep_fork_type_distinction(data_train)
    data_valid = keep_fork_type_distinction(data_valid)
    data_test = keep_fork_type_distinction(data_test)
}

debug_data <- prepare.save.datasets(data_train_balanced0, "train_signals_balanced_0xempty", p)
debug_data <- prepare.save.datasets(data_train_balanced3, "train_signals_balanced_3xempty", p)
debug_data <- prepare.save.datasets(data_train, "train_signals_nonbalanced", p)
debug_data <- prepare.save.datasets(data_valid, "valid_signals", p)
debug_data <- prepare.save.datasets(data_test, "test_signals", p)