#!/usr/bin/env Rscript

library("optparse")
# Option parser setup
option_list <- list(
  make_option(c("-i", "--path_data"), type="character", help="Path annotated data (.rds)"),
  # make_option(c("-a", "--path_annotation"), type="character", help="Path annotation"),
  make_option(c("-m", "--smoothing_mode"), type="character", help="Type of smoothing to perform"),
  make_option(c("-w", "--smoothing_window"), type="integer", default=NA, help="Size smoothing window"),
  make_option(c("-o", "--path_output"), type="character", help="Path output")
)

# Parse options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

path_data <- opt$path_data
smoothing_window <- opt$smoothing_window
smoothing_mode <- opt$smoothing_mode
path_output <- opt$path_output

library("tibble")
library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")

if(smoothing_mode %in% c("by_bin","by_bin_med","by_bin_med2")){
	if(is.na(smoothing_window)){
		stop("Missing --smoothing_window")
	}
}

smoothing.raw.signal <- function(data, smoothing_mode, smoothing_window=NA){
	if(smoothing_mode == "by_bin"){
		results <- data %>%
			mutate(x=round(positions/smoothing_window)*smoothing_window) %>%
			group_by(mod_type, x) %>%
			summarise(y=round(mean(mod_prob/255, na.rm=TRUE), digits=2), .groups='drop')
	}else if(smoothing_mode == "by_bin_hd"){
		results <- data %>%
			mutate(x=round(positions/smoothing_window)*smoothing_window) %>%
			group_by(mod_type, x) %>%
			summarise(y=round(mean(mod_prob/255, na.rm=TRUE), digits=3), .groups='drop')
	}else if(smoothing_mode == "by_bin_centered"){
		results <- data %>%
			mutate(x=round(positions/smoothing_window)*smoothing_window + round(smoothing_window/2)) %>%
			group_by(mod_type, x) %>%
			summarise(y=round(mean(mod_prob/255, na.rm=TRUE), digits=2), .groups='drop')
	}else if(smoothing_mode == "by_bin_med"){
		results <- data %>%
			mutate(x=round(positions/smoothing_window)*smoothing_window) %>%
			group_by(mod_type, x) %>%
			summarise(y=round(median(mod_prob/255, na.rm=TRUE), digits=2), .groups='drop')
	}else if(smoothing_mode == "by_bin_med2"){
		results <- data %>%
			mutate(x=round(positions/smoothing_window)*smoothing_window) %>%
			group_by(mod_type, x) %>%
			summarise(y=nth(mod_prob/255, round(n()/2)), .groups='drop') # Do not handle NA
	}

	return(results)
}

annotate.signal <- function(signal_smooth, signal_annotation){
	signal_annotated <- signal_smooth %>%
		arrange(x) %>%
		mutate(
			fork_right = 0,
			fork_left = 0,
			fork_none = 0
		)

	# Loop through each annotation and update the signal data
	for (i in seq_along(signal_annotation$x0)) {
		x0 <- signal_annotation$x0[i]
		x1 <- signal_annotation$x1[i]
		x2 <- signal_annotation$x2[i]
		fork_dir <- signal_annotation$fork_dir[i]

		if(!is.na(x0)){
			signal_annotated <- signal_annotated %>%
				mutate(
					fork_right = if_else(fork_dir == "right" & x >= x0 & x <= x2, 1, fork_right, 0),
					fork_left = if_else(fork_dir == "left" & x >= x2 & x <= x0, 1, fork_left, 0),
					fork_none = if_else(fork_dir == "unknown" & x >= x0 & x <= x1, 1, fork_none, 0), # Likely not useful because filtered out
					fork_none = if_else(fork_dir == is.na(fork_dir) & x >= x0 & x <= x2, 1, fork_none, 0) # Shouldn't be needed at this point
				)
		}
	}

	# Set default value for no forks and set to integer
	signal_annotated <- signal_annotated %>%
		mutate(fork_none=if_else(fork_right + fork_left == 0, 1, 0)) %>%
  	mutate(across(starts_with("fork_"), as.integer))

	# Remove blindspot from annotation
	size_blindspot <- 15
	signal_annotated <- signal_annotated[((size_blindspot + 1):(nrow(signal_annotated) - size_blindspot)),]

	return(signal_annotated)
}

draw.annotated.signal <- function(processed_chunks_data, path_output){
	pdf(gsub(".rds", ".pdf", path_output), width=15)
	for(i in seq(1, nrow(processed_chunks_data))){
		tmp <- processed_chunks_data[i,] %>%
			unnest(signal_annotated) %>%
			mutate(type=ifelse(fork_right>0,2,ifelse(fork_left>0,3,1))) # Define integer matching left, right, no fork.
		
		gp <- ggplot(tmp) +
			geom_point(aes(x=x, y=y, col=as.factor(type))) +
			theme_bw() +
			scale_color_manual(values=c("1"=1,"2"=2,"3"=3)) +
			scale_y_continuous(limits=c(0,1)) +
			labs(title=paste0(unique(tmp$read_id)), x="Genomic position", y="Annotated signal")
		print(gp)
	}
	dev.off()
}

processed_chunks_data <- readRDS(path_data) %>%
	mutate(signal_smooth=map(signal_raw, ~ smoothing.raw.signal(.x, smoothing_mode, smoothing_window))) %>%
	dplyr::select(-signal_raw) %>%
	dplyr::filter(map_int(signal_smooth, nrow) >= 10000 / smoothing_window) %>% # Keep signal longer than 10 kbp
	mutate(signal_annotated=map2(signal_smooth, annotation, annotate.signal)) %>%
	dplyr::select(-annotation, -signal_smooth)

saveRDS(processed_chunks_data, file=path_output)

draw.annotated.signal(processed_chunks_data, path_output)
