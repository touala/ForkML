#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
	make_option(c("-s", "--sample"), type = "character", help = "Sample name (e.g., VR53)"),
	make_option(c("-i", "--input_dir"), type = "character", help = "Input directory"),
	make_option(c("-o", "--output_dir"), type = "character", help = "Output directory for results"),
	make_option("--legacy", action = "store_true", default = FALSE, help = "Run in legacy mode"),
	make_option(c("-d", "--pulse_duration"), type = "integer", default = 30, help = "Time between pulses (in minutes) [default: %default]"),
	make_option(c("-p", "--plot"), action = "store_true", default = FALSE, help = "Whether to produce PDF plots [default: %default]"),
	make_option(c("-b", "--bin_size"), type = "integer", default = 40000, help = "Bin size for read length analysis [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

library(foreach)
library(tibble)
library(dplyr)
library(gtools)
library(ggplot2)
library(knitr)
library(tidyr)
library(purrr)
library(readr)

load.forkml.results <- function(sample_name, input_directory){
	res <- foreach(current_file = list.files(path=gsub("<sample>", sample_name, input_directory), pattern="tracks_[^fs]+", full.names=TRUE, recursive=TRUE), .combine = rbind) %do% {
		# print(current_file)
		forkml_track <- read.table(current_file, header = TRUE, stringsAsFactors = TRUE, sep="\t")
		forkml_track$sample <- sample_name

		return(as_tibble(forkml_track))
	}
	res$sample <- as.factor(res$sample)

	return(res)
}

load.forkml.signals <- function(sample_name, input_directory){
	res <- foreach(current_file = list.files(path=gsub("<sample>", sample_name, input_directory), pattern="signals_.*.tsv.gz", full.names=TRUE, recursive=TRUE), .combine = rbind) %do% {
		# print(current_file)
		forkml_signal <- read.table(current_file, header = TRUE, stringsAsFactors = TRUE, sep="\t") %>%
			group_by(read_id, type) %>%
			# mutate(nb_bins=max(x)) %>%
			mutate(nb_bins=diff(range(x)) + median(diff(x))) %>%
			group_by(read_id, type, nb_bins) %>%
			nest()
		forkml_signal$sample <- sample_name

		return(as_tibble(forkml_signal))
	}
	res$sample <- as.factor(res$sample)

	return(res)
}

filter.out.suspicious.signal <- function(forkml_raw, path_score, score_th=0.4, capped_score_th=0.2){
	list_read_id_good <- read.signal.score(path_score, capped_score_th) %>%
		mutate(mod_type=toupper(mod_type)) %>%
		dplyr::filter(grepl("B", mod_type)) %>%
		dplyr::filter(capped_signal_score >= score_th) %>%
		mutate(read_id=paste0("read_", read_id)) %>%
		distinct(read_id) %>%
		pull(read_id)

	forkml_raw <- forkml_raw %>% dplyr::filter(read_id %in% list_read_id_good)

	return(forkml_raw)
}

read.signal.score <- function(path_file, threshold_score){
	signal_score <- read.table(path_file, header=TRUE, sep='\t', stringsAsFactors=TRUE) %>%
		as_tibble() %>%
		rowwise() %>%
		mutate(capped_signal_score=ifelse(signal_score < threshold_score, threshold_score, signal_score))
	signal_score$chromosome <- factor(signal_score$chromosome, levels=mixedsort(unique(as.character(signal_score$chromosome))))

	return(signal_score)
}

refine.forkml.prediction <- function(information_raw, advanced=FALSE){
	min_iso_fork_size <- 2500
	max_merging_dist <- 1500
	max_merging_size <- 5000
	max_dist_border <- 500

	information_raw_refined <- information_raw %>%
		group_by(read_id, read_chrm, read_strand, read_start, read_end, read_inf, read_sup, sample) %>%
		arrange(pred_inf) %>%
		mutate(detected_length=pred_sup-pred_inf, dist_next=lead(pred_inf) - pred_sup, dist_prev=pred_inf - lag(pred_sup)) %>%
		mutate(dist_next=ifelse(is.na(dist_next), max_merging_dist+1, dist_next), dist_prev=ifelse(is.na(dist_prev), max_merging_dist+1, dist_prev)) %>%
		mutate(is_short=detected_length < min_iso_fork_size, is_close_next=dist_next < max_merging_dist, is_close_prev=dist_prev < max_merging_dist) %>%
		mutate(is_codir_next=direction==lead(direction), is_codir_prev=direction==lag(direction)) %>%
		mutate(is_codir_next=ifelse(is.na(is_codir_next), FALSE, is_codir_next), is_codir_prev=ifelse(is.na(is_codir_prev), FALSE, is_codir_prev)) %>%
		mutate(is_border=((pred_inf - read_inf) < max_dist_border) | ((read_sup - pred_sup) < max_dist_border)) %>%
		# mutate(is_mergeable=(is_close_next & is_codir_next) | (is_close_prev & is_codir_prev)) %>%
		mutate(is_mergeable=(is_close_next & is_codir_next & ( detected_length < max_merging_size | lead(detected_length) < max_merging_size )) | (is_close_prev & is_codir_prev & ( detected_length < max_merging_size | lag(detected_length) < max_merging_size ))) %>%
		mutate(merging_grp=if_else(is_mergeable, if_else(lag(is_mergeable) & is_close_prev & is_codir_prev, 0, 1), 1)) %>%
		mutate(merging_grp=if_else(is.na(merging_grp), 1, merging_grp)) %>%
		mutate(merging_grp=cumsum(merging_grp)) %>%
		mutate(is_discardable=(is_short & !is_close_next & !is_mergeable & !is_border) & (is_short & !is_close_prev & !is_mergeable & !is_border)) %>% # print(width=Inf)
		dplyr::filter(!is_discardable) %>%
		group_by(sample, read_id, read_chrm, read_strand, read_start, read_end, read_inf, read_sup, merging_grp) %>%
		summarize(pred_inf=first(pred_inf), pred_sup=last(pred_sup), direction=first(direction), .groups="drop") %>%
		group_by(sample, read_id, read_chrm, read_strand, read_start, read_end, read_inf, read_sup) %>%
		mutate(is_codir_next=direction==lead(direction), is_codir_prev=direction==lag(direction)) %>%
		mutate(is_codir_next=ifelse(is.na(is_codir_next), FALSE, is_codir_next), is_codir_prev=ifelse(is.na(is_codir_prev), FALSE, is_codir_prev)) %>%
		mutate(pred_mid_inf=ifelse(is_codir_prev, pred_inf - (pred_inf - lag(pred_sup))/2, pred_inf), pred_mid_sup=ifelse(is_codir_next, pred_sup + (lead(pred_inf) - pred_sup)/2, pred_sup)) %>%
		dplyr::select(-is_codir_next, -is_codir_prev) %>%
		ungroup()

	if(advanced){
		min_iso_fork_size <- 1000 # Update

		information_raw_refined <- information_raw_refined %>%
			group_by(read_id, read_chrm, read_strand, read_start, read_end, read_inf, read_sup, sample) %>%
			arrange(pred_inf) %>%
			mutate(detected_length=pred_sup-pred_inf, dist_next=lead(pred_inf) - pred_sup, dist_prev=pred_inf - lag(pred_sup)) %>%
			mutate(dist_next=ifelse(is.na(dist_next), max_merging_dist+1, dist_next), dist_prev=ifelse(is.na(dist_prev), max_merging_dist+1, dist_prev)) %>%
			mutate(is_short=detected_length < min_iso_fork_size, is_close_next=dist_next < max_merging_dist, is_close_prev=dist_prev < max_merging_dist)  %>%
			mutate(is_codir_next=direction==lead(direction), is_codir_prev=direction==lag(direction)) %>%
			mutate(is_codir_next=ifelse(is.na(is_codir_next), FALSE, is_codir_next), is_codir_prev=ifelse(is.na(is_codir_prev), FALSE, is_codir_prev)) %>%
			mutate(is_mergeable=(is_close_next & is_codir_next & ( detected_length < max_merging_size | lead(detected_length) < max_merging_size )) | (is_close_prev & is_codir_prev & ( detected_length < max_merging_size | lag(detected_length) < max_merging_size ))) %>%
			mutate(is_discardable=(is_short & !is_mergeable) & (is_short & !is_mergeable)) %>%

			mutate(direction=if_else(is_discardable, direction * -1, direction)) %>%

			mutate(detected_length=pred_sup-pred_inf, dist_next=lead(pred_inf) - pred_sup, dist_prev=pred_inf - lag(pred_sup)) %>%
			mutate(dist_next=ifelse(is.na(dist_next), max_merging_dist+1, dist_next), dist_prev=ifelse(is.na(dist_prev), max_merging_dist+1, dist_prev)) %>%
			mutate(is_short=detected_length < min_iso_fork_size, is_close_next=dist_next < max_merging_dist, is_close_prev=dist_prev < max_merging_dist)  %>%
			mutate(is_codir_next=direction==lead(direction), is_codir_prev=direction==lag(direction)) %>%
			mutate(is_codir_next=ifelse(is.na(is_codir_next), FALSE, is_codir_next), is_codir_prev=ifelse(is.na(is_codir_prev), FALSE, is_codir_prev)) %>%
			mutate(is_border=((pred_inf - read_inf) < max_dist_border) | ((read_sup - pred_sup) < max_dist_border)) %>%
			mutate(is_mergeable=(is_close_next & is_codir_next & ( detected_length < max_merging_size | lead(detected_length) < max_merging_size )) | (is_close_prev & is_codir_prev & ( detected_length < max_merging_size | lag(detected_length) < max_merging_size ))) %>%
			mutate(merging_grp=if_else(is_mergeable, if_else(lag(is_mergeable) & is_close_prev & is_codir_prev, 0, 1), 1)) %>%
			mutate(merging_grp=if_else(is.na(merging_grp), 1, merging_grp)) %>%
			mutate(merging_grp=cumsum(merging_grp)) %>%

			mutate(is_discardable=(is_short & !is_mergeable) & (is_short & !is_mergeable)) %>% #print(width=Inf)
			dplyr::filter(!is_discardable) %>%

			group_by(sample, read_id, read_chrm, read_strand, read_start, read_end, read_inf, read_sup, merging_grp) %>%
			summarize(pred_inf=first(pred_inf), pred_sup=last(pred_sup), direction=first(direction), .groups="drop") %>%
			group_by(sample, read_id, read_chrm, read_strand, read_start, read_end, read_inf, read_sup) %>%
			mutate(is_codir_next=direction==lead(direction), is_codir_prev=direction==lag(direction)) %>%
			mutate(is_codir_next=ifelse(is.na(is_codir_next), FALSE, is_codir_next), is_codir_prev=ifelse(is.na(is_codir_prev), FALSE, is_codir_prev)) %>%
			mutate(pred_mid_inf=ifelse(is_codir_prev, pred_inf - (pred_inf - lag(pred_sup))/2, pred_inf), pred_mid_sup=ifelse(is_codir_next, pred_sup + (lead(pred_inf) - pred_sup)/2, pred_sup)) %>%
			dplyr::select(-is_codir_next, -is_codir_prev) %>%
			ungroup()
	}

	return(information_raw_refined)
}

detect.speeds.forkml <- function(forkml_annotation, pulse_duration, debug=FALSE){
	min_distance_initiation <- 2000 # TODO maybe process events to polish up
	maximum_fork_distance <- 200000
	amplitude_threshold <- 0.1
	border_distance_threshold <- 2000
	border_distance_relative_threshold <- 0.01
	read_length_threshold <- 60000
	relative_position_threshold <- 0.20

	if(!"person" %in% colnames(forkml_annotation)){
		forkml_annotation <- forkml_annotation %>% mutate(person="none")
	}
	information_speed <- forkml_annotation %>%
		rowwise() %>%
		mutate(read_length=read_end - read_start) %>%
		mutate(fork_dir=ifelse(direction == 1, "right", ifelse(direction == -1, "left", "unknown"))) %>%
		mutate(x0=ifelse(fork_dir=="right", pred_inf, pred_sup), x2=ifelse(fork_dir=="right", pred_sup, pred_inf)) %>%
		mutate(xrs=min(x0, x2), yrs=amplitude_threshold, y0=0, y2=0) %>% # TODO temporary
		mutate(dist_x0_rs=abs(x0-read_start), dist_x2_rs=abs(x2-read_start), min_dist_rs=min(dist_x0_rs, dist_x2_rs)) %>%
		mutate(perc_rs=min(dist_x0_rs/read_length, dist_x2_rs/read_length)) %>%
		mutate(xre=max(x0, x2), yre=amplitude_threshold) %>% # TODO temporary
		mutate(dist_x0_re=abs(x0-read_end), dist_x2_re=abs(x2-read_end), min_dist_re=min(dist_x0_re, dist_x2_re)) %>%
		mutate(perc_re=min(dist_x0_re/read_length, dist_x2_re/read_length)) %>%
		mutate(on_lb=ifelse(read_length < read_length_threshold, 
			ifelse(min_dist_rs < border_distance_threshold & yrs >= amplitude_threshold, TRUE, FALSE), 
			ifelse((min_dist_rs < border_distance_threshold | perc_rs < border_distance_relative_threshold) & yrs >= amplitude_threshold, TRUE, FALSE))) %>%
		mutate(on_rb=ifelse(read_length < read_length_threshold,
			ifelse(min_dist_re < border_distance_threshold && yre >= amplitude_threshold, TRUE, FALSE),
			ifelse((min_dist_re < border_distance_threshold || perc_re < border_distance_relative_threshold) && yre >= amplitude_threshold, TRUE, FALSE))) %>%
		mutate(lpos=min(pred_inf, pred_sup), rpos=max(pred_inf, pred_sup)) %>%
		group_by(sample, person, read_chrm, read_id) %>%
		arrange(sample, person, read_chrm, read_id, pred_inf) %>%
		mutate(dist_next_fork=lead(lpos) - rpos, dist_prev_fork=lpos - lag(rpos)) %>%
		mutate(step=ifelse(lag(fork_dir)==fork_dir & dist_prev_fork < maximum_fork_distance,0,1)) %>%
		mutate(step=ifelse(is.na(step),1,step)) %>%
		mutate(grp=cumsum(step)) %>%

		arrange(read_chrm, read_id, grp, fork_dir) %>%
		group_by(sample, person, read_chrm, read_id, grp, fork_dir) %>%
		mutate(fork_speed=abs(x0-lead(x0))/pulse_duration, x0_first_fork=x0, x0_second_fork=lead(x0)) %>%
		group_by(sample, person, read_chrm, read_id) %>%
		mutate(fork_speed=ifelse(fork_dir=="unknown", 
			NA, #"A", 
			ifelse(fork_dir=="left" & !is.na(lead(fork_dir)), 
				ifelse(fork_dir=="left" & lead(fork_dir)=="left" & !is.na(lead(fork_dir,2)),
					ifelse(fork_dir=="left" & lead(fork_dir)=="left" & lead(fork_dir,2)=="right" & abs(lead(x0)-lead(x0,2)) < min_distance_initiation,
						NA, #"B",
						fork_speed
					),
					ifelse(lead(on_rb),NA,fork_speed)
				), 
				ifelse(fork_dir=="right" & !is.na(lag(fork_dir)), 
					ifelse(fork_dir=="right" & lag(fork_dir)=="left" & abs(x0-lag(x0)) < min_distance_initiation, 
						NA, #"C", 
						fork_speed
					),
					ifelse(on_lb,NA,fork_speed)
				)
			)
		)) %>%
		rowwise() %>%
		mutate(fork_type=ifelse(read_strand=="+",ifelse(direction==1,"leading","lagging"),ifelse(direction==1,"lagging","leading"))) # New

	information_speed <- information_speed %>%
		group_by(sample, person, read_chrm, read_id, grp, fork_dir) %>%
		arrange(sample, person, read_chrm, read_id, grp, fork_dir, pred_inf) %>%
		mutate(pulse_rel_idx=row_number(), pulse_rel_idx=ifelse(direction==1, pulse_rel_idx, max(pulse_rel_idx) - pulse_rel_idx + 1)) # New

	if(!debug){
		information_speed <- information_speed %>% dplyr::filter(!is.na(fork_speed))
	}

	return(information_speed %>% ungroup() %>% dplyr::select(-c(merging_grp, x0, x2, y0, y2, dist_x0_rs, dist_x2_rs, min_dist_rs, perc_rs, xre, yre, dist_x0_re, dist_x2_re, min_dist_re, perc_re, on_lb, on_rb, lpos, rpos, dist_next_fork, dist_prev_fork, step, grp, pulse_rel_idx)))
}

detect.events.forkml <- function(forkml_annotation){
	if(!"person" %in% colnames(forkml_annotation)){
		forkml_annotation <- forkml_annotation %>% mutate(person="none")
	}
	if(!"sample" %in% colnames(forkml_annotation)){
		forkml_annotation <- forkml_annotation %>% mutate(sample="none")
	}
	if(!"haplotype" %in% colnames(forkml_annotation)){
		forkml_annotation <- forkml_annotation %>% mutate(haplotype="No_tag")
	}
	information_event <- forkml_annotation %>%
		# dplyr::filter(read_chrm=="chr18") %>%
		mutate(read_length=read_end - (read_start - 1)) %>%
		dplyr::select(sample, haplotype, read_id, read_chrm, read_strand, read_length, person, direction, pred_inf, pred_sup) %>%
		group_by(sample, haplotype, read_id, read_chrm, read_strand, read_length, person) %>%
		arrange(sample, haplotype, read_id, read_chrm, read_strand, read_length, person, pred_inf) %>%
		mutate(direction=as.integer(direction)) %>% # Forces format
		mutate(next_event=lead(direction) != direction, next_event=ifelse(is.na(next_event), FALSE, next_event)) %>%
		dplyr::filter(any(next_event)) %>% # Keep reads with events
		mutate(
			type_event=ifelse(next_event,ifelse(direction>lead(direction),"Ter","Ini"),NA), 
			pos_event=ifelse(next_event, (pred_sup + lead(pred_inf))/2,NA)
		) %>%
		mutate(pos_event=as.integer(round(pos_event))) %>% # Forces format
		dplyr::filter(!is.na(type_event)) %>%
		ungroup() %>%
		dplyr::select(-c(direction, next_event, person))
	information_event$read_chrm <- factor(information_event$read_chrm, levels=mixedsort(unique(as.character(information_event$read_chrm))))

	return(information_event %>% dplyr::select(-haplotype))
}

extract.input.stats <- function(path_stats_file){
	fields_of_interest <- c("raw total sequences", "reads mapped", "average length", "bases mapped (cigar)")

	input_stats <- tibble(raw = readLines(path_stats_file)) %>%
		dplyr::filter(grepl(":", raw)) %>%
		mutate(
			key = sub(":.*", "", raw) %>% trimws(),
			value = sub(".*:\t", "", raw) %>% sub("\t.*$", "", .) %>% parse_number()
		) %>%
		dplyr::select(key, value) %>%
		dplyr::filter(key %in% fields_of_interest)

	return(input_stats)
}

save.datasets <- function(forkml_refined, forkml_speeds, forkml_events, input_directory, sample_name){
	if (!dir.exists(input_directory)) {
		dir.create(input_directory, recursive = TRUE)
	}

	write.table(forkml_refined, file=paste0(input_directory,"/ForkML.forks.",sample_name,".tsv"), quote=FALSE, sep="\t", row.names=FALSE)
	write.table(forkml_speeds, file=paste0(input_directory,"/ForkML.speeds.",sample_name,".tsv"), quote=FALSE, sep="\t", row.names=FALSE)
	write.table(forkml_events, file=paste0(input_directory,"/ForkML.events.",sample_name,".tsv"), quote=FALSE, sep="\t", row.names=FALSE)
}

analyze.forkml.speed <- function(forkml_refined, forkml_speeds, forkml_events, forkml_signal, input_stats, bin_size = 40000, min_read_length = 90000, plot = TRUE, output_pdf=NULL){
    
    if (!is.null(output_pdf)) {
        output_txt <- gsub(".pdf$", ".txt", gsub("ForkML.speeds_distribution.","ForkML.basic_results_summary.",output_pdf))
        sink(output_txt, split = TRUE)  # Duplicate output to file and terminal
        on.exit(sink(), add = TRUE)     # Ensure it stops when function exits
    }

	cat("🔍 ForkML Speed Analysis\n")
	cat("-------------------------\n\n")

	# Number of reads with forks
	list_read_id_w_forks <- forkml_refined %>%
		distinct(read_id) %>%
		pull(read_id) %>%
		as.character()

	nb_reads_with_forks <- length(list_read_id_w_forks)

	# Number of forks
	nb_fitted_forks <- forkml_refined %>% nrow()

	input_detection_stats <- forkml_signal %>%
		dplyr::filter(read_id %in% list_read_id_w_forks) %>%
		ungroup() %>%
		summarize(n=n(), yield=sum(nb_bins), mean=round(mean(nb_bins), 0), median=round(median(nb_bins), 0), sd=round(sd(nb_bins), 0))

	cat("🧮 Input Summary:\n")
	cat(sprintf("  • Input yield (Mb)             : %d\n", round(input_stats$value[3]/1000000, 0)))
	cat(sprintf("  • Number of reads              : %d\n", input_stats$value[1]))
	cat(sprintf("  • Number of mapped reads       : %d\n", input_stats$value[2]))
	cat(sprintf("  • Average read length          : %d\n", input_stats$value[4]))
	cat("\n")

	cat("🧮 Fork Detection Summary:\n")
	cat(sprintf("  • Number of forks              : %d\n", nb_fitted_forks))
	cat(sprintf("  • Input yield w/ forks (Mb)    : %d (%.2f%%)\n", round(input_detection_stats$yield/1000000, 0), round((input_detection_stats$yield * 100)/input_stats$value[3], 2)))
	cat(sprintf("  • Number of reads w/ forks     : %d (%.2f%%)\n", nb_reads_with_forks, round((nb_reads_with_forks * 100)/input_stats$value[2], 2)))
	cat(sprintf("  • Read length w/ forks:\n"))
	cat(sprintf("    • Mean                       : %d\n", input_detection_stats$mean))
	cat(sprintf("    • Median                     : %d\n", input_detection_stats$median))
	cat(sprintf("    • SD                         : %d\n", input_detection_stats$sd))
	cat("\n\n")

	# Summary statistics
	summary_stats <- forkml_speeds %>%
		filter(!is.na(fork_speed)) %>%
		summarise(Mean = round(mean(fork_speed), 1), Median = round(median(fork_speed), 1), SD = round(sd(fork_speed), 1), N = n())
	summary_stats_long <- forkml_speeds %>%
		filter(read_length >= min_read_length) %>%
		filter(!is.na(fork_speed)) %>%
		summarise(Mean = round(mean(fork_speed), 1), Median = round(median(fork_speed), 1), SD = round(sd(fork_speed), 1), N = n())

	cat("📈 Fork Speed (overall - biased):")
	print(kable(summary_stats, format = "simple"))
	cat("\n\n")
	cat(paste0("📈 Fork Speed (>= ",min_read_length,"):"))
	print(kable(summary_stats_long, format = "simple"))
	cat("\n\n")

	# Binning
	res1_binned <- forkml_speeds %>%
		mutate(read_length_bin = floor(read_length / bin_size) * bin_size)

	fork_speed_by_bin <- res1_binned %>%
		group_by(read_length_bin) %>%
		summarise(Mean = round(mean(fork_speed), 1), Median = round(median(fork_speed), 1), Count = n()) %>%
		arrange(read_length_bin) %>%
		rename(`Read Length Bin`=read_length_bin)

	cat("📊 Fork Speed by Read Length Bin:")
	print(kable(fork_speed_by_bin, format = "simple"))
	cat("\n")

	# Event types
	nb_initiations <- forkml_events %>% filter(type_event == "Ini") %>% nrow()
	nb_terminations <- forkml_events %>% filter(type_event == "Ter") %>% nrow()

	cat("🧬 Event Counts:\n")
	cat(sprintf("  • Initiations (Ini) : %d\n", nb_initiations))
	cat(sprintf("  • Terminations (Ter): %d\n", nb_terminations))
	cat("\n")

	# Optional plot
	if (plot) {
		pdf(output_pdf, width=6, height=4)
		gp <- ggplot(res1_binned, aes(x = as.factor(read_length_bin), y = fork_speed)) +
			geom_boxplot(outliers = FALSE) +
			geom_jitter(height = 0, width = 0.3, pch = 16, size = 0.2, col="grey") +
			stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "red", color = "red") +
			theme_bw() +
			labs(title = paste0("Fork speed distribution - ",levels(forkml_signal$sample)), subtitle="By Read Length Bin", x = "Read Length Bin (bp)", y = "Fork Speed (bp/min)") +
			theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
			scale_x_discrete(labels = function(x) scales::comma(as.numeric(x)))
		print(gp)
		gp <- ggplot(res1_binned %>% rowwise() %>% mutate(capped_fork_speed=min(fork_speed, 3000)), aes(x = as.factor(read_length_bin), y = capped_fork_speed)) +
			geom_boxplot(outliers = FALSE) +
			geom_jitter(height = 0, width = 0.3, pch = 16, size = 0.2, col="grey") +
			stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "red", color = "red") +
			theme_bw() +
			labs(title = paste0("Fork speed distribution - ",levels(forkml_signal$sample)), subtitle="By Read Length Bin", x = "Read Length Bin (bp)", y = "Fork Speed (capped at 3000 bp/min)") +
			theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
			scale_x_discrete(labels = function(x) scales::comma(as.numeric(x)))
		print(gp)
		dev.off()
	}
}

plot.forkml.example <- function(forkml_signal, forkml_raw, forkml_refined, forkml_speeds, forkml_events, output_pdf = "reads_plot_with_arrows.pdf", min_read_length = 90000) {
	forkml_signal <- forkml_signal %>%
		dplyr::filter(nb_bins >= min_read_length)

	list_read_id_w_speed <- forkml_speeds %>%
		mutate(median_fork_speed=median(fork_speed), q1_fork_speed=quantile(fork_speed, 0.25), q3_fork_speed=quantile(fork_speed, 0.75)) %>%
		dplyr::filter((fork_speed > q1_fork_speed) & (fork_speed < q3_fork_speed)) %>%
		distinct(read_id) %>%
		pull(read_id) %>%
		as.character()

	forkml_signal <- forkml_signal %>%
		dplyr::filter(read_id %in% list_read_id_w_speed)

	set.seed(101)
	forkml_signal <- forkml_signal[sample(1:nrow(forkml_signal)), ]

	pdf(output_pdf, width = 12, height = 4)
	walk(1:min(50, nrow(forkml_signal)), function(i) {
		row <- forkml_signal[i, ]
		data_tbl <- row$data[[1]]
		read_id <- as.character(row$read_id)
		read_length <- diff(range(data_tbl$x), na.rm = TRUE)

		# Match annotations for this read_id
		ann_refined <- forkml_refined %>% filter(read_id == !!read_id)
		ann_speed <- forkml_speeds %>% filter(read_id == !!read_id)
		ann_event <- forkml_events %>% filter(read_id == !!read_id)

		# Build plot
		p <- ggplot(data_tbl, aes(x = x, y = y)) +
			geom_line() +
			theme_bw() +
			labs( title = paste0(read_id, " | Length: ", read_length), x = "x", y = "y")

		# Add arrows if annotation available
		if(!is.null(ann_refined)){
			if (nrow(ann_refined) > 0) {
				p <- p +
					geom_segment(data=ann_refined %>% dplyr::filter(direction==1), aes(x=pred_inf, xend=pred_sup, y=0.9, yend=0.9, color=as.factor(direction)), arrow=arrow(length=unit(0.15, "inches")), inherit.aes=FALSE, linewidth=0.8) +
					geom_segment(data=ann_refined %>% dplyr::filter(direction==-1), aes(x=pred_sup, xend=pred_inf, y=0.9, yend=0.9, color=as.factor(direction)), arrow=arrow(length=unit(0.15, "inches")), inherit.aes=FALSE, linewidth=0.8)
			}

		}

		# Fork speeds
		if(!is.null(ann_speed)){
			if (nrow(ann_speed) > 0) {
				p <- p +
					geom_text(data=ann_speed %>% dplyr::filter(direction==1), aes(x=(x0_first_fork + x0_second_fork)/2, y=0.95, label=round(fork_speed, 0)), col=3) +
					geom_text(data=ann_speed %>% dplyr::filter(direction==-1), aes(x=(x0_first_fork + x0_second_fork)/2, y=0.95, label=round(fork_speed, 0)), col=3)
			}
		}

		if(!is.null(ann_event)){
			if (nrow(ann_event) > 0) {
				p <- p +
					geom_point(data=ann_event, aes(x=pos_event, y=0.85, col=type_event, pch=type_event), size=2) +
					scale_colour_manual(values=c("1"="blue", "-1"="red", "Ini"="#23D508","Ter"="#BC0000"), name="Features") +
					scale_shape_manual(values=c("Ini"=17, "Ter"=15))
			}else{
				p <- p +
					scale_color_manual(values = c("1"="blue", "-1"="red"), name="Relative direction")
			}
		}

		print(p)
	})

	dev.off()
}

compare.speed.referential <- function(forkml_speeds_genome, forkml_speeds_read, output_pdf){
	forkml_speeds_matched <- forkml_speeds_genome %>%
		group_by(read_id) %>%
		arrange(x0) %>%
		dplyr::mutate(idx=row_number()) %>%
		ungroup() %>%
		select(read_id, read_chrm, read_strand, idx, fork_speed, x0_first_fork, x0_second_fork) %>%
		full_join(forkml_speeds_read %>%
				group_by(read_id) %>%
				arrange(x0) %>%
				mutate(idx=row_number()) %>%
				mutate(idx=if_else(read_strand=="-", abs(idx - max(idx)) + 1, idx)) %>%
				ungroup() %>%
				select(read_id, read_chrm, read_strand, idx, fork_speed, x0_first_fork, x0_second_fork)
			, by=join_by(read_id, read_chrm, read_strand, idx), suffix=c(".genome",".read")) # TODO bad macthes by strands ?

	pdf(output_pdf)
	gp <- ggplot(forkml_speeds_matched) +
		geom_point(aes(fork_speed.genome, fork_speed.read, col=read_strand)) +
		theme_bw()
	print(gp)
	gp <- ggplot(forkml_speeds_matched) +
		geom_point(aes(log2(fork_speed.read * fork_speed.genome)/2, log2(fork_speed.read/fork_speed.genome), col=read_strand)) +
		theme_bw()
	print(gp)
	dev.off()
}

sample_name <- opt$sample
input_directory <- opt$input_dir
output_dir <- opt$output_dir
is_legacy <- opt$legacy
pulse_duration <- opt$pulse_duration
bin_size <- opt$bin_size
plot_enabled <- opt$plot

# Load raw detection
cat("Loading raw data\n")
forkml_raw <- load.forkml.results(sample_name, input_directory)

if(is_legacy){ # Swap coordinate & Properly define borders
	# Filter out reads with
	cat("Filter out suspicious reads\n")
	forkml_raw <- filter.out.suspicious.signal(forkml_raw, paste0(input_directory, "/scores_", sample_name, ".txt"), 0.4)

	cat("Refine raw detection\n")
	forkml_forks_genome <- refine.forkml.prediction(forkml_raw %>% mutate(pred_inf=pred_inf_bin, pred_sup=pred_sup_bin), TRUE)
	cat("Estimate fork speeds\n")
	forkml_speeds_genome <- detect.speeds.forkml(forkml_forks_genome, pulse_duration)
	cat("Localize replication events\n")
	forkml_events_genome <- detect.events.forkml(forkml_forks_genome)

	# To save
	forkml_forks <- forkml_forks_genome
	forkml_speeds <- forkml_speeds_genome
	forkml_events <- forkml_events_genome

	# To plot
	forkml_forks_plot <- forkml_forks_genome
	forkml_speeds_plot <- forkml_speeds_genome
	forkml_events_plot <- forkml_events_genome
}else{
	cat("Refine raw detection\n")
	forkml_forks_genome <- refine.forkml.prediction(forkml_raw %>% mutate(read_inf=read_start, read_sup=read_end, direction=if_else(read_strand=="-", direction * -1, direction)), TRUE)
	
	forkml_forks_read <- refine.forkml.prediction(forkml_raw %>% mutate(read_start=read_inf, read_end=read_sup, pred_inf=pred_inf_bin, pred_sup=pred_sup_bin), TRUE)
	cat("Estimate fork speeds\n")
	forkml_speeds_read <- detect.speeds.forkml(forkml_forks_read, pulse_duration) %>%
		mutate(fork_type=if_else(read_strand=="-", if_else(fork_type=="leading", "lagging", "leading"), fork_type))
	forkml_speeds_genome <- detect.speeds.forkml(forkml_forks_genome, pulse_duration) %>%
		mutate(fork_type=if_else(read_strand=="-", if_else(fork_type=="leading", "lagging", "leading"), fork_type))
	cat("Localize replication events\n")
	forkml_events_read <- detect.events.forkml(forkml_forks_read)
	forkml_events_genome <- detect.events.forkml(forkml_forks_genome)

	# To save
	forkml_forks <- forkml_forks_genome
	forkml_speeds <- forkml_speeds_genome
	forkml_events <- forkml_events_genome

	# To plot
	forkml_forks_plot <- forkml_forks_read
	forkml_speeds_plot <- forkml_speeds_read
	forkml_events_plot <- forkml_events_read

	# Save detection alternate referential
	save.datasets(forkml_forks_plot, forkml_speeds_plot, forkml_events_plot, paste0(input_directory,"/alt_referential"), sample_name)
}

# Save detection
save.datasets(forkml_forks, forkml_speeds, forkml_events, input_directory, sample_name)

# Load binned signals
forkml_signal <- load.forkml.signals(sample_name, input_directory)

# Load detection statistics
input_stats <- extract.input.stats(paste0(input_directory, "/inputStats_", sample_name, ".txt"))

# Generate basic fork speed analysis
analyze.forkml.speed(forkml_forks, forkml_speeds, forkml_events, forkml_signal, input_stats, bin_size, 90000, plot_enabled, paste0(output_dir, "ForkML.speeds_distribution.",sample_name,".pdf"))

# Generate
plot.forkml.example(forkml_signal, forkml_raw, forkml_forks_plot, forkml_speeds_plot, forkml_events_plot, paste0(output_dir,"ForkML.typical_signals_examples.",unique(forkml_forks_genome$sample),".pdf"), 0)
