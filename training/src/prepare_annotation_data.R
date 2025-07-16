#!/usr/bin/env Rscript

library("optparse")

# Option parser setup
option_list <- list(
  make_option(c("-a", "--path_annotations"), type="character", help="Path annotations"),
  make_option(c("-o", "--path_output"), type="character", help="Path output directory"),
  make_option(c("-s", "--seed"), type="integer", default=101, help="Random seed"),
  make_option(c("--noise_control"), type = "character", default = NULL, help = "Path to list read_id with noise")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

path_annotations <- opt$path_annotations
path_output <- opt$path_output
noise_control <- opt$noise_control

pulse_duration <- 30

library("tibble")
library("dplyr")
library("tidyr")
library("purrr")

process.annotations <- function(annotations, pulse_duration, data_type="doublepulse"){
    min_distance_initiation <- 2000
    maximum_fork_distance <- 200000
    amplitude_threshold <- 0.1
    border_distance_threshold <- 2000
    border_distance_relative_threshold <- 0.01
    read_length_threshold <- 60000
    relative_position_threshold <- 0.20

    annotations <- annotations %>%
        rowwise() %>%
        mutate(read_length=read_end - read_start) %>%

        mutate(xrs=min(x0, x2), yrs=ifelse(x0 <= x2, y0, y2)) %>%
        mutate(dist_x0_rs=abs(x0-read_start), dist_x2_rs=abs(x2-read_start), min_dist_rs=min(dist_x0_rs, dist_x2_rs)) %>%
        mutate(perc_x0_rs=dist_x0_rs/read_length, perc_x2_rs=dist_x2_rs/read_length) %>%
        mutate(on_lb=ifelse(read_length < read_length_threshold, ifelse(min_dist_rs < border_distance_threshold & yrs > amplitude_threshold, TRUE, FALSE), ifelse((min_dist_rs < border_distance_threshold | min(perc_x0_rs, perc_x2_rs) < border_distance_relative_threshold) & yrs > amplitude_threshold,TRUE,FALSE))) %>%

        mutate(xre=max(x0, x2), yre=ifelse(x0 >= x2, y0, y2)) %>%
        mutate(dist_x0_re=abs(x0-read_end), dist_x2_re=abs(x2-read_end), min_dist_re=min(dist_x0_re, dist_x2_re)) %>%
        mutate(perc_x0_re=dist_x0_re/read_length, perc_x2_re=dist_x2_re/read_length) %>%
        mutate(on_rb=ifelse(read_length < read_length_threshold, ifelse(min_dist_re < border_distance_threshold && yre > amplitude_threshold,TRUE,FALSE), ifelse((min_dist_re < border_distance_threshold || min(perc_x0_re, perc_x2_re) < border_distance_relative_threshold) && yre > amplitude_threshold,TRUE,FALSE))) %>%

        mutate(dist_x1_x0=abs(x1-x0), dist_x1_x2=abs(x1-x2), dist_x0_x2=abs(x0-x2)) %>%
        mutate(perc_x1_x0=dist_x1_x0/read_length, perc_x1_x2=dist_x1_x2/read_length, perc_x0_x2=dist_x0_x2/read_length) %>%
        mutate(near_ini=ifelse(read_length < read_length_threshold, ifelse(y0 > amplitude_threshold && on_lb==FALSE && on_rb==FALSE, TRUE, FALSE), ifelse(y0 > amplitude_threshold && on_lb==FALSE && on_rb==FALSE, TRUE, FALSE))) %>%
        mutate(near_ter=ifelse(read_length < read_length_threshold, ifelse(y2 > amplitude_threshold && on_lb==FALSE && on_rb==FALSE, TRUE, FALSE), ifelse(y2 > amplitude_threshold && on_lb==FALSE && on_rb==FALSE, TRUE, FALSE))) %>%
        mutate(is_unknown=ifelse((x2 > x0 & x2 < x1 & (dist_x0_x2/dist_x1_x0) > relative_position_threshold) | (x2 > x1 & x2 < x0 & (dist_x0_x2/dist_x1_x0) > relative_position_threshold), TRUE, FALSE)) %>%
        mutate(fork_dir=as.factor(ifelse(is_unknown,"unknown",ifelse(x0 < x2,"right","left")))) %>%
        mutate(near_ter=ifelse(is_unknown, FALSE, near_ter))

    res <- annotations %>%
        # dplyr::filter(!is.na(x0) & !is.na(x1) & !is.na(x2)) %>%
        group_by(read_chrm, read_id, person) %>%
        arrange(xrs) %>%
        dplyr::select(person, read_chrm, read_id, x0, x1, x2, y0, y1, y2, on_lb, on_rb, near_ini, near_ter, is_unknown, fork_dir) %>%
        rowwise() %>%
        mutate(lpos=min(c(x0, x1, x2)), rpos=max(c(x0, x1, x2))) %>%
        group_by(read_chrm, read_id, person) %>%
        mutate(dist_next_fork=lead(lpos) - rpos, dist_prev_fork=lpos - lag(rpos)) %>%
        mutate(step=ifelse(lag(fork_dir)==fork_dir & dist_prev_fork < maximum_fork_distance,0,1)) %>%
        mutate(step=ifelse(is.na(step),1,step)) %>%
        mutate(grp=cumsum(step))

    if(data_type=="doublepulse"){
        res <- res %>%
            arrange(read_chrm, read_id, person, grp, fork_dir) %>%
            group_by(read_chrm, read_id, person, grp, fork_dir) %>%
            mutate(fork_speed=abs(x0-lead(x0))/pulse_duration) %>%
            group_by(read_chrm, read_id, person) %>%
            mutate(x0=ifelse(fork_dir!="unknown" & lead(fork_dir)!="unknown" & fork_dir!=lead(fork_dir) & !is.na(fork_dir) & !is.na(lead(fork_dir)) & y0 > amplitude_threshold & lead(y0) > amplitude_threshold & abs(x0-lead(x0)) < min_distance_initiation, x0 + round((lead(x0) - x0)/2), x0)) %>% # UNTESTED
            mutate(x0=ifelse(fork_dir!="unknown" & lag(fork_dir)!="unknown" & fork_dir!=lag(fork_dir) & !is.na(fork_dir) & !is.na(lag(fork_dir)) & y0 > amplitude_threshold & lag(y0) > amplitude_threshold & abs(x0-lag(x0)) < min_distance_initiation, lag(x0) + 1, x0)) %>% # UNTESTED
            mutate(fork_speed=ifelse(fork_dir=="unknown", 
                NA, #"A", 
                ifelse(fork_dir=="left" & !is.na(lead(fork_dir)), 
                    ifelse(fork_dir=="left" & lead(fork_dir)=="left" & !is.na(lead(fork_dir,2)),
                        ifelse(fork_dir=="left" & lead(fork_dir)=="left" & lead(fork_dir,2)=="right" & abs(lead(x0)-lead(x0,2)) < min_distance_initiation,
                            NA, #"B",
                            ifelse(fork_dir=="left" & lead(fork_dir)=="left" & lead(y0) > amplitude_threshold, NA, fork_speed)
                        ),
                        ifelse(lead(on_rb),NA,fork_speed)
                    ), 
                    ifelse(fork_dir=="right" & !is.na(lag(fork_dir)), 
                        ifelse(fork_dir=="right" & lag(fork_dir)=="left" & abs(x0-lag(x0)) < min_distance_initiation, 
                            NA, #"C", 
                            fork_speed
                        ),
                        ifelse(on_lb,NA,ifelse(y0 > amplitude_threshold, NA, fork_speed))
                    )
                )
            )) %>% group_by(read_chrm, read_id, person, grp, fork_dir)
        }else{
            res <- res %>%
                mutate(fork_speed=abs(x0-x1)/pulse_duration) %>%
                mutate(fork_speed=ifelse(fork_dir=="unknown", NA, fork_speed)) %>%
                mutate(fork_speed=ifelse(fork_dir=="left" & on_rb, NA, fork_speed)) %>%
                mutate(fork_speed=ifelse(fork_dir=="right" & on_lb, NA, fork_speed)) %>%
                mutate(fork_speed=ifelse(y0 > amplitude_threshold, NA, fork_speed)) %>%
                group_by(read_chrm, read_id, person, grp, fork_dir)
        }

    return(list(annotations=annotations, res=res))
}

filter.annotation <- function(x, list_annotation_to_conserve){
    return(inner_join(x, list_annotation_to_conserve %>% dplyr::select(-has_unknown), by=c("read_id","person")))
}

print(paste0("Searching for ",path_annotations))
if(file.exists(path_annotations)){
    if(file.info(path_annotations)$isdir){
        manual_annotation_all <- list.files(path=path_annotations, pattern="*.rds", full.name=TRUE, recursive=TRUE) %>%
            map_dfr(readRDS)
    }
}else{
    # Load all files
    list_files <- unlist(map(strsplit(path_annotations, ",")[[1]], Sys.glob))
    # print(list_files)
    manual_annotation_all <- bind_rows(map(list_files, readRDS))
}
manual_annotation_all <- process.annotations(manual_annotation_all, pulse_duration)

# Check if annotation not fully broken
print(paste0("Number of annotated reads (with something): ",length(unique(manual_annotation_all$res$read_id))))

# Subsample one annotation per read only
set.seed(opt$seed)
list_annotation_to_conserve <- manual_annotation_all$annotations %>%
    group_by(read_id, person) %>%
    summarize(nb=n(), has_unknown=any(fork_dir=="unknown", na.rm=TRUE), .groups="drop") %>%
    group_by(read_id) %>%
    dplyr::filter(!any(has_unknown)) %>% # Remove read_id with dubious signal
    sample_n(1)

manual_annotation_selected <- manual_annotation_all
manual_annotation_selected$res <- filter.annotation(manual_annotation_selected$res, list_annotation_to_conserve)
manual_annotation_selected$annotations <- filter.annotation(manual_annotation_selected$annotations, list_annotation_to_conserve)
list_annotation_to_conserve <- list_annotation_to_conserve %>% dplyr::select(read_id) # Simplify

# Handle noise dataset
if(!is.null(noise_control)){
    # Read the list of read_ids from a text file
    list_read_id_noise <- readLines(noise_control)
 print(head(list_read_id_noise))
    # Get the column names and types from an existing dataset
    template <- manual_annotation_selected$res %>% ungroup() %>% slice(0)
 print(template)
    # Create a new tibble with the same columns, filling read_id and leaving others as NA
    new_annotations <- template %>%
        add_row(!!!set_names(rep(list(NA), length(template)), names(template))) %>%
        slice(rep(1, length(list_read_id_noise))) %>%
        mutate(read_id=list_read_id_noise)
 print(new_annotations)
    # Append to main datasets
    manual_annotation_selected$res <- bind_rows(manual_annotation_selected$res, new_annotations)
    list_annotation_to_conserve <- rbind(list_annotation_to_conserve, tibble(read_id=list_read_id_noise))
}

saveRDS(manual_annotation_selected, file=path_output)
write.table(list_annotation_to_conserve, file=gsub("annotation.rds", "annotation.read_id.txt", path_output), quote=FALSE, col.names=FALSE, row.names=FALSE)
