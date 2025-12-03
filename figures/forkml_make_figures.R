# Rscript --vanilla this_script.R <path_to_data>

options(repos=c(CRAN="https://cloud.r-project.org"))
conda_lib <- .libPaths()[grepl("\\.conda", .libPaths())][1]
if (!is.na(conda_lib)) .libPaths(conda_lib)

devtools::install_github("teunbrand/ggarrow")

library("dplyr")
library("tibble")
library("gtools")
library("GenomicRanges")
library("tidyr")
library("ggplot2")
library("ggarrow")
library("patchwork")
library("ggrastr")
library("ggpubr")
library("scales")
library("ggbeeswarm")


load.gloe.data <- function(path_data){
    gloe_seqx <- read.table(path_data, header=FALSE, stringsAsFactors=TRUE) %>% as_tibble()
    colnames(gloe_seqx) <- c("chr","start","nb_w","nb_c","rfd")
    gloe_seqx <- gloe_seqx %>% arrange(chr, start)
    gloe_seqx$chr <- factor(gloe_seqx$chr, levels=mixedsort(unique(as.character(gloe_seqx$chr))))

    return(gloe_seqx)    
}

compute.iz.content <- function(forkml_event_annotation, IZ_wavelet){
    gr_events <- GRanges(
        seqnames = forkml_event_annotation$read_chrm,
        ranges = IRanges(
            start = forkml_event_annotation$pos_event_better2,
            end = forkml_event_annotation$pos_event_better2
        ),
        type_event = forkml_event_annotation$type_event
    )

    gr_IZ <- GRanges(
        seqnames = IZ_wavelet$chr,
        ranges = IRanges(
            start = IZ_wavelet$start,
            end = IZ_wavelet$end
        )
    )

    overlaps <- findOverlaps(gr_events, gr_IZ)

    # Annotate events
    forkml_event_annotation$is_IZ <- FALSE
    forkml_event_annotation$is_IZ[queryHits(overlaps)] <- TRUE

    res <- forkml_event_annotation %>%
        group_by(sample, type_event, is_IZ) %>%
        summarize(n=n(), .groups="drop") %>%
        pivot_wider(names_from=is_IZ, values_from=n) %>%
        mutate(ratio=`TRUE`/(`TRUE`+`FALSE`)) %>%
        arrange(type_event, sample)

    return(res)
}

rebinning.replication.timing <- function(replication_timing, bin_size=100000){
    binned_expanded_timing <- replication_timing %>%
        rowwise() %>%
        mutate(read_chrm=chrm) %>%
        mutate(bin_start=floor(start / bin_size) * bin_size, bin_end=ceiling(end / bin_size) * bin_size - 1) %>%
        arrange(read_chrm, bin_start) %>%
        group_by(read_chrm, start, end) %>%
        summarise(bin=list(seq(bin_start, bin_end, by=bin_size)), start=list(start), end=list(end), weighted_average=list(weighted_average), .groups='drop') %>%
        unnest(cols=c(bin, start, end, weighted_average)) %>%
        rowwise() %>%  # Use rowwise to calculate overlap for each row individually
        mutate(bin_end=bin + bin_size) %>%
        mutate(overlap_start=max(start, bin), overlap_end=min(end, bin_end), overlap_length=max(0, overlap_end - overlap_start)) %>%
        mutate(relative_weight=weighted_average * overlap_length) %>%
        group_by(read_chrm, bin) %>%
        summarize(coverage_bin=sum(overlap_length), binned_average=sum(relative_weight)/sum(overlap_length), .groups='drop')

    binned_expanded_timing$read_chrm <- factor(binned_expanded_timing$read_chrm, levels=mixedsort(unique(as.character(binned_expanded_timing$read_chrm))))

    return(binned_expanded_timing)
}

#  ______                              _                
#  | ___ \                            | |               
#  | |_/ /_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___ 
#  |  __/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
#  | | | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
#  \_|  \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
#                                                       
#                                                       

args <- commandArgs(trailingOnly=TRUE)
print(args)

main_path <- args[1] # path dataset directory

cmb_name <- "VRComplete"
nb_pulse <- 2
length_threshold <- 90000 # Default
length_threshold_hela <- length_threshold

df_length_threshold <- data.frame(
    sample=c("WT_rep1", "WT_rep2", "WT_rep3", "WT_rep4", "HU_rep1", "HU_rep2", "APH_rep1", "APH_rep2", "WT_rep4_R10", "HeLa_rep1", "HeLa_rep2", "A549", "BJ", "JEG3", "LoVo"),
    min_read_length=c(length_threshold, length_threshold, length_threshold, length_threshold, length_threshold, length_threshold, length_threshold, length_threshold, length_threshold, length_threshold_hela, length_threshold_hela, 70000, 100000, 80000, 90000)
)

sample_clean_name <- data.frame(
    original=c("VR34", "VR36", "VR37", "VR53", "VR38", "VR40", "VR39", "VR41", "VR54", "FP7", "VR57", "FP10", "FP12"),
    simple=c(paste0("WT_rep",seq(1,4)), paste0("HU_rep",seq(1,2)), paste0("APH_rep",seq(1,2)), "WT_rep4_R10", "HeLa_rep0_R9", paste0("HeLa_rep",seq(0,2))),
    pore_type=c("R9", "R9", "R9", "R9", "R9", "R9", "R9", "R9", "R10", "R9", "R10", "R10", "R10"),
    clean=c(paste0("Untreated (",seq(1,4),")"), paste0("HU 100 µM (",seq(1,2),")"), paste0("APH 100 nM (",seq(1,2),")"), "Untreated (4) R10", "HeLa Untreated (0) R9", "HeLa Untreated (0)", "HeLa Untreated (1)", "HeLa Untreated (2)"),
    mycol=c("#A7CAAF", "#88CF98", "#42C05E", "#1B9437", "#F78C41", "#ED6506", "#D3618C", "#F35A94", "#1B9437", "#5A9BD5", "#5A9BD5", "#6D95C9", "#3E6FB0")
)
sample_bonus_clean_name <- data.frame(
    original=c("FP20", "FP19", "FP21", "FP22"),
    simple=c("BJ", "A549", "JEG3", "LoVo"),
    pore_type=c("R10", "R10", "R10", "R10"),
    clean=c("BJ-hTERT Untreated", "A549 Untreated", "JEG-3 Untreated", "LoVo Untreated"),
    mycol=c("#BF743D", "#1BB3C1", "#DDB5F7", "#F0C44D")
)
sample_clean_name <- rbind(sample_clean_name, sample_bonus_clean_name)
sample_clean_name$clean <- factor(sample_clean_name$clean, levels=sample_clean_name$clean)
sample_clean_name$simple <- factor(sample_clean_name$simple, levels=sample_clean_name$simple)

# For chromatin
clean_name <- c("A1"="Active chromatin,\nhigh transcription", "A2"="Active chromatin,\nweak transcription", "B0"="Poised\nheterochromatin", "B1"="Facultative\nheterochromatin", "B4"="Constitutive\nheterochromatin", "Unknown"="#808080")
pretty_colors1 <- c("A1"="#85C1E9", "A2"="#D6EAF8", "B0"="#FCEb4E", "B1"="#FFD700", "B4"="#FFB133", "Unknown"="#808080")
pretty_colors2 <- pretty_colors1
names(pretty_colors2) <- clean_name

bin_size <- 30000 # For problematic regions

mypal <- c("Raw"="#AEAEAE", "100 bp average"="#808080", "Gaussian smoothing (2500 bp)"="#000000", "Rightward fork"="#F8766D", "Leftward fork"="#619CFF", "Initiation event"="#23D508", "Termination event"="#BC0000")
# mypal2 <- c("1"="solid", "-1"="solid")
mypal2 <- c("Initiation event"=17, "Termination event"=15)


#  ______                     _                   _     
#  | ___ \                   (_)                 | |    
#  | |_/ /__ ___      __  ___ _  __ _ _ __   __ _| |___ 
#  |    // _` \ \ /\ / / / __| |/ _` | '_ \ / _` | / __|
#  | |\ \ (_| |\ V  V /  \__ \ | (_| | | | | (_| | \__ \
#  \_| \_\__,_| \_/\_/   |___/_|\__, |_| |_|\__,_|_|___/
#                                __/ |                  
#                               |___/                   

example_amp <- readRDS(file=paste0(main_path, "/example_amp.RDS"))
example_ipd <- readRDS(file=paste0(main_path, "/example_ipd.RDS"))
example_HCT116_UT_R9_rep1 <- readRDS(file=paste0(main_path, "/example_reads.HCT116_UT_R9_rep1.RDS"))
example_HCT116_UT_R10_rep4 <- readRDS(file=paste0(main_path, "/example_reads.HCT116_UT_R10_rep4.RDS"))
example_HCT116_HU_R9_rep1 <- readRDS(file=paste0(main_path, "/example_reads.HCT116_HU_R9_rep1.RDS"))
example_HCT116_APH_R9_rep1 <- readRDS(file=paste0(main_path, "/example_reads.HCT116_APH_R9_rep1.RDS"))
example_HeLa_UT_R10_rep1 <- readRDS(file=paste0(main_path, "/example_reads.HeLa_UT_R10_rep1.RDS"))
example_bonus <- readRDS(file=paste0(main_path, "/example_reads.bonus.RDS"))

#   _____  _     _____ _____                      
#  |  __ \| |   |  _  |  ___|                     
#  | |  \/| |   | | | | |__ ______ ___  ___  __ _ 
#  | | __ | |   | | | |  __|______/ __|/ _ \/ _` |
#  | |_\ \| |___\ \_/ / |___      \__ \  __/ (_| |
#   \____/\_____/\___/\____/      |___/\___|\__, |
#                                              | |
#                                              |_|

gloe_seqx2 <- load.gloe.data(paste0(main_path, "/all_SRR1105855x_trimmed_filtered.txt"))

#  ______         _   ___  ___ _           _      _            _   _                 
#  |  ___|       | |  |  \/  || |         | |    | |          | | (_)                
#  | |_ ___  _ __| | _| .  . || |       __| | ___| |_ ___  ___| |_ _  ___  _ __  ___ 
#  |  _/ _ \| '__| |/ / |\/| || |      / _` |/ _ \ __/ _ \/ __| __| |/ _ \| '_ \/ __|
#  | || (_) | |  |   <| |  | || |____ | (_| |  __/ ||  __/ (__| |_| | (_) | | | \__ \
#  \_| \___/|_|  |_|\_\_|  |_/\_____/  \__,_|\___|\__\___|\___|\__|_|\___/|_| |_|___/
#                                                                                    
#                                                                                    

forkml_fork_annotation_hct116 <- readRDS(file=paste0(main_path, "/forkml_fork_annotation_hct116.RDS"))
forkml_speed_annotation_hct116 <- readRDS(file=paste0(main_path, "/forkml_speed_annotation_hct116.RDS"))
forkml_event_annotation_hct116 <- readRDS(file=paste0(main_path, "/forkml_event_annotation_hct116.RDS"))
forkml_fork_annotation_r9_ctrl <- readRDS(file=paste0(main_path, "/forkml_fork_annotation_r9_ctrl.RDS"))
forkml_speed_annotation_r9_ctrl <- readRDS(file=paste0(main_path, "/forkml_speed_annotation_r9_ctrl.RDS"))
forkml_event_annotation_r9_ctrl <- readRDS(file=paste0(main_path, "/forkml_event_annotation_r9_ctrl.RDS"))
forkml_fork_annotation_r10 <- readRDS(file=paste0(main_path, "/forkml_fork_annotation_r10.RDS"))
forkml_speed_annotation_r10 <- readRDS(file=paste0(main_path, "/forkml_speed_annotation_r10.RDS"))
forkml_event_annotation_r10 <- readRDS(file=paste0(main_path, "/forkml_event_annotation_r10.RDS"))
forkml_fork_annotation_r10_ctrl <- readRDS(file=paste0(main_path, "/forkml_fork_annotation_r10_ctrl.RDS"))
forkml_speed_annotation_r10_ctrl <- readRDS(file=paste0(main_path, "/forkml_speed_annotation_r10_ctrl.RDS"))
forkml_event_annotation_r10_ctrl <- readRDS(file=paste0(main_path, "/forkml_event_annotation_r10_ctrl.RDS"))
forkml_fork_annotation_bonus <- readRDS(file=paste0(main_path, "/forkml_fork_annotation_bonus.RDS"))
forkml_speed_annotation_bonus <- readRDS(file=paste0(main_path, "/forkml_speed_annotation_bonus.RDS"))
forkml_event_annotation_bonus <- readRDS(file=paste0(main_path, "/forkml_event_annotation_bonus.RDS"))

forkml_event_annotation_random_marked <- readRDS(file=paste0(main_path, "/forkml_event_annotation_random_marked.RDS"))
replication_timing_v1 <- read.table(paste0(main_path, "/timing_HR_allbin_liftover_v1.bed"), header=TRUE, stringsAsFactors=TRUE)
all_chromatine_stats <- readRDS(file=paste0(main_path, "/all_chromatine_stats.RDS"))
IZ_wavelet <- read.table(paste0(main_path, "/GLOE_HCT116_waveplus_v1.bedgraph"), skip=1, col.names=c("chr","start","end","DelaRFD")) %>% as_tibble()

df_speed_manual <- readRDS(file=paste0(main_path, "/df_speed_manual.RDS"))
df_speed_auto <- readRDS(file=paste0(main_path, "/df_speed_auto.RDS"))
df_speed_manual_r10 <- readRDS(file=paste0(main_path, "/df_speed_manual_r10.RDS"))
df_speed_manual_hela <- readRDS(file=paste0(main_path, "/df_speed_manual_hela.RDS"))
df_speed_manual_bonus <- readRDS(file=paste0(main_path, "/df_speed_manual_bonus.RDS"))

repartition_WT <- compute.iz.content(forkml_event_annotation_hct116 %>% dplyr::filter(grepl("Untreated", full_name) & !is_outlier) %>% mutate(sample=full_name), IZ_wavelet) %>% dplyr::rename(non_IZ = `FALSE`, IZ = `TRUE`)
repartition_WT_all <- compute.iz.content(forkml_event_annotation_hct116 %>% dplyr::filter(grepl("Untreated", full_name) & !is_outlier) %>% mutate(sample="Untreated - All"), IZ_wavelet) %>% dplyr::rename(non_IZ = `FALSE`, IZ = `TRUE`)
repartition_random <- compute.iz.content(forkml_event_annotation_random_marked %>% mutate(sample="Random"), IZ_wavelet) %>% dplyr::rename(non_IZ = `FALSE`, IZ = `TRUE`)
repartition_matrix <- rbind(repartition_WT_all %>% dplyr::filter(type_event=="Ini"), repartition_random) %>% dplyr::select(c("non_IZ","IZ")) %>% as.data.frame()
rownames(repartition_matrix) <- c("Untreated","Random")


#  ______ _                           
#  |  ___(_)                          
#  | |_   _  __ _ _   _ _ __ ___  ___ 
#  |  _| | |/ _` | | | | '__/ _ \/ __|
#  | |   | | (_| | |_| | | |  __/\__ \
#  \_|   |_|\__, |\__,_|_|  \___||___/
#            __/ |                    
#           |___/                     

plot_example <- function(example_reads, forkml_fork_annotation, forkml_speed_annotation, forkml_event_annotation, cor_factor_point=0.7, cor_factor_line=0.7, opt_title=NULL, make_raster=TRUE, forced_arrow=TRUE, potluck=FALSE){
    simple_labels <- example_reads %>%
        group_by(clean_read_id) %>%
        mutate(read_length=n()) %>%
        dplyr::filter(!is.na(avg)) %>%
        mutate(xlab=min(positions), read_length_w_data=diff(range(positions)) + 1) %>%
        distinct(clean_read_id, chrom, read_length, read_length_w_data, xlab) %>%
        mutate(full_lab=as.character(clean_read_id), lab=paste0("#", substr(full_lab, nchar(full_lab), nchar(full_lab)), ", ",chrom," (",round(read_length_w_data/1000,0)," kb)"))        
    if(potluck){
        simple_labels <- simple_labels %>%
            mutate(full_lab=as.character(clean_read_id), lab=paste0(sub(" .*", "", full_lab), " #", substr(full_lab, nchar(full_lab), nchar(full_lab)), ", ",chrom," (",round(read_length_w_data/1000,0)," kb)"))
    }

    forkml_fork_annotation <- forkml_fork_annotation %>%
        dplyr::filter(read_id %in% unique(example_reads$read_id)) %>%
        left_join(example_reads %>% distinct(clean_read_id, read_id), by=join_by(read_id)) %>%
        mutate(direction=as.factor(ifelse(direction==1, "Rightward fork", "Leftward fork"))) %>%
        mutate(xstart=ifelse(direction=="Rightward fork", pred_inf, pred_sup)) %>%
        mutate(xend=ifelse(direction=="Rightward fork", pred_sup, pred_inf)) %>%
        mutate(clean_read_id = droplevels(clean_read_id)) # Fix plotting issue

    forkml_speed_annotation <- forkml_speed_annotation %>%
        dplyr::filter(read_id %in% unique(example_reads$read_id)) %>%
        left_join(example_reads %>% distinct(clean_read_id, read_id), by=join_by(read_id)) %>%
        mutate(direction=as.factor(ifelse(fork_dir=="right", "Rightward fork", "Leftward fork"))) %>%
        mutate(pos_speed=(x0_first_fork + x0_second_fork)/2) %>%
        mutate(pos_left=ifelse(direction=="Rightward fork",x0_first_fork,x0_second_fork), pos_right=ifelse(direction=="Rightward fork",x0_second_fork,x0_first_fork))

    forkml_event_annotation <- forkml_event_annotation %>%
        dplyr::filter(read_id %in% unique(example_reads$read_id)) %>%
        left_join(example_reads %>% distinct(clean_read_id, read_id), by=join_by(read_id)) %>%
        mutate(type_event=ifelse(type_event=="Ini", "Initiation event","Termination event"))

    label_size <- 5/.pt
    linewidth_target <- 0.5
    y_speeds <- 0.950
    half_b_height <- 0.015

    if(is.null(make_raster)){
        gp <- ggplot(example_reads)
        custom_legend_size <- c(2.5, 2.5)
        custom_legend_pch <- c(17, 15)
    }else if(make_raster){
        gp <- ggplot(example_reads) +
            rasterise(geom_point(aes(x=positions, y=avg, col="100 bp average"), pch=16, alpha=0.5, size=0.05*cor_factor_point), dev="cairo", dpi=600)
        custom_legend_size <- c(1, 2.5, 2.5)
        custom_legend_pch <- c(16, 17, 15)
    }else{
        gp <- ggplot(example_reads) +
            geom_point(aes(x=positions, y=avg, col="100 bp average"), pch=16, alpha=0.5, size=0.05*cor_factor_point)        
        custom_legend_size <- c(1, 2.5, 2.5)
        custom_legend_pch <- c(16, 17, 15)
    }

    gp <- gp +
        geom_line(aes(x=positions, y=smooth_signal_norm, col="Gaussian smoothing (2500 bp)"), linewidth=linewidth_target*cor_factor_line) +

        # geom_segment(data=forkml_fork_annotation, aes(x=xstart, xend=xend, y=0.9, yend=0.9, col=direction), linewidth=0.5, lineend="round", linejoin="round", arrow=arrow(length=unit(0.05, "inches"))) +
        geom_arrow_segment(data=forkml_fork_annotation, aes(x=xstart, xend=xend, y=0.9, yend=0.9, color=direction), linewidth=0.3, lineend="round", linejoin="round", arrow_head=arrow_head_wings(offset=30, inset=40), length_head=2.6, force_arrow=forced_arrow) + # , arrow_angle=20, arrow_length=0.03

        geom_segment(data=forkml_speed_annotation, aes(x=pos_left, xend=pos_right, y=y_speeds, yend=y_speeds, col=direction), linewidth=0.12, lineend="round", linejoin="round", show.legend=FALSE) +
        geom_segment(data=forkml_speed_annotation, aes(x=pos_left, xend=pos_left, y=y_speeds + half_b_height, yend=y_speeds - half_b_height, col=direction), linewidth=0.2, lineend="round", linejoin="round", show.legend=FALSE) +
        geom_segment(data=forkml_speed_annotation, aes(x=pos_right, xend=pos_right, y=y_speeds + half_b_height, yend=y_speeds - half_b_height, col=direction), linewidth=0.2, lineend="round", linejoin="round", show.legend=FALSE) +
        geom_text(data=forkml_speed_annotation, aes(x=pos_speed, y=y_speeds + half_b_height*4, label=paste0(round(fork_speed, 0)," bp/min"), col=direction), size=label_size, show.legend=FALSE) +

        geom_point(data=forkml_event_annotation, aes(x=pos_event_better2, y=0.85, col=as.factor(type_event), pch=as.factor(type_event)), size=1) +

        geom_text(data=simple_labels, aes(x=xlab, y=1.10, label=lab), col=1, hjust=0, size=label_size) +
        facet_wrap(vars(clean_read_id), ncol=2, scales="free_x") +
        theme_bw() +
        scale_color_manual(name=NULL, values=mypal, breaks=names(mypal)) +
        scale_shape_manual(name=NULL, values=mypal2) +
        scale_x_continuous(expand=expansion(mult=c(0.05, 0.04)), labels=function(x) scales::comma(x / 1000), name="Genomic position (kb)") +
        scale_y_continuous(expand=expansion(mult=c(0.05, 0.05)), name="Predicted BrdU content") +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_text(size = 5),
            axis.title = element_text(size = 7),
            legend.title = element_text(size = 7, face = "bold"),
            legend.position = "bottom",
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-5, -5, -5, -5),
            legend.text = element_text(size = 6)
        ) +
        guides(pch="none", color=guide_legend(nrow=1, override.aes=list(alpha=1, linewidth=1, size=custom_legend_size, pch=custom_legend_pch)))

    if(!is.null(opt_title)){
        gp <- gp +
            labs(title=opt_title) +
            theme(plot.title=element_text(size=7))
    }

    return(gp)
}

make.length.plots <- function(forkml_speed_annotation, individual_sample){
    bin_size <- 10000
    max_read_length1 <- 300000
    max_read_length2 <- 150000
    th_metric <- 5

    data <- forkml_speed_annotation %>% dplyr::filter(grepl(individual_sample, full_name)) %>%
        mutate(capped_read_length=ifelse(read_length > max_read_length1, max_read_length1/1000, read_length/1000)) %>% # Convert to kb
        mutate(read_length_group=cut_interval(capped_read_length, length=bin_size/1000, dig.lab=10))

    levels(data$read_length_group) <- gsub(sprintf("%d", max_read_length1/1000), "Inf", levels(data$read_length_group))
    levels(data$read_length_group) <- gsub("Inf\\]", "Inf)", levels(data$read_length_group))

    data_sample <- forkml_speed_annotation %>%
        mutate(capped_read_length=ifelse(read_length > max_read_length2, max_read_length2/1000, read_length/1000)) %>% # Convert to kb
        mutate(read_length_group=cut_interval(capped_read_length, length=bin_size/1000, dig.lab=10))

    summary_data_sample <- data_sample %>%
        group_by(full_name, read_length_group) %>%
        summarize(nb_fork_speed=n(), median_fork_speed=median(fork_speed, na.rm=TRUE), sd_fork_speed=sd(fork_speed, na.rm=TRUE), .groups="drop")

    summary_data <- data %>%
        group_by(full_name, read_length_group) %>%
        summarize(nb_fork_speed=n(), median_fork_speed=median(fork_speed, na.rm=TRUE), sd_fork_speed=sd(fork_speed, na.rm=TRUE), .groups="drop") %>%
        dplyr::select(colnames(summary_data_sample))

    summary_data_sample_type <- data_sample %>%
        group_by(full_name, read_length_group) %>%
        summarize(nb_fork_speed=n(), median_fork_speed=median(fork_speed, na.rm=TRUE), sd_fork_speed=sd(fork_speed, na.rm=TRUE), .groups="drop")

    summary_data_all <- rbind(summary_data_sample) %>%
        group_by(full_name) %>%
        mutate(delta_median=abs(median_fork_speed - lag(median_fork_speed))) %>%
        mutate(perc_delta_median=(delta_median/lag(median_fork_speed))*100)

    df <- summary_data_all %>%
        dplyr::select(full_name, read_length_group, perc_delta_median) %>%
        pivot_longer(cols=starts_with("perc_delta"), names_to="metrics", values_to="value")
    df_th1 <- df %>%
        group_by(full_name, metrics) %>%
        dplyr::filter(abs(value) <= th_metric) %>%
        dplyr::filter(as.integer(read_length_group)==min(as.integer(read_length_group))) 

    all_levels <- levels(df_th1$read_length_group)
    largest_used_level <- tail(sort(unique(as.character(df_th1$read_length_group))), 1)
    idx <- match(largest_used_level, all_levels)
    final_threshold <- if (idx < length(all_levels)) all_levels[idx + 1] else NA

    col_vector <- c("All"="#000000", setNames(sample_clean_name$mycol, sample_clean_name$clean))

    gp1 <- ggplot(data, aes(x=read_length_group, y=fork_speed, fill=full_name)) +
        geom_boxplot(fill=NA, outlier.shape=NA, show.legend=FALSE) +
        stat_summary(fun=mean, geom="point", size=1, color="red", show.legend=FALSE) +  # Add median points
        geom_text(data=summary_data, aes(x=read_length_group, y=1, label=nb_fork_speed), size=2) +
        labs(x=paste0("Read length groups (kb)"), y=paste0("Fork speed in ",as.character(unique(data$full_name))," (bp/min)")) +
        scale_fill_manual(values=col_vector) +
        scale_y_continuous(limits=c(0, 2000)) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        theme(
            plot.title = element_text(size = 8, face = "bold"),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8)
        )

    gp2 <- ggplot(summary_data_sample_type, aes(x=read_length_group, y=median_fork_speed, color=full_name)) +
        geom_line(aes(group=full_name), linewidth=0.3) + 
        geom_point(size=1) + 
        labs(x=paste0("Read length groups (capped at ",max_read_length2/1000," kb)"), y="Median fork speed (bp/min)", color="Samples") +
        scale_y_continuous(limits=c(0, 2000)) +
        scale_color_manual(values=col_vector) +
        theme_bw() + 
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        theme(
            plot.title = element_text(size = 8, face = "bold"),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.title = element_text(size = 8, face = "bold"),
            legend.text = element_text(size = 8)
        )

    gp3 <- ggplot(df, aes(x=read_length_group, y=value, color=as.character(full_name), group=full_name)) +
        geom_point(size=1) +
        geom_line(linewidth=0.5) +
        geom_hline(yintercept=5, lty="dashed", linewidth=0.3) +
        scale_color_manual(values=col_vector) +
        scale_y_continuous(position = "right") +
        labs(x=paste0("Read length groups (capped at ",max_read_length2/1000," kb)"), y="Absolute median fork speed change (%)", color="Samples") +
        theme_bw() + 
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        theme(legend.position="none") +
        theme(
            plot.title = element_text(size = 8, face = "bold"),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8)
        )

    return(list(gp1=gp1, gp2=gp2, gp3=gp3))
}

plot.relationship.forkspeed.raw_timing.paper <- function(forkml_speed_annotation, replication_timing){
    localized_forkml_speed_annotation_all <- forkml_speed_annotation %>%
        # ungroup() %>%
        dplyr::select(read_chrm, read_strand, read_id, fork_speed, x0_first_fork, x0_second_fork) %>%
        mutate(fork_speed_position=(x0_first_fork + x0_second_fork)/2)

    # Creating GRanges for fork speeds
    fork_speed_ranges <- GRanges(
        seqnames=localized_forkml_speed_annotation_all$read_chrm,
        ranges=IRanges(start=localized_forkml_speed_annotation_all$fork_speed_position, end=localized_forkml_speed_annotation_all$fork_speed_position),
        fork_speed=localized_forkml_speed_annotation_all$fork_speed
    )

    # Creating GRanges for replication timing
    timing_ranges <- GRanges(
      seqnames = replication_timing$chrm,
      ranges = IRanges(start = replication_timing$start,end = replication_timing$end
      ),
      weighted_average = replication_timing$weighted_average
    )

    # Finding overlaps
    overlaps <- findOverlaps(fork_speed_ranges, timing_ranges)

    # Extracting data based on overlaps
    result <- data.frame(
        fork_speed = fork_speed_ranges$fork_speed[queryHits(overlaps)],
        weighted_average = timing_ranges$weighted_average[subjectHits(overlaps)],
        fork_speed_position = start(fork_speed_ranges[queryHits(overlaps)])
    )

    result_no_extremum <- result %>%
        ungroup() %>%
        dplyr::filter(between(weighted_average, quantile(weighted_average, 0.01, na.rm=TRUE), quantile(weighted_average, 0.99, na.rm=TRUE)))

    gp <- ggplot(result, aes(x = weighted_average, y = fork_speed)) +
        geom_point(pch=46) +
        geom_smooth(data=result_no_extremum, aes(x=weighted_average, y=fork_speed), method="gam", se=TRUE, color="red") +
        labs(x="Replication timing (early to late)", y="Fork speed (bp/min)") +
        scale_x_continuous(limits=c(1,16)) +
        coord_cartesian(xlim=c(1,16), ylim=c(0, 2500)) +
        theme_bw()

    return(gp)
}

plot.fork.speed.chromatine.timing.paper <- function(forkml_speed_annotation_subset, binned_expanded_timing, all_chromatine_stats){
    # Add chromatine information on timing
    gr_speeds <- GRanges(
        seqnames = forkml_speed_annotation_subset$read_chrm,
        ranges = IRanges(
            start = (forkml_speed_annotation_subset$x0_first_fork + forkml_speed_annotation_subset$x0_second_fork)/2,
            end = (forkml_speed_annotation_subset$x0_first_fork + forkml_speed_annotation_subset$x0_second_fork)/2
        )
    )

    gr_chromatine <- GRanges(
        seqnames = all_chromatine_stats$chrom,
        ranges = IRanges(
            start = all_chromatine_stats$chromStart,
            end = all_chromatine_stats$chromEnd
        )
    )

    overlaps <- findOverlaps(gr_speeds, gr_chromatine)

    # Filter overlaps: keep only those with unique queryHits
    query_counts <- table(queryHits(overlaps))
    multiple_query <- as.integer(names(query_counts[query_counts > 1]))
    overlaps <- overlaps[!queryHits(overlaps) %in% multiple_query]

    # Annotate speeds
    forkml_speed_annotation_subset$chromatine <- "Unknown"
    forkml_speed_annotation_subset$chromatine[queryHits(overlaps)] <- as.character(all_chromatine_stats$name[subjectHits(overlaps)])
    forkml_speed_annotation_subset$chromatine_name <- "Unknown"
    forkml_speed_annotation_subset$chromatine_name[queryHits(overlaps)] <- as.character(all_chromatine_stats$clean_name[subjectHits(overlaps)])

    gr_timing <- GRanges(
        seqnames = binned_expanded_timing$read_chrm,
        ranges = IRanges(
            start = binned_expanded_timing$bin + 1,
            end = binned_expanded_timing$bin + 100000
        )
    )

    overlaps <- findOverlaps(gr_speeds, gr_timing)

    # Filter overlaps: keep only those with unique queryHits
    query_counts <- table(queryHits(overlaps))
    multiple_query <- as.integer(names(query_counts[query_counts > 1]))
    overlaps <- overlaps[!queryHits(overlaps) %in% multiple_query]

    # Annotate speeds
    forkml_speed_annotation_subset$timing <- NA
    forkml_speed_annotation_subset$timing[queryHits(overlaps)] <- binned_expanded_timing$binned_average[subjectHits(overlaps)]

    forkml_speed_annotation_subset <- forkml_speed_annotation_subset %>% dplyr::filter(chromatine!="Unknown")
    # forkml_speed_annotation_subset$chromatine_name <- factor(forkml_speed_annotation_subset$chromatine_name, levels=unique(forkml_speed_annotation_subset$chromatine_name))

    gp1 <- ggplot(forkml_speed_annotation_subset, aes(x=timing, y=fork_speed, col=chromatine, fill=chromatine)) +
        geom_smooth(method="gam", se=TRUE, alpha=0.4) +
        scale_color_manual(values=pretty_colors1, breaks=levels(all_chromatine_stats$name)) +
        scale_fill_manual(values=pretty_colors1, breaks=levels(all_chromatine_stats$name)) +
        scale_x_continuous(limits=c(1,16)) +
        guides(fill="none", col=guide_legend(nrow=2, byrow=TRUE, override.aes=list(fill="white", alpha=1))) + # fill="none", 
        labs(x="Replication timing (early to late)", y="Fork speed (bp/min)", col=NULL) +
        theme_bw() +
        theme(
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 8),

            legend.position = c(0.05, 0.95),
            legend.justification = c(0, 1),
            legend.box.background = element_rect(color = "black", fill = "white"),
            legend.margin = margin(1,3,1,3),
            legend.text = element_text(size = 6),
            legend.key = element_rect(fill = "white", color = NA),
            legend.key.height = unit(3, "mm"),
            legend.key.width = unit(3, "mm"),
            legend.key.spacing.y = unit(2, 'mm')
        )

    n_labels <- forkml_speed_annotation_subset %>%
        group_by(chromatine, chromatine_name) %>%
        summarize(n=n(), .groups="drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("n = ",n), n))

    median_labels <- forkml_speed_annotation_subset %>%
        group_by(chromatine, chromatine_name) %>%
        summarize(median_fork_speed=round(median(fork_speed, na.rm=TRUE),0), .groups="drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("Med. = ",median_fork_speed), median_fork_speed))

    label_size <- 5/.pt
    gp2 <- ggplot(forkml_speed_annotation_subset, aes(x=chromatine, y=fork_speed, fill=chromatine_name)) +
        geom_violin() +
        geom_boxplot(width=0.2, outlier.shape=NA, show.legend=FALSE) +
        stat_summary(fun="mean", geom="point", color="red", show.legend=FALSE, size=0.7) +
        geom_text(data=n_labels, aes(x=chromatine, y=20, label=lab), size=label_size) +
        geom_text(data=median_labels, aes(x=chromatine, y=-100, label=lab), size=label_size) +
        scale_fill_manual(values=pretty_colors2, breaks=levels(all_chromatine_stats$clean_name)) +
        labs(x="Chromatin subcompartments", y="Fork speed (bp/min)", fill=NULL) +
        coord_cartesian(ylim=c(-100, 2500)) +
        theme_bw() +
        theme(
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 8),
            legend.position="right",
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6),
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-5, -5, -5, -5),
            legend.key.height = unit(5, "mm"),
            legend.key.width = unit(5, "mm"),
            legend.key.spacing.y = unit(2, 'mm')
        ) +
        guides(fill=guide_legend(ncol=1, byrow=TRUE))
        
    return(list(gp1=gp1, gp2=gp2))
}

plot.initiation.bias <- function(repartition_matrix, repartition_WT, repartition_WT_all, repartition_random){
    region_names <- c("Initiation Zone (IZ)", "Non-IZ")

    # Compute p-value
    pval <- chisq.test(repartition_matrix)$p.value
    pval_label <- ifelse(pval < 2.2e-16, "Pearson's Chi-squared\np < 2.2e-16", paste0("p = ", sprintf("%.2g", pval)))

    # Prepare plot data
    df_ratio <- rbind(repartition_WT, repartition_WT_all, repartition_random) %>% 
        dplyr::filter(type_event == "Ini") %>%
        left_join(sample_clean_name, by = join_by("sample" == "original")) %>%
        mutate(clean = ifelse(is.na(clean), sample, as.character(clean)))

    # Preserve input order
    df_ratio$clean <- factor(df_ratio$clean, levels = unique(df_ratio$clean))

    # Define comparison for annotation
    df_pval <- data.frame(group1 = "Untreated - All", group2 = "Random", y.position = 1 + 0.05, p = pval, label = pval_label)

    df_ratio <- rbind(df_ratio %>% mutate(stats=region_names[1]), df_ratio %>% mutate(stats=region_names[2], ratio=1-ratio))
    df_ratio$stats <- factor(df_ratio$stats, levels = rev(unique(df_ratio$stats)))

    label_size <- 5/.pt
    y_offset_top <- 0.05
    y_offset_bot <- 0.04
    gp <- ggplot(df_ratio, aes(x = clean, y = ratio)) +
        geom_bar(aes(fill = stats), width = 0.7, stat="identity", position="fill") +
        geom_text(aes(y=ifelse(stats==region_names[1], 0 + y_offset_bot, 1 - y_offset_top), label=ifelse(stats==region_names[1], IZ, non_IZ)), size=label_size) + # ratio + y_offset_top
        scale_fill_manual(values=setNames(c("#C78383", "#E8BABA"), region_names), breaks=region_names) +
        stat_pvalue_manual(df_pval, label = "label", tip.length = 0.02, size = label_size, vjust=-0.1) +
        coord_cartesian(ylim = c(0, 1.2)) +
        scale_y_continuous(breaks=seq(0, 1, 0.2), labels = function(x) scales::number(x * 100, accuracy = 1)) +
        labs(x = NULL, y = "Initiation events (%)", fill="Genomic region") +
        theme_bw() +
        theme(
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6),
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-5, -5, -5, -5),
            legend.key.height = unit(5, "mm"),
            legend.key.width = unit(5, "mm"),
            legend.key.spacing.y = unit(2, 'mm')
        )

    return(gp)
}



plot_Fig1 <- function(forkml_speed_annotation, example_reads, forkml_fork_annotation, forkml_event_annotation, length_threshold, make_raster=TRUE){
    gp1 <- plot_example(example_reads, forkml_fork_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), forkml_speed_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), forkml_event_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), 0.7, 0.7, NULL, make_raster)

    forkml_speed_annotation <- forkml_speed_annotation %>%
        dplyr::filter(read_length > length_threshold)
    
    n_labels <- forkml_speed_annotation %>%
        group_by(full_name) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("n = ", n), n))

    median_labels <- forkml_speed_annotation %>%
        group_by(full_name) %>%
        summarise(median_fork_speed = round(median(fork_speed, na.rm = TRUE), 0), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("Med. = ", median_fork_speed), median_fork_speed))

    label_size <- 5/.pt
    gp2 <- ggplot(forkml_speed_annotation, aes(x=full_name, y=fork_speed)) +
        geom_violin(aes(fill=full_name), linewidth=0.4, scale="width", show.legend=TRUE) +
        geom_boxplot(aes(fill=full_name), linewidth=0.2, width=0.2, outlier.shape=NA, show.legend=FALSE) +
        stat_summary(fun="mean", geom="point", color="red", size=0.7) +
        geom_text(data = n_labels, aes(x = full_name, y = 0, label = lab), col=1, size = label_size) +
        geom_text(data = median_labels, aes(x = full_name, y = -200, label = lab), col=1, size = label_size) +
        scale_fill_manual(values=forkml_speed_annotation %>% distinct(full_name, mycol) %>% { setNames(.$mycol, .$full_name) }, breaks=levels(forkml_speed_annotation$full_name)) +
        coord_cartesian(ylim=c(-200, 2500)) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(), # element_text(angle=45, hjust=1),
            axis.text = element_text(size = 5),
            axis.title = element_text(size = 7),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(), # element_text(size = 7),

            legend.position = "bottom",
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 5),

            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-5, -5, -5, -5),
            legend.key.height = unit(5, "mm"),
            legend.key.width = unit(5, "mm"),
            legend.key.spacing.y = unit(2, 'mm')
        ) +
        labs(y="Fork speed (bp/min)", fill=NULL)

    blank <- ggplot() + theme_void()
    layout <- (blank / gp1  / (blank + gp2) / plot_spacer()) + plot_layout(heights = c(1, 6, 2.5, 2)) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size = 7))
    ggsave("Fig1.pdf", layout, width = 180, height = 210, units = "mm")

    # return(gp1)
}
plot_Fig1(forkml_speed_annotation_hct116 %>% filter(!is_outlier), example_HCT116_UT_R9_rep1, forkml_fork_annotation_hct116, forkml_event_annotation_hct116, length_threshold, FALSE)

plot_Fig2 <- function(forkml_speed_annotation, replication_timing, all_chromatine_stats, repartition_matrix, repartition_WT, repartition_WT_all, repartition_random){
    gp1 <- plot.relationship.forkspeed.raw_timing.paper(forkml_speed_annotation, replication_timing) +
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
        )

    binned_expanded_timing <- rebinning.replication.timing(replication_timing, 100000)
    list_gp <- plot.fork.speed.chromatine.timing.paper(forkml_speed_annotation, binned_expanded_timing, all_chromatine_stats)

    gp4 <- plot.initiation.bias(repartition_matrix, repartition_WT, repartition_WT_all, repartition_random) +
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6)
        )

    layout_design <- "
        AB
        CD
        FF
        FF
    "
    layout <- gp1 + list_gp$gp2 + free(list_gp$gp1) + gp4 +  plot_spacer() +
        plot_layout(design=layout_design) +
        plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size = 7))
    ggsave("Fig2.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_Fig2(forkml_speed_annotation_hct116 %>% dplyr::filter(read_length > length_threshold & grepl("Untreated", full_name) & !is_outlier), replication_timing_v1, all_chromatine_stats, repartition_matrix, repartition_WT, repartition_WT_all, repartition_random)

plot_Fig3 <- function(forkml_speed_annotation, df_length_threshold){
    sample_order <- unique(forkml_speed_annotation$full_name)
    forkml_speed_annotation <- forkml_speed_annotation %>%
        left_join(df_length_threshold, by=join_by(sample)) %>%
        dplyr::filter(read_length > min_read_length & !is_outlier)

    forkml_speed_annotation <- forkml_speed_annotation %>%
        mutate(full_name=factor(full_name, levels=sample_order))

    n_labels <- forkml_speed_annotation %>%
        group_by(full_name) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("n = ", n), n))

    median_labels <- forkml_speed_annotation %>%
        group_by(full_name) %>%
        summarise(median_fork_speed = round(median(fork_speed, na.rm = TRUE), 0), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("Med. = ", median_fork_speed), median_fork_speed))

    label_size <- 5/.pt
    gp1 <- ggplot(forkml_speed_annotation, aes(x=full_name, y=fork_speed)) +
        geom_violin(aes(fill=full_name), linewidth=0.4, scale="width", show.legend=TRUE) +
        geom_boxplot(aes(fill=full_name), linewidth=0.2, width=0.2, outlier.shape=NA, show.legend=FALSE) +
        stat_summary(fun="mean", geom="point", color="red", size=0.7) +
        geom_text(data = n_labels, aes(x = full_name, y = 0, label = lab), col=1, size = label_size) +
        geom_text(data = median_labels, aes(x = full_name, y = -200, label = lab), col=1, size = label_size) +
        scale_fill_manual(values=forkml_speed_annotation %>% distinct(full_name, mycol) %>% { setNames(.$mycol, .$full_name) }, breaks=levels(forkml_speed_annotation$full_name)) +
        coord_cartesian(ylim=c(-200, 2500)) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle=45, hjust=1),
            axis.text = element_text(size = 5),
            axis.title = element_text(size = 7),
            axis.title.x = element_blank(),

            legend.position = "none"
        ) +
        labs(y="Fork speed (bp/min)", fill=NULL)

    blank <- ggplot() + theme_void()
    layout <- ((gp1 + blank) / plot_spacer()) + plot_layout(heights = c(2.5, 9))
    ggsave("Fig3.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_Fig3(rbind(forkml_speed_annotation_r10 %>% dplyr::filter(full_name %in% c("Untreated (4) R10","HeLa Untreated (1)")) %>% mutate(full_name=gsub(" Untreated \\(1\\)","",gsub("Untreated \\(4\\) R10","HCT116",as.character(full_name)))), forkml_speed_annotation_bonus %>% mutate(full_name=gsub(" Untreated","",as.character(full_name)))), df_length_threshold)

plot_FigS1 <- function(example_amp){
    df_offsetting <- data.frame(
        read_id=c("3e24582d-2e0b-4af6-a3b5-975a883452c5","5a05d420-5f20-4589-95cb-3de632211b60","ee3f62a1-6961-4777-82b5-a65ba3b98cd9","d4626891-abe5-4f95-a11a-05925a90a61b","82a24dcc-450e-4625-b3e2-ebbafd8bf7eb"),
        centering=c(0, 0, 0, 0, 0),
        direction=c(1, 1, 1, 1, 1)
    ) %>% mutate(read_id=paste0("read_", read_id))

    levels(example_amp$clean_read_id) <- gsub("uM","µM",levels(example_amp$clean_read_id))

    example_amp_fixed <- example_amp %>%
        left_join(df_offsetting, by=join_by(read_id==read_id)) %>%
        mutate(cen_pos=rel_pos - centering) %>%
        mutate(cen_pos=ifelse(direction == -1, cen_pos * -1, cen_pos))
    example_amp_fixed$clean_read_id <- factor(example_amp_fixed$clean_read_id, levels=c("No BrdU", "0.1 µM BrdU", "1 µM BrdU", "10 µM BrdU", "100 µM BrdU"))

    linewidth_target <- 0.4
    gp1 <- ggplot(example_amp_fixed) +
        rasterise(geom_point(aes(x=cen_pos, y=avg, col="100 bp average"), pch=16, alpha=0.5, size=0.05), dev="cairo", dpi=600) +
        geom_line(aes(x=cen_pos, y=smooth_signal_norm, col="Gaussian smoothing (2500 bp)"), linewidth=linewidth_target) +
        facet_grid(clean_read_id ~ ., space="free_x") +
        theme_bw() +
        scale_color_manual(name="BrdU signal:", values=mypal) +
        scale_x_continuous(limits=c(-3000, 62000), expand=expansion(mult=c(0.01, 0), add=c(0, 1))) +
        scale_y_continuous(expand=expansion(mult=c(0.05, 0.05))) +
        labs(x="Relative position within reads (bp)", y="Predicted BrdU content") +
        theme(legend.position="bottom") +
        guides(color = guide_legend(override.aes = list(linewidth = 1, size = 2)))

    blank <- plot_spacer()
    layout <- (gp1 / blank) + plot_layout(heights = c(10, 5))
    ggsave("FigS1.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS1(example_amp)

plot_FigS2 <- function(example_ipd){
    df_offsetting <- data.frame(
        read_id=c("52a6ec83-d600-4d83-b7aa-dc8f44b5769c", "a68fb01d-4640-44e2-9492-1eb14dd4761a", "ad1e9410-04e1-4253-8071-ef119391c041", "51a19d54-5ab5-4b26-a6e3-f80320115f25"),
        centering=c(18500, 15800, 77500, 78500),
        direction=c(1, 1, -1, -1)
    ) %>% mutate(read_id=paste0("read_", read_id))

    example_ipd_fixed <- example_ipd %>%
        left_join(df_offsetting, by=join_by(read_id==read_id)) %>%
        mutate(cen_pos=rel_pos - centering) %>%
        mutate(cen_pos=ifelse(direction == -1, cen_pos * -1, cen_pos)) %>%
        mutate(clean_read_id=gsub("interpulse","inter-pulse",clean_read_id))

    linewidth_target <- 0.4
    gp1 <- ggplot(example_ipd_fixed) +
        rasterise(geom_point(aes(x=cen_pos, y=avg, col="100 bp average"), pch=16, alpha=0.5, size=0.05), dev="cairo", dpi=600) +
        geom_line(aes(x=cen_pos, y=smooth_signal_norm, col="Gaussian smoothing (2500 bp)"), linewidth=linewidth_target) +
        facet_grid(clean_read_id ~ ., space="free_x") +
        theme_bw() +
        scale_color_manual(name="BrdU signal:", values=mypal) +
        scale_x_continuous(limits=c(-3000, 62000), expand=expansion(mult=c(0.01, 0), add=c(0, 1))) +
        scale_y_continuous(expand=expansion(mult=c(0.05, 0.05))) +
        labs(x="Position centered on first pulse (bp)", y="Predicted BrdU content") +
        theme(legend.position="bottom") +
        guides(color = guide_legend(override.aes = list(linewidth = 1, size = 2)))

    blank <- plot_spacer()
    layout <- (gp1 / blank) + plot_layout(heights = c(10, 5))
    ggsave("FigS2.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS2(example_ipd)

plot_FigS4 <- function(forkml_fork_annotation, rfd_data, rfd_bin_size, nb_subgroups){
    merged_data <- rfd_data %>%
        mutate(n=nb_w + nb_c, rel_n=(n*100)/mean(n)) %>%
        dplyr::filter(rel_n > 20 & rel_n < 200)

    gr_merged_data <- GRanges(
      seqnames=merged_data$chr,
      ranges=IRanges(start=merged_data$start, end=merged_data$start + rfd_bin_size -1),
      rfd=merged_data$rfd
    )

    gr_all_forks <- GRanges(
      seqnames=forkml_fork_annotation$read_chrm,
      ranges=IRanges(start=forkml_fork_annotation$pred_inf, end=forkml_fork_annotation$pred_sup),
      strand=forkml_fork_annotation$read_strand
    )

    overlaps <- findOverlaps(gr_all_forks, gr_merged_data)

    forkml_fork_annotation_res <- forkml_fork_annotation[queryHits(overlaps),]
    forkml_fork_annotation_res$rfd <- mcols(gr_merged_data)$rfd[subjectHits(overlaps)]
    forkml_fork_annotation_res$width <- width(pintersect(gr_all_forks[queryHits(overlaps)], gr_merged_data[subjectHits(overlaps)]))

    res <- forkml_fork_annotation_res %>%
        group_by(full_name, mycol, read_id, pred_inf, pred_sup, read_chrm, read_strand, direction) %>%
        summarize(avg_rfd=mean(rfd), w_avg_rfd=sum(rfd*width)/sum(width), .groups="drop") %>%
        mutate(avg_rfd_bin=cut(avg_rfd, breaks=seq(-1, 1, length.out=21)), w_avg_rfd_bin=cut(w_avg_rfd, breaks=seq(-1, 1, length.out=21)))

    res <- res %>% # Select weighted avg.
        mutate(rfd_bin=w_avg_rfd_bin) %>% dplyr::select(-avg_rfd_bin, -w_avg_rfd_bin)

    set.seed(101)
    res_sum <- res %>%
        ungroup() %>%
        dplyr::select(full_name, mycol, read_id, pred_inf, pred_sup, read_chrm, read_strand, rfd_bin, direction) %>%
        group_by(full_name, mycol, rfd_bin) %>%
        mutate(rdm_id = sample(1:nb_subgroups, size = n(), replace = TRUE)) %>%  # Generate random IDs within each rfd_bin group
        group_by(full_name, mycol, rfd_bin, rdm_id) %>%
        summarize(n=n(), score=mean(direction, na.rm = TRUE), .groups="drop") %>%
        dplyr::filter(n >= 10 & !is.na(rfd_bin))

    gp <- ggplot(res_sum, aes(x=rfd_bin, y=score, fill=mycol)) +
        geom_boxplot(outlier.shape=NA) +
        geom_jitter(height=0, width=0.3, size=0.4, col="grey") +
        geom_hline(yintercept=0) +
        facet_wrap(vars(full_name)) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        scale_y_continuous(limits=c(-1, 1)) +
        scale_x_discrete(drop = FALSE) +
        scale_fill_identity() +
        labs(x="10 kbp genomic bins from GLOE-seq RFD", y="Binned RFD from double BrdU pulse-labelling", fill="Sample")

    blank <- plot_spacer()
    layout <- (gp / blank) + plot_layout(heights = c(10, 5))
    ggsave("FigS4.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS4(forkml_fork_annotation_hct116 %>% dplyr::filter(grepl("Untreated", full_name) & !is_outlier), gloe_seqx2, 10000, 5)

plot_FigS6 <- function(df_speed_manual, df_speed_auto){
    df_all <- rbind(df_speed_manual, df_speed_auto)

    n_labels <- df_all %>%
        group_by(full_name, data_type) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("n = ", n), n))

    median_labels <- df_all %>%
        group_by(full_name, data_type) %>%
        summarise(median_fork_speed = round(median(fork_speed, na.rm = TRUE), 0), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("Med. = ", median_fork_speed), median_fork_speed))

    label_size <- 6/.pt
    plot_offsets <- 0.93
    gp1 <- ggplot(df_all, aes(x=full_name, y=fork_speed, group=interaction(full_name, data_type))) +
        geom_violin(aes(lty=data_type, fill=mycol), linewidth=0.4, position=position_dodge(width=plot_offsets), scale="width") +
        geom_boxplot(aes(fill=mycol), linewidth=0.2, width=0.2, outlier.shape=NA, position=position_dodge(width=plot_offsets), show.legend=FALSE) + # lty=data_type, 
        geom_text(data = n_labels, aes(x = full_name, y = 0, label = lab), size = label_size, position=position_dodge(width=plot_offsets)) +
        geom_text(data = median_labels, aes(x = full_name, y = -200,  label = lab), size = label_size, position=position_dodge(width=plot_offsets)) +
        stat_summary(fun="mean", geom="point", color="red", position=position_dodge(width=plot_offsets)) +
        scale_fill_identity() +
        coord_cartesian(ylim=c(-200, 2500)) +
        labs(x="Samples", y="Fork speed (bp/min)", lty="Annotation type") +
        theme_bw() +
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7),
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6),
            legend.position=c(0.95, 0.95),
            legend.justification=c(1, 1),
            legend.box.background=element_rect(color="black", fill=NA)
        )

    blank <- plot_spacer()
    layout <- (gp1 / blank / blank)
    ggsave("FigS6.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS6(df_speed_manual, df_speed_auto)

plot_FigS7 <- function(forkml_speed_annotation, individual_sample){
    plots <- make.length.plots(forkml_speed_annotation, individual_sample)

    blank <- plot_spacer()
    layout <- ((plots$gp1) / (plots$gp2 + plots$gp3) / blank) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size = 7))
    ggsave("FigS7.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS7(forkml_speed_annotation_hct116 %>% dplyr::filter(!is_outlier & !is.na(full_name)), "Untreated \\(1\\)")

plot_FigS8 <- function(example_reads_a, example_reads_b, forkml_fork_annotation, forkml_speed_annotation, forkml_event_annotation, forced_arrow){
    gp1 <- plot_example(example_reads_a, forkml_fork_annotation %>% dplyr::filter(sample %in% unique(example_reads_a$simple)), forkml_speed_annotation %>% dplyr::filter(sample %in% unique(example_reads_a$simple)), forkml_event_annotation %>% dplyr::filter(sample %in% unique(example_reads_a$simple)), 0.1, 0.7, "HCT116 HU 100 µM (1)", TRUE, forced_arrow)
    gp2 <- plot_example(example_reads_b, forkml_fork_annotation %>% dplyr::filter(sample %in% unique(example_reads_b$simple)), forkml_speed_annotation %>% dplyr::filter(sample %in% unique(example_reads_b$simple)), forkml_event_annotation %>% dplyr::filter(sample %in% unique(example_reads_b$simple)), 0.1, 0.7, "HCT116 APH 100 nM (1)")

    layout <- (gp1 / gp2 / plot_spacer()) + plot_layout(heights = c(5, 5, 1)) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size = 7))
    ggsave("FigS8.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS8(example_HCT116_HU_R9_rep1, example_HCT116_APH_R9_rep1, forkml_fork_annotation_hct116, forkml_speed_annotation_hct116, forkml_event_annotation_hct116, FALSE)

plot_FigS9 <- function(forkml_speed_annotation){
    df_fork_speed <- forkml_speed_annotation %>%
        mutate(fork_type=ifelse(fork_type=="leading","Leading","Lagging")) %>%
        mutate(fork_type=factor(fork_type, levels = c("Leading", "Lagging")))
        df_fork_speed

    n_labels <- df_fork_speed %>%
        group_by(full_name, fork_type, mycol) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("n = ", n), n))

    median_labels <- df_fork_speed %>%
        group_by(full_name, fork_type, mycol) %>%
        summarise(median_fork_speed = round(median(fork_speed, na.rm = TRUE), 0), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("Med. = ", median_fork_speed), median_fork_speed))

    label_size <- 5/.pt
    plot_offsets <- 0.93
    gp1 <- ggplot(df_fork_speed, aes(x=full_name, y=fork_speed, group=interaction(full_name, fork_type), fill=mycol)) +
        geom_violin(aes(lty=fork_type), linewidth=0.4, position=position_dodge(width=plot_offsets), scale="width") +
        geom_boxplot(linewidth=0.2, width=0.2, outlier.shape=NA, position=position_dodge(width=plot_offsets), show.legend=FALSE) +
        stat_summary(fun="mean", geom="point", color="red", size=0.2, position=position_dodge(width=plot_offsets)) +
        geom_text(data=n_labels, aes(x=full_name, y=0, label=lab), size=label_size, position=position_dodge(width=plot_offsets)) +
        geom_text(data=median_labels, aes(x=full_name, y=-200, label=lab), size=label_size, position=position_dodge(width=plot_offsets)) +
        scale_fill_identity() +
        scale_linetype_manual(values=c("Leading"=1, "Lagging"=2), breaks=c("Leading", "Lagging")) +
        coord_cartesian(ylim=c(-200, 2500)) +
        theme_bw() +
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6),
            legend.position="bottom"
        ) +
        guides(lty=guide_legend(override.aes=list(col=c("#030BFC", "#D99DFC")))) +
        labs(x=NULL, y="Fork speed (bp/min)", lty="Strand")

    df_fork_speed_sum <- df_fork_speed %>%
        group_by(full_name, fork_type, mycol) %>%
        summarize(mean_fork_speed=mean(fork_speed), .groups="drop") %>%
        mutate(sample_type=gsub(" \\([0-9]\\)","",full_name))  %>%
        mutate(sample_type=factor(sample_type, levels=unique(sample_type)))

    df_pval <- data.frame(sample_type=unique(df_fork_speed_sum$sample_type), group1="Leading", group2="Lagging", y.position=1500, label="0.37\n(n.s.)")

    df_fork_speed_sum_stat <- df_fork_speed_sum %>%
        group_by(sample_type, fork_type) %>% summarize(mean_fork_speed=mean(mean_fork_speed), .groups="drop")

    gp2 <- ggplot(df_fork_speed_sum, aes(x=sample_type, y=mean_fork_speed, col=fork_type)) +  
        geom_beeswarm(aes(group=fork_type), pch=21, show.legend=FALSE, size=0.25, dodge.width=plot_offsets, cex=4) +
        geom_point(data=df_fork_speed_sum_stat, aes(group=fork_type), color="red", pch=95, size=3, position=position_dodge(width=plot_offsets)) +
        geom_text(data=data.frame(sample_type="HU 100 µM", y=2000, label="p-value", fork_type="Leading"), aes(x=sample_type, y=y, label=label, col="black"), size=label_size, show.legend=FALSE) +
        stat_pvalue_manual(df_pval, x="sample_type", label="label", size=label_size) +
        scale_color_manual(values=c("Leading"="#030BFC","Lagging"="#D99DFC","black"=1)) +
        coord_cartesian(ylim=c(-200, 2500)) +
        theme_bw() +
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(x=NULL, y=NULL, lty=NULL)

    layout_design <- "
        AAAAAAAAB
        CCCCCCCCC
        CCCCCCCCC
        CCCCCCCCC
    "
    layout <- gp1 + gp2 + plot_spacer() +
        plot_layout(design=layout_design) +
        plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size = 7))
    ggsave("FigS9.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS9(forkml_speed_annotation_hct116 %>% dplyr::filter(read_length > length_threshold & !is_outlier))

plot_FigS10 <- function(example_reads, forkml_fork_annotation, forkml_speed_annotation, forkml_event_annotation, df_speed_manual_r10, forkml_speed_annotation_hct116, length_threshold){
    gp1 <- plot_example(example_reads, forkml_fork_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), forkml_speed_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), forkml_event_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), 0.1, 0.7, "Untreated (4) R10")

    df_manual <- df_speed_manual_r10 %>% dplyr::filter(read_length > length_threshold & grepl("good$", subset_name)) %>% dplyr::select(full_name, mycol, data_type, read_id, fork_speed)
    df_auto_r10 <- forkml_speed_annotation %>% mutate(read_id=gsub("read_","",read_id)) %>% dplyr::filter(read_id %in% as.character(unique(df_manual$read_id))) %>% dplyr::select(full_name, mycol, data_type, read_id, fork_speed)
    df_auto_r10_all <- forkml_speed_annotation %>% dplyr::filter(read_length > length_threshold & !grepl("\\(0\\)", as.character(full_name)) & !is_outlier) %>% dplyr::select(full_name, mycol, data_type, read_id, fork_speed)
    df_auto_r9_all <- forkml_speed_annotation_hct116 %>% dplyr::filter(read_length > length_threshold & !is_outlier) %>% dplyr::select(full_name, mycol, data_type, read_id, fork_speed)

    df_all <- rbind(
        rbind(df_manual, df_auto_r10) %>% 
            mutate(full_name = paste0(full_name, " - subset")), 
        rbind(df_auto_r10_all, df_auto_r9_all) %>% 
            mutate(full_name = paste0(full_name, " - full dataset"))
    )
    df_all$data_type <- factor(df_all$data_type, levels=unique(df_all$data_type), ordered = TRUE)
    df_all$full_name <- factor(df_all$full_name, levels=unique(df_all$full_name), ordered = TRUE)

    n_labels <- df_all %>%
        group_by(full_name, data_type) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("n = ", n), n))

    median_labels <- df_all %>%
        group_by(full_name, data_type) %>%
        summarise(median_fork_speed = round(median(fork_speed, na.rm = TRUE), 0), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("Med. = ", median_fork_speed), median_fork_speed))

    plot_offsets <- 0.93
    label_size <- 7/.pt
    gp2 <- ggplot(df_all, aes(x=full_name, y=fork_speed, lty=data_type)) +
        geom_violin(aes(fill=mycol), linewidth=0.4, position=position_dodge(width=plot_offsets), scale="width") +
        geom_boxplot(aes(fill=mycol), linewidth=0.2, width=0.2, outlier.shape=NA, position=position_dodge(width=plot_offsets), show.legend=FALSE) +
        stat_summary(fun="mean", geom="point", color="red", position=position_dodge(width=plot_offsets), size=0.2) +
        geom_text(data = n_labels, aes(x = full_name, y = 0, label = lab), size = label_size, position=position_dodge(width=plot_offsets)) +
        geom_text(data = median_labels, aes(x = full_name, y = -200, label = lab), size = label_size, position=position_dodge(width=plot_offsets)) +
        scale_fill_identity() +
        scale_linetype_manual(values=c("Manual"="dashed", "Automated"="solid")) +
        coord_cartesian(ylim=c(-200, 2500)) +
        labs(x=NULL, y="Fork speed (bp/min)", lty="Annotation type") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 5),
            axis.title = element_text(size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6)
        )

    top_height <- 3
    layout_design <- c(
        area(t=1, b=top_height, l=1, r=7),  # gp1
        area(t=top_height + 1, b=top_height + 2, l=1, r=4),  # gp2
        area(t=top_height + 3, b=top_height + 3, l=1, r=7)
    )
    layout <- wrap_plots(gp1, gp2, design = layout_design) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size = 7))
    ggsave("FigS10.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS10(example_HCT116_UT_R10_rep4, forkml_fork_annotation_r10 %>% dplyr::filter(!grepl("HeLa", sample)), forkml_speed_annotation_r10 %>% dplyr::filter(!grepl("HeLa", sample)), forkml_event_annotation_r10 %>% dplyr::filter(!grepl("HeLa", sample)), df_speed_manual_r10, forkml_speed_annotation_hct116 %>% dplyr::filter(full_name=="Untreated (4)") %>% mutate(full_name="Untreated (4) R9"), length_threshold)

plot_FigS11 <- function(example_reads, forkml_fork_annotation, forkml_speed_annotation, forkml_event_annotation, df_manual_hela, length_threshold_hela){
    gp1 <- plot_example(example_reads, forkml_fork_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), forkml_speed_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), forkml_event_annotation %>% dplyr::filter(sample %in% unique(example_reads$simple)), 0.1, 0.7, "HeLa Untreated (1)")

    df_manual_hela <- df_manual_hela %>% dplyr::filter(read_length > length_threshold_hela & grepl("all$", subset_name)) %>% dplyr::select(full_name, mycol, data_type, read_id, fork_speed) 
    df_auto_hela <- forkml_speed_annotation %>% dplyr::mutate(read_id=gsub("read_","",read_id)) %>% dplyr::filter(read_id %in% as.character(unique(df_manual_hela$read_id))) %>% dplyr::select(full_name, mycol, data_type, read_id, fork_speed)
    df_auto_hela_all <- forkml_speed_annotation %>% dplyr::filter(read_length > length_threshold & !is_outlier) %>% dplyr::select(full_name, mycol, data_type, read_id, fork_speed)

    df_all_hela <- rbind(
        rbind(df_manual_hela, df_auto_hela) %>% 
            mutate(full_name = paste0(full_name, " - subset")), 
        rbind(df_auto_hela_all) %>% 
            mutate(full_name = paste0(full_name, " - full dataset"))
    )
    df_all_hela$data_type <- factor(df_all_hela$data_type, levels=unique(df_all_hela$data_type), ordered = TRUE)
    df_all_hela$full_name <- gsub(" R10", "", df_all_hela$full_name)
    df_all_hela$full_name <- factor(df_all_hela$full_name, levels=unique(df_all_hela$full_name), ordered = TRUE)

    n_labels <- df_all_hela %>%
        group_by(full_name, mycol, data_type) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("n = ", n), n))

    median_labels <- df_all_hela %>%
        group_by(full_name, mycol, data_type) %>%
        summarise(median_fork_speed = round(median(fork_speed, na.rm = TRUE), 0), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("Med. = ", median_fork_speed), median_fork_speed))

    plot_offsets <- 0.93
    label_size <- 7/.pt
    gp2 <- ggplot(df_all_hela, aes(x=full_name, y=fork_speed, lty=data_type)) +
        geom_violin(aes(fill=mycol), linewidth=0.4, position=position_dodge(width=plot_offsets), scale="width") +
        geom_boxplot(aes(fill=mycol), linewidth=0.2, width=0.2, outlier.shape=NA, position=position_dodge(width=plot_offsets), show.legend=FALSE) +
        stat_summary(fun="mean", geom="point", color="red", position=position_dodge(width=plot_offsets), size=0.2) +

        geom_text(data = n_labels, aes(x = full_name, y = 0, label = lab), size = label_size, position=position_dodge(width=plot_offsets)) +
        geom_text(data = median_labels, aes(x = full_name, y = -200, label = lab), size = label_size, position=position_dodge(width=plot_offsets)) +

        scale_fill_identity() +
        scale_linetype_manual(values=c("Manual"="dashed", "Automated"="solid")) +
        coord_cartesian(ylim=c(-200, 2500)) +
        labs(x=NULL, y="Fork speed (bp/min)", lty="Annotation type") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 5),
            axis.title = element_text(size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6)
        )

    top_height <- 3
    layout_design <- c(
        area(t=1, b=top_height, l=1, r=7),  # gp1
        area(t=top_height + 1, b=top_height + 2, l=1, r=4),  # gp2
        area(t=top_height + 3, b=top_height + 3, l=1, r=7)
    )
    layout <- wrap_plots(gp1, gp2, design = layout_design) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size = 7))
    ggsave("FigS11.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS11(example_HeLa_UT_R10_rep1, forkml_fork_annotation_r10 %>% dplyr::filter(grepl("HeLa", sample)), forkml_speed_annotation_r10 %>% dplyr::filter(grepl("HeLa", sample)), forkml_event_annotation_r10 %>% dplyr::filter(grepl("HeLa", sample)), df_speed_manual_hela, length_threshold_hela)

plot_FigS12 <- function(example_reads, forkml_fork_annotation_bonus, forkml_speed_annotation_bonus, forkml_event_annotation_bonus, df_speed_manual_bonus, df_length_threshold){
    gp1 <- plot_example(example_reads, forkml_fork_annotation_bonus %>% dplyr::filter(!is_outlier & sample %in% unique(example_reads$simple)), forkml_speed_annotation_bonus %>% dplyr::filter(!is_outlier & sample %in% unique(example_reads$simple)), forkml_event_annotation_bonus %>% dplyr::filter(!is_outlier & sample %in% unique(example_reads$simple)), 0.1, 0.7, "Additional cell lines", potluck=TRUE)

    df_manual_bonus <- df_speed_manual_bonus %>%
        left_join(df_length_threshold, by=join_by(sample)) %>%
        dplyr::filter(read_length > min_read_length & grepl("all$", subset_name)) %>%
        dplyr::select(full_name, mycol, data_type, read_id, fork_speed) 
    df_auto_bonus <- forkml_speed_annotation_bonus %>%
        dplyr::mutate(read_id=gsub("read_","",read_id)) %>%
        dplyr::filter(read_id %in% as.character(unique(df_manual_bonus$read_id))) %>%
        dplyr::select(full_name, mycol, data_type, read_id, fork_speed)
    df_auto_bonus_all <- forkml_speed_annotation_bonus %>%
        left_join(df_length_threshold, by=join_by(sample)) %>%
        dplyr::filter(read_length > min_read_length & !is_outlier) %>%
        dplyr::select(full_name, mycol, data_type, read_id, fork_speed)

    df_all_bonus <- rbind(
        rbind(df_manual_bonus, df_auto_bonus) %>% 
            mutate(full_name = paste0(full_name, " - subset")), 
        rbind(df_auto_bonus_all) %>% 
            mutate(full_name = paste0(full_name, " - full dataset"))
    )
    df_all_bonus$data_type <- factor(df_all_bonus$data_type, levels=unique(df_all_bonus$data_type), ordered = TRUE)
    df_all_bonus$full_name <- gsub(" R10", "", df_all_bonus$full_name)
    df_all_bonus$full_name <- factor(df_all_bonus$full_name, levels=unique(df_all_bonus$full_name), ordered = TRUE)

    base_order <- levels(df_speed_manual_bonus$full_name)
    lvl_long <- levels(df_all_bonus$full_name)
    new_levels <- unlist(lapply(base_order, function(b) {grep(paste0("^", b), lvl_long, value = TRUE)}))
    df_all_bonus$full_name <- factor(df_all_bonus$full_name, levels = new_levels)

    n_labels <- df_all_bonus %>%
        group_by(full_name, mycol, data_type) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("n = ", n), n))

    median_labels <- df_all_bonus %>%
        group_by(full_name, mycol, data_type) %>%
        summarise(median_fork_speed = round(median(fork_speed, na.rm = TRUE), 0), .groups = "drop") %>%
        mutate(lab=ifelse(row_number()==1, paste0("Med. = ", median_fork_speed), median_fork_speed))

    plot_offsets <- 0.93
    label_size <- 5/.pt
    gp2 <- ggplot(df_all_bonus, aes(x=full_name, y=fork_speed, lty=data_type)) +
        geom_violin(aes(fill=mycol), linewidth=0.4, position=position_dodge(width=plot_offsets), scale="width") +
        geom_boxplot(aes(fill=mycol), linewidth=0.2, width=0.2, outlier.shape=NA, position=position_dodge(width=plot_offsets), show.legend=FALSE) +
        stat_summary(fun="mean", geom="point", color="red", position=position_dodge(width=plot_offsets), size=0.2) +

        geom_text(data = n_labels, aes(x = full_name, y = 0, label = lab), size = label_size, position=position_dodge(width=plot_offsets)) +
        geom_text(data = median_labels, aes(x = full_name, y = -200, label = lab), size = label_size, position=position_dodge(width=plot_offsets)) +

        scale_fill_identity() +
        scale_linetype_manual(values=c("Manual"="dashed", "Automated"="solid")) +
        coord_cartesian(ylim=c(-200, 2500)) +
        labs(x=NULL, y="Fork speed (bp/min)", lty="Annotation type") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 5),
            axis.title = element_text(size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6)
        )

    top_height <- 5
    layout_design <- c(
        area(t=1, b=top_height, l=1, r=7),  # gp1
        area(t=top_height + 1, b=top_height + 2, l=1, r=7),  # gp2
        area(t=top_height + 3, b=top_height + 3, l=1, r=7)
    )
    layout <- wrap_plots(gp1, gp2, design = layout_design) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', size = 7))
    ggsave("FigS12.pdf", layout, width = 180, height = 210, units = "mm")
}
plot_FigS12(example_bonus, forkml_fork_annotation_bonus, forkml_speed_annotation_bonus, forkml_event_annotation_bonus, df_speed_manual_bonus, df_length_threshold)


