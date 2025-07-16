#!/usr/bin/env Rscript
# version: 1.1.2

library(shiny)
library(shinyWidgets)

# suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(purrr)))

suppressWarnings(suppressMessages(library(scales)))
library(markdown)
suppressWarnings(suppressMessages(library(filesstrings)))
library(rdrop2)
suppressWarnings(suppressMessages(library(formattable)))
suppressWarnings(suppressMessages(library(DT)))
suppressWarnings(suppressMessages(library(glue)))
`%+%` = paste0

################### Parameters ###################
nb_max_reads = 50
nb_tracks = 6  # Do not change. If modified, the UI in trackInput function must be modified accordingly.
fork_colors = c("red"="#f44336", "navy"="#076cb2", "goldenrod"="#f5b342")

# Platform - "dropbox" or "local"
data_storage = "local"

# Local parameters
if(data_storage == "local"){
    file_samples = getwd() %+% '/samples.csv'
    dir_samples = getwd() %+% '/'
    dir_dt_tmp = getwd() %+% '/reads_annotating/'
    dir_dt_out = getwd() %+% '/reads_done/'
    dir_tks_out = getwd() %+% '/tracks_done/'

    # Create structure if needed
    if(!dir.exists(dir_dt_tmp)){dir.create(dir_dt_tmp, recursive=TRUE)}
    if(!dir.exists(dir_dt_out)){dir.create(dir_dt_out, recursive=TRUE)}
    if(!dir.exists(dir_tks_out)){dir.create(dir_tks_out, recursive=TRUE)}
}else if(data_storage == "dropbox"){
    file_samples = 'samples.csv'
    dir_samples ='dataset/'                 # must end with a '/'
    dir_dt_out = 'dataset_annotated/'     # must end with a '/'
    dir_dt_tmp = 'dataset_annotating/'    # must end with a '/'
    dir_tks_out = 'tracks_annotated/'     # must end with a '/'
    droptoken = readRDS('droptoken.rds')
}

################### Preprocessing ###################
# List all read files and extract sample names
samples = str_match(list.files(path=dir_samples, pattern="data_"), "data_(.+)")[,2]

# Read list of prepared datasets
samples_data = read_csv(file_samples, col_names=TRUE, col_types="cccccc") %>%
    filter(sample %in% samples) %>% 
    mutate(label = case_when(
        !is.na(notes) ~ sample %+%" - "%+% strain %+%", "%+% pulse_duration %+%", "%+% pulse_conc %+%", "%+% notes,
        is.na(notes) ~ sample %+%" - "%+% strain %+%", "%+% pulse_duration %+%", "%+% pulse_conc,
    ))
sample_choices = samples_data$sample
names(sample_choices) = samples_data$label

vect_persons = unique(samples_data$sample) %>% sort()
names(vect_persons) = vect_persons

################### Functions ###################
load_data_local = function(dir_dt_in, dir_dt_tmp, file) {
    # load data for local use
    file.move(dir_dt_in%+%file, dir_dt_tmp)
    readRDS(dir_dt_tmp%+%file) %>% ungroup()
}

load_data_dropbox = function (dir_dt_in, dir_dt_tmp, file) {
    # download & read dataset
    is_downloaded = drop_download(dir_dt_in%+%file, dtoken=droptoken, overwrite=TRUE)
    drop_move(from_path=dir_dt_in%+%file, to_path=dir_dt_tmp%+%file, dtoken=droptoken)
    readRDS(file)
}

save_data_local = function(data, dir_tks_out, person, dir_dt_tmp, dir_dt_out, file) {
    # filename
    input_name = str_match(file, "(.+).rds$")[2]
    filename = input_name %+% "_" %+% format(Sys.time(),"%F_%T") %+% '_' %+% person %+% '.rds'
    
    # save results of the read
    saveRDS(data, dir_tks_out %+%'/'%+% filename)

    # move dataset to "annotated"
    file.move(dir_dt_tmp%+%file, dir_dt_out)
}

save_data_dropbox = function(data, dir_tks_out, person, dir_dt_tmp, dir_dt_out, file) {
    # filename
    input_name = str_match(file, "(.+).rds$")[2]
    filename = input_name %+% "_" %+% format(Sys.time(),"%F_%T") %+% '_' %+% person %+% '.rds'
    # save results of the read & upload on dropbox
    saveRDS(data, filename)
    drop_upload(filename, path=dir_tks_out, dtoken=droptoken)
    # move dataset to "annotated"
    drop_move(from_path=dir_dt_tmp%+%file, to_path=dir_dt_out%+%file, dtoken=droptoken)
}

collect_data = function(rv_table_tracks, person) {
    rv_list <- reactiveValuesToList(rv_table_tracks)

    tracks_table_all <- map(rv_list, ~ .x()) %>% 
        bind_rows() %>%
        mutate(person=person, date_time=format(Sys.time(), "%F %T %z"))

    return(tracks_table_all)
}

plot_init = function() {
    out <- ggplot(NULL) +
        theme_classic() +
        theme(
            legend.position="none",
            axis.title=element_blank(),
            axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.x=element_blank(), axis.ticks.x=element_line(color="#888888") # TODO maybe remove tick too
        ) +
        scale_color_identity()

    return(out)
}

plot_signal <- function(read_data, is_mirrored) {
    ticks_step <- 2000
    signal <- read_data$signal[[1]] %>% arrange(x)

    if("annotation" %in% colnames(read_data)){
        if(!is.null(read_data$annotation[[1]])){ # Maybe dup
            read_inf <- min(read_data$read_inf, read_data$annotation[[1]]$pred_inf, read_data$annotation[[1]]$pred_sup, na.rm=TRUE)
            read_sup <- max(read_data$read_sup, read_data$annotation[[1]]$pred_inf, read_data$annotation[[1]]$pred_sup, na.rm=TRUE)
        }else{
            read_inf <- read_data$read_inf
            read_sup <- read_data$read_sup   
        }
    }else{
        read_inf <- read_data$read_inf
        read_sup <- read_data$read_sup   
    }

    ticks_start <- round(read_inf/ticks_step) * ticks_step
    ticks_end <- round(read_sup/ticks_step) * ticks_step

    if("sample" %in% colnames(signal)){ # Currently only for rDNA
        min_size_window <- 80000
        out <- list(
            geom_point(data=signal %>% dplyr::filter(!is.na(avg)), mapping=aes(x, avg, group=sample_col, color=sample_col), size=0.5, alpha=0.2),
            geom_line(data=signal, mapping=aes(x, signal, color=sample_col)),
            scale_y_continuous(limits=c(0, 1), expand=c(0.02,0), sec.axis=dup_axis())
        )
    }else{
        min_size_window <- 140000
        out <- list(
            geom_point(data=signal %>% dplyr::filter(!is.na(avg)), mapping=aes(x, avg), color="#000000", pch=46),
            geom_line(data=signal, mapping=aes(x, y), color="#19CF22"),
            scale_y_continuous(limits=c(0, 1), expand=c(0.02,0), sec.axis=dup_axis())
        )        
    }

    if(is_mirrored){
        out <- c(out, scale_x_reverse(
            limits=c(read_sup, min(read_sup - min_size_window, read_inf)), # Min plot length is around 65kb , na.rm=TRUE
            expand=c(0,0),
            breaks=seq(ticks_end, read_inf, by=-ticks_step),
            labels=comma,
            sec.axis=dup_axis()
        ))
    }else{
        out <- c(out, scale_x_continuous(
            limits=c(read_inf, max(read_inf + min_size_window, read_sup)), # Min plot length is around 65kb
            expand=c(0,0),
            breaks=seq(ticks_start, read_sup, by=ticks_step),
            labels=comma,
            sec.axis=dup_axis()
        ))
    }

    return(out)
}

plot_ranges <- function(fork_data, fork_colors, is_mirrored){
    fork_data <- fork_data %>% dplyr::filter(!is.na(x2))

    if(nrow(fork_data)==0){
        out <- list()
    }else{
        fork_data <- fork_data %>%
            rowwise() %>%
            mutate(right_fork=ifelse(is_mirrored, x0>x2, x0<x2), unknown=ifelse((x1<x2 & x2<x0) || (x0<x2 & x2<x1), TRUE, FALSE), arrow_start=x0, arrow_end=ifelse(unknown, x1, x2), rec_start=min(x0,x1,x2), rec_end=max(x0,x1,x2), mid_pos=ifelse(unknown, x2, x1), fork_color=ifelse(unknown, fork_colors[3], ifelse(right_fork, fork_colors[1], fork_colors[2])))
        # print(fork_data)

        out <- list(
            geom_rect(
                data=fork_data,
                mapping=aes(xmin=rec_start, xmax=rec_end, ymin=0, ymax=1, color=fork_color, fill=fork_color),
                alpha=0
            ),
            geom_point(
                data=fork_data,
                mapping=aes(x=mid_pos, y=0.8, color=fork_color)
            ),
            geom_text(
                data=fork_data,
                mapping=aes(x=mid_pos, y=0.9, label=fork_id, color=fork_color, size=0.6)
            ))
        
        known_forks <- fork_data %>% dplyr::filter(!unknown)
        unknown_forks <- fork_data %>% dplyr::filter(unknown)
        if(nrow(known_forks)>0){
            out <- c(out, geom_segment(
                data=known_forks,
                mapping=aes(x=arrow_start, xend=arrow_end, y=0.8, yend=0.8, color=fork_color),
                size=0.8,
                arrow=arrow(length=unit(2, "mm")),
                lineend="butt", linejoin="bevel"
            ))            
        }
        if(nrow(unknown_forks)>0){
            out <- c(out, geom_segment(
                data=unknown_forks %>% dplyr::filter(unknown),
                mapping=aes(x=arrow_start, xend=arrow_end, y=0.8, yend=0.8, color=fork_color),
                size=0.8,
                arrow=arrow(length=unit(2, "mm"), ends="both"),
                lineend="butt", linejoin="bevel"
            ))            
        }
    }
    return(out)
}

show_info_page = function() {
    showModal(modalDialog(
      easyClose=TRUE, size="l",
      includeMarkdown("ABOUT.md")
    ))
}

show_thanks_page = function() {
    showModal(modalDialog(
      title="Thanks", easyClose=FALSE, size="m",
      "Thank you for your participation. You can now close the window, or refresh the page to continue annotating reads.",
      footer=NULL
    ))
}

################### Shiny Modules ###################
# read module -------------------------------------------------
my_readUI <- function(id){
    ns <- NS(id)
    tagList(
        plotOutput(outputId=ns("signal"), click=ns("plot_click"), height=400),
        wellPanel(
            actionButton(inputId=ns("toggle_tracks"), label="Toggle tracks", icon=icon("arrows-alt-v")),
            conditionalPanel(
                condition="input.toggle_tracks % 2 == 1", ns=ns,
                br(),
                verbatimTextOutput(ns("clicked")), # TODO replace by table output
                dataTableOutput(outputId=ns("table"))
            )
        )
    )
}

my_read <- function(input, output, session, read_data, session_id) {
    cat(paste0("Rendering read: ",read_data$read_id,"\n"))
    
    if(runif(1) > 0.5){
        is_mirrored <- TRUE
    }else{
        is_mirrored <- FALSE
    }

    # read plot
    g <- plot_init()
    sig <- plot_signal(read_data, is_mirrored)
    graph_forks <- list()

    read_inf <- read_data$read_inf
    read_length <- read_data$read_sup - read_data$read_inf

    output$clicked <- renderPrint({
        req(input$plot_click)
        paste("Clicked coordinates - x:", round(input$plot_click$x, 0), "y:", round(input$plot_click$y, 2))
    })
    
    init_track <- data.frame(read_id=read_data$read_id, sample=read_data$sample, read_chrm=read_data$read_chrm, read_strand=read_data$read_strand, read_start=read_data$read_start, read_end=read_data$read_end, fork_id=1, x0=NA, x1=NA, x2=NA, y0=NA, y1=NA, y2=NA)
    init_track <- cbind(init_track, data.frame(
            delete = glue("<button rowid=1 onclick='Shiny.setInputValue(\"{session_id}-removeRow\",this.getAttribute(\"rowid\"))'>Delete</button>"),
            rowid = 1,
            stringsAsFactors = FALSE
        )
    )
    empty_track <- init_track
    default_buttons <- NULL
    if("annotation" %in% colnames(read_data)){
        read_ann_data <- read_data$annotation[[1]]
        if(!is.null(read_ann_data)){
            list_tracks <- seq_along(read_ann_data$direction)
            auto_track <- data.frame(read_id=read_data$read_id, sample=read_data$sample, read_chrm=read_data$read_chrm, read_strand=read_data$read_strand, read_start=read_data$read_start, read_end=read_data$read_end, fork_id=list_tracks, x0=ifelse(read_ann_data$direction==1, read_ann_data$pred_inf, read_ann_data$pred_sup), x1=(read_ann_data$pred_inf + read_ann_data$pred_sup)/2, x2=ifelse(read_ann_data$direction==1, read_ann_data$pred_sup, read_ann_data$pred_inf), y0=1.1, y1=1.1, y2=1.1)
            default_buttons <- sapply(list_tracks, function(i) {
                glue::glue("<button rowid={i} onclick='Shiny.setInputValue(\"{session_id}-removeRow\", this.getAttribute(\"rowid\"))'>Delete</button>")
            })
            auto_track <- cbind(auto_track, data.frame(
                    delete = default_buttons,
                    rowid = list_tracks,
                    stringsAsFactors = FALSE
                )
            )
            init_track <- auto_track
        }
    }
    colnames(init_track)[ncol(init_track)-1] <- " "
    colnames(empty_track)[ncol(empty_track)-1] <- " "

    table_tracks <- reactiveValues(data=init_track)
    table_tracks_ready <- reactiveValues(data=init_track)
    observeEvent(input$plot_click, {
        nb_missing_value <- sum(is.na(table_tracks$data[nrow(table_tracks$data),]))
        nb_forks <- nrow(table_tracks$data)
        fork_id <- max(table_tracks$data$fork_id) # TODO issue here depending when line is removed

        # print(nb_missing_value)
        if(nb_missing_value==0){
            new_idx_fork <- fork_id + 1
            nb_forks <- nb_forks + 1

            empty_track[1, ncol(empty_track)-1] <- gsub("button rowid=[0-9]+ ", paste0("button rowid=",new_idx_fork," "), empty_track[1, ncol(empty_track)-1])
            empty_track$rowid[1] <- new_idx_fork

            table_tracks$data <- bind_rows(table_tracks$data, empty_track) # Maybe needs same names
            table_tracks$data$fork_id[nb_forks] <- new_idx_fork
        }else if(nb_missing_value==7){
            table_tracks$data$fork_id[nb_forks] <- fork_id
        }

        nb_missing_value <- sum(is.na(table_tracks$data[nrow(table_tracks$data),]))
        idx_fork <- nrow(table_tracks$data)
        idx_coordinate <- abs((nb_missing_value/2) - 3)

        table_tracks$data[[paste0("x", idx_coordinate)]][idx_fork] <- round(input$plot_click$x, 0)
        table_tracks$data[[paste0("y", idx_coordinate)]][idx_fork] <- round(input$plot_click$y, 2)

        nb_missing_value <- sum(is.na(table_tracks$data[nrow(table_tracks$data),]))
        if(nb_missing_value==0){ # Block not clean
            table_tracks_ready$data <- table_tracks$data
        }
    })

    observeEvent(input$removeRow, {
        removeRow <- as.integer(input$removeRow)
        table_tracks$data <- table_tracks$data[-which(table_tracks$data$rowid == removeRow),]
        table_tracks_ready$data <- table_tracks$data
        if(nrow(table_tracks$data)==0){
            table_tracks$data <- empty_track
        }
    })

    table_tracks_render <- reactive({
        table_tracks$data %>%
            dplyr::select(-read_id, -read_chrm, -read_strand, -read_start, -read_end)
    })

    table_tracks_results <- reactive({
        table_tracks$data %>%
            dplyr::select(-" ", -rowid)
    })

    output$table <- renderDataTable(
        datatable(
            table_tracks_render(),
            escape=FALSE, # Allow HTML in cells
            selection="none",
            options=list(columnDefs = list(list(targets=ncol(table_tracks_render()), visible=FALSE))) # Hide "delete" column
        )
    )

    graph_forks <- reactive({plot_ranges(table_tracks_ready$data, fork_colors, is_mirrored)})
    graph <- reactive({g + graph_forks() + sig})
    output$signal <- renderPlot({graph()}, res=100)
    
    return(table_tracks_results)
}

################### Shiny application ###################
ui <- fluidPage(
    # for correcting modal dialog that keep adding 13px on the right padding. 
    # https://stackoverflow.com/questions/32862394/bootstrap-modals-keep-adding-padding-right-to-body-after-closed
    tags$head(tags$style(HTML("body:not(.modal-open){padding-right: 0px !important;}"))), 
    
    # Header
    fluidRow(
        column(10, offset=1, h1("ForkML annotation", align="center", style="font-family:helvetica,verdana")),
        column(1,  offset=0, br(), div(
            actionButton(inputId="info", label="About", class="btn-secondary"), 
            style="float:right"
        ))
    ),

    # Sample selection
    conditionalPanel(
        condition="input.load == 0",
        column(4, offset=4, wellPanel(
            p("Select the sample you want to annotate"),
            selectInput(inputId="sample", label="Sample", choices=sample_choices),
            actionButton(inputId="load", label="Load data", class="btn-secondary", icon=icon('upload'))
        )) # ,
    ),
    
    # Reads
    uiOutput('reads'),
    
    # Footer/save data
    conditionalPanel(
        condition="input.load >= 1",
        column(4, offset=4, wellPanel(
            p("Please select your name in the list to sign your annotations"),
            selectInput(inputId="person", label="Name", choices=c("", vect_persons)),
            actionButton(inputId="save", label="Save Annotations", class="btn-primary", icon=icon('check'))
        )) # ,
    ),
    
    # Blank space
    div(style="margin-bottom:800px;")
)

# Server function ---------------------------------------------
server = function(input, output, session) {
    # Info button -> Info page
    observeEvent(input$info, show_info_page())

    # Reactive list for storing file and directory names when loading
    sample_path <- reactiveValues(dir_dt_in=NULL, file=NULL)

    # Load button -> Load data + "Loading" modal dialog
    observeEvent(input$load, once=TRUE, {
        # Show modal loading...
        showModal(modalDialog(easyClose=FALSE, size="m", footer=NULL, h3("Loading...")))

        # Define sample
        dir_dt_in <- dir_samples %+% "data_" %+% input$sample %+% "/reads_toannotate/"

        cat("Loading data.\n")
        if(data_storage == "local"){
            file <- list.files(dir_dt_in)[1]
            reads_data <- load_data_local(dir_dt_in, dir_dt_tmp, file)
        }else if(data_storage == "dropbox"){
            file <- drop_dir(dir_dt_in, dtoken=droptoken)[[1, "name"]]
            reads_data <- load_data_dropbox(dir_dt_in, dir_dt_tmp, file)
        }
        cat(paste0(file,"\n"))
        sample_path$file <- file
        sample_path$dir_dt_in <- dir_dt_in

        # Removing rare reads with long query but short mapping
        reads_data <- reads_data %>% dplyr::filter(!is.infinite(read_inf) & !is.infinite(read_sup))

        # nb_reads <- nrow(reads_data %>% head(1)) # For testing
        nb_reads <- nrow(reads_data)
        if(nb_reads > nb_max_reads){
            cat("Too many reads in input files, review data preparation.\n")
        }

        cat("Call all read modules.\n")
        rv_table_tracks <- reactiveValues()
        for(i in 1:min(nb_max_reads, nb_reads)){
            id <- 'read'%+%i
            rv_table_tracks[[id]] <- callModule(my_read, id, read_data=reads_data %>% slice(i), session_id=id) # tweak by read_id
        }
        
        # RenderUI : render all reads
        cat("Render all reads.\n")
        list_readUI <- tagList()
        for(i in 1:min(nb_max_reads, nb_reads)){
            if(i==1){
                list_readUI <- tagList(list_readUI, tagList(my_readUI('read'%+%i)))
                insertUI(selector='#reads', where='beforeEnd', ui=tagList(my_readUI('read'%+%i)), immediate=TRUE)
            }else{
                insertUI(selector='#reads', where='beforeEnd', ui=tagList(my_readUI('read'%+%i)), immediate=TRUE)
            }
        }

        # remove modal loading...
        removeModal()
        
        cat("Waiting for input.\n")

        # Save button -> Save data + "Thanks" modal dialog
        observeEvent(input$save, once=TRUE, {
            # collect data
            tracks_table_all <- collect_data(rv_table_tracks, input$person)
            
            # save data
            showModal(modalDialog(easyClose=FALSE, size="m", footer=NULL, h3("Saving data...")))
            if(data_storage=="local"){
                save_data_local(tracks_table_all, dir_tks_out, input$person, dir_dt_tmp, dir_dt_out, file)
            }else if(data_storage=="dropbox"){
                save_data_dropbox(tracks_table_all, dir_tks_out, input$person, dir_dt_tmp, dir_dt_out, file)
            }
            removeModal()
            
            # "Thanks" modal dialog
            show_thanks_page()
        })
    })
    
    # Move input data back to dataset if it was not saved
    session$onSessionEnded(function(){
        if (isolate(input$load)>=1) {
            if (data_storage=="local" & !isolate(as.logical(input$save))) {
                file = isolate(sample_path$file)
                dir_dt_in = isolate(sample_path$dir_dt_in)
                file.move(dir_dt_tmp %+% file, dir_dt_in)
            }else if(data_storage=="dropbox" & !isolate(as.logical(input$save))) {
                file = isolate(sample_path$file)
                dir_dt_in = isolate(sample_path$dir_dt_in)
                drop_move(from_path=dir_dt_tmp%+%file, to_path=dir_dt_in%+%file, dtoken=droptoken)
            }
        }
    })
}

# Start Shiny app
shinyApp(ui=ui, server=server, options=list("port"=9999))
