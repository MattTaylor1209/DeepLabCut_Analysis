# app.R
# Shiny app for DeepLabCut filtered.csv batch processing + plots + stats

library(shiny)
library(shinyFiles)
library(shinythemes)
library(tidyverse)
library(stringr)
library(slider)
library(DT)
library(multcomp)
library(broom)
library(colourpicker)
library(plotly)

source("R/helpers_dlc.R")

options(shiny.maxRequestSize = 200 * 1024^2)  # 200 MB


# ----------------------------
# Shiny app
# ----------------------------

# To add
# Angles/headswings
# toggle to invert y axis
# Time moving
# Distance traveled
# Rolling?
# Speed decay over time - are they getting tired?
# Pause duration

ui <- fluidPage(
  theme = shinytheme("simplex"),
  titlePanel("DeepLabCut CSV Explorer (filtered.csv)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$h4("1) Get data"),
      
      radioButtons(
        "data_source",
        "Data source",
        choices = c(
          "Local folder (desktop)" = "folder",
          "Upload filtered.csv files (shinyapps.io)" = "upload"
        ),
        selected = "folder"
      ),
      
      conditionalPanel(
        condition = "input.data_source == 'folder'",
        shinyDirButton("dir", "Select folder", "Choose folder containing *filtered.csv"),
        verbatimTextOutput("folder_txt")
      ),
      
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput(
          "upload_files",
          "Upload one or more *filtered.csv* files",
          multiple = TRUE,
          accept = c(".csv")
        )
      ),
      
      tags$hr(),
      
      
      tags$h4("2) Processing settings"),
      textInput("pattern", "File pattern", value = "filtered.csv"),
      
      radioButtons(
        "dlc_format",
        "DLC CSV format",
        choices = c(
          "Multi-animal (4 header rows: scorer/individuals/bodyparts/coords)" = "multi",
          "Single-animal (flat header)" = "single",
          "Auto-detect per file" = "auto"
        ),
        selected = "multi"
      ),
      selectizeInput(
        "bodyparts_keep",
        "Bodyparts to keep",
        choices = NULL,
        selected = "mid",
        multiple = TRUE,
        options = list(plugins = list("remove_button"))
      ),
      numericInput("likelihood_min", "Min likelihood (optional)", value = NA, min = 0, max = 1, step = 0.01),
      numericInput("threshold_px", "Movement threshold (px/frame)", value = 2, min = 0, step = 0.1),
      numericInput("window_n", "Stillness window (frames)", value = 10, min = 1, step = 1),
      checkboxInput("compute_angle", "Compute angle (requires exactly 3 bodyparts)", value = TRUE),
      uiOutput("angle_vertex_ui"),
      
      tags$h4("Jump filtering"),
      checkboxInput("drop_big_jumps", "Drop big jumps in plots", value = TRUE),
      numericInput("max_jump_px", "Max allowed step distance (px)", value = 10, min = 0, step = 5),
      checkboxInput("use_robust_jump", "Use robust per-track threshold (median * multiplier)", value = FALSE),
      numericInput("robust_mult", "Robust multiplier", value = 10, min = 1, step = 1),
      checkboxInput("drop_bad_excursions", "Drop short excursions between big jumps", value = TRUE),
      numericInput("max_excursion_frames", "Max excursion length (frames)", value = 10, min = 1, step = 1),
      
      tags$h4("Speed QC"),
      numericInput(
        "max_speed",
        "Max speed to keep (px/frame)",
        value = 20,
        min = 0,
        step = 1
      ),
      checkboxInput("drop_over_max_speed", "Drop speeds above max in plots + summary", value = TRUE),
      
      
      tags$hr(),
      
      
      actionButton("run", "Load & process", class = "btn-primary"),
      
      tags$hr(),
      
      tags$h4("3) Plot settings"),
      numericInput("x_min", "x min", value = 0),
      numericInput("x_max", "x max", value = 1280),
      numericInput("y_min", "y min", value = 0),
      numericInput("y_max", "y max", value = 960),
      checkboxInput("flip_y", "Flip Y axis (image coordinates)", value = TRUE),
      
      tags$h4("Trajectory styling"),
      colourInput("moving_col", "Moving colour", value = "#11F011"),
      colourInput("still_col",  "Still colour",  value = "#FA05EE"),
      numericInput("line_width", "Line width", value = 0.6, min = 0.1, step = 0.1),
      selectInput("legend_pos", "Legend position",
                  choices = c("top","bottom","left","right","none"),
                  selected = "top"),
      
      tags$h4("Text sizes"),
      numericInput("axis_title_size", "Axis title size", value = 13, min = 6, step = 1),
      numericInput("axis_text_size",  "Axis text size",  value = 11, min = 6, step = 1),
      numericInput("strip_text_size", "Facet strip size", value = 13, min = 6, step = 1),
      
      tags$br(), tags$br(),
      uiOutput("file_picker_ui")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Files",
                 tags$h4("Detected files"),
                 DTOutput("files_table")
        ),
        tabPanel("Trajectory",
                 plotlyOutput("traj_plotly", height = 1000),
                 downloadButton("download_traj", "Download trajectory plot (PNG)")
        ),
        tabPanel("Trajectory (moving only) NOT IN USE",
                 plotOutput("traj_moving_plot", height = 700),
                 downloadButton("download_traj_moving", "Download moving-only plot (PNG)")
        ),
        tabPanel("Speed trace",
                 plotOutput("speed_plot", height = 700),
                 downloadButton("download_speed", "Download speed plot (PNG)")
        ),
        tabPanel(
          "Angle",
          plotOutput("angle_plot", height = 700),
          downloadButton("download_angle", "Download angle plot (PNG)")
        ),
        tabPanel("Summary",
                 tags$h4("Mean moving speed per individual (mid by default)"),
                 downloadButton("download_summary", "Download summary CSV"),
                 DTOutput("summary_table")
        ),
        tabPanel("Dunnett vs reference",
                 uiOutput("ref_group_ui"),
                 checkboxInput("use_median_dunnett", "Use median for stats", value = FALSE),
                 actionButton("run_dunnett", "Run Dunnett"),
                 tags$br(), tags$br(),
                 DTOutput("dunnett_table")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  volumes <- c(Home = path.expand("~"), "R()" = R.home(), getVolumes()())
  shinyDirChoose(input, "dir", roots = volumes, session = session)
  
  chosen_dir <- reactiveVal(NULL)
  
  uploaded_paths <- reactive({
    req(input$upload_files)
    input$upload_files$datapath   # temp file paths on the server
  })
  
  uploaded_names <- reactive({
    req(input$upload_files)
    input$upload_files$name
  })
  
  
  observeEvent(input$dir, {
    dir <- parseDirPath(volumes, input$dir)
    if (length(dir) && nzchar(dir)) chosen_dir(dir)
  })
  
  output$folder_txt <- renderText({
    d <- chosen_dir()
    if (is.null(d)) "No folder selected yet."
    else d
  })
  
  files_df <- reactive({
    src <- input$data_source
    
    if (src == "upload") {
      req(input$upload_files)
      tibble(
        file_name = uploaded_names(),
        file_path = uploaded_paths()
      ) %>% arrange(file_name)
      
    } else {
      d <- chosen_dir()
      req(d)
      pat <- input$pattern
      paths <- list.files(d, pattern = pat, full.names = TRUE)
      tibble(
        file_name = basename(paths),
        file_path = paths
      ) %>% arrange(file_name)
    }
  })
  
  
  bodyparts_available <- reactive({
    df <- files_df()
    req(nrow(df) > 0)
    
    # union of bodyparts across all files (only reads header lines, so it's fast)
    bps <- unique(unlist(lapply(df$file_path, get_dlc_bodyparts, format = input$dlc_format)))
    sort(bps)
  })
  
  observeEvent(bodyparts_available(), {
    choices <- bodyparts_available()
    
    # keep current selection if possible; otherwise default to "mid" if present
    current <- isolate(input$bodyparts_keep)
    if (is.null(current) || length(current) == 0) {
      current <- if ("mid" %in% choices) "mid" else choices[1]
    } else {
      current <- intersect(current, choices)
      if (length(current) == 0) current <- if ("mid" %in% choices) "mid" else choices[1]
    }
    
    updateSelectizeInput(session, "bodyparts_keep",
                         choices = choices,
                         selected = current,
                         server = TRUE)
  }, ignoreInit = FALSE)
  
  
  output$files_table <- renderDT({
    req(files_df())
    datatable(files_df(), options = list(pageLength = 10))
  })
  
  output$file_picker_ui <- renderUI({
    df <- files_df()
    req(nrow(df) > 0)
    selectInput("selected_file", "Select file to plot",
                choices = df$file_name, selected = df$file_name[[1]])
  })
  
  processed_all <- reactiveVal(NULL)
  
  output$angle_vertex_ui <- renderUI({
    req(input$bodyparts_keep)
    if (length(input$bodyparts_keep) != 3) return(NULL)
    selectInput("angle_vertex", "Vertex bodypart (point B)", choices = input$bodyparts_keep,
                selected = input$bodyparts_keep[2])
  })
  
  
  observeEvent(input$run, {
    df <- files_df()
    req(nrow(df) > 0)
    
    bodyparts_keep <- input$bodyparts_keep
    if (is.null(bodyparts_keep) || length(bodyparts_keep) == 0) bodyparts_keep <- NULL
    
    likelihood_min <- if (is.na(input$likelihood_min)) NULL else input$likelihood_min
    
    shiny::withProgress(message = "Processing files...", value = 0, {
      out <- purrr::map2_dfr(df$file_path, df$file_name, ~{
        shiny::incProgress(1 / nrow(df))
        
        process_one_file(
          path = .x,
          file_name = .y,   # ✅ preserves original upload name
          dlc_format = input$dlc_format,
          bodyparts_keep = if (length(bodyparts_keep)) bodyparts_keep else NULL,
          likelihood_min = likelihood_min,
          threshold_px   = input$threshold_px,
          window_n       = input$window_n,
          max_jump_px     = input$max_jump_px,
          use_robust_jump = input$use_robust_jump,
          robust_mult     = input$robust_mult,
          max_excursion_frames = input$max_excursion_frames,
          compute_angle   = isTRUE(input$compute_angle),
          angle_vertex    = input$angle_vertex
        )
      })
      
      processed_all(out)
    })
    
    showNotification(
      paste0("Done. Processed ", n_distinct(processed_all()$file_name), " file(s)."),
      type = "message",
      duration = 4
    )
  }, ignoreInit = TRUE)
  
  
  processed_selected <- reactive({
    dat <- processed_all()
    req(dat, input$selected_file)
    dat %>% filter(file_name == input$selected_file)
  })
  
  summary_df <- reactive({
    dat <- processed_all()
    req(dat)
    
    dat2 <- dat
    if ("mid" %in% unique(dat2$bodypart)) dat2 <- dat2 %>% filter(bodypart == "mid")
    
    dat2 <- apply_speed_cap(
      dat2,
      max_speed = input$max_speed,
      enabled = input$drop_over_max_speed
    )
    
    dat2 %>%
      filter(moving %in% TRUE, !(is_big_jump %in% TRUE)) %>%
      group_by(group, file_name, individual) %>%
      summarise(mean_speed = mean(speed, na.rm = TRUE), 
                median_speed = median(speed, na.rm = TRUE),
                .groups = "drop")
  })
  
  
  output$summary_table <- renderDT({
    req(summary_df())
    datatable(summary_df(), options = list(pageLength = 10))
  })
  
  output$ref_group_ui <- renderUI({
    df <- summary_df()
    req(df)
    grps <- sort(unique(df$group))
    selectInput("ref_group", "Reference group", choices = grps, selected = grps[[1]])
  })
  
  dunnett_res <- eventReactive(input$run_dunnett, {
    df <- summary_df()
    req(df, input$ref_group)
    
    df <- df %>%
      mutate(group = factor(group, levels = c(input$ref_group, setdiff(sort(unique(group)), input$ref_group))))
    
    if(input$use_median_dunnett == TRUE) {
      aov_fit <- aov(median_speed ~ group, data = df)
    }
    else {
      aov_fit <- aov(mean_speed ~ group, data = df)
    }
    gl <- multcomp::glht(aov_fit, linfct = multcomp::mcp(group = "Dunnett"))
    broom::tidy(gl)
  }, ignoreInit = TRUE)
  
  output$dunnett_table <- renderDT({
    dat <- dunnett_res()
    req(dat)
    datatable(dat, options = list(pageLength = 10))
  })
  
  output$traj_plotly <- renderPlotly({
    df <- processed_selected()
    req(nrow(df) > 0)
    
    x_range <- c(input$x_min, input$x_max)
    y_range <- if (isTRUE(input$flip_y)) c(input$y_max, input$y_min) else c(input$y_min, input$y_max)
    
    
    ttl <- paste0("Trajectory (moving segments coloured): ", unique(df$title))
    
    g <- plot_trajectory_coloured_segments(
      df_long = df,
      xlim = x_range, ylim = y_range,
      ttl = ttl,
      moving_col = input$moving_col,
      still_col  = input$still_col,
      line_width = input$line_width,
      axis_title_size = input$axis_title_size,
      axis_text_size  = input$axis_text_size,
      strip_text_size = input$strip_text_size,
      legend_position = if (input$legend_pos == "none") "none" else input$legend_pos,
      drop_big_jumps = input$drop_big_jumps
    )
    
    p <- ggplotly(g, dynamicTicks = TRUE) %>%
      layout(dragmode = "zoom") %>%
      config(
        displaylogo = FALSE,
        scrollZoom = TRUE,
        modeBarButtonsToAdd = c("zoom2d","pan2d","select2d","lasso2d","resetScale2d")
      )
    
    set_all_subplot_ranges(p, x_range = x_range, y_range = y_range)
  })
  
  
  output$traj_moving_plot <- renderPlot({
    df <- processed_selected()
    req(nrow(df) > 0)
    dfm <- df %>% filter(moving %in% TRUE, !(is_big_jump %in% TRUE))
    validate(need(nrow(dfm) > 0, "No moving frames after threshold/window settings."))
    xlim <- c(input$x_min, input$x_max)
    ylim <- c(input$y_min, input$y_max)
    ttl <- paste0("Trajectory (moving): ", unique(df$title))
    plot_trajectory(dfm, xlim = xlim, ylim = ylim, ttl = ttl)
  })
  
  output$speed_plot <- renderPlot({
    df <- processed_selected()
    req(nrow(df) > 0)
    
    df <- apply_speed_cap(
      df,
      max_speed = input$max_speed,
      enabled = input$drop_over_max_speed
    )
    
    ttl <- paste0("Speed trace: ", unique(df$title))
    plot_speed_trace(df, ttl = ttl)
  })
  
  
  output$download_summary <- downloadHandler(
    filename = function() paste0("summary_mean_moving_speed_", Sys.Date(), ".csv"),
    content = function(file) {
      readr::write_csv(summary_df(), file)
    }
  )
  
  output$download_traj <- downloadHandler(
    filename = function() paste0("trajectory_", tools::file_path_sans_ext(input$selected_file), ".png"),
    content = function(file) {
      df <- processed_selected()
      xlim <- c(input$x_min, input$x_max)
      ylim <- c(input$y_min, input$y_max)
      
      p <- plot_trajectory_coloured_segments(
        df_long = df,
        xlim = xlim, ylim = ylim,
        ttl = paste0("Trajectory (moving segments coloured): ", unique(df$title)),
        moving_col = input$moving_col,
        still_col  = input$still_col,
        line_width = input$line_width,
        axis_title_size = input$axis_title_size,
        axis_text_size  = input$axis_text_size,
        strip_text_size = input$strip_text_size,
        legend_position = if (input$legend_pos == "none") "none" else input$legend_pos,
        drop_big_jumps = input$drop_big_jumps
      )
      
      ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_traj_moving <- downloadHandler(
    filename = function() paste0("trajectory_moving_", tools::file_path_sans_ext(input$selected_file), ".png"),
    content = function(file) {
      df <- processed_selected() %>% filter(moving %in% TRUE, !(is_big_jump %in% TRUE))
      xlim <- c(input$x_min, input$x_max)
      ylim <- c(input$y_min, input$y_max)
      p <- plot_trajectory(df, xlim = xlim, ylim = ylim,
                           ttl = paste0("Trajectory (moving): ", unique(df$title)))
      ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_speed <- downloadHandler(
    filename = function() paste0("speed_", tools::file_path_sans_ext(input$selected_file), ".png"),
    content = function(file) {
      df <- processed_selected()
      req(nrow(df) > 0)
      
      df <- apply_speed_cap(
        df,
        max_speed = input$max_speed,
        enabled = input$drop_over_max_speed
      )
      
      p <- plot_speed_trace(df, ttl = paste0("Speed trace: ", unique(df$title)))
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )
  
  output$angle_plot <- renderPlot({
    req(input$compute_angle)
    df <- processed_selected()
    req(nrow(df) > 0)
    
    validate(need("angle" %in% names(df), "Angle column not found (did you re-run processing?)"))
    df_ang <- df %>% filter(!is.na(angle))
    
    validate(need(nrow(df_ang) > 0, "No angle values available (check bodyparts/likelihood)."))
    
    ttl <- paste0("Angle over time: ", unique(df$title))
    plot_angle_trace(df_ang, ttl = ttl)
  })
  
  
  output$download_angle <- downloadHandler(
    filename = function() paste0("angle_", tools::file_path_sans_ext(input$selected_file), ".png"),
    content = function(file) {
      req(input$compute_angle)
      df <- processed_selected() %>% filter(!is.na(angle))
      ttl <- paste0("Angle over time: ", unique(df$title))
      p <- plot_angle_trace(df, ttl = ttl)
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )
  
}

shinyApp(ui, server)
