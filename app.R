# app.R
# ============================================================================
# Shiny app for DeepLabCut filtered.csv batch processing, plotting, and stats
#
# Features:
#   - Batch import from local folder or file upload
#   - Auto-detection of DLC CSV formats (multi, single, flat)
#   - Interactive trajectory plots (plotly) with moving/still colouring
#   - Speed traces, body angle plots, head-swing detection
#   - Summary statistics with optional unit conversion
#   - Dunnett's test against a reference group
#
# Dependencies loaded below. Helper functions are in R/helpers_dlc.R.
# ============================================================================


# --- Load packages ---
library(shiny)
library(shinyFiles)
library(shinythemes)
library(tidyverse)     # includes dplyr, ggplot2, tidyr, purrr, readr, stringr, etc.
library(slider)        # for sliding window calculations
library(DT)            # interactive data tables
library(multcomp)      # Dunnett's test via glht()
library(broom)         # tidy model output
library(colourpicker)  # colour picker input widgets
library(plotly)        # interactive plots

# --- Load helper functions ---
source("R/helpers_dlc.R")

# --- App-level settings ---
options(shiny.maxRequestSize = 200 * 1024^2)  # allow uploads up to 200 MB


# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  theme = shinytheme("simplex"),
  titlePanel("DeepLabCut CSV Explorer (filtered.csv)"),
  
  sidebarLayout(
    
    # ------------------------------------------------------------------
    # Sidebar: data source, processing settings, plot settings
    # ------------------------------------------------------------------
    sidebarPanel(
      width = 3,
      
      # --- Section 1: Data source ---
      tags$h4("1) Get data"),
      
      radioButtons(
        "data_source", "Data source",
        choices = c(
          "Local folder (desktop)"                   = "folder",
          "Upload filtered.csv files (shinyapps.io)" = "upload"
        ),
        selected = "folder"
      ),
      
      # Show folder picker when using local data
      conditionalPanel(
        condition = "input.data_source == 'folder'",
        shinyDirButton("dir", "Select folder", "Choose folder containing *filtered.csv"),
        verbatimTextOutput("folder_txt")
      ),
      
      # Show file upload widget when deploying to shinyapps.io
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput(
          "upload_files",
          "Upload one or more *filtered.csv* files",
          multiple = TRUE,
          accept   = c(".csv")
        )
      ),
      
      tags$hr(),
      
      # --- Section 2: Processing settings ---
      tags$h4("2) Processing settings"),
      
      textInput("pattern", "File pattern", value = "filtered.csv"),
      
      radioButtons(
        "dlc_format", "DLC CSV format",
        choices = c(
          "Auto-detect"                                                    = "auto",
          "Multi-animal (4 header rows: scorer/individuals/bodyparts/coords)" = "multi",
          "Single-animal (3 header rows: scorer/bodyparts/coords)"         = "single",
          "Flat header (single row)"                                       = "flat"
        ),
        selected = "auto"
      ),
      
      selectizeInput(
        "bodyparts_keep", "Bodyparts to keep",
        choices  = NULL,
        selected = "mid",
        multiple = TRUE,
        options  = list(plugins = list("remove_button"))
      ),
      
      numericInput("likelihood_min", "Min likelihood (optional)",
                   value = NA, min = 0, max = 1, step = 0.01),
      numericInput("threshold_px", "Movement threshold (px/frame)",
                   value = 2, min = 0, step = 0.1),
      numericInput("window_n", "Stillness window (frames)",
                   value = 10, min = 1, step = 1),
      checkboxInput("compute_angle",
                    "Compute angle (requires exactly 3 bodyparts)", value = TRUE),
      uiOutput("angle_vertex_ui"),
      
      # Jump filtering settings
      tags$h4("Jump filtering"),
      checkboxInput("drop_big_jumps", "Drop big jumps in plots", value = TRUE),
      numericInput("max_jump_px", "Max allowed step distance (px)",
                   value = 10, min = 0, step = 5),
      checkboxInput("use_robust_jump",
                    "Use robust per-track threshold (median * multiplier)",
                    value = FALSE),
      numericInput("robust_mult", "Robust multiplier",
                   value = 10, min = 1, step = 1),
      checkboxInput("drop_bad_excursions",
                    "Drop short excursions between big jumps", value = TRUE),
      numericInput("max_excursion_frames", "Max excursion length (frames)",
                   value = 10, min = 1, step = 1),
      
      # Speed QC settings
      tags$h4("Speed QC"),
      numericInput("max_speed", "Max speed to keep (px/frame)",
                   value = 20, min = 0, step = 1),
      checkboxInput("drop_over_max_speed",
                    "Drop speeds above max in plots + summary", value = TRUE),
      
      tags$hr(),
      
      # --- Run button ---
      actionButton("run", "Load & process", class = "btn-primary"),
      
      tags$hr(),
      
      # --- Section 3: Plot settings ---
      tags$h4("3) Plot settings"),
      numericInput("x_min", "x min", value = 0),
      numericInput("x_max", "x max", value = 1280),
      numericInput("y_min", "y min", value = 0),
      numericInput("y_max", "y max", value = 960),
      checkboxInput("flip_y", "Flip Y axis (image coordinates)", value = TRUE),
      
      # Trajectory styling
      tags$h4("Trajectory styling"),
      colourInput("moving_col", "Moving colour", value = "#11F011"),
      colourInput("still_col",  "Still colour",  value = "#FA05EE"),
      numericInput("line_width", "Line width", value = 0.6, min = 0.1, step = 0.1),
      selectInput("legend_pos", "Legend position",
                  choices  = c("top", "bottom", "left", "right", "none"),
                  selected = "top"),
      
      # Text size controls
      tags$h4("Text sizes"),
      numericInput("axis_title_size", "Axis title size", value = 13, min = 6, step = 1),
      numericInput("axis_text_size",  "Axis text size",  value = 11, min = 6, step = 1),
      numericInput("strip_text_size", "Facet strip size", value = 13, min = 6, step = 1),
      
      tags$br(), tags$br(),
      
      # File picker for individual file plots (dynamically generated)
      uiOutput("file_picker_ui")
    ),
    
    # ------------------------------------------------------------------
    # Main panel: tabbed output
    # ------------------------------------------------------------------
    mainPanel(
      tabsetPanel(
        # Tab 1: File listing
        tabPanel("Files",
                 tags$h4("Detected files"),
                 DTOutput("files_table")
        ),
        
        # Tab 2: Interactive trajectory plot (plotly)
        tabPanel("Trajectory",
                 plotlyOutput("traj_plotly", height = 1000),
                 downloadButton("download_traj", "Download trajectory plot (PNG)")
        ),
        
        # Tab 3: Moving-only trajectory (static ggplot, currently not in active use)
        tabPanel("Trajectory (moving only) NOT IN USE",
                 plotOutput("traj_moving_plot", height = 700),
                 downloadButton("download_traj_moving", "Download moving-only plot (PNG)")
        ),
        
        # Tab 4: Speed trace
        tabPanel("Speed trace",
                 plotOutput("speed_plot", height = 700),
                 downloadButton("download_speed", "Download speed plot (PNG)")
        ),
        
        # Tab 5: Body angle over time
        tabPanel("Angle",
                 plotOutput("angle_plot", height = 700),
                 downloadButton("download_angle", "Download angle plot (PNG)")
        ),
        
        # Tab 6: Summary statistics
        tabPanel("Summary",
                 tags$h4("Mean moving speed per individual (mid by default)"),
                 checkboxInput("convertunits", "Convert pixel units to mm?", value = FALSE),
                 conditionalPanel(
                   condition = "input.convertunits == true",
                   numericInput("scalefactor", "Pixel scale factor (mm per pixel)",
                                min = 0.01, max = 10, step = 0.05, value = 0.1)
                 ),
                 checkboxInput("convertspeed", "Convert speed units to per second?",
                               value = FALSE),
                 conditionalPanel(
                   condition = "input.convertspeed == true",
                   numericInput("framerate", "Video frame rate",
                                min = 1, max = 240, step = 1, value = 5)
                 ),
                 downloadButton("download_summary", "Download summary CSV"),
                 DTOutput("summary_table")
        ),
        
        # Tab 7: Dunnett's test
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


# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {
  
  # --- File system volumes for shinyFiles folder picker ---
  volumes <- c(Home = path.expand("~"), "R()" = R.home(), getVolumes()())
  shinyDirChoose(input, "dir", roots = volumes, session = session)
  
  
  # ========================================================================
  # REACTIVE VALUES
  # ========================================================================
  
  # Store the user-selected directory path
  chosen_dir <- reactiveVal(NULL)
  
  # Store the fully processed data (all files combined)
  processed_all <- reactiveVal(NULL)
  
  
  # ========================================================================
  # DATA SOURCE: folder or upload
  # ========================================================================
  
  # Update chosen_dir when the user picks a folder via shinyFiles
  observeEvent(input$dir, {
    dir <- parseDirPath(volumes, input$dir)
    if (length(dir) && nzchar(dir)) chosen_dir(dir)
  })
  
  # Display the selected folder path
  output$folder_txt <- renderText({
    d <- chosen_dir()
    if (is.null(d)) "No folder selected yet." else as.character(d)
  })
  
  # Build a tibble of file paths + names from either the folder or uploaded files
  files_df <- reactive({
    src <- input$data_source
    
    if (src == "upload") {
      # --- Upload mode: use temp paths from fileInput ---
      req(input$upload_files)
      tibble(
        file_name = input$upload_files$name,
        file_path = input$upload_files$datapath
      ) %>%
        arrange(file_name)
      
    } else {
      # --- Folder mode: scan directory for matching files ---
      d   <- chosen_dir()
      req(d)
      pat   <- input$pattern
      paths <- list.files(d, pattern = pat, full.names = TRUE)
      
      # warn the user if no matching files were found
      validate(need(
        length(paths) > 0,
        paste0("No files matching '", pat, "' found in: ", d)
      ))
      
      tibble(
        file_name = basename(paths),
        file_path = paths
      ) %>%
        arrange(file_name)
    }
  })
  
  
  # ========================================================================
  # BODYPART DETECTION (from file headers)
  # ========================================================================
  
  # Scan headers of all detected files to find the union of bodypart names.
  # Uses map() + reduce() instead of lapply() + unlist().
  bodyparts_available <- reactive({
    df <- files_df()
    req(nrow(df) > 0)
    
    fmt <- input$dlc_format
    
    # Read bodypart names from each file header (fast, no data rows read)
    map(df$file_path, get_dlc_bodyparts, format = fmt) %>%
      reduce(union) %>%     # union across all files
      sort()
  })
  
  # Update the bodypart selector whenever available bodyparts change
  observeEvent(bodyparts_available(), {
    choices <- bodyparts_available()
    
    # Try to preserve the current selection; default to "mid" if present
    current <- isolate(input$bodyparts_keep)
    
    if (is.null(current) || length(current) == 0) {
      # Nothing selected yet — pick "mid" if available, else the first bodypart
      mid_idx <- which(tolower(choices) == "mid")
      current <- if (length(mid_idx)) choices[mid_idx[1]] else choices[1]
    } else {
      # Keep only bodyparts that still exist in the new file set
      current <- intersect(current, choices)
      if (length(current) == 0) {
        mid_idx <- which(tolower(choices) == "mid")
        current <- if (length(mid_idx)) choices[mid_idx[1]] else choices[1]
      }
    }
    
    updateSelectizeInput(
      session, "bodyparts_keep",
      choices  = choices,
      selected = current,
      server   = TRUE
    )
  }, ignoreInit = FALSE)
  
  
  # ========================================================================
  # FILE TABLE AND FILE PICKER
  # ========================================================================
  
  output$files_table <- renderDT({
    df <- files_df()
    validate(need(nrow(df) > 0, "No files detected. Select a folder or upload files."))
    datatable(df, options = list(pageLength = 10))
  })
  
  # Dynamically generate a dropdown to pick which file to plot
  output$file_picker_ui <- renderUI({
    df <- files_df()
    req(nrow(df) > 0)
    selectInput("selected_file", "Select file to plot",
                choices  = df$file_name,
                selected = df$file_name[[1]])
  })
  
  # Dynamically generate the angle vertex selector (only shown with 3 bodyparts)
  output$angle_vertex_ui <- renderUI({
    req(input$bodyparts_keep)
    if (length(input$bodyparts_keep) != 3) return(NULL)
    selectInput("angle_vertex", "Vertex bodypart (point B)",
                choices  = input$bodyparts_keep,
                selected = input$bodyparts_keep[2])
  })
  
  
  # ========================================================================
  # MAIN PROCESSING: Load & Process button
  # ========================================================================
  
  observeEvent(input$run, {
    df <- files_df()
    req(nrow(df) > 0)
    
    # Gather processing parameters from UI inputs
    bodyparts_keep <- input$bodyparts_keep
    if (is.null(bodyparts_keep) || length(bodyparts_keep) == 0) bodyparts_keep <- NULL
    
    likelihood_min <- if (is.na(input$likelihood_min)) NULL else input$likelihood_min
    
    # Wrap process_one_file with possibly() so one bad file doesn't crash the batch
    safe_process <- possibly(process_one_file, otherwise = NULL)
    
    shiny::withProgress(message = "Processing files...", value = 0, {
      
      # Process each file in parallel (paths + names), collecting results as a list
      results <- map2(df$file_path, df$file_name, ~ {
        shiny::incProgress(1 / nrow(df))
        
        safe_process(
          path                 = .x,
          file_name            = .y,
          bodyparts_keep       = bodyparts_keep,
          likelihood_min       = likelihood_min,
          threshold_px         = input$threshold_px,
          window_n             = input$window_n,
          max_jump_px          = input$max_jump_px,
          use_robust_jump      = input$use_robust_jump,
          robust_mult          = input$robust_mult,
          max_excursion_frames = input$max_excursion_frames,
          compute_angle        = isTRUE(input$compute_angle),
          angle_vertex         = input$angle_vertex,
          dlc_format           = input$dlc_format
        )
      })
    })
    
    # --- Post-processing: separate successes from failures ---
    
    # Identify which files returned NULL (failed)
    failed_names <- df$file_name[map_lgl(results, is.null)]
    
    # Bind successful results, dropping NULLs
    out <- compact(results) %>% list_rbind()
    
    # Notify user about failures
    if (length(failed_names) > 0) {
      showNotification(
        paste0(length(failed_names), " file(s) failed: ",
               paste(failed_names, collapse = ", ")),
        type     = "warning",
        duration = 10
      )
    }
    
    # Store results and notify on success
    if (!is.null(out) && nrow(out) > 0) {
      processed_all(out)
      showNotification(
        paste0("Done. Processed ", n_distinct(out$file_name),
               " of ", nrow(df), " file(s)."),
        type     = "message",
        duration = 4
      )
    } else {
      showNotification(
        "No files were processed successfully. Check the R console for details.",
        type = "error"
      )
    }
  }, ignoreInit = TRUE)
  
  
  # ========================================================================
  # REACTIVE: selected file's data (for individual plots)
  # ========================================================================
  
  processed_selected <- reactive({
    dat <- processed_all()
    req(dat, input$selected_file)
    
    selected <- dat %>% filter(file_name == input$selected_file)
    
    validate(need(
      nrow(selected) > 0,
      paste0("No data for '", input$selected_file,
             "'. Try re-processing or selecting a different file.")
    ))
    
    selected
  })
  
  
  # ========================================================================
  # SUMMARY TABLE
  # ========================================================================
  
  summary_df <- reactive({
    dat <- processed_all()
    req(dat)
    
    # Default to "mid" bodypart for summary if available
    dat2 <- dat
    if ("mid" %in% unique(dat2$bodypart)) {
      dat2 <- filter(dat2, bodypart == "mid")
    }
    
    # Apply optional speed cap
    dat2 <- apply_speed_cap(
      dat2,
      max_speed = input$max_speed,
      enabled   = input$drop_over_max_speed
    )
    
    # Calculate percent time moving (from all frames, including still ones)
    pct_moving <- dat2 %>%
      group_by(group, file_name, individual) %>%
      summarise(
        percent_time_moving = mean(moving, na.rm = TRUE) * 100,
        .groups = "drop"
      )
    
    # Calculate speed statistics (only from moving, non-jump frames)
    speed_stats <- dat2 %>%
      filter(moving == TRUE, !is_big_jump) %>%
      group_by(group, file_name, individual) %>%
      summarise(
        mean_moving_speed    = mean(speed, na.rm = TRUE),
        median_moving_speed  = median(speed, na.rm = TRUE),
        total_distance_moved = sum(speed, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Combine speed stats with percent moving
    result <- left_join(
      speed_stats, pct_moving,
      by = c("group", "file_name", "individual")
    )
    
    # Apply optional unit conversions (pixels -> mm, per-frame -> per-second)
    scale_factor <- if (isTRUE(input$convertunits)) input$scalefactor else 1
    fps_factor   <- if (isTRUE(input$convertspeed)) input$framerate   else 1
    
    result %>%
      mutate(
        mean_moving_speed    = mean_moving_speed    * scale_factor * fps_factor,
        median_moving_speed  = median_moving_speed  * scale_factor * fps_factor,
        # Distance is not a rate, so no fps conversion
        total_distance_moved = total_distance_moved * scale_factor
      )
  })
  
  output$summary_table <- renderDT({
    df <- summary_df()
    validate(need(nrow(df) > 0, "No summary data available. Process files first."))
    datatable(df, options = list(pageLength = 10))
  })
  
  # Download handler for the summary CSV
  output$download_summary <- downloadHandler(
    filename = function() paste0("summary_mean_moving_speed_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(summary_df(), file)
  )
  
  
  # ========================================================================
  # DUNNETT'S TEST
  # ========================================================================
  
  # Dynamic UI: dropdown to select the reference group
  output$ref_group_ui <- renderUI({
    df <- summary_df()
    req(df)
    grps <- sort(unique(df$group))
    selectInput("ref_group", "Reference group",
                choices = grps, selected = grps[[1]])
  })
  
  dunnett_res <- eventReactive(input$run_dunnett, {
    df <- summary_df()
    req(df, input$ref_group)
    
    # need at least 2 groups for a Dunnett test
    validate(need(
      n_distinct(df$group) >= 2,
      "Dunnett's test requires at least 2 groups."
    ))
    
    # Set the reference group as the first factor level
    df <- df %>%
      mutate(group = factor(
        group,
        levels = c(input$ref_group, setdiff(sort(unique(group)), input$ref_group))
      ))
    
    # Choose between mean and median speed for the ANOVA
    # NOTE: the column names are mean_moving_speed / median_moving_speed
    speed_col <- if (isTRUE(input$use_median_dunnett)) {
      "median_moving_speed"
    } else {
      "mean_moving_speed"
    }
    
    # Run one-way ANOVA then Dunnett's post-hoc test
    formula <- as.formula(paste(speed_col, "~ group"))
    
    tryCatch(
      {
        aov_fit <- aov(formula, data = df)
        gl      <- multcomp::glht(aov_fit, linfct = multcomp::mcp(group = "Dunnett"))
        broom::tidy(gl)
      },
      error = function(e) {
        showNotification(
          paste("Dunnett's test failed:", e$message),
          type = "error", duration = 10
        )
        NULL
      }
    )
  }, ignoreInit = TRUE)
  
  output$dunnett_table <- renderDT({
    dat <- dunnett_res()
    validate(need(!is.null(dat), "No results. Run the Dunnett test first."))
    datatable(dat, options = list(pageLength = 10))
  })
  
  
  # ========================================================================
  # TRAJECTORY PLOT (interactive plotly)
  # ========================================================================
  
  output$traj_plotly <- renderPlotly({
    df <- processed_selected()
    
    # Build axis ranges, optionally flipping Y for image coordinates
    x_range <- c(input$x_min, input$x_max)
    y_range <- if (isTRUE(input$flip_y)) {
      c(input$y_max, input$y_min)
    } else {
      c(input$y_min, input$y_max)
    }
    
    ttl <- paste0("Trajectory (moving segments coloured): ", unique(df$title))
    
    # Build the ggplot
    g <- plot_trajectory_coloured_segments(
      df_long         = df,
      xlim            = x_range,
      ylim            = y_range,
      ttl             = ttl,
      moving_col      = input$moving_col,
      still_col       = input$still_col,
      line_width      = input$line_width,
      axis_title_size = input$axis_title_size,
      axis_text_size  = input$axis_text_size,
      strip_text_size = input$strip_text_size,
      legend_position = input$legend_pos,
      drop_big_jumps  = input$drop_big_jumps
    )
    
    # Convert to plotly with zoom and pan enabled
    p <- ggplotly(g, dynamicTicks = TRUE) %>%
      layout(dragmode = "zoom") %>%
      config(
        displaylogo = FALSE,
        scrollZoom  = TRUE,
        modeBarButtonsToAdd = c("zoom2d", "pan2d", "select2d",
                                "lasso2d", "resetScale2d")
      )
    
    # Ensure all facet subplots share the same axis ranges
    set_all_subplot_ranges(p, x_range = x_range, y_range = y_range)
  })
  
  
  # ========================================================================
  # TRAJECTORY PLOT: moving only (static, not currently in active use)
  # ========================================================================
  
  output$traj_moving_plot <- renderPlot({
    df  <- processed_selected()
    dfm <- df %>% filter(moving %in% TRUE, !(is_big_jump %in% TRUE))
    
    validate(need(
      nrow(dfm) > 0,
      "No moving frames after threshold/window settings."
    ))
    
    xlim <- c(input$x_min, input$x_max)
    ylim <- c(input$y_min, input$y_max)
    ttl  <- paste0("Trajectory (moving): ", unique(df$title))
    
    plot_trajectory(dfm, xlim = xlim, ylim = ylim, ttl = ttl)
  })
  
  
  # ========================================================================
  # SPEED TRACE PLOT
  # ========================================================================
  
  output$speed_plot <- renderPlot({
    df <- processed_selected()
    
    df <- apply_speed_cap(
      df,
      max_speed = input$max_speed,
      enabled   = input$drop_over_max_speed
    )
    
    validate(need(nrow(df) > 0, "No data to plot after speed cap applied."))
    
    ttl <- paste0("Speed trace: ", unique(df$title))
    plot_speed_trace(df, ttl = ttl)
  })
  
  
  # ========================================================================
  # ANGLE PLOT
  # ========================================================================
  
  output$angle_plot <- renderPlot({
    req(input$compute_angle)
    df <- processed_selected()
    
    validate(need(
      "angle" %in% names(df),
      "Angle column not found. Ensure exactly 3 bodyparts are selected and re-process."
    ))
    
    df_ang <- filter(df, !is.na(angle))
    
    validate(need(
      nrow(df_ang) > 0,
      "No angle values available (check bodyparts and likelihood settings)."
    ))
    
    ttl <- paste0("Angle over time: ", unique(df$title))
    plot_angle_trace(df_ang, ttl = ttl)
  })
  
  
  # ========================================================================
  # DOWNLOAD HANDLERS
  # ========================================================================
  
  # Helper to build a ggplot for the trajectory (reused in render + download)
  build_traj_plot <- function(df) {
    xlim <- c(input$x_min, input$x_max)
    ylim <- c(input$y_min, input$y_max)
    
    plot_trajectory_coloured_segments(
      df_long         = df,
      xlim            = xlim,
      ylim            = ylim,
      ttl             = paste0("Trajectory (moving segments coloured): ", unique(df$title)),
      moving_col      = input$moving_col,
      still_col       = input$still_col,
      line_width      = input$line_width,
      axis_title_size = input$axis_title_size,
      axis_text_size  = input$axis_text_size,
      strip_text_size = input$strip_text_size,
      legend_position = input$legend_pos,
      drop_big_jumps  = input$drop_big_jumps
    )
  }
  
  output$download_traj <- downloadHandler(
    filename = function() {
      paste0("trajectory_", tools::file_path_sans_ext(input$selected_file), ".png")
    },
    content = function(file) {
      p <- build_traj_plot(processed_selected())
      ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_traj_moving <- downloadHandler(
    filename = function() {
      paste0("trajectory_moving_", tools::file_path_sans_ext(input$selected_file), ".png")
    },
    content = function(file) {
      df <- processed_selected() %>%
        filter(moving %in% TRUE, !(is_big_jump %in% TRUE))
      
      xlim <- c(input$x_min, input$x_max)
      ylim <- c(input$y_min, input$y_max)
      ttl  <- paste0("Trajectory (moving): ", unique(df$title))
      
      p <- plot_trajectory(df, xlim = xlim, ylim = ylim, ttl = ttl)
      ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_speed <- downloadHandler(
    filename = function() {
      paste0("speed_", tools::file_path_sans_ext(input$selected_file), ".png")
    },
    content = function(file) {
      df <- processed_selected()
      req(nrow(df) > 0)
      
      df <- apply_speed_cap(
        df,
        max_speed = input$max_speed,
        enabled   = input$drop_over_max_speed
      )
      
      ttl <- paste0("Speed trace: ", unique(df$title))
      p   <- plot_speed_trace(df, ttl = ttl)
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )
  
  output$download_angle <- downloadHandler(
    filename = function() {
      paste0("angle_", tools::file_path_sans_ext(input$selected_file), ".png")
    },
    content = function(file) {
      req(input$compute_angle)
      df  <- processed_selected() %>% filter(!is.na(angle))
      ttl <- paste0("Angle over time: ", unique(df$title))
      p   <- plot_angle_trace(df, ttl = ttl)
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )
}


# ============================================================================
# RUN
# ============================================================================

shinyApp(ui, server)
