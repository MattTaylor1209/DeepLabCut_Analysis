# app.R
# Shiny app for DeepLabCut filtered.csv batch processing + plots + stats

library(shiny)
library(shinyFiles)
library(tidyverse)
library(stringr)
library(slider)
library(DT)
library(multcomp)
library(broom)
library(colourpicker)
library(plotly)

# ----------------------------
# Helper functions (pipeline)
# ----------------------------

read_dlc_filtered_csv <- function(path) {
  header_lines <- readLines(path, n = 4, warn = FALSE)
  if (length(header_lines) < 4) return(NULL)
  
  scorer      <- str_split(header_lines[1], ",")[[1]]
  individuals <- str_split(header_lines[2], ",")[[1]]
  bodyparts   <- str_split(header_lines[3], ",")[[1]]
  coords      <- str_split(header_lines[4], ",")[[1]]
  
  # sanity: first column labels
  if (!tolower(individuals[1]) %in% c("individuals", "individual")) return(NULL)
  if (!tolower(bodyparts[1])   %in% c("bodyparts", "bodypart"))     return(NULL)
  if (!tolower(coords[1])      %in% c("coords", "coord"))           return(NULL)
  
  # build column names
  multi_names <- paste(individuals[-1], bodyparts[-1], coords[-1], sep = "_")
  all_colnames <- c("frame", multi_names)
  
  df_raw <- readr::read_csv(
    path,
    skip = 4,               # <-- IMPORTANT (skip coords too)
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  # If readr created fewer/more cols than expected, pad/truncate safely
  if (ncol(df_raw) < length(all_colnames)) {
    # pad missing columns with NA
    for (i in (ncol(df_raw) + 1):length(all_colnames)) df_raw[[i]] <- NA
  } else if (ncol(df_raw) > length(all_colnames)) {
    df_raw <- df_raw[, seq_len(length(all_colnames))]
  }
  
  colnames(df_raw) <- all_colnames
  
  df_raw %>%
    mutate(frame = suppressWarnings(as.integer(frame))) %>%
    mutate(across(-frame, ~ suppressWarnings(as.numeric(.x))))
}


get_dlc_bodyparts <- function(path) {
  hdr <- readLines(path, n = 3, warn = FALSE)
  if (length(hdr) < 3) return(character(0))
  
  bp <- stringr::str_split(hdr[3], ",")[[1]]
  bp <- bp[-1] # drop the "bodyparts" label
  bp <- bp[bp != ""]
  unique(bp)
}


tidy_dlc <- function(df_raw, bodyparts_keep = NULL, likelihood_min = NULL) {
  df_long <- df_raw %>%
    pivot_longer(-frame, names_to = "key", values_to = "value") %>%
    tidyr::extract(
      key,
      into = c("individual", "bodypart", "coord"),
      regex = "^(.*)_(.*)_(x|y|likelihood)$",
      remove = TRUE
    ) %>%
    pivot_wider(names_from = coord, values_from = value) %>%
    drop_na(x, y)
  
  if (!is.null(bodyparts_keep)) df_long <- df_long %>% filter(bodypart %in% bodyparts_keep)
  if (!is.null(likelihood_min) && "likelihood" %in% names(df_long)) {
    df_long <- df_long %>% filter(is.na(likelihood) | likelihood >= likelihood_min)
  }
  df_long
}

add_speed <- function(df_long) {
  df_long %>%
    group_by(individual, bodypart) %>%
    arrange(frame, .by_group = TRUE) %>%
    mutate(
      dx = x - lag(x),
      dy = y - lag(y),
      dframe = frame - lag(frame),
      dist = sqrt(dx^2 + dy^2),
      speed = dist / dframe
    ) %>%
    ungroup()
}

filter_big_jumps <- function(df_long,
                             max_jump_px = 50,
                             use_robust_jump = FALSE,
                             robust_mult = 10) {
  df_long %>%
    group_by(individual, bodypart) %>%
    mutate(
      med_dist = median(dist, na.rm = TRUE),
      thr = if (use_robust_jump) pmax(max_jump_px, robust_mult * med_dist) else max_jump_px,
      is_big_jump = dist > thr
    ) %>%
    ungroup()
}

flag_short_excursions <- function(df_long, max_excursion_frames = 10) {
  # Needs: columns x, y, frame, individual, bodypart, and 'thr' from filter_big_jumps()
  df_long %>%
    group_by(individual, bodypart) %>%
    arrange(frame, .by_group = TRUE) %>%
    mutate(
      xend = lead(x),
      yend = lead(y),
      seg_dist = sqrt((xend - x)^2 + (yend - y)^2),
      
      # big jump on the segment i -> i+1
      big_jump_next = seg_dist > thr,
      
      # region increments *after* a big jump segment
      region = cumsum(dplyr::lag(big_jump_next, default = FALSE)),
      
      max_region = max(region, na.rm = TRUE),
      
      # candidate excursion regions are those after the first jump AND before the last jump
      is_excursion_region = region >= 1 & region < max_region,
      
      region_n = ave(region, region, FUN = length),
      
      # mark short excursions (whole region)
      bad_excursion = is_excursion_region & (region_n <= max_excursion_frames)
    ) %>%
    ungroup()
}


add_movement_flag <- function(df_long, threshold_px = 2, window_n = 5) {
  df_long %>%
    group_by(individual, bodypart) %>%
    arrange(frame, .by_group = TRUE) %>%
    mutate(
      low_speed = speed < threshold_px,
      still_count = slider::slide_int(
        .x = low_speed,
        .f = sum,
        .before = window_n - 1,
        .complete = TRUE
      ),
      moving = dplyr::case_when(
        is.na(still_count) ~ NA,
        still_count == window_n ~ FALSE,
        TRUE ~ TRUE
      )
    ) %>%
    ungroup()
}

extract_group <- function(file_name) {
  grp <- stringr::str_extract(file_name, "^.*(?=_[0-9]+DLC)")
  if (is.na(grp)) {
    grp <- stringr::str_extract(file_name, "^.*(?=DLC)")
  }
  grp
}

extract_title <- function(file_name) {
  extract_group(file_name)
}

process_one_file <- function(path,
                             bodyparts_keep = c("mid"),
                             likelihood_min = NULL,
                             threshold_px = 2,
                             window_n = 5,
                             max_jump_px = 50,
                             use_robust_jump = FALSE,
                             robust_mult = 10,
                             max_excursion_frames = 10) {
  
  fn <- basename(path)
  df_raw <- read_dlc_filtered_csv(path)
  if (is.null(df_raw)) return(NULL)
  
  df_raw %>%
    tidy_dlc(bodyparts_keep = bodyparts_keep, likelihood_min = likelihood_min) %>%
    add_speed() %>%
    filter_big_jumps(
      max_jump_px = max_jump_px,
      use_robust_jump = use_robust_jump,
      robust_mult = robust_mult
    ) %>%
    flag_short_excursions(max_excursion_frames = max_excursion_frames) %>% 
    add_movement_flag(threshold_px = threshold_px, window_n = window_n) %>%
    mutate(
      file_path = path,
      file_name = fn,
      title = extract_title(fn),
      group = extract_group(fn)
    )
}

plot_trajectory <- function(df_long, xlim = c(0, 1280), ylim = c(0, 960), ttl = "") {
  ggplot(df_long, aes(x = x, y = y, colour = bodypart)) +
    geom_path(linewidth = 0.4, alpha = 0.85) +
    facet_wrap(~ individual, ncol = 2) +
    coord_fixed(xlim = xlim, ylim = ylim) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    ) +
    labs(title = ttl, x = "x (px)", y = "y (px)", colour = "Bodypart")
}

plot_speed_trace <- function(df_long, ttl = "") {
  ggplot(df_long, aes(x = frame, y = speed)) +
    geom_line(alpha = 0.6) +
    facet_grid(individual ~ bodypart, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = ttl, x = "Frame", y = "Speed (px/frame)")
}

make_segments <- function(df_long, drop_big_jumps = TRUE) {
  df_long %>%
    arrange(individual, bodypart, frame) %>%
    group_by(individual, bodypart) %>%
    mutate(
      xend = lead(x),
      yend = lead(y),
      moving_seg = moving,
      big_jump_seg = is_big_jump %in% TRUE
    ) %>%
    ungroup() %>%
    filter(!is.na(xend), !is.na(yend)) %>%
    { if (drop_big_jumps) dplyr::filter(., !big_jump_seg) else . }
}

plot_trajectory_coloured_segments <- function(df_long,
                                              xlim = c(0, 1280),
                                              ylim = c(0, 960),
                                              ttl = "",
                                              moving_col = "#11F011",
                                              still_col  = "#FA05EE",
                                              line_width = 0.6,
                                              axis_title_size = 13,
                                              axis_text_size  = 11,
                                              strip_text_size = 13,
                                              legend_position = "top",
                                              drop_big_jumps = TRUE) {
  
  seg <- make_segments(df_long, drop_big_jumps = drop_big_jumps)
  
  seg <- seg %>%
    mutate(
      moving_seg = dplyr::case_when(
        is.na(moving_seg) ~ "Unknown",
        moving_seg ~ "Moving",
        TRUE ~ "Still"
      ),
      moving_seg = factor(moving_seg, levels = c("Moving", "Still", "Unknown"))
    )
  
  ggplot(seg, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment(aes(colour = moving_seg), linewidth = line_width, alpha = 0.9) +
    facet_wrap(~ individual, ncol = 2) +
    coord_fixed(xlim = xlim, ylim = ylim) +
    scale_colour_manual(
      values = c(Moving = moving_col, Still = still_col, Unknown = "black"),
      drop = FALSE
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = legend_position,
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = strip_text_size),
      axis.title = element_text(size = axis_title_size),
      axis.text  = element_text(size = axis_text_size),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = ttl, x = "x (px)", y = "y (px)", colour = NULL)
}

# ----------------------------
# Shiny app
# ----------------------------

ui <- fluidPage(
  titlePanel("DeepLabCut CSV Explorer (filtered.csv)"),
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("1) Choose folder"),
      shinyDirButton("dir", "Select folder", "Choose folder containing *filtered.csv"),
      verbatimTextOutput("folder_txt"),
      tags$hr(),
      
      tags$h4("2) Processing settings"),
      textInput("pattern", "File pattern", value = "filtered.csv"),
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
      numericInput("window_n", "Stillness window (frames)", value = 5, min = 1, step = 1),
      
      tags$h4("Jump filtering"),
      checkboxInput("drop_big_jumps", "Drop big jumps in plots", value = TRUE),
      numericInput("max_jump_px", "Max allowed step distance (px)", value = 50, min = 0, step = 5),
      checkboxInput("use_robust_jump", "Use robust per-track threshold (median * multiplier)", value = FALSE),
      numericInput("robust_mult", "Robust multiplier", value = 10, min = 1, step = 1),
      checkboxInput("drop_bad_excursions", "Drop short excursions between big jumps", value = TRUE),
      numericInput("max_excursion_frames", "Max excursion length (frames)", value = 10, min = 1, step = 1),
    
      tags$hr(),
      
      tags$h4("3) Plot settings"),
      numericInput("x_min", "x min", value = 0),
      numericInput("x_max", "x max", value = 1280),
      numericInput("y_min", "y min", value = 0),
      numericInput("y_max", "y max", value = 960),
      
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
      
      tags$hr(),
      actionButton("run", "Load & process", class = "btn-primary"),
      tags$br(), tags$br(),
      uiOutput("file_picker_ui"),
      tags$hr(),
      downloadButton("download_summary", "Download summary CSV")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Files",
                 tags$h4("Detected files"),
                 DTOutput("files_table")
        ),
        tabPanel("Trajectory",
                 plotlyOutput("traj_plotly", height = 700),
                 downloadButton("download_traj", "Download trajectory plot (PNG)")
        ),
        tabPanel("Trajectory (moving only)",
                 plotOutput("traj_moving_plot", height = 700),
                 downloadButton("download_traj_moving", "Download moving-only plot (PNG)")
        ),
        tabPanel("Speed trace",
                 plotOutput("speed_plot", height = 700),
                 downloadButton("download_speed", "Download speed plot (PNG)")
        ),
        tabPanel("Summary",
                 tags$h4("Mean moving speed per individual (mid by default)"),
                 DTOutput("summary_table")
        ),
        tabPanel("Dunnett vs reference",
                 uiOutput("ref_group_ui"),
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
    d <- chosen_dir()
    req(d)
    pat <- input$pattern
    paths <- list.files(d, pattern = pat, full.names = TRUE)
    tibble(
      file_name = basename(paths),
      file_path = paths
    ) %>% arrange(file_name)
  })
  
  bodyparts_available <- reactive({
    df <- files_df()
    req(nrow(df) > 0)
    
    # union of bodyparts across all files (only reads header lines, so it's fast)
    bps <- unique(unlist(lapply(df$file_path, get_dlc_bodyparts)))
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
  
  observeEvent(input$run, {
    df <- files_df()
    req(nrow(df) > 0)
    
    bodyparts_keep <- input$bodyparts_keep
    if (is.null(bodyparts_keep) || length(bodyparts_keep) == 0) bodyparts_keep <- NULL
    
    likelihood_min <- if (is.na(input$likelihood_min)) NULL else input$likelihood_min
    
    shiny::withProgress(message = "Processing files...", value = 0, {
      out <- purrr::map_dfr(df$file_path, ~{
        shiny::incProgress(1 / nrow(df))
        
        process_one_file(
          .x,
          bodyparts_keep = if (length(bodyparts_keep)) bodyparts_keep else NULL,
          likelihood_min = likelihood_min,
          threshold_px   = input$threshold_px,
          window_n       = input$window_n,
          
          # ✅ PASS YOUR JUMP FILTER UI VALUES
          max_jump_px     = input$max_jump_px,
          use_robust_jump = input$use_robust_jump,
          robust_mult     = input$robust_mult,
          max_excursion_frames = input$max_excursion_frames
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
    
    dat2 %>%
      filter(moving %in% TRUE, !(is_big_jump %in% TRUE)) %>%
      group_by(group, file_name, individual) %>%
      summarise(mean_speed = mean(speed, na.rm = TRUE), .groups = "drop")
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
    
    aov_fit <- aov(mean_speed ~ group, data = df)
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
    
    xlim <- c(input$x_min, input$x_max)
    ylim <- c(input$y_min, input$y_max)
    ttl  <- paste0("Trajectory (moving segments coloured): ", unique(df$title))
    
    p <- plot_trajectory_coloured_segments(
      df_long = df,
      xlim = xlim, ylim = ylim,
      ttl = ttl,
      moving_col = input$moving_col,
      still_col  = input$still_col,
      line_width = input$line_width,
      axis_title_size = input$axis_title_size,
      axis_text_size  = input$axis_text_size,
      strip_text_size = input$strip_text_size,
      legend_position = if (input$legend_pos == "none") "none" else input$legend_pos,
      drop_big_jumps = input$drop_big_jumps  # ✅ checkbox now controls this
    )
    
    ggplotly(p, dynamicTicks = TRUE) %>%
      layout(dragmode = "zoom") %>%
      config(
        displaylogo = FALSE,
        scrollZoom = TRUE,
        modeBarButtonsToAdd = c("zoom2d","pan2d","select2d","lasso2d","resetScale2d")
      )
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
      p <- plot_speed_trace(df, ttl = paste0("Speed trace: ", unique(df$title)))
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )
}

shinyApp(ui, server)
