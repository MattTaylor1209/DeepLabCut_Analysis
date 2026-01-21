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

# ----------------------------
# Helper functions (pipeline)
# ----------------------------

read_dlc_filtered_csv <- function(path) {
  header_lines <- readLines(path, n = 3)
  scorer      <- str_split(header_lines[1], ",")[[1]]
  individuals <- str_split(header_lines[2], ",")[[1]]
  bodyparts   <- str_split(header_lines[3], ",")[[1]]
  
  coords <- rep(c("x", "y", "likelihood"), length.out = length(individuals) - 1)
  multi_names <- paste(individuals[-1], bodyparts[-1], coords, sep = "_")
  all_colnames <- c("frame", multi_names)
  
  df_raw <- readr::read_csv(
    path,
    skip = 3,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  if (ncol(df_raw) != length(all_colnames)) return(NULL)
  
  colnames(df_raw) <- all_colnames
  df_raw %>%
    mutate(frame = as.numeric(frame)) %>%
    mutate(across(-frame, ~ suppressWarnings(as.numeric(.x))))
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
  # capture everything before _<digits>DLC
  grp <- stringr::str_extract(file_name, "^.*(?=_[0-9]+DLC)")
  if (is.na(grp)) {
    # fallback if pattern not found: use everything before "DLC" (old behaviour)
    grp <- stringr::str_extract(file_name, "^.*(?=DLC)")
  }
  grp
}

extract_title <- function(file_name) {
  # If you want "title" to match group, just reuse it:
  extract_group(file_name)
}


process_one_file <- function(path,
                             bodyparts_keep = c("mid"),
                             likelihood_min = NULL,
                             threshold_px = 2,
                             window_n = 5) {
  fn <- basename(path)
  df_raw <- read_dlc_filtered_csv(path)
  if (is.null(df_raw)) return(NULL)
  
  df_raw %>%
    tidy_dlc(bodyparts_keep = bodyparts_keep, likelihood_min = likelihood_min) %>%
    add_speed() %>%
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
      textInput("bodyparts", "Bodyparts to keep (comma-separated)", value = "mid"),
      numericInput("likelihood_min", "Min likelihood (optional)", value = NA, min = 0, max = 1, step = 0.01),
      numericInput("threshold_px", "Movement threshold (px/frame)", value = 2, min = 0, step = 0.1),
      numericInput("window_n", "Stillness window (frames)", value = 5, min = 1, step = 1),
      tags$hr(),
      
      tags$h4("3) Plot settings"),
      numericInput("x_min", "x min", value = 0),
      numericInput("x_max", "x max", value = 1280),
      numericInput("y_min", "y min", value = 0),
      numericInput("y_max", "y max", value = 960),
      
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
                 plotOutput("traj_plot", height = 700),
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
                 tags$p("Tip: The reference group is the FIRST level of the factor. You can relevel by editing 'Group order' below."),
                 textInput("group_order", "Group order (comma-separated; first = reference)", value = ""),
                 actionButton("apply_levels", "Apply group order"),
                 tags$br(), tags$br(),
                 DTOutput("dunnett_table")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Folder selection (local machine)
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
  
  # Scan files
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
  
  output$files_table <- renderDT({
    req(files_df())
    datatable(files_df(), options = list(pageLength = 10))
  })
  
  # UI for choosing which file to plot
  output$file_picker_ui <- renderUI({
    df <- files_df()
    req(nrow(df) > 0)
    selectInput("selected_file", "Select file to plot", choices = df$file_name, selected = df$file_name[[1]])
  })
  
  # Process all files when user clicks "Load & process"
  processed_all <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    df <- files_df()
    req(nrow(df) > 0)
    
    bodyparts_keep <- str_split(input$bodyparts, ",")[[1]] %>%
      str_trim() %>% discard(~ .x == "")
    
    likelihood_min <- if (is.na(input$likelihood_min)) NULL else input$likelihood_min
    
    shiny::withProgress(message = "Processing files...", value = 0, {
      out <- purrr::map_dfr(df$file_path, ~{
        shiny::incProgress(1 / nrow(df))
        process_one_file(
          .x,
          bodyparts_keep = if (length(bodyparts_keep)) bodyparts_keep else NULL,
          likelihood_min = likelihood_min,
          threshold_px = input$threshold_px,
          window_n = input$window_n
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
  
  
  # Subset for chosen file
  processed_selected <- reactive({
    dat <- processed_all()
    req(dat)
    req(input$selected_file)
    dat %>% filter(file_name == input$selected_file)
  })
  
  # Summary table (mean moving speed per individual)
  summary_df <- reactive({
    dat <- processed_all()
    req(dat)
    # if "mid" exists, default to it; otherwise keep all bodyparts and let user filter later
    dat2 <- dat
    if ("mid" %in% unique(dat2$bodypart)) dat2 <- dat2 %>% filter(bodypart == "mid")
    
    dat2 %>%
      filter(moving %in% TRUE) %>%
      group_by(group, file_name, individual) %>%
      summarise(mean_speed = mean(speed, na.rm = TRUE), .groups = "drop")
  })
  
  output$summary_table <- renderDT({
    req(summary_df())
    datatable(summary_df(), options = list(pageLength = 10))
  })
  
  # Allow user-defined group order for Dunnett
  group_levels <- reactiveVal(NULL)
  
  observeEvent(input$apply_levels, {
    s <- input$group_order
    if (!nzchar(s)) {
      group_levels(NULL)
    } else {
      lev <- str_split(s, ",")[[1]] %>% str_trim() %>% discard(~ .x == "")
      group_levels(lev)
    }
  })
  
  dunnett_res <- reactive({
    df <- summary_df()
    req(df)
    
    # apply custom level order if provided
    lev <- group_levels()
    if (!is.null(lev) && length(lev) > 1) {
      df <- df %>% mutate(group = factor(group, levels = lev))
    } else {
      df <- df %>% mutate(group = factor(group))
    }
    
    # Need >=2 groups for Dunnett
    if (nlevels(df$group) < 2) return(tibble(note = "Need at least 2 groups."))
    
    aov_fit <- aov(mean_speed ~ group, data = df)
    gl <- multcomp::glht(aov_fit, linfct = multcomp::mcp(group = "Dunnett"))
    broom::tidy(gl)
  })
  
  output$dunnett_table <- renderDT({
    dat <- dunnett_res()
    req(dat)
    datatable(dat, options = list(pageLength = 10))
  })
  
  # ---- Plots ----
  
  output$traj_plot <- renderPlot({
    df <- processed_selected()
    req(nrow(df) > 0)
    xlim <- c(input$x_min, input$x_max)
    ylim <- c(input$y_min, input$y_max)
    ttl <- paste0("Trajectory: ", unique(df$title))
    plot_trajectory(df, xlim = xlim, ylim = ylim, ttl = ttl)
  })
  
  output$traj_moving_plot <- renderPlot({
    df <- processed_selected()
    req(nrow(df) > 0)
    dfm <- df %>% filter(moving %in% TRUE)
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
  
  # ---- Downloads ----
  
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
      p <- plot_trajectory(df, xlim = xlim, ylim = ylim, ttl = paste0("Trajectory: ", unique(df$title)))
      ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_traj_moving <- downloadHandler(
    filename = function() paste0("trajectory_moving_", tools::file_path_sans_ext(input$selected_file), ".png"),
    content = function(file) {
      df <- processed_selected() %>% filter(moving %in% TRUE)
      xlim <- c(input$x_min, input$x_max)
      ylim <- c(input$y_min, input$y_max)
      p <- plot_trajectory(df, xlim = xlim, ylim = ylim, ttl = paste0("Trajectory (moving): ", unique(df$title)))
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
