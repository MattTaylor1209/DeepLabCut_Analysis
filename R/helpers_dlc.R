# R/helpers_dlc.R
# Helper functions for DLC filtered.csv processing + plotting
# (expects the required packages are loaded by app.R)

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


# Calculate angle between 3 body parts

angle_deg <- function(A, B, C) {
  # A, B, C are numeric vectors length 2: c(x, y)
  BA <- A - B
  BC <- C - B
  
  # protect against zero-length vectors
  if (any(!is.finite(c(BA, BC))) || (sum(BA^2) == 0) || (sum(BC^2) == 0)) {
    return(NA_real_)
  }
  
  # angle via atan2(|cross|, dot)  -> stable, returns [0, pi]
  cross <- BA[1] * BC[2] - BA[2] * BC[1]
  dot   <- BA[1] * BC[1] + BA[2] * BC[2]
  
  ang <- atan2(cross, dot)                # [-pi, pi]
  ang_deg <- ang * 180 / pi               # [-180, 180]
  
  if (ang_deg < 0) ang_deg <- ang_deg + 360  # [0, 360)
  ang_deg
}


compute_angles <- function(df_long, parts3, vertex = parts3[2]) {
  stopifnot(length(parts3) == 3, vertex %in% parts3)
  
  others <- setdiff(parts3, vertex)
  A_name <- others[1]
  C_name <- others[2]
  B_name <- vertex
  
  df_wide <- df_long %>%
    filter(bodypart %in% parts3) %>%
    select(frame, individual, bodypart, x, y) %>%
    tidyr::pivot_wider(
      names_from = bodypart,
      values_from = c(x, y),
      names_glue = "{bodypart}_{.value}"
    )
  
  needed <- c(
    paste0(A_name, c("_x","_y")),
    paste0(B_name, c("_x","_y")),
    paste0(C_name, c("_x","_y"))
  )
  if (!all(needed %in% names(df_wide))) {
    return(df_wide %>% transmute(frame, individual, angle = NA_real_))
  }
  
  ax <- paste0(A_name, "_x"); ay <- paste0(A_name, "_y")
  bx <- paste0(B_name, "_x"); by <- paste0(B_name, "_y")
  cx <- paste0(C_name, "_x"); cy <- paste0(C_name, "_y")
  
  df_wide %>%
    mutate(
      angle = purrr::pmap_dbl(
        list(.data[[ax]], .data[[ay]],
             .data[[bx]], .data[[by]],
             .data[[cx]], .data[[cy]]),
        ~ angle_deg(c(..1, ..2), c(..3, ..4), c(..5, ..6))
      )
    ) %>%
    select(frame, individual, angle)
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

apply_speed_cap <- function(df_long, max_speed = NULL, enabled = TRUE) {
  if (!enabled) return(df_long)
  if (is.null(max_speed) || is.na(max_speed)) return(df_long)
  df_long %>% dplyr::filter(is.na(speed) | speed <= max_speed)
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
                             max_excursion_frames = 10,
                             compute_angle = FALSE,
                             angle_vertex = NULL) {
  
  fn <- basename(path)
  df_raw <- read_dlc_filtered_csv(path)
  if (is.null(df_raw)) return(NULL)
  
  # 1) df_long created here (THIS is what angles should use)
  df_long <- tidy_dlc(df_raw, bodyparts_keep = bodyparts_keep, likelihood_min = likelihood_min)
  
  # 2) Optionally compute angle from df_long (NOT df_raw)
  if (isTRUE(compute_angle) && !is.null(bodyparts_keep) && length(bodyparts_keep) == 3) {
    
    vertex <- if (!is.null(angle_vertex) && angle_vertex %in% bodyparts_keep) {
      angle_vertex
    } else {
      bodyparts_keep[2]  # default: middle selection is vertex
    }
    
    angle_tbl <- compute_angles(df_long, parts3 = bodyparts_keep, vertex = vertex)
    
    # join angle back onto df_long so it flows through the rest of your pipeline
    df_long <- df_long %>%
      left_join(angle_tbl, by = c("frame", "individual"))
    
  } else {
    df_long <- df_long %>% mutate(angle = NA_real_)
  }
  
  # 3) Continue pipeline as normal
  df_long %>%
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


plot_angle_trace <- function(df_long, ttl = "") {
  ggplot(df_long, aes(x = frame, y = angle)) +
    geom_line(alpha = 0.6) +
    geom_hline(yintercept = 180, colour = "red3", linetype = "dashed")+
    facet_wrap(~ individual, ncol = 2, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = ttl, x = "Frame", y = "Angle (degrees)")
}
