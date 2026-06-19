# R/helpers_dlc.R
# ============================================================================
# Helper functions for DeepLabCut / SLEAP tracking CSV processing + plotting
#
# This file contains the full pipeline for reading, tidying, processing, and
# plotting animal tracking data exported from DeepLabCut (filtered.csv) or
# SLEAP (analysis.csv). Functions are organised into sections:
#
#   1.  DLC CSV format detection and reading
#   1b. SLEAP CSV reading (parallel reader/tidier for SLEAP's analysis.csv)
#   2.  Data tidying (wide -> long)
#   3.  Speed, jump, excursion, and movement calculations
#   4.  Angle computation and head-swing detection
#   5.  File metadata extraction
#   6.  Main processing pipeline (process_one_file)
#   7.  Plotting functions
#   8.  Plotly helpers
#
# DLC and SLEAP each get their own reader + tidier (Sections 1/1b/2), but both
# converge on the same long-format shape (frame, individual, bodypart, x, y,
# likelihood) so every later step is completely format-agnostic.
#
# Expects tidyverse, slider, and plotly to be loaded by app.R.
# ============================================================================


# ============================================================================
# 1. DLC CSV FORMAT DETECTION AND READING
# ============================================================================

#' Detect the tracking CSV format by inspecting the first few header rows.
#'
#' Recognises DLC's three export formats plus SLEAP's analysis.csv export:
#'   - "multi"  : 4-row header (scorer / individuals / bodyparts / coords)
#'   - "single" : 3-row header (scorer / bodyparts / coords, no individuals)
#'   - "flat"   : 1-row header (already flattened column names)
#'   - "sleap"  : SLEAP analysis.csv export (1-row header starting "track")
#'   - "unknown": file is empty or unreadable
#'
#' @param path File path to the CSV.
#' @return Character string: "multi", "single", "flat", "sleap", or "unknown".
detect_dlc_format <- function(path) {
  hdr <- readLines(path, n = 4, warn = FALSE)
  
  # Empty file guard
  if (length(hdr) < 1) return("unknown")
  
  first_cell <- tolower(strsplit(hdr[1], ",", fixed = TRUE)[[1]][1])
  
  # SLEAP's analysis.csv export always starts with a "track" column, which
  # is distinct enough from DLC's "scorer" / pre-flattened headers to detect
  # up front, before checking for the DLC multi-row header signature
  if (identical(first_cell, "track")) return("sleap")
  
  # Check whether the first cell is "scorer" (DLC multi-row header signature)
  if (!identical(first_cell, "scorer")) {
    # No DLC multi-index header -- treat as a pre-flattened single-row header
    return("flat")
  }
  
  # Determine whether there is an "individuals" row (multi) or just bodyparts (single)
  if (length(hdr) >= 2) {
    second_cell <- tolower(strsplit(hdr[2], ",", fixed = TRUE)[[1]][1])
    if (second_cell %in% c("individuals", "individual")) return("multi")
    if (second_cell %in% c("bodyparts", "bodypart"))     return("single")
  }
  
  # Safest fallback for DLC-style files that start with 'scorer' but don't
  # match expected patterns
  "multi"
}


#' Sanitize a DLC header token for safe use as part of a column name.
#'
#' Replaces whitespace with underscores and strips non-alphanumeric characters
#' (except '.', '-', '_') so that bodypart and individual names are safe to use
#' as R column name fragments.
#'
#' @param x Character vector of header tokens.
#' @return Character vector of sanitized tokens.
sanitize_dlc_token <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s+", "_")
  x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
  x
}


#' Normalize a raw data frame to match expected column names.
#'
#' Handles mismatches between the number of data columns and the number of
#' expected column names (padding with NA or trimming). Coerces the frame
#' column to integer and all other columns to numeric, and drops rows where
#' frame is NA (e.g. trailing empty rows from Excel re-saves).
#'
#' @param df_raw  Data frame read from CSV (no header names yet).
#' @param all_colnames Character vector of expected column names.
#' @return Tibble with correct column names and types.
normalize_df_cols <- function(df_raw, all_colnames) {
  
  # Pad missing columns with NA or trim extra columns
  if (ncol(df_raw) < length(all_colnames)) {
    for (i in (ncol(df_raw) + 1):length(all_colnames)) df_raw[[i]] <- NA
  } else if (ncol(df_raw) > length(all_colnames)) {
    df_raw <- df_raw[, seq_len(length(all_colnames))]
  }
  
  colnames(df_raw) <- all_colnames
  
  df_raw %>%
    # Coerce frame to integer -- non-numeric values (empty rows) become NA
    mutate(frame = suppressWarnings(as.integer(frame))) %>%
    # drop rows with NA frame (trailing empty rows from Excel, etc.)
    filter(!is.na(frame)) %>%
    # Coerce all measurement columns to numeric
    mutate(across(-frame, ~ suppressWarnings(as.numeric(.x))))
}


#' Read a DeepLabCut filtered CSV file into a wide data frame.
#'
#' Supports three DLC export formats (multi-animal, single-animal, flat).
#' Column names are normalised to the pattern: individual_bodypart_coord
#' (e.g. "animal1_Mid_x", "animal1_Mid_likelihood").
#'
#' @param path   File path to the CSV.
#' @param format One of "auto", "multi", "single", "flat".
#' @return Tibble with columns: frame, plus individual_bodypart_coord columns.
#'         Returns NULL if the file cannot be parsed.
read_dlc_filtered_csv <- function(path, format = c("auto", "multi", "single", "flat")) {
  format <- match.arg(format)
  
  # Auto-detect from file headers if not specified
  if (format == "auto") {
    format <- detect_dlc_format(path)
    if (format == "unknown") {
      warning("Could not detect DLC format for: ", basename(path))
      return(NULL)
    }
  }
  
  # --- Multi-animal format: 4-row header (scorer / individuals / bodyparts / coords) ---
  if (format == "multi") {
    header_lines <- readLines(path, n = 4, warn = FALSE)
    if (length(header_lines) < 4) return(NULL)
    
    # Parse each header row into a character vector
    individuals <- str_split(header_lines[2], ",")[[1]]
    bodyparts   <- str_split(header_lines[3], ",")[[1]]
    coords      <- str_split(header_lines[4], ",")[[1]]
    
    # Validate header labels
    if (!tolower(individuals[1]) %in% c("individuals", "individual")) return(NULL)
    if (!tolower(bodyparts[1])   %in% c("bodyparts", "bodypart"))     return(NULL)
    if (!tolower(coords[1])      %in% c("coords", "coord"))          return(NULL)
    
    # Build column names: individual_bodypart_coord (skip the label in position 1)
    n <- min(length(individuals), length(bodyparts), length(coords)) - 1
    if (n <= 0) return(NULL)
    
    indiv <- sanitize_dlc_token(individuals[2:(n + 1)])
    bp    <- sanitize_dlc_token(bodyparts[2:(n + 1)])
    cd    <- sanitize_dlc_token(coords[2:(n + 1)])
    
    all_colnames <- c("frame", paste(indiv, bp, cd, sep = "_"))
    
    # Read data rows (skip the 4 header lines)
    df_raw <- readr::read_csv(
      path, skip = 4, col_names = FALSE,
      show_col_types = FALSE, progress = FALSE
    )
    
    return(normalize_df_cols(df_raw, all_colnames))
  }
  
  # --- Single-animal format: 3-row header (scorer / bodyparts / coords) ---
  if (format == "single") {
    header_lines <- readLines(path, n = 3, warn = FALSE)
    if (length(header_lines) < 3) return(NULL)
    
    bodyparts <- str_split(header_lines[2], ",")[[1]]
    coords    <- str_split(header_lines[3], ",")[[1]]
    
    if (!tolower(bodyparts[1]) %in% c("bodyparts", "bodypart")) return(NULL)
    if (!tolower(coords[1])    %in% c("coords", "coord"))      return(NULL)
    
    n <- min(length(bodyparts), length(coords)) - 1
    if (n <= 0) return(NULL)
    
    # Single-animal: assign a default individual name
    indiv <- rep("animal1", n)
    bp    <- sanitize_dlc_token(bodyparts[2:(n + 1)])
    cd    <- sanitize_dlc_token(coords[2:(n + 1)])
    
    all_colnames <- c("frame", paste(indiv, bp, cd, sep = "_"))
    
    df_raw <- readr::read_csv(
      path, skip = 3, col_names = FALSE,
      show_col_types = FALSE, progress = FALSE
    )
    
    return(normalize_df_cols(df_raw, all_colnames))
  }
  
  # --- Flat format: single-row header (already flattened) ---
  if (format == "flat") {
    df_raw <- readr::read_csv(
      path, col_names = TRUE,
      show_col_types = FALSE, progress = FALSE
    )
    
    if (ncol(df_raw) < 4) return(NULL)
    
    # First column is always the frame index (may be named "X1" by readr)
    nm    <- names(df_raw)
    nm[1] <- "frame"
    
    # Rename columns to the individual_bodypart_coord pattern expected by tidy_dlc()
    # Uses map_chr instead of a for loop for cleaner functional style
    new_names <- map_chr(seq_along(nm), function(i) {
      if (i == 1) return("frame")
      
      col <- nm[i]
      
      # Already in individual_bodypart_coord format
      if (str_detect(col, "^(.*)_(.*)_(x|y|likelihood)$")) return(col)
      
      # bodypart_coord format (common single-animal flat export) -- prepend "animal1"
      m <- str_match(col, "^(.*)_(x|y|likelihood)$")
      if (!is.na(m[1, 1])) {
        return(paste("animal1", sanitize_dlc_token(m[1, 2]), m[1, 3], sep = "_"))
      }
      
      # Unrecognised pattern -- prefix with animal1 and sanitize
      paste0("animal1_", sanitize_dlc_token(col))
    })
    
    names(df_raw) <- new_names
    
    df_raw %>%
      mutate(frame = suppressWarnings(as.integer(frame))) %>%
      filter(!is.na(frame)) %>%
      mutate(across(-frame, ~ suppressWarnings(as.numeric(.x))))
  }
}


#' Extract unique bodypart names from a DLC CSV file header.
#'
#' Only reads the header lines (no data), so it's fast enough to call on every
#' file in a batch to build a union of available bodyparts.
#'
#' @param path   File path to the CSV.
#' @param format One of "auto", "multi", "single", "flat", "sleap".
#' @return Character vector of unique sanitized bodypart names.
get_dlc_bodyparts <- function(path, format = c("auto", "multi", "single", "flat", "sleap")) {
  format <- match.arg(format)
  if (format == "auto") format <- detect_dlc_format(path)
  
  # SLEAP has its own header layout (bodyparts encoded as "{bodypart}.x" etc.
  # column suffixes rather than dedicated header rows), so it gets its own helper
  if (format == "sleap") return(get_sleap_bodyparts(path))
  
  # Multi-animal: bodyparts are in header row 3
  if (format == "multi") {
    hdr <- readLines(path, n = 3, warn = FALSE)
    if (length(hdr) < 3) return(character(0))
    bp <- str_split(hdr[3], ",")[[1]][-1]  # drop the "bodyparts" label
    return(unique(sanitize_dlc_token(bp[bp != ""])))
  }
  
  # Single-animal: bodyparts are in header row 2
  if (format == "single") {
    hdr <- readLines(path, n = 2, warn = FALSE)
    if (length(hdr) < 2) return(character(0))
    bp <- str_split(hdr[2], ",")[[1]][-1]
    return(unique(sanitize_dlc_token(bp[bp != ""])))
  }
  
  # Flat: parse bodypart names from column headers using map_chr
  nm <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  if (length(nm) <= 1) return(character(0))
  
  # Drop the frame column, then extract the bodypart portion of each name
  bps <- map_chr(nm[-1], function(col) {
    if (str_detect(col, "^(.*)_(.*)_(x|y|likelihood)$")) {
      # individual_bodypart_coord -> extract bodypart (second capture group)
      return(str_replace(col, "^(.*)_(.*)_(x|y|likelihood)$", "\\2"))
    }
    # bodypart_coord -> strip the _coord suffix
    str_replace(col, "_(x|y|likelihood)$", "")
  })
  
  unique(sanitize_dlc_token(bps))
}


# ============================================================================
# 1b. SLEAP SUPPORT
# ============================================================================
#
# SLEAP's "analysis.csv" export (one of the formats produced by sleap-convert
# or File > Export Analysis CSV in the SLEAP GUI) has a single header row and
# is already "one row per individual (track) per frame":
#
#   track,frame_idx,instance.score,Head.score,Head.x,Head.y,Mid.score,...
#
# This is structurally different from DLC's wide export (one column per
# measurement, individual/bodypart encoded into the column name), so it gets
# its own reader + tidier rather than being shoehorned through
# read_dlc_filtered_csv()/tidy_dlc(). Both tidiers converge on the same long
# format (frame, individual, bodypart, x, y, likelihood), so every function
# downstream of tidying (add_speed(), filter_big_jumps(), plotting, etc.) is
# completely unaware of which format the data originally came from.

#' Read a SLEAP analysis.csv export into a wide-per-individual data frame.
#'
#' @param path File path to the SLEAP analysis CSV.
#' @return Tibble with columns: frame, individual, plus "{bodypart}.score" /
#'         ".x" / ".y" columns for each tracked bodypart. Returns NULL if the
#'         expected "track"/"frame_idx" columns aren't present.
read_sleap_csv <- function(path) {
  df_raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  
  if (!all(c("track", "frame_idx") %in% names(df_raw))) return(NULL)
  
  df_raw %>%
    rename(individual = track, frame = frame_idx) %>%
    mutate(frame = suppressWarnings(as.integer(frame))) %>%
    # Drop detections SLEAP couldn't assign a track identity to -- there's no
    # individual to attribute them to, analogous to DLC's trailing NA-frame rows
    filter(!is.na(individual), individual != "", !is.na(frame)) %>%
    # Real-world SLEAP exports can contain a rare duplicate row for the same
    # (track, frame) with every score column blank but identical x/y -- an
    # export artefact, not a second detection. instance.score is otherwise
    # not used downstream, so drop it here once it's served this purpose.
    filter(!is.na(instance.score)) %>%
    select(-instance.score)
}

#' Extract unique bodypart names from a SLEAP analysis.csv header.
#'
#' Only reads the header row (no data), mirroring get_dlc_bodyparts()'s
#' fast header-only scanning so it's cheap to call across a whole batch.
#'
#' @param path File path to the SLEAP CSV.
#' @return Character vector of unique sanitized bodypart names.
get_sleap_bodyparts <- function(path) {
  nm <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  
  # Bodypart columns look like "Head.score" / "Head.x" / "Head.y" -- capture
  # the part before the suffix, excluding the per-instance "instance.score"
  bp <- str_match(nm, "^(.*)\\.(score|x|y)$")[, 2]
  bp <- bp[!is.na(bp) & bp != "instance"]
  
  unique(sanitize_dlc_token(bp))
}

#' Convert a SLEAP wide-per-individual data frame to tidy long format.
#'
#' Deliberately mirrors tidy_dlc()'s output shape -- frame, individual,
#' bodypart, x, y, likelihood -- so the rest of the pipeline is completely
#' format-agnostic. SLEAP's "{bodypart}.score" (a per-node confidence score)
#' is renamed to "likelihood" to match DLC's terminology, since both serve
#' the same role (e.g. for the likelihood_min filter below).
#'
#' @param df_raw         Wide tibble from read_sleap_csv().
#' @param bodyparts_keep Optional character vector of bodyparts to retain.
#' @param likelihood_min Optional minimum likelihood (SLEAP node score) threshold.
#' @return Long-format tibble: frame, individual, bodypart, x, y, likelihood.
tidy_sleap <- function(df_raw, bodyparts_keep = NULL, likelihood_min = NULL) {
  
  df_long <- df_raw %>%
    # Pivot all "{bodypart}.{score,x,y}" columns into long format in one step,
    # capturing bodypart and coord directly from the column name
    pivot_longer(
      cols          = matches("\\.(score|x|y)$"),
      names_to      = c("bodypart", "coord"),
      names_pattern = "^(.*)\\.(score|x|y)$",
      values_to     = "value"
    ) %>%
    mutate(
      bodypart = sanitize_dlc_token(bodypart),
      coord    = if_else(coord == "score", "likelihood", coord)
    ) %>%
    # Spread coord back into separate x / y / likelihood columns
    pivot_wider(names_from = coord, values_from = value) %>%
    drop_na(x, y)
  
  # Same defensive check as tidy_dlc(): list-columns mean duplicate
  # frame/individual/bodypart combinations slipped through pivot_wider
  if (is.list(df_long$x) || is.list(df_long$y)) {
    stop(
      "pivot_wider produced list-columns for x/y coordinates. ",
      "This usually means the CSV contains duplicate frame/individual/bodypart ",
      "combinations. Check the raw CSV for duplicates.",
      call. = FALSE
    )
  }
  
  if (!is.null(bodyparts_keep)) {
    df_long <- filter(df_long, bodypart %in% bodyparts_keep)
  }
  
  if (!is.null(likelihood_min) && "likelihood" %in% names(df_long)) {
    df_long <- filter(df_long, is.na(likelihood) | likelihood >= likelihood_min)
  }
  
  df_long
}


# ============================================================================
# 2. DATA TIDYING (WIDE -> LONG)
# ============================================================================

#' Convert a wide DLC data frame to tidy long format.
#'
#' Pivots from one-column-per-measurement (animal1_Mid_x, animal1_Mid_y, ...)
#' to one-row-per-bodypart-per-frame with columns: frame, individual, bodypart,
#' x, y, likelihood.
#'
#' @param df_raw         Wide data frame from read_dlc_filtered_csv().
#' @param bodyparts_keep Optional character vector of bodyparts to retain.
#' @param likelihood_min Optional minimum likelihood threshold.
#' @return Long-format tibble.
tidy_dlc <- function(df_raw, bodyparts_keep = NULL, likelihood_min = NULL) {
  
  df_long <- df_raw %>%
    # Pivot all measurement columns into key-value pairs
    pivot_longer(-frame, names_to = "key", values_to = "value") %>%
    # Split the composite key into individual, bodypart, and coordinate
    tidyr::extract(
      key,
      into = c("individual", "bodypart", "coord"),
      regex = "^(.*)_(.*)_(x|y|likelihood)$",
      remove = TRUE
    ) %>%
    # Spread coordinates back into separate columns (x, y, likelihood)
    pivot_wider(names_from = coord, values_from = value) %>%
    # Drop rows where either coordinate is missing
    drop_na(x, y)
  
  # check that pivot_wider didn't produce list-columns
  # (happens when duplicate frame/individual/bodypart combinations exist,
  if (is.list(df_long$x) || is.list(df_long$y)) {
    stop(
      "pivot_wider produced list-columns for x/y coordinates. ",
      "This usually means the CSV contains duplicate frame/individual/bodypart ",
      "combinations (e.g. trailing empty rows from Excel). ",
      "Check the raw CSV for duplicates.",
      call. = FALSE
    )
  }
  
  # Filter to requested bodyparts if specified
  if (!is.null(bodyparts_keep)) {
    df_long <- filter(df_long, bodypart %in% bodyparts_keep)
  }
  
  # Apply likelihood threshold if specified
  if (!is.null(likelihood_min) && "likelihood" %in% names(df_long)) {
    df_long <- filter(df_long, is.na(likelihood) | likelihood >= likelihood_min)
  }
  
  df_long
}


# ============================================================================
# 3. SPEED, JUMP, EXCURSION, AND MOVEMENT CALCULATIONS
# ============================================================================

#' Calculate frame-to-frame speed and add it to the long data.
#'
#' For each individual x bodypart track, computes the Euclidean distance moved
#' between consecutive frames and divides by the frame gap (handles dropped
#' frames correctly).
#'
#' @param df_long Long-format tibble with columns: frame, individual, bodypart, x, y.
#' @return Same tibble with added columns: dx, dy, dframe, dist, speed.
add_speed <- function(df_long) {
  
  # verify that x and y are numeric before attempting arithmetic
  stopifnot(
    "x column must be numeric in add_speed()" = is.numeric(df_long$x),
    "y column must be numeric in add_speed()" = is.numeric(df_long$y)
  )
  
  df_long %>%
    group_by(individual, bodypart) %>%
    arrange(frame, .by_group = TRUE) %>%
    mutate(
      # Frame-to-frame displacement in each axis
      dx     = x - lag(x),
      dy     = y - lag(y),
      # Number of frames elapsed (usually 1, but accounts for dropped frames)
      dframe = frame - lag(frame),
      # Euclidean distance moved
      dist   = sqrt(dx^2 + dy^2),
      # Speed = distance / frames elapsed
      speed  = dist / dframe
    ) %>%
    ungroup()
}


#' Flag frames with implausibly large jumps (tracking artefacts).
#'
#' Adds a threshold column (thr) and a boolean flag (is_big_jump) to mark
#' frames where the step distance exceeds either a fixed pixel threshold or
#' a robust per-track threshold (median distance x multiplier), whichever
#' is larger.
#'
#' @param df_long        Long-format tibble with a 'dist' column from add_speed().
#' @param max_jump_px    Fixed maximum allowed step distance in pixels.
#' @param use_robust_jump If TRUE, also consider a per-track adaptive threshold.
#' @param robust_mult    Multiplier for the median distance (used when use_robust_jump = TRUE).
#' @return Same tibble with added columns: med_dist, thr, is_big_jump.
filter_big_jumps <- function(df_long,
                             max_jump_px     = 50,
                             use_robust_jump = FALSE,
                             robust_mult     = 10) {
  df_long %>%
    group_by(individual, bodypart) %>%
    mutate(
      # Per-track median step distance (for robust thresholding)
      med_dist = median(dist, na.rm = TRUE),
      # Threshold: either fixed, or the larger of fixed and robust
      thr = if (use_robust_jump) {
        pmax(max_jump_px, robust_mult * med_dist)
      } else {
        max_jump_px
      },
      # Flag frames that exceed the threshold
      is_big_jump = dist > thr
    ) %>%
    ungroup()
}


#' Flag short excursions -- brief tracking jumps to a wrong location and back.
#'
#' Identifies regions between two big jumps that are shorter than
#' max_excursion_frames. These are likely artefacts where the tracker briefly
#' locked onto the wrong feature and then snapped back.
#'
#' @param df_long              Long-format tibble with 'thr' column from filter_big_jumps().
#' @param max_excursion_frames Maximum number of frames for a region to be
#'                             considered a short excursion.
#' @return Same tibble with added columns: xend, yend, seg_dist, big_jump_next,
#'         region, max_region, is_excursion_region, region_n, bad_excursion.
flag_short_excursions <- function(df_long, max_excursion_frames = 10) {
  
  df_long %>%
    group_by(individual, bodypart) %>%
    arrange(frame, .by_group = TRUE) %>%
    mutate(
      # Look-ahead coordinates for segment distance calculation
      xend = lead(x),
      yend = lead(y),
      seg_dist = sqrt((xend - x)^2 + (yend - y)^2),
      
      # Does the segment from this frame to the next exceed the jump threshold?
      big_jump_next = seg_dist > thr,
      
      # Assign each frame to a region (increments after each big jump)
      region     = cumsum(lag(big_jump_next, default = FALSE)),
      max_region = max(region, na.rm = TRUE),
      
      # Excursion candidates: regions between the first and last big jump
      is_excursion_region = region >= 1 & region < max_region,
      
      # Count frames in each region
      region_n = ave(region, region, FUN = length),
      
      # Mark as bad excursion if the region is short enough
      bad_excursion = is_excursion_region & (region_n <= max_excursion_frames)
    ) %>%
    ungroup()
}


#' Classify each frame as "moving" or "still" using a sliding window.
#'
#' A frame is classified as "still" if every frame in a trailing window of
#' window_n frames has speed below threshold_px. Also flags frames near the
#' start or end of each track (edge_flag).
#'
#' @param df_long      Long-format tibble with a 'speed' column.
#' @param threshold_px Speed below this (px/frame) counts as "not moving".
#' @param window_n     Number of consecutive slow frames required to be "still".
#' @param edge_n       Number of frames to flag at track start/end.
#' @return Same tibble with added columns: low_speed, still_count, moving, edge_flag.
add_movement_flag <- function(df_long, threshold_px = 2, window_n = 5, edge_n = 5) {
  
  df_long %>%
    group_by(individual, bodypart) %>%
    arrange(frame, .by_group = TRUE) %>%
    mutate(
      # Is speed below the movement threshold?
      low_speed = speed < threshold_px,
      
      # Count how many of the last window_n frames were slow
      still_count = slider::slide_int(
        .x        = low_speed,
        .f        = sum,
        .before   = window_n - 1,
        .complete = TRUE
      ),
      
      # Classify: still only if the entire window was slow
      moving = case_when(
        is.na(still_count)      ~ NA,
        still_count == window_n ~ FALSE,
        TRUE                    ~ TRUE
      ),
      
      # Flag frames near the start or end of the track
      .row_in_track = row_number(),
      .n_in_track   = n(),
      
      edge_flag = case_when(
        .n_in_track <= 2 * edge_n            ~ "START",
        .row_in_track <= edge_n              ~ "START",
        .row_in_track > .n_in_track - edge_n ~ "END",
        TRUE                                 ~ NA_character_
      ),
      edge_flag = factor(edge_flag, levels = c("START", "END"))
    ) %>%
    # Clean up temporary columns
    select(-.row_in_track, -.n_in_track) %>%
    ungroup()
}


#' Apply a maximum speed cap, removing frames above the threshold.
#'
#' @param df_long   Long-format tibble with a 'speed' column.
#' @param max_speed Maximum allowed speed (px/frame). NULL or NA to skip.
#' @param enabled   If FALSE, return data unchanged.
#' @return Filtered tibble.
apply_speed_cap <- function(df_long, max_speed = NULL, enabled = TRUE) {
  if (!enabled) return(df_long)
  if (is.null(max_speed) || is.na(max_speed)) return(df_long)
  
  df_long %>% filter(is.na(speed) | speed <= max_speed)
}


# ============================================================================
# 4. ANGLE COMPUTATION AND HEAD-SWING DETECTION
# ============================================================================

#' Compute the angle at vertex B formed by points A-B-C, in degrees.
#'
#' Uses atan2 for numerical stability. Returns a value in [0, 360) where
#' 180 degrees means the three points are collinear.
#'
#' @param A Numeric vector of length 2: c(x, y) for point A.
#' @param B Numeric vector of length 2: c(x, y) for the vertex.
#' @param C Numeric vector of length 2: c(x, y) for point C.
#' @return Angle in degrees, or NA if any input is invalid / zero-length.
angle_deg <- function(A, B, C) {
  
  BA <- A - B
  BC <- C - B
  
  # Guard against non-finite values or zero-length vectors
  if (any(!is.finite(c(BA, BC))) || (sum(BA^2) == 0) || (sum(BC^2) == 0)) {
    return(NA_real_)
  }
  
  # Cross product and dot product
  cross <- BA[1] * BC[2] - BA[2] * BC[1]
  dot   <- BA[1] * BC[1] + BA[2] * BC[2]
  
  # atan2 gives angle in [-pi, pi]; convert to [0, 360)
  ang_deg <- atan2(cross, dot) * 180 / pi
  if (ang_deg < 0) ang_deg <- ang_deg + 360
  
  ang_deg
}


#' Compute body angle for each frame given exactly 3 bodyparts.
#'
#' Pivots the long data to wide (one row per frame), then applies angle_deg()
#' row-wise. Frames where any of the 3 bodyparts are flagged as big jumps or
#' bad excursions get NA angles.
#'
#' @param df_long                Long-format tibble with x, y, is_big_jump, bad_excursion.
#' @param parts3                 Character vector of exactly 3 bodypart names.
#' @param vertex                 Which bodypart is the vertex (point B). Defaults to parts3[2].
#' @param drop_big_jump_frames   Set angle to NA if any bodypart has a big jump on that frame.
#' @param drop_bad_excursions    Set angle to NA if any bodypart is in a bad excursion.
#' @return Tibble with columns: frame, individual, angle.
compute_angles <- function(df_long,
                           parts3,
                           vertex                = parts3[2],
                           drop_big_jump_frames = TRUE,
                           drop_bad_excursions  = TRUE) {
  
  # validate inputs
  stopifnot(
    "compute_angles() requires exactly 3 bodyparts" = length(parts3) == 3,
    "vertex must be one of the 3 bodyparts"          = vertex %in% parts3
  )
  
  # Identify the three named points: A (arm1), B (vertex), C (arm2)
  others <- setdiff(parts3, vertex)
  A_name <- others[1]
  B_name <- vertex
  C_name <- others[2]
  
  # Pivot to wide: one row per frame x individual, with columns for each bodypart
  df_wide <- df_long %>%
    filter(bodypart %in% parts3) %>%
    select(frame, individual, bodypart, x, y, is_big_jump, bad_excursion) %>%
    pivot_wider(
      names_from  = bodypart,
      values_from = c(x, y, is_big_jump, bad_excursion),
      names_glue  = "{bodypart}_{.value}"
    )
  
  # Build column name references for each point
  ax <- paste0(A_name, "_x"); ay <- paste0(A_name, "_y")
  bx <- paste0(B_name, "_x"); by <- paste0(B_name, "_y")
  cx <- paste0(C_name, "_x"); cy <- paste0(C_name, "_y")
  
  # Column names for the big-jump and bad-excursion flags across all 3 bodyparts
  bj_cols <- paste0(parts3, "_is_big_jump")
  be_cols <- paste0(parts3, "_bad_excursion")
  
  df_wide <- df_wide %>%
    mutate(
      # Check whether ANY of the 3 bodyparts are flagged on this frame
      any_big_jump = if (drop_big_jump_frames) {
        rowSums(across(all_of(bj_cols), ~ .x %in% TRUE), na.rm = TRUE) > 0
      } else {
        FALSE
      },
      any_excursion = if (drop_bad_excursions) {
        rowSums(across(all_of(be_cols), ~ .x %in% TRUE), na.rm = TRUE) > 0
      } else {
        FALSE
      },
      
      # Compute angle, setting to NA for flagged frames
      angle = if_else(
        any_big_jump | any_excursion,
        NA_real_,
        pmap_dbl(
          list(.data[[ax]], .data[[ay]],
               .data[[bx]], .data[[by]],
               .data[[cx]], .data[[cy]]),
          ~ angle_deg(c(..1, ..2), c(..3, ..4), c(..5, ..6))
        )
      )
    )
  
  df_wide %>% select(frame, individual, angle)
}


#' Flag head swings based on angle deviation from 180 degrees (straight).
#'
#' A head swing is defined as any frame where the body angle deviates from
#' 180 degrees by more than swing_deg degrees.
#'
#' @param df_long   Long-format tibble (must contain an 'angle' column).
#' @param swing_deg Deviation threshold in degrees (default 40).
#' @return Same tibble with added boolean column 'head_swing'.
flag_head_swing <- function(df_long, swing_deg = 40) {
  
  # If no angle column exists, return unchanged (angle computation may have been skipped)
  if (!"angle" %in% names(df_long)) return(df_long)
  
  centre <- 180
  low    <- centre - swing_deg
  high   <- centre + swing_deg
  
  df_long %>%
    mutate(head_swing = !is.na(angle) & (angle < low | angle > high))
}


# ============================================================================
# 5. FILE METADATA EXTRACTION
# ============================================================================

#' Extract the experimental group name from a tracking CSV filename.
#'
#' Expects DLC-style filenames like "GroupName_123456789DLC_..._filtered.csv"
#' and extracts everything before the numeric ID + "DLC" portion. SLEAP
#' filenames have no such marker, so as a fallback we strip a trailing
#' "_analysis" suffix (SLEAP's analysis.csv export convention); if that still
#' doesn't match, the full filename (minus extension) is used so that group
#' is never NA -- an NA group would silently collapse files together in the
#' Summary/Dunnett tabs rather than raising an obvious error.
#'
#' NOTE: this SLEAP fallback is a guess based on one example filename. If your
#' SLEAP files encode group identity differently, this regex will need
#' adjusting -- let me know your naming convention and I can tailor it.
#'
#' @param file_name The basename of the CSV file.
#' @return Character string of the group name (never NA).
extract_group <- function(file_name) {
  # Try the more specific pattern first: group_<digits>DLC
  grp <- str_extract(file_name, "^.*(?=_[0-9]+DLC)")
  
  # Fall back to: group + DLC (no underscore + digits)
  if (is.na(grp)) {
    grp <- str_extract(file_name, "^.*(?=DLC)")
  }
  
  # SLEAP fallback: strip a trailing "_analysis" suffix
  if (is.na(grp)) {
    grp <- str_extract(file_name, "^.*(?=_[Aa]nalysis)")
  }
  
  # Last resort: filename without its extension (guarantees a non-NA result)
  if (is.na(grp)) {
    grp <- tools::file_path_sans_ext(file_name)
  }
  
  grp
}


#' Extract a display title from the filename (currently same as group name).
#'
#' @param file_name The basename of the CSV file.
#' @return Character string for use as a plot title.
extract_title <- function(file_name) {
  extract_group(file_name)
}


# ============================================================================
# 6. MAIN PROCESSING PIPELINE
# ============================================================================

#' Process a single tracking CSV (DLC or SLEAP) through the full analysis pipeline.
#'
#' Reads the file, tidies to long format, computes speed, flags jumps and
#' excursions, classifies movement, and optionally computes body angles.
#' Attaches file metadata (path, name, group, title) to every row.
#'
#' @param path                 File path to the CSV.
#' @param file_name            Display name for the file (defaults to basename).
#' @param bodyparts_keep       Character vector of bodyparts to retain.
#' @param likelihood_min       Minimum likelihood threshold (NULL to skip).
#' @param threshold_px         Movement speed threshold (px/frame).
#' @param window_n             Stillness window size (frames).
#' @param max_jump_px          Fixed max step distance for jump detection (px).
#' @param use_robust_jump      Use adaptive per-track jump threshold?
#' @param robust_mult          Multiplier for robust threshold.
#' @param max_excursion_frames Max frames for a region to be flagged as excursion.
#' @param compute_angle        Compute body angles? (requires exactly 3 bodyparts).
#' @param angle_vertex         Which bodypart is the angle vertex.
#' @param dlc_format           CSV format: "auto", "multi", "single", "flat", or "sleap".
#' @return Long-format tibble with all computed columns, or NULL on failure.
process_one_file <- function(path,
                             file_name            = NULL,
                             bodyparts_keep       = c("mid"),
                             likelihood_min       = NULL,
                             threshold_px         = 2,
                             window_n             = 5,
                             max_jump_px          = 50,
                             use_robust_jump      = FALSE,
                             robust_mult          = 10,
                             max_excursion_frames = 10,
                             compute_angle        = FALSE,
                             angle_vertex         = NULL,
                             dlc_format           = c("auto", "multi", "single", "flat", "sleap")) {
  
  dlc_format <- match.arg(dlc_format)
  
  # Resolve "auto" once, up front, so the rest of this function only ever
  # deals with a concrete format
  if (dlc_format == "auto") dlc_format <- detect_dlc_format(path)
  
  # Use the provided display name, or fall back to the file basename
  fn <- if (!is.null(file_name) && nzchar(file_name)) file_name else basename(path)
  message("[process] Starting: ", fn)
  
  # --- Step 1: Read raw CSV ---
  # SLEAP's export shape is different enough from DLC's that it gets its own
  # reader + tidier pair (read_sleap_csv()/tidy_sleap()); both converge on the
  # same long-format shape, so every step from here on is format-agnostic
  if (dlc_format == "sleap") {
    df_raw  <- read_sleap_csv(path)
    tidy_fn <- tidy_sleap
  } else {
    df_raw  <- read_dlc_filtered_csv(path, format = dlc_format)
    tidy_fn <- tidy_dlc
  }
  
  if (is.null(df_raw)) {
    warning("Could not read file: ", fn, call. = FALSE)
    return(NULL)
  }
  
  #  check that we have data rows after reading
  if (nrow(df_raw) == 0) {
    warning("File contains no data rows after cleaning: ", fn, call. = FALSE)
    return(NULL)
  }
  
  message("[process] ", fn, ": ", nrow(df_raw), " data rows, ",
          "frame range ", min(df_raw$frame), " to ", max(df_raw$frame))
  
  # --- Step 2: Tidy to long format, compute speed, flag artefacts ---
  df_long <- df_raw %>%
    tidy_fn(bodyparts_keep = bodyparts_keep, likelihood_min = likelihood_min) %>%
    add_speed() %>%
    filter_big_jumps(
      max_jump_px     = max_jump_px,
      use_robust_jump = use_robust_jump,
      robust_mult     = robust_mult
    ) %>%
    flag_short_excursions(max_excursion_frames = max_excursion_frames) %>%
    add_movement_flag(threshold_px = threshold_px, window_n = window_n)
  
  # check we still have data after all filtering
  if (nrow(df_long) == 0) {
    warning("No rows remain after filtering for file: ", fn, call. = FALSE)
    return(NULL)
  }
  
  # --- Step 3: Compute body angle (optional, requires exactly 3 bodyparts) ---
  if (isTRUE(compute_angle) && !is.null(bodyparts_keep) && length(bodyparts_keep) == 3) {
    
    # Determine the vertex bodypart
    vtx <- if (!is.null(angle_vertex) && angle_vertex %in% bodyparts_keep) {
      angle_vertex
    } else {
      bodyparts_keep[2]
    }
    
    ang <- compute_angles(
      df_long,
      parts3                = bodyparts_keep,
      vertex                = vtx,
      drop_big_jump_frames = TRUE,
      drop_bad_excursions  = TRUE
    )
    
    df_long <- df_long %>%
      left_join(ang, by = c("frame", "individual")) %>%
      flag_head_swing(swing_deg = 40)
  }
  
  # --- Step 4: Attach file metadata ---
  df_long %>%
    mutate(
      file_path = path,
      file_name = fn,
      title     = extract_title(fn),
      group     = extract_group(fn)
    )
}


# ============================================================================
# 7. PLOTTING FUNCTIONS
# ============================================================================

#' Plot basic bodypart trajectory coloured by bodypart.
#'
#' @param df_long Long-format tibble with x, y, bodypart, individual columns.
#' @param xlim    x-axis limits (pixels).
#' @param ylim    y-axis limits (pixels).
#' @param ttl     Plot title string.
#' @return A ggplot object.
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


#' Plot speed over time, faceted by individual and bodypart.
#'
#' @param df_long Long-format tibble with frame, speed, individual, bodypart columns.
#' @param ttl     Plot title string.
#' @return A ggplot object.
plot_speed_trace <- function(df_long, ttl = "") {
  
  ggplot(df_long, aes(x = frame, y = speed)) +
    geom_line(alpha = 0.6) +
    facet_grid(individual ~ bodypart, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = ttl, x = "Frame", y = "Speed (px/frame)")
}


#' Build segment data for trajectory plots coloured by moving/still status.
#'
#' Creates segment endpoints (x -> xend, y -> yend) and flags whether each
#' segment crosses a big jump boundary.
#'
#' @param df_long        Long-format tibble with x, y, moving, thr columns.
#' @param drop_big_jumps If TRUE, remove segments that exceed the jump threshold.
#' @return Tibble with segment columns added.
make_segments <- function(df_long, drop_big_jumps = TRUE) {
  
  seg <- df_long %>%
    arrange(individual, bodypart, frame) %>%
    group_by(individual, bodypart) %>%
    mutate(
      # Segment endpoint: the next frame's coordinates
      xend = lead(x),
      yend = lead(y),
      
      # Segment distance (for jump detection on the drawn segment)
      seg_dist = sqrt((xend - x)^2 + (yend - y)^2),
      
      # Carry forward the movement status for colouring
      moving_seg = moving,
      
      # Flag if this segment itself is a big jump
      big_jump_seg = seg_dist > thr
    ) %>%
    ungroup() %>%
    # Drop incomplete segments (last frame in each track has no endpoint)
    filter(!is.na(xend), !is.na(yend))
  
  if (drop_big_jumps) seg <- filter(seg, !big_jump_seg)
  
  seg
}


#' Plot trajectory segments coloured by moving vs still status.
#'
#' @param df_long         Long-format tibble.
#' @param xlim            x-axis limits.
#' @param ylim            y-axis limits.
#' @param ttl             Plot title.
#' @param moving_col      Colour for moving segments.
#' @param still_col       Colour for still segments.
#' @param line_width      Segment line width.
#' @param axis_title_size Font size for axis titles.
#' @param axis_text_size  Font size for axis tick labels.
#' @param strip_text_size Font size for facet strip labels.
#' @param legend_position Legend position string.
#' @param drop_big_jumps  Remove big-jump segments?
#' @return A ggplot object.
plot_trajectory_coloured_segments <- function(df_long,
                                              xlim            = c(0, 1280),
                                              ylim            = c(0, 960),
                                              ttl             = "",
                                              moving_col      = "#11F011",
                                              still_col       = "#FA05EE",
                                              line_width      = 0.6,
                                              axis_title_size = 13,
                                              axis_text_size  = 11,
                                              strip_text_size = 13,
                                              legend_position = "top",
                                              drop_big_jumps  = TRUE) {
  
  seg <- make_segments(df_long, drop_big_jumps = drop_big_jumps)
  
  # Convert the boolean moving flag to a labelled factor for the legend
  seg <- seg %>%
    mutate(
      moving_seg = case_when(
        is.na(moving_seg) ~ "Start",
        moving_seg         ~ "Moving",
        TRUE               ~ "Still"
      ),
      moving_seg = factor(moving_seg, levels = c("Moving", "Still", "Start"))
    )
  
  ggplot(seg, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment(
      aes(colour = moving_seg, linetype = bodypart),
      linewidth = line_width,
      alpha     = 0.9
    ) +
    facet_wrap(~ individual, ncol = 2) +
    coord_fixed(xlim = xlim, ylim = ylim) +
    scale_colour_manual(
      values = c(Moving = moving_col, Still = still_col, Start = "black"),
      drop   = FALSE
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = legend_position,
      legend.title    = element_text(face = "bold"),
      strip.text      = element_text(face = "bold", size = strip_text_size),
      axis.title      = element_text(size = axis_title_size),
      axis.text       = element_text(size = axis_text_size),
      plot.title      = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = ttl, x = "x (px)", y = "y (px)", colour = NULL)
}


#' Plot body angle over time with head-swing highlighting.
#'
#' Draws a shaded band around 180 degrees (straight body) and colours frames
#' flagged as head swings in red.
#'
#' @param df_long Long-format tibble with frame, angle, head_swing, individual columns.
#' @param ttl     Plot title string.
#' @return A ggplot object.
plot_angle_trace <- function(df_long, ttl = "") {
  
  ggplot(df_long, aes(x = frame, y = angle)) +
    # Shaded band around 180 degrees (the "straight body" zone)
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = 140,  ymax = 220,
      alpha = 0.15
    ) +
    # Reference line at 180 degrees
    geom_hline(yintercept = 180, colour = "red3", linetype = "dashed") +
    # Points coloured by head-swing status
    geom_point(
      aes(group = interaction(individual), colour = head_swing),
      alpha = 0.7
    ) +
    # Connecting line
    geom_line(alpha = 0.6, linetype = "longdash") +
    scale_colour_manual(
      values   = c(`TRUE` = "red", `FALSE` = "black"),
      na.value = "black"
    ) +
    guides(colour = "none") +
    facet_wrap(~ individual, ncol = 2, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = ttl, x = "Frame", y = "Angle (degrees)")
}


# ============================================================================
# 8. PLOTLY HELPERS
# ============================================================================

#' Set consistent axis ranges across all subplots in a plotly object.
#'
#' After converting a faceted ggplot to plotly, each facet gets its own axis
#' (xaxis, xaxis2, etc.). This function overrides all of them to use the
#' same range, ensuring consistent zoom across facets.
#'
#' @param p       A plotly object (from ggplotly()).
#' @param x_range Numeric vector of length 2 for x-axis range, or NULL.
#' @param y_range Numeric vector of length 2 for y-axis range, or NULL.
#' @return The modified plotly object.
set_all_subplot_ranges <- function(p, x_range = NULL, y_range = NULL) {
  
  b <- plotly::plotly_build(p)
  lay_names <- names(b$x$layout)
  
  # Find all x-axis and y-axis names (xaxis, xaxis2, xaxis3, ...)
  x_axes <- grep("^xaxis", lay_names, value = TRUE)
  y_axes <- grep("^yaxis", lay_names, value = TRUE)
  
  # Override each axis with the specified range using walk() instead of for loop
  if (!is.null(x_range)) {
    walk(x_axes, ~ {
      b$x$layout[[.x]]$range     <<- x_range
      b$x$layout[[.x]]$autorange <<- FALSE
    })
  }
  
  if (!is.null(y_range)) {
    walk(y_axes, ~ {
      b$x$layout[[.x]]$range     <<- y_range
      b$x$layout[[.x]]$autorange <<- FALSE
    })
  }
  
  b
}