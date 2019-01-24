
# This function create dataframe with features for single user
features.extract <- function(df) {
  
  #df: input dataframe
  
  require(tidyverse)
  require(zoo)
  
  features <- list()
  
  ### 1. Mean group ###
  
  # RAW DATA
  #
  # (1)-1.1 raw_mean_ax
  features$raw_mean_ax <- zoo(df$ax) %>%
    rollapply(width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round
  
  # (2)-1.2 raw_mean_ay
  features$raw_mean_ay <- zoo(df$ay) %>%
    rollapply(width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round
  
  # (3)-1.3 raw_mean_az
  features$raw_mean_az <- zoo(df$az) %>%
    rollapply(width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round
  
  # FFT DATA
  # 
  # (4)-1.4 fft_mean_ax
  features$fft_mean_ax <- df$ax %>% diff %>% fft %>% Mod %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round
  
  # (5)-1.5 fft_mean_ay
  features$fft_mean_ay <- df$ay %>% diff %>% fft %>% Mod %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round
  
  # (6)-1.6 fft_mean_az
  features$fft_mean_az <- df$az %>% diff %>% fft %>% Mod %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round
  
  ### 2. Median group ###
  
  # RAW DATA
  #
  # (7)-2.1 raw_median_ax
  features$raw_median_ax <- df$ax %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = median) %>% as.numeric %>% round
  
  # (8)-2.2 raw_median_ay
  features$raw_median_ay <- df$ay %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = median) %>% as.numeric %>% round  
  
  # (9)-2.3 raw_median_az
  features$raw_median_az <- df$az %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = median) %>% as.numeric %>% round  
  
  # FFT DATA
  #
  # (10)-2.4 fft_median_ax
  features$fft_median_ax <- df$ax %>% diff %>% fft %>% Mod %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = median) %>% as.numeric %>% round
  
  # (11)-2.5 fft_median_ay
  features$fft_median_ay <- df$ay %>% diff %>% fft %>% Mod %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = median) %>% as.numeric %>% round
  
  # (12)-2.6 fft_median_az
  features$fft_median_az <- df$az %>% diff %>% fft %>% Mod %>% zoo %>% 
    rollapply(width = sample_N, by = sample_N * overlap, FUN = median) %>% as.numeric %>% round
  
  ### 3. Magnitude
  
  # (13)-3.1 raw_mag
  features$raw_mag <- list(df$ax, df$ay, df$az) %>% 
    pmap_dbl(function(ax, ay, az) sqrt(ax^2 + ay^2 + az^2)) %>%  
    rollapply(width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round
  
  # (14)-3.2 fft_mag
  features$fft_mag <- list(df$ax, df$ay, df$az) %>%
    map(~ .x %>% diff %>% fft %>% Mod) %>%                       # make fft from raw data
    pmap_dbl(function(ax, ay, az) sqrt(ax^2 + ay^2 + az^2)) %>%  # make magnitude for one widow
    rollapply(width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round

  ### 4. Cross-correlation (both for time and frequancy domain)
  ### 
  ### RAW DATA
  # (15)-4.1 corr_xy
  features$corr_xy <- features$raw_mean_ax / features$raw_mean_ay
  
  # (16)-4.2 corr_zy
  features$corr_zy <- features$raw_mean_az / features$raw_mean_ay
  
  ### FFT DATA
  ### 
  # (17)-4.3 corr_xy_fft
  features$corr_xy_fft <- features$fft_mean_ax / features$fft_mean_ay
  
  # (18)-4.4 corr_zy_fft
  features$corr_zy_fft <- features$fft_mean_az / features$fft_mean_ay
  
  ### 5. Peak Count (only for time domain)
  
  # (19)-5.1 mean_peak_ax
  # (20)-5.2 mean_peak_ay  
  # (21)-5.3 mean_peak_az
  
  # Set names for peaks features
  mean_peak <- c("mean_peak_ax", "mean_peak_ay", "mean_peak_az")
  
  # Calculate tree features all at once
  map(list(df$ax, df$ay, df$az),
      function(x) { rollapply(x, width = sample_N, by = sample_N * overlap,
                             FUN = function(z) {
                               id <- 1:length(z)
                               peaks <- argmax(x = id, y = z,
                                               w = peak.detection.window,
                                               span = peak.detection.span)
                               return(length(peaks$x) / sample_N)}
                             ) %>% as.numeric}) %>%
    set_names(mean_peak) %>% append(features) -> features
  
  ### 6. Distance between Peaks (time domain only)
  #
  # (22)-6.1 peak_avg_dist_ax
  # (23)-6.2 peak_avg_dist_ay
  # (24)-6.3 peak_avg_dist_az
  # 
  # Set names for peaks features
  dist_peak <- c("peak_avg_dist_ax", "peak_avg_dist_ay", "peak_avg_dist_az")
  
  map(list(df$ax, df$ay, df$az),
      function(x) { rollapply(x, width = sample_N, by = sample_N * overlap,
                              FUN = function(z) {
                                id <- 1:length(z)
                                peaks <- argmax(x = id, y = z,
                                                w = peak.detection.window,
                                                span = peak.detection.span)
                                return(mean(diff(peaks$i)))}
                              ) %>% as.numeric}) %>% 
    set_names(dist_peak) %>% append(features) -> features
  
  ### 7. Standard deviation of mean distance of peaks diffs (time domain only)
  ### 
  # (25)-7.1 peak_sd_dist_ax
  # (26)-7.2 peak_sd_dist_ay
  # (27)-7.3 peak_sd_dist_az
  # 
  # Set names for peaks features
  sd_peak <- c("peak_sd_dist_ax", "peak_sd_dist_ay", "peak_sd_dist_az")
  
  map(list(df$ax, df$ay, df$az),
      function(x) { rollapply(x, width = sample_N, by = sample_N * overlap,
                              FUN = function(z) {
                                id <- 1:length(z)
                                peaks <- argmax(x = id, y = z,
                                                w = peak.detection.window,
                                                span = peak.detection.span)
                                return(sd(diff(peaks$i)) / sample_N)}
      ) %>% as.numeric}) %>% 
    set_names(sd_peak) %>% append(features) -> features
  
  ### 8. Spectral Centroid (both time and frequancy domain)
  ### 
  # (28)-8.1 spectral_centr_ax
  # (29)-8.2 spectral_centr_ay
  # (30)-8.3 spectral_centr_az
  # 
  # Set names for peaks features
  spec_centr <- c("spectral_centr_ax", "spectral_centr_ay", "spectral_centr_az")
  
  map(list(df$ax, df$ay, df$az),
      function(x) {
        # Calculate FFT. After applaying diff() input array has n-1 length!
        mod.fft <- Mod(fft(diff(x))) %>% round
        
        # Align the vector along the length by adding the median as the last element
        mod.fft[length(mod.fft) + 1] <- mod.fft[median(mod.fft)]
        
        # Calculate spectral centroid through all windows
        rollapply(zoo(x * mod.fft), width = sample_N, by = sample_N * overlap, FUN = mean) %>% as.numeric %>% round
  }) %>% set_names(spec_centr) %>% append(features) -> features
  
  ### 9. Average Difference from Mean  (time domain only)
  ### 
  ### (31)-9.1 avg_diff_mean_ax
  ### (32)-9.2 avg_diff_mean_ay
  ### (33)-9.3 avg_diff_mean_az
  ### 
  ### Set names for peaks features
  avg_diff <- c("avg_diff_mean_ax", "avg_diff_mean_ay", "avg_diff_mean_az")
  
  map(list(df$ax, df$ay, df$az), 
      function(z) { rollapply(z, width = sample_N, by = sample_N * overlap,
                            FUN = function(x) mean(abs(x - mean(x)))) %>% as.numeric %>% round}) %>% 
    set_names(avg_diff) %>% append(features) -> features
  
  return(features)
}

##### Auxiliary functions for peak detections 

# Function find peaks by given parameters
argmax <- function(x, y, w = 1, ...) {
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x = x[i.max], i = i.max, y.hat = y.smooth)
}


# === TEST ZONE ===

## For testing run user1 data from features_extraction_single_user.Rmd

##########################
### TEST features.extract 
### 
#Parameters for peak.detection function:
peak.detection.window <- 3 # w
peak.detection.span <- 0.1 # span

f.e.test <- features.extract(user1)
map(f.e.test, summary)

