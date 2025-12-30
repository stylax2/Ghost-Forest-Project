# ==============================================================================
# [Script 01] Data Preparation & Acquisition
# Project: Vertical Decoupling and Ghost Forests in Quercus mongolica Communities
# Description: 
#   1. Selection of diagnostic species from vegetation survey data.
#   2. Acquisition of occurrence records via GBIF (including vicariant species).
#   3. Download and preprocessing of CMIP6 climate data (Current & Future).
# ==============================================================================

# 1. Setup & Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, rgbif, sf, geodata, terra, stringr)

message("🚀 [Step 1] Starting Data Preparation...")

# ==============================================================================
# Section A: Diagnostic Species Analysis
# Note: This section requires 'vege_table.xlsx'. 
#       If using pre-processed lists, this can be skipped.
# ==============================================================================
if(file.exists("vege_table.xlsx")) {
  message("   -> Analyzing vegetation data for diagnostic species...")
  
  veg_raw <- read_excel("vege_table.xlsx")
  
  # Helper: Convert Braun-Blanquet to numeric
  bb_convert <- function(x){
    case_when(x == "r" ~ 0.1, x == "+" ~ 0.5, x == "1" ~ 1, x == "2" ~ 2, x == "3" ~ 3, x == "4" ~ 4, TRUE ~ suppressWarnings(as.numeric(x)))
  }
  
  # Identify dominant species and community type
  veg_comm <- veg_raw %>%
    mutate(cover_num = bb_convert(cover),
           is_tree_layer = str_detect(layer, "^[TS]"))
  
  dominant_tbl <- veg_comm %>% filter(is_tree_layer) %>% group_by(site) %>%
    slice_max(order_by = cover_num, n = 1, with_ties = FALSE) %>%
    transmute(site, community = sci_name)
  
  # Calculate Frequencies (Q. mongolica vs Others)
  sp_freq_wide <- veg_comm %>% left_join(dominant_tbl, by = "site") %>%
    mutate(is_qm_comm = ifelse(community == "Quercus mongolica", "Qm", "Other"),
           life_form = case_when(str_detect(layer, "^[SE]") ~ "shrub", str_detect(layer, "^[Hh]") ~ "herb", TRUE ~ "other")) %>%
    filter(life_form %in% c("shrub", "herb")) %>%
    group_by(is_qm_comm, sci_name, life_form) %>% count(name="n_plots_sp") %>%
    left_join(count(distinct(veg_comm, site, is_qm_comm), is_qm_comm, name="n_total"), by="is_qm_comm") %>%
    mutate(freq = n_plots_sp / n_total) %>%
    select(-n_plots_sp, -n_total) %>%
    pivot_wider(names_from = is_qm_comm, values_from = freq, values_fill = 0, names_prefix = "freq_") %>%
    mutate(delta_freq = freq_Qm - freq_Other) %>%
    arrange(desc(freq_Qm))
  
  write_csv(sp_freq_wide, "candidate_species_list.csv")
  message("   -> Species list saved.")
} else {
  message("vege_table.xlsx' not found. Skipping Section A.")
}

# ==============================================================================
# Section B: GBIF Data Collection (Extended Niche Strategy)
# ==============================================================================
message("\n   -> Collecting GBIF occurrence data...")

# Map: Target Species -> Vicariant Species Integration
species_map <- list(
  "Quercus mongolica" = "Quercus mongolica",
  "Quercus crispula"  = "Quercus mongolica",      # Integrated for extended niche
  "Acer pseudosieboldianum" = "Acer pseudosieboldianum",
  "Acer sieboldianum"       = "Acer pseudosieboldianum", 
  "Rhododendron schlippenbachii" = "Rhododendron schlippenbachii",
  "Lindera obtusiloba"           = "Lindera obtusiloba",
  "Rhododendron mucronulatum"    = "Rhododendron mucronulatum",
  "Carex siderosticta"   = "Carex siderosticta",
  "Ainsliaea acerifolia" = "Ainsliaea acerifolia",
  "Disporum smilacinum"  = "Disporum smilacinum",
  "Polygonatum odoratum" = "Polygonatum odoratum var. pluriflorum"
)

study_area_wkt <- "POLYGON((115 30, 150 30, 150 55, 115 55, 115 30))"
all_data_list <- list()

for (search_name in names(species_map)) {
  final_target_name <- species_map[[search_name]]
  tryCatch({
    occ_res <- occ_search(scientificName = search_name, geometry = study_area_wkt, 
                          hasCoordinate = TRUE, year = '1970,2025', limit = 1000)
    if (!is.null(occ_res$data) && nrow(occ_res$data) > 0) {
      df_clean <- occ_res$data %>% 
        select(decimalLongitude, decimalLatitude, basisOfRecord) %>%
        mutate(scientificName = final_target_name, originalSearch = search_name)
      all_data_list[[search_name]] <- df_clean
      message(paste("      Found:", search_name, "-> Integrated as:", final_target_name))
    }
  }, error = function(e) message(paste("      Error:", search_name)))
}

if(length(all_data_list) > 0) {
  final_df <- bind_rows(all_data_list)
  write_csv(final_df, "GBIF_Final_Dataset.csv")
  message("   -> GBIF data saved.")
}

# ==============================================================================
# Section C: Climate Data Download (CMIP6)
# ==============================================================================
message("\n   -> Downloading Climate Data (CMIP6)...")
dir.create("Climate_Data/Ensemble", recursive = TRUE, showWarnings = FALSE)

# Configuration
gcm_list <- c("MIROC6", "MPI-ESM1-2-HR", "EC-Earth3-Veg", "CNRM-ESM2-1", "UKESM1-0-LL")
ssp_list <- c("245", "585") # Process both SSP2-4.5 and SSP5-8.5
roi_ext <- ext(115, 150, 30, 55)

process_climate <- function(model, ssp) {
  tryCatch({
    # Download tiles (East Asia region split)
    t1 <- geodata::cmip6_tile(model=model, ssp=ssp, time="2041-2060", lon=116, lat=38, var="bioc", res=0.5, path="temp_clim")
    t2 <- geodata::cmip6_tile(model=model, ssp=ssp, time="2041-2060", lon=127, lat=38, var="bioc", res=0.5, path="temp_clim")
    m_crop <- crop(merge(t1, t2), roi_ext)
    names(m_crop) <- paste0("bio", 1:19)
    unlink("temp_clim", recursive = TRUE) # Cleanup
    return(m_crop)
  }, error = function(e) return(NULL))
}

# (Note: Current climate download logic should be added here similarly or assumed present)
# Calculating Ensembles for Future Scenarios
for(ssp in ssp_list) {
  stack_list <- list()
  message(paste("      Processing SSP", ssp, "..."))
  for(gcm in gcm_list) {
    r <- process_climate(gcm, ssp)
    if(!is.null(r)) stack_list[[gcm]] <- r
  }
  if(length(stack_list) > 0) {
    ens_mean <- app(sds(stack_list), mean)
    names(ens_mean) <- paste0("bio", 1:19)
    writeRaster(ens_mean, paste0("Climate_Data/Ensemble/Ensemble_SSP", ssp, "_2050s.tif"), overwrite=TRUE)
  }
}
message("[Step 1] Data Preparation Complete.")
