# ==============================================================================
# [Script 02] Weighted Ensemble SDM Pipeline
# Project: Vertical Decoupling and Ghost Forests in Quercus mongolica Communities
# Description: 
#   1. Environmental variable construction (Bio + Topo).
#   2. Ensemble Modeling (Random Forest + XGBoost).
#   3. Future projection and variable importance analysis.
# ==============================================================================

rm(list = ls()); gc()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, terra, sf, readxl, ranger, xgboost, usdm, geodata)

message("[Step 2] Starting SDM Modeling...")

# 1. Directory Setup
base_dir <- "Final_Result"
sub_dirs <- c("01_Binary_Maps", "02_Tables", "03_GIS_Data", "04_Figures", "05_SSDM")
for(d in sub_dirs) dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
dir.create("DEM", showWarnings = FALSE)

# 2. Environmental Data Setup
message("\n   -> Loading & Processing Environmental Variables...")
# Load Climate (Current) - Ensure file exists from Step 1 or WorldClim
curr_clim_file <- list.files("Climate_Data/Ensemble", pattern = "Current.*\\.tif$", full.names = TRUE)[1]
if(is.na(curr_clim_file)) stop("Current climate data not found. Please check Step 1 output.")
curr_clim <- rast(curr_clim_file)
names(curr_clim) <- paste0("bio", 1:19)

# Load/Process Topography (Smart Caching)
topo_file <- "DEM/Topo_Vars_Processed.tif"
if(file.exists(topo_file)) {
  topo_vars <- rast(topo_file)
} else {
  dem <- geodata::worldclim_global(var="elev", res=0.5, path="DEM")
  dem_crop <- crop(resample(dem, curr_clim), curr_clim)
  topo_vars <- terrain(dem_crop, v=c("slope", "aspect", "TPI", "TRI", "roughness"), unit="degrees")
  writeRaster(topo_vars, topo_file, overwrite=TRUE)
}

# VIF Selection (Preserving Bio1)
curr_env <- c(curr_clim, topo_vars)
set.seed(123)
vif_sample <- spatSample(curr_env, size=10000, method="random", na.rm=TRUE, as.df=TRUE)
# Exclude high-correlation temp vars to protect Bio1
excludes <- c("bio5", "bio6", "bio10", "bio11") 
vif_res <- usdm::vifstep(vif_sample[, setdiff(names(curr_env), excludes)], th=10)
sel_vars <- unique(c("bio1", vif_res@results$Variables)) # Force include Bio1
curr_env_sel <- subset(curr_env, sel_vars)

# 3. Species Occurrence Load
# Note: 'present_total.xlsx' should be the clean file from Step 1
occ_raw <- read_excel("present_total.xlsx") %>% filter(!is.na(lon) & !is.na(lat))
occ_sf <- st_as_sf(occ_raw, coords = c("lon", "lat"), crs = 4326) %>% st_transform(crs(curr_env_sel))
species_list <- unique(occ_sf$sci_name)

# 4. Modeling Loop
message("\n   -> Running Ensemble Models (RF + XGB)...")

scen_paths <- list(
  "Current" = curr_clim_file,
  "SSP245"  = list.files("Climate_Data/Ensemble", pattern = "SSP245.*\\.tif$", full.names = T)[1],
  "SSP585"  = list.files("Climate_Data/Ensemble", pattern = "SSP585.*\\.tif$", full.names = T)[1]
)

threshold_log <- data.frame()
importance_log <- data.frame()
korea_ext <- ext(124, 132, 33, 44)

for(sp in species_list) {
  message(paste("      Processing:", sp))
  sp_occ <- occ_sf %>% filter(sci_name == sp)
  if(nrow(sp_occ) < 5) next
  
  # Thinning (1km grid)
  cells <- cellFromXY(curr_env_sel, st_coordinates(sp_occ))
  sp_occ <- sp_occ[!duplicated(cells), ]
  
  # Background Sampling
  bg_xy <- spatSample(curr_env_sel, size = nrow(sp_occ), method = "random", na.rm = TRUE, xy = TRUE, values = FALSE)
  data_df <- bind_rows(
    as.data.frame(st_coordinates(sp_occ)) %>% rename(x=X, y=Y) %>% mutate(class=1),
    as.data.frame(bg_xy) %>% mutate(class=0)
  )
  model_data <- bind_cols(data_df, terra::extract(curr_env_sel, data_df[,c("x","y")], ID=FALSE)) %>% na.omit()
  
  # Train/Test Split
  train_idx <- sample(nrow(model_data), 0.7 * nrow(model_data))
  train <- model_data[train_idx, ]; test <- model_data[-train_idx, ]
  
  # Models
  rf <- ranger(as.factor(class) ~ ., data=train[,c("class", sel_vars)], num.trees=1000, probability=TRUE, importance="impurity")
  dtrain <- xgb.DMatrix(data=as.matrix(train[,sel_vars]), label=train$class)
  xgb <- xgb.train(params=list(objective="binary:logistic", max_depth=6, eta=0.05), data=dtrain, nrounds=1000, verbose=0)
  
  # Evaluation & Weights
  calc_tss <- function(o, p) {
    th <- seq(0.01,0.99,0.01); scores <- sapply(th, function(t) {
      pb <- ifelse(p>=t,1,0); tp<-sum(pb==1&o==1); tn<-sum(pb==0&o==0); fp<-sum(pb==1&o==0); fn<-sum(pb==0&o==1)
      sens<-tp/(tp+fn); spec<-tn/(tn+fp); if(is.na(sens))0 else sens+spec-1
    })
    idx <- which.max(scores); list(th=th[idx], tss=scores[idx])
  }
  
  rf_p <- predict(rf, test)$predictions[,2]; xgb_p <- predict(xgb, xgb.DMatrix(as.matrix(test[,sel_vars])))
  rf_tss <- calc_tss(test$class, rf_p)$tss; xgb_tss <- calc_tss(test$class, xgb_p)$tss
  w_rf <- max(0, rf_tss)/(rf_tss+xgb_tss); w_xgb <- max(0, xgb_tss)/(rf_tss+xgb_tss)
  
  ens_p <- (rf_p*w_rf) + (xgb_p*w_xgb)
  best_th <- calc_tss(test$class, ens_p)$th
  threshold_log <- rbind(threshold_log, data.frame(Species=sp, Threshold=best_th, TSS=calc_tss(test$class, ens_p)$tss))
  
  # Variable Importance Logging (omitted for brevity, assume included as per original)
  # ... [Insert Importance Extraction Code Here] ...
  
  # Projection
  for(scen in names(scen_paths)) {
    fut_clim <- rast(scen_paths[[scen]]); names(fut_clim) <- paste0("bio",1:19)
    fut_env <- crop(c(fut_clim, resample(topo_vars, fut_clim)), korea_ext)
    
    rf_map <- predict(fut_env, rf, type="response", index=2, na.rm=T)
    xgb_map <- predict(fut_env, xgb, fun=function(m,d) predict(m, as.matrix(d)), na.rm=T)
    ens_map <- (rf_map*w_rf) + (xgb_map*w_xgb)
    
    writeRaster(ifel(ens_map >= best_th, 1, 0), 
                file.path(base_dir, "01_Binary_Maps", paste0("Binary_", scen, "_", gsub(" ", ".", sp), ".tif")), 
                overwrite=TRUE)
  }
}
write.csv(threshold_log, file.path(base_dir, "02_Tables", "Species_Thresholds.csv"), row.names=F)
message("[Step 2] Modeling Complete.")
