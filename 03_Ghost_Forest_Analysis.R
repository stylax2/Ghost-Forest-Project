# ==============================================================================
# [Script 03] Ghost Forest Analysis & Visualization
# Project: Vertical Decoupling and Ghost Forests in Quercus mongolica Communities
# Description: 
#   1. Calculate Species Richness (S-SDM) for Canopy & Understory.
#   2. Compute 'Ghost Forest Risk Index' based on reference baseline.
#   3. Generate statistical summaries and Figure 5.
# ==============================================================================

rm(list = ls()); gc()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(terra, tidyverse, patchwork, ggpubr)

message("[Step 3] Ghost Forest Analysis...")

base_dir <- "Final_Result"
dir_bin <- file.path(base_dir, "01_Binary_Maps")
dir_ssdm <- file.path(base_dir, "05_SSDM")
dir_out <- file.path(base_dir, "06_Ghost_Forest")
if(!dir.exists(dir_out)) dir.create(dir_out, recursive=TRUE)

scenarios <- c("Current", "SSP245", "SSP585")

# 1. Calculate SSDM (Stacked Species Distribution Models)
message("\n   -> Calculating S-SDM (Species Richness)...")
for(scen in scenarios) {
  files <- list.files(dir_bin, pattern = paste0("Binary_", scen, "_.*\\.tif$"), full.names = T)
  if(length(files) == 0) next
  
  rich_map <- sum(rast(files), na.rm=TRUE)
  # Separate Understory: Total - Canopy(Quercus)
  qm_file <- files[str_detect(files, "Quercus.*mongolica")]
  if(length(qm_file) > 0) {
    und_rich <- rich_map - rast(qm_file)
    und_rich <- ifel(und_rich < 0, 0, und_rich)
    writeRaster(und_rich, file.path(dir_ssdm, paste0("SSDM_Understory_", scen, ".tif")), overwrite=TRUE)
  }
}

# 2. Ghost Forest Classification
message("\n   -> Classifying Ghost Forests...")
# Define Baseline (Reference Condition)
r_qm_curr <- rast(list.files(dir_bin, pattern="Current_Quercus.*mongolica", full.names=T))
r_und_curr <- rast(file.path(dir_ssdm, "SSDM_Understory_Current.tif"))
S_ref <- as.numeric(global(mask(r_und_curr, r_qm_curr, maskvalues=0), "mean", na.rm=T))
message(paste("      Reference Baseline (S_ref):", round(S_ref, 2)))

stats_list <- list()
for(scen in scenarios) {
  f_qm <- list.files(dir_bin, pattern=paste0(scen, "_Quercus.*mongolica"), full.names=T)
  f_und <- file.path(dir_ssdm, paste0("SSDM_Understory_", scen, ".tif"))
  
  if(length(f_qm)==0) next
  
  # Logic: 1=Collapsed, 2=Ghost (Und < 50% ref), 3=Healthy
  r_qm <- rast(f_qm); r_und <- rast(f_und)
  ghost_map <- ifel(r_qm == 0, 1, ifel(r_und < (S_ref * 0.5), 2, 3))
  writeRaster(ghost_map, file.path(dir_out, paste0("Ghost_Map_", scen, ".tif")), overwrite=TRUE)
  
  # Stats
  df <- as.data.frame(freq(ghost_map)) %>% 
    mutate(Scenario=scen, Category=case_when(value==1~"Collapsed", value==2~"Ghost Forest", value==3~"Healthy Forest")) %>%
    filter(!is.na(Category)) %>% select(Scenario, Category, Pixel_Count=count)
  stats_list[[scen]] <- df
}

final_stats <- bind_rows(stats_list)
write.csv(final_stats, file.path(dir_out, "Ghost_Forest_Stats.csv"), row.names=F)

# 3. Figure 5 Visualization
message("\n   -> Generating Figure 5...")
df_summary <- final_stats %>% filter(Category != "Collapsed") %>%
  group_by(Scenario) %>% mutate(Total=sum(Pixel_Count), Ratio=Pixel_Count/Total*100) %>%
  mutate(Scenario = factor(Scenario, levels=c("Current","SSP245","SSP585"), labels=c("Current","SSP2-4.5","SSP5-8.5")))

# (A) Stacked Bar
p1 <- ggplot(df_summary, aes(x=Scenario, y=Pixel_Count, fill=Category)) +
  geom_bar(stat="identity", width=0.6) +
  geom_text(aes(label=paste0(round(Ratio,1),"%")), position=position_stack(vjust=0.5), color="white", fontface="bold") +
  scale_fill_manual(values=c("Ghost Forest"="#D62728", "Healthy Forest"="#2CA02C")) +
  labs(title="(A) Degradation of Forest Composition", y="Area (Pixels)") + theme_classic()

# (B) Decoupling Gap
df_line <- df_summary %>% select(Scenario, Category, Pixel_Count) %>%
  pivot_wider(names_from=Category, values_from=Pixel_Count) %>%
  mutate(Species_Range=`Healthy Forest`+`Ghost Forest`, Community_Range=`Healthy Forest`)

p2 <- ggplot(df_line, aes(x=Scenario)) +
  geom_ribbon(aes(ymin=Community_Range, ymax=Species_Range, group=1), fill="#D62728", alpha=0.2) +
  geom_line(aes(y=Species_Range, group=1, color="Species Range"), size=1.2) +
  geom_line(aes(y=Community_Range, group=1, color="Community Range"), size=1.2) +
  scale_color_manual(values=c("Species Range"="black", "Community Range"="#2CA02C")) +
  labs(title="(B) Decoupling of Species vs. Community", y=NULL) + theme_classic() + theme(axis.text.y=element_blank())

ggsave(file.path(base_dir, "04_Figures", "Figure_5.tiff"), p1+p2, width=12, height=6, dpi=600)
message("[Step 3] Analysis Complete.")
