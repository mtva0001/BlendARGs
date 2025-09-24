# BlendARGs project - R scripts
# 16/10/2024

#Load packages
library(readxl)
library(ggOceanMaps)
library(elementalist)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggh4x)
library(viridis)
library(dplyr)
library(gggenomes)
library(tidyr)
library(reshape2)
library(forcats)
library(ggforce)
library(stringr)
library(rstatix)
library(gsw)


##### World map #####
# Marine data
metaM=read_xlsx("/Users/matev/Documents/Research/Chalmers/BlendARGs/Geotraces_filtered_metadata.xlsx") #for the coordinates
hgtsM=read.table("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/HGT_summary.txt", header = T) #simple summary file of HGT counts
sum(hgtsM$HGT_count)
gene_counts=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/gene_counts.csv")
hgtsM <- hgtsM %>%
  left_join(gene_counts, by = c("Sample" = "SampleID"))

hgtsM$HGT_activity=hgtsM$HGT_genes_unique/hgtsM$GeneCount*100

merged_tableM <- merge(hgtsM, metaM, by.x = "Sample", by.y = "group", all = F)
merged_tableM=merged_tableM[,c("Sample", "HGT_activity", "Latitude and Longitude")]
hgt_tableM=unique(merged_tableM)


split_coordinates <- strsplit(as.character(hgt_tableM$`Latitude and Longitude`), " ")

latitude <- numeric(length(split_coordinates))
latitude_direction <- character(length(split_coordinates))
longitude <- numeric(length(split_coordinates))
longitude_direction <- character(length(split_coordinates))

for (i in seq_along(split_coordinates)) {
  split_coordinate <- split_coordinates[[i]]
  latitude[i] <- as.numeric(split_coordinate[1])
  latitude_direction[i] <- split_coordinate[2]
  longitude[i] <- as.numeric(split_coordinate[3])
  longitude_direction[i] <- split_coordinate[4]
  
  # Convert latitude and longitude to negative if direction is South or West
  if (latitude_direction[i] == "S") {
    latitude[i] <- -latitude[i]
  }
  if (longitude_direction[i] == "W") {
    longitude[i] <- -longitude[i]
  }
}

hgt_tableM$lat <- latitude
hgt_tableM$long <- longitude

# Freshwater data
meta=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/meta_lakes.csv")
hgts=read.table("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/HGT_summary.txt", header = T)
sum(hgts$HGT_count)
gene_counts=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/gene_counts.csv")
hgts <- hgts %>%
  left_join(gene_counts, by = c("Sample" = "SampleID"))

hgts$HGT_activity=hgts$HGT_genes_unique/hgts$GeneCount*100

merged_table <- merge(hgts, meta, by.x = "Sample", by.y = "group", all = F)
merged_table=merged_table[,c("Sample", "HGT_activity", "geographic.location..longitude.", "geographic.location..latitude.", "environment..feature.")]
hgt_table=unique(merged_table)

names(hgt_table)[names(hgt_table) == "geographic.location..longitude."] = 'long'
names(hgt_table)[names(hgt_table) == "geographic.location..latitude."] = 'lat'

hgt_table$lat=as.numeric(hgt_table$lat)
hgt_table$long=as.numeric(hgt_table$long)


# Base map
base_map=basemap(c(-180,180,-70, 80), crs = 4326, legends = T, bathymetry = T, 
                 land.col = "#fffbf7", bathy.style = "rcb")


size_legend=range(c(hgt_table$HGT_activity, hgt_tableM$HGT_activity))


HGT_world_map=base_map+
  geom_point(data = hgt_table, 
             aes(x = long, y = lat, shape="Freshwater",
                 colour = ifelse(HGT_activity == 0, "0", "non-zero"),
                 size=HGT_activity)) +
  geom_point(data = hgt_tableM, 
             aes(x = long, y = lat, shape="Marine",
                 colour = ifelse(HGT_activity == 0, "0", "non-zero"), 
                 size=HGT_activity)) +
  scale_size_continuous(range = c(1, 5), limits = size_legend)+
  scale_color_manual(values=c("#f1c232", "#ec008a"), labels=c("absent", "present"))+
  scale_shape_manual(name="Environment", labels=c("Freshwater", "Marine"), values=c(16, 15))+
  guides(size = guide_legend(title = "HGT activity (%)", order = 3, title.position="top", title.hjust = 0.5,
                             override.aes = list(fill=NA)),
         shape = guide_legend(title = "Environment", order = 1, nrow=2, title.position="top", title.hjust = 0.5,
                              override.aes = list(size = 5, fill=NA)),
         color = guide_legend(title = "HGT across depths", order = 2, nrow=2, title.position="top", title.hjust = 0.5,
                              override.aes = list(size = 5, fill=NA)))+
  theme_linedraw(12)+
  theme(legend.position = "top",
        legend.direction = "vertical",
        legend.key = element_rect(fill='white', color='white'),
        plot.margin = unit(c(0,  # Top margin
                             0.5,  # Right margin
                             0,  # Bottom margin
                             0.5), "inches"),
        legend.box.margin = margin(0, 0, 0, 0))

HGT_world_map

base_map_lake=basemap(c(12,28,59,69), crs = 4326, legends = T, bathymetry = T, 
                      land.col = "#fffbf7")

zoom_map=base_map_lake+
  geom_point(data = hgt_table, 
             aes(x = long, y = lat, shape="Freshwater",
                 colour = ifelse(HGT_activity == 0, "0", "non-zero"),
                 size=HGT_activity), alpha=0.8, position = position_jitter(width = 0.4, height = 0.4)) +
  geom_point(data = hgt_tableM, 
             aes(x = long, y = lat, shape="Marine",
                 colour = ifelse(HGT_activity == 0, "0", "non-zero"), 
                 size=HGT_activity)) +
  scale_size_continuous(range = c(1, 5), limits = size_legend)+
  scale_color_manual(values=c("#f1c232", "#ec008a"), labels=c("absent", "present"))+
  scale_shape_manual(name="Environment", labels=c("Freshwater", "Marine"), values=c(16, 15))+
  guides(size = guide_legend(title = "HGT activity (%)", order = 3, title.position="top", title.hjust = 0.5,
                             override.aes = list(fill=NA)),
         shape = guide_legend(title = "Environment", order = 1, nrow=2, title.position="top", title.hjust = 0.5,
                              override.aes = list(size = 5, fill=NA)),
         color = guide_legend(title = "HGT across depths", order = 2, nrow=2, title.position="top", title.hjust = 0.5,
                              override.aes = list(size = 5, fill=NA)))+
  theme_linedraw(12)+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.direction = "vertical",
        legend.key = element_rect(fill='white', color='white'),
        plot.margin = unit(c(0,  # Top margin
                             0,  # Right margin
                             0,  # Bottom margin
                             0), "inches"),
        legend.box.margin = margin(0, 0, 0, 0))
zoom_map


HGT_world_map+inset_element(zoom_map, left = 0.5, bottom = 0.25, right = 0.95, top = 0.67, align_to = "full")

quartz.save("BlendARGs_worldmap.pdf", type="pdf", dpi=800, width = 9, height = 8)
quartz.save("BlendARGs_worldmap_small.pdf", type="pdf", dpi=800, width = 8, height = 7)

#Stat for HGT activity#
hgts$Source=c("Freshwater")
hgtsM$Source=c("Marine")
merged_HGT=rbind(hgts,hgtsM)
write.csv(merged_HGT, "/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/Figure1.csv")

shapiro.test(merged_HGT$HGT_activity[merged_HGT$Source=="Marine"])
shapiro.test(merged_HGT$HGT_activity[merged_HGT$Source=="Freshwater"])
#Not normally distributed based on Shapiro tests:
wilcox_test <- wilcox_test(merged_HGT, HGT_activity ~ Source)
wilcox_effect <- wilcox_effsize(merged_HGT, HGT_activity ~ Source)
wilcox_test
wilcox_effect

#within freshwater sites:
wilcox_test <- wilcox_test(hgt_table, HGT_activity ~ environment..feature.)
wilcox_effect <- wilcox_effsize(hgt_table, HGT_activity ~ environment..feature.)
wilcox_test
wilcox_effect





##### Stratification check #####
#### Marine 
hydrodata=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/stratificationcheck/merged_hydrodata.csv") #data collected from https://geotraces.webodv.awi.de/ 
head(metaM)
head(hydrodata)

metaM2 <- metaM %>%
  rename(
    Cruise = `GEOTRACES section`,
    Station = `Cruise Station`
  )


# Filter hydrodata by matching Cruise and Station with metaM
hydrodata <- hydrodata %>%
  mutate(
    # Clean Station values
    Station = str_remove(Station, "^Station\\s*"),
    Station = case_when(
      Station == "7019" ~ "19",     # keep as character
      TRUE ~ Station
    ),
    # Apply special replacement
    Station = case_when(
      Cruise == "GP13" & Station == "19" & Operator.s.Cruise.Name == "SS2011-1 - GP13 Leg1" ~ "GT19",
      TRUE ~ Station
    )
  ) %>%
  # Remove the other unwanted case
  filter(!(Cruise == "GA02" & Station == "14" & Operator.s.Cruise.Name == "PE319"))


hydrodata <- hydrodata %>%
  mutate(Station = as.character(Station))


# Filter by Cruise and Station
filtered_hydrodata <- hydrodata %>%
  semi_join(metaM2, by = c("Cruise", "Station"))


filtered_hydrodata <- filtered_hydrodata %>%
  mutate(
    Cruise = as.character(Cruise),
    Station = as.character(Station)
  ) %>%
  group_by(Cruise, Station) %>%
  arrange(`Operator.s.Cruise.Name`, .by_group = TRUE) %>%
  mutate(
    # Number distinct operator names encountered up to current row
    op_index = match(`Operator.s.Cruise.Name`, unique(`Operator.s.Cruise.Name`)),
    # Build group label: base + optional J
    group = paste0(Cruise, "_", Station, ifelse(op_index > 1, "J", ""))
  ) %>%
  ungroup() %>%
  select(-op_index)




#Calculating density for each sampling location

# Function to compute σ0 directly from T, S, depth
add_sigma0 <- function(df) {
  SA <- gsw_SA_from_SP(df$CTDSAL, df$DEPTH..m., 0, 0)   # use (0,0) as dummy
  CT <- gsw_CT_from_t(SA, df$CTDTMP..deg.C., df$DEPTH..m.)
  rho0 <- gsw_rho(SA, CT, 0)  # potential density (kg/m^3) referenced to 0 dbar
  df$sigma0 <- rho0 - 1000    # density anomaly
  df
}

# Function to compute MLD by density threshold
compute_mld <- function(depth, sigma0, d_sigma = 0.03) {
  # Remove missing values
  valid <- !is.na(depth) & !is.na(sigma0)
  depth <- depth[valid]
  sigma0 <- sigma0[valid]
  
  # If not enough points, return 0
  if (length(depth) < 2) return(0)
  
  # Order by depth
  ord <- order(depth)
  depth <- depth[ord]
  sigma0 <- sigma0[ord]
  
  # Shallowest depth as reference
  z_ref <- depth[1]
  sigma_ref <- sigma0[1]
  target <- sigma_ref + d_sigma
  
  # Find first index where sigma0 exceeds target
  idx <- which(sigma0 >= target)
  if (length(idx) == 0) return(0)  # fully mixed
  
  i <- idx[1]
  if (i == 1) return(z_ref)
  
  # Linear interpolation
  z0 <- depth[i-1]; z1 <- depth[i]
  s0 <- sigma0[i-1]; s1 <- sigma0[i]
  
  # Prevent division by zero
  if (is.na(s0) | is.na(s1) | s1 == s0) return(z1)
  
  frac <- (target - s0) / (s1 - s0)
  return(z0 + frac * (z1 - z0))
}


df_sigma <- filtered_hydrodata %>%
  group_by(group) %>%
  group_modify(~ add_sigma0(.x)) %>%
  ungroup()

mld_results <- df_sigma %>%
  group_by(group) %>%
  summarise(MLD = compute_mld(DEPTH..m., sigma0), .groups = "drop") %>%
  # Ensure all original groups are included
  right_join(
    filtered_hydrodata %>% distinct(group),
    by = "group"
  ) %>%
  mutate(MLD = ifelse(is.na(MLD), 0, MLD))

# Map MLD to hgtsM based on Sample -> group
hgtsM$MLD <- mld_results$MLD[match(hgtsM$Sample, mld_results$group)]


ggplot(hgtsM, aes(x = MLD, y = HGT_activity)) +
  geom_point(size = 3, alpha = 0.7) +  # points colored by Sample
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +  # regression line
  theme_linedraw(base_size = 12) +  # clean theme
  labs(
    title = "Relationship between HGT Activity and Mixed Layer Depth",
    x = "Mixed Layer Depth (MLD)",
    y = "HGT Activity",
    color = "Sample"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )
write.csv(mld_results, "/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/stratificationcheck/MLD_results_Marine.csv", row.names = FALSE)

min(mld_results$MLD)
max(mld_results$MLD)
median(mld_results$MLD)

hist(metaM$`Depth (m)`)
median(metaM$`Depth (m)`)
min(metaM$`Depth (m)`)
max(metaM$`Depth (m)`)

##### Taxonomic composition #####
tax_lake = read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/gtdbtk_summary_metachip.tsv", sep = "\t")
tax_m = read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/gtdbtk_summary_metachip.tsv", sep = "\t")

tax_lake <- tax_lake %>%
  separate(classification, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(
    Kingdom = gsub("d__", "", Kingdom),
    Phylum = gsub("p__", "", Phylum),
    Class = gsub("c__", "", Class),
    Order = gsub("o__", "", Order),
    Family = gsub("f__", "", Family),
    Genus = gsub("g__", "", Genus),
    Species = gsub("s__", "", Species)
  )
tax_lake$Phylum[is.na(tax_lake$Phylum) | tax_lake$Phylum == ""] <- "unclassified Bacteria"

tax_m <- tax_m %>%
  separate(classification, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(
    Kingdom = gsub("d__", "", Kingdom),
    Phylum = gsub("p__", "", Phylum),
    Class = gsub("c__", "", Class),
    Order = gsub("o__", "", Order),
    Family = gsub("f__", "", Family),
    Genus = gsub("g__", "", Genus),
    Species = gsub("s__", "", Species)
  )
tax_m$Phylum[is.na(tax_m$Phylum) | tax_m$Phylum == ""] <- "unclassified Bacteria"


lakes_tax <- tax_lake %>%
  mutate(group = sub("_.*", "", user_genome))
marine_tax <- tax_m %>%
  mutate(group = sub("_.*", "", user_genome))


lakes_tax_meta= merge(meta, lakes_tax, by.x="sample", by.y="group")
marine_tax_meta= merge(metaM, marine_tax, by.x="Sample_name", by.y="group")

## Loading sequencing depth data
bin_depth_L=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/bin_depthsummary_Freshwater.tsv", header = T, sep = '\t')
bin_depth_L$bin <- gsub("\\.", "_", gsub("^MEGAHIT-MetaBAT2-|\\.fa$", "", bin_depth_L$bin))
colnames(bin_depth_L) <- gsub("^Depth\\.", "", colnames(bin_depth_L))
colnames(bin_depth_L) <- gsub("\\.", "-", colnames(bin_depth_L))
bin_depth_L=bin_depth_L[,-1]
bin_depth_L=bin_depth_L[,-159:-209]
bin_depth_c_L <- bin_depth_L %>%
  pivot_longer(
    cols = -bin,   # All columns except 'bin' (assumed to be the identifier)
    names_to = "Sample",  # Convert column names (Depth samples) to a new column 'Sample'
    values_to = "Depth"    # Store corresponding values in 'Depth' column
  ) %>%
  filter(!is.na(Depth))

bin_depth_c_L <- bin_depth_c_L %>%
  mutate(
    bin_prefix = str_remove(bin, "(_[^_]+$|-[^_]+$)")  # Remove the last underscore or hyphen part
  ) %>%
  filter(bin_prefix == Sample)


lakes_tax_meta_depth <- bin_depth_c_L %>%
  right_join(lakes_tax_meta, by = c("bin" = "user_genome"))

lakes_tax_meta_depth <- lakes_tax_meta_depth %>%
  filter(!is.na(Depth)) 
range(lakes_tax_meta_depth$Depth)


bin_depth_M=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/bin_depthsummary_Marine.tsv", header = T, sep = '\t')
bin_depth_M$bin <- gsub("\\.", "_", gsub("^MEGAHIT-MetaBAT2-|\\.fa$", "", bin_depth_M$bin))
colnames(bin_depth_M) <- gsub("^Depth\\.", "", colnames(bin_depth_M))
colnames(bin_depth_M) <- gsub("\\.", "-", colnames(bin_depth_M))
bin_depth_M=bin_depth_M[,-1]
bin_depth_M=bin_depth_M[,-198:-248]
bin_depth_c_M <- bin_depth_M %>%
  pivot_longer(
    cols = -bin,   # All columns except 'bin' (assumed to be the identifier)
    names_to = "Sample",  # Convert column names (Depth samples) to a new column 'Sample'
    values_to = "Depth"    # Store corresponding values in 'Depth' column
  ) %>%
  filter(!is.na(Depth))

bin_depth_c_M <- bin_depth_c_M %>%
  filter(str_extract(bin, "^[^_]+") == Sample)

marine_tax_meta_depth <- bin_depth_c_M %>%
  right_join(marine_tax_meta, by = c("bin" = "user_genome"))



lakes_tax_meta_depth$Class[is.na(lakes_tax_meta_depth$Class) | lakes_tax_meta_depth$Class == ""] <- 
  paste("unclassified")
marine_tax_meta_depth$Class[is.na(marine_tax_meta_depth$Class) | marine_tax_meta_depth$Class == ""] <- 
  paste("unclassified")


int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}


size_legend_tax=range(c(lakes_tax_meta_depth$Depth, marine_tax_meta_depth$Depth))

lakes_tax_meta_depth$Phylum=factor(lakes_tax_meta_depth$Phylum)
marine_tax_meta_depth$Phylum=factor(marine_tax_meta_depth$Phylum)
# 1. Get a unified and sorted list of Class levels across both datasets
phylum_order <- sort(unique(c(levels(lakes_tax_meta_depth$Phylum), 
                              levels(marine_tax_meta_depth$Phylum))))
# 2. Sort Class levels within each Phylum
all_classes <- lakes_tax_meta_depth %>% 
  select(Phylum, Class) %>% 
  bind_rows(marine_tax_meta_depth %>% select(Phylum, Class)) %>%
  distinct() %>%
  arrange(match(Phylum, phylum_order), Class) %>%  # Sort by Phylum first, then Class
  pull(Class) %>%
  unique()
# 3. Convert Class factor levels in both datasets
lakes_tax_meta_depth$Class <- factor(lakes_tax_meta_depth$Class, levels = all_classes)
marine_tax_meta_depth$Class <- factor(marine_tax_meta_depth$Class, levels = all_classes)
# 4. Generate color mapping
all_class_colors <- viridis::viridis(length(all_classes), option = 'turbo', end = 0.9)
names(all_class_colors) <- all_classes

all_class_colors["unclassified"] <- "black"


lakes_tax_meta_depth$Phylum <- fct_relevel(lakes_tax_meta_depth$Phylum, "unclassified Bacteria", after = Inf)
lakes_tax_meta_depth$Class <- fct_relevel(lakes_tax_meta_depth$Class, "unclassified", after = Inf)

marine_tax_meta_depth$Phylum <- fct_relevel(marine_tax_meta_depth$Phylum, "unclassified Bacteria", after = Inf)
marine_tax_meta_depth$Class <- fct_relevel(marine_tax_meta_depth$Class, "unclassified", after = Inf)


set.seed(123)
Lakes_tax_plot=ggplot(lakes_tax_meta_depth, aes(x = geographic.location..depth., y=Phylum , size = Depth, fill=Class, color=Class)) +
  geom_point(alpha=0.8, shape=21, stroke=1, position = position_dodge(.6), show.legend = T)+
  theme_linedraw(base_size = 12) +
  coord_flip()+
  scale_x_reverse(breaks = int_breaks)+
  #facet_nested_wrap(~Tax, scales = "free", nrow = 8, ncol = 8)+
  labs(x = "Depth (m)", title = "Freshwater environments") +
  theme(legend.position = "right",
        plot.title = element_text(face='bold'),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        strip.text = element_text(face="bold", size=12, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='lightgrey'),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = all_class_colors, drop=F) +  # Apply modified colors
  scale_colour_manual(values = all_class_colors, drop=F) +  # Ensure outline color matches
  scale_size_continuous(range=c(1,10), breaks=c(1,10,100,1000, 2000), limits=size_legend_tax)+
  guides(size=guide_legend(title="Genome bin depth"), 
         colour='none',
         fill=guide_legend(override.aes = list(size=3), ncol=3, title='Class'))


Marine_tax_plot=ggplot(marine_tax_meta_depth, aes(x = `Depth (m)`, y=Phylum , size = Depth, fill=Class, color=Class)) +
  geom_point(alpha=0.8, shape=21, stroke=1, position = position_dodge(.6), show.legend = T)+
  theme_linedraw(base_size = 12) +
  coord_flip()+
  scale_x_reverse(breaks = int_breaks)+
  #facet_nested_wrap(~Tax, scales = "free", nrow = 8, ncol = 8)+
  labs(x = "Depth (m)", title = "Marine environments") +
  theme(legend.position = "right",
        plot.title = element_text(face='bold'),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        strip.text = element_text(face="bold", size=12, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='lightgrey'),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = all_class_colors, drop=F) +  # Apply modified colors
  scale_colour_manual(values = all_class_colors, drop=F) +  # Ensure outline color matches
  scale_size_continuous(range=c(1,10), breaks=c(1, 10,100,1000, 2000), limits=size_legend_tax)+
  guides(size=guide_legend(title="Genome bin depth"), 
         colour='none',
         fill=guide_legend(override.aes = list(size=3), ncol=3, title='Class'))

ggarrange(Lakes_tax_plot, Marine_tax_plot, common.legend = T, ncol = 2, widths = c(0.59, 0.41), legend = "right")
ggsave("Taxonomy_depth.pdf", dpi=600, width = 16, height = 8, bg='white', device = cairo_pdf)


#Lake Erken
erken_tax_meta_depth <- lakes_tax_meta_depth %>% filter(group == "Erken")
Erken_tax_plot=ggplot(erken_tax_meta_depth, aes(x = geographic.location..depth., y=Phylum , size = Depth, fill=Class, color=Class)) +
  geom_point(alpha=0.8, shape=21, stroke=1, position = position_dodge(.6), show.legend = T)+
  theme_linedraw(base_size = 12) +
  coord_flip()+
  scale_x_reverse(breaks = int_breaks)+
  labs(x = "Depth (m)", title = "Lake Erken") +
  theme(legend.position = "right",
        plot.title = element_text(face='bold'),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        strip.text = element_text(face="bold", size=12, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='lightgrey'),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = all_class_colors, drop=T) +  # Apply modified colors
  scale_colour_manual(values = all_class_colors, drop=T) +  # Ensure outline color matches
  scale_size_continuous(range=c(1,10), breaks=c(1,10,30,50))+
  guides(size=guide_legend(title="Genome bin depth"), 
         colour='none',
         fill=guide_legend(override.aes = list(size=3), ncol=2, title='Class'))

Erken_tax_plot
ggsave("Taxonomy_depth_Erken.pdf", dpi=600, width = 8, height = 6, bg='white', device = cairo_pdf)


##### Presence of stressors #####
dn_m=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/summary_DeepNOG_allbins_COG.csv", header = T)
dn_l=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/summary_DeepNOG_allbins_COG.csv", header = T)

dn_m <- dn_m %>% filter(Category_code == 'V')
dn_l <- dn_l %>% filter(Category_code == 'V')

dn_m <- dn_m %>%
  mutate(group = sub("_.*", "", ID))
dn_l <- dn_l %>%
  mutate(group = sub("_.*", "", ID))

dn_m=merge(dn_m, metaM, by.x="group", by.y='Sample_name')  
dn_l=merge(dn_l, meta, by.x="group", by.y='sample')  


legend_stress=union(sort(dn_m$group.y),sort(dn_l$group.y))

dn_l$group.y=factor(dn_l$group.y, levels = legend_stress)
dn_m$group.y=factor(dn_m$group.y, levels = legend_stress)


pond <- subset(dn_l, environment..feature. == "pond")
length(unique(pond$group.y))
lake <- subset(dn_l, environment..feature. == "lake")
length(unique(lake$group.y))

freshw_label=c("lake" = "lake (n=21)",
               "pond" = "pond (n=5)",
               "reservoir" = "reservoir (n=1)")

atlantic <- subset(dn_m, geo_loc_name == "Atlantic Ocean")
length(unique(atlantic$group.y))
pacific <- subset(dn_m, geo_loc_name == "Pacific Ocean")
length(unique(pacific$group.y))

marine_label=c("Atlantic Ocean" = "Atlantic Ocean (n=36)",
               "Pacific Ocean" = "Pacific Ocean (n=5)")

dn_depth_l=ggplot(dn_l, aes(x = geographic.location..depth., , fill = group.y)) +
  geom_density(aes(y=after_stat(count)),linewidth=0.1, alpha=0.5, show.legend = TRUE)+
  theme_linedraw(base_size = 12) +
  coord_flip()+
  facet_wrap(~environment..feature., scales = "free_y", labeller = as_labeller(freshw_label))+
  scale_x_reverse()+
  labs(x = "Depth (m)", y = "Genes mapped to defense mechanisms", title = "Freshwater environments",
       tag = expression(bold(a))) +
  theme(legend.position = "right",
        legend.key.size = unit(0.3, 'cm'),
        legend.key.spacing.y = unit(0.1, 'cm'),
        strip.text = element_text(face="bold", size=12),
        panel.grid = element_line(color='#f2f2f2')) +
  scale_fill_viridis_d(option = 'turbo', drop=F)+
  guides(fill=guide_legend(ncol=2, title="Sample location",by.row=T, override.aes = list(size=2)))


dn_depth_m=ggplot(dn_m, aes(x = `Depth (m)`, fill = group.y)) +
  geom_density(aes(y=after_stat(count)),linewidth=0.1, alpha=0.5, show.legend = TRUE)+
  theme_linedraw(base_size = 12) +
  coord_flip()+
  facet_wrap(~geo_loc_name, scales = "fixed", labeller = as_labeller(marine_label))+
  scale_x_reverse()+
  labs(x = "Depth (m)", y = "Genes mapped to defense mechanisms", title = "Marine environments",
       tag = expression(bold(b))) +
  theme(legend.position = "right",
        legend.key.size = unit(0.3, 'cm'),
        legend.key.spacing.y = unit(0.1, 'cm'),
        strip.text = element_text(face="bold", size=12),
        panel.grid = element_line(color='#f2f2f2')) +
  scale_fill_viridis_d(option = 'turbo', drop=F)+
  guides(fill=guide_legend(ncol=2, title="Sample location", by.row=T, override.aes = list(size=2)))


ggarrange(dn_depth_l, dn_depth_m, ncol = 2, common.legend = T, widths = c(0.59, 0.41), legend = "right")
quartz.save("Stress_Depth.pdf", type="pdf", dpi=800, width = 15, height = 6)

dn_l=dn_l[,c("group", "ID", "COG_ID", "Value","Category_code","Annotation","Category","group.y","geographic.location..depth.")]
dn_m=dn_m[,c("group", "ID", "COG_ID", "Value","Category_code","Annotation","Category","group.y","Depth (m)")]

names(dn_l) <- names(dn_m) 
write.csv(rbind(dn_l, dn_m),"/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/Figure_stress.csv")


####### HGT-mediated functions #####
hgtf_m=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/summary_DeepNOG_DorR_COG.csv", header = T)
hgtf_l=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/summary_DeepNOG_DorR_COG.csv", header = T)
hgtf_l <- hgtf_l[!duplicated(hgtf_l$ID), ]

hgtf_m <- hgtf_m %>%
  mutate(group = sub("_.*", "", ID))
hgtf_l <- hgtf_l %>%
  mutate(group = sub("_.*", "", ID))


hgtf_l <- hgtf_l %>%
  mutate(tax_bin = sub("_[^_]+$", "", ID))
hgtf_m <- hgtf_m %>%
  mutate(tax_bin = sub("_[^_]+$", "", ID))


hgtf_m=merge(hgtf_m, metaM, by.x="group", by.y='Sample_name')  
hgtf_l=merge(hgtf_l, meta, by.x="group", by.y='sample')  


hgtf_l_tax = merge(hgtf_l, lakes_tax_meta_depth, by.x="tax_bin", by.y="bin")
hgtf_m_tax = merge(hgtf_m, marine_tax_meta_depth, by.x="tax_bin", by.y="bin")

hgtf_l_tax=hgtf_l_tax[,-11]
counts_df_l <- hgtf_l_tax %>%
  group_by(Category_code, environment..feature..x, geographic.location..depth..x, DorR, Phylum, Class) %>%
  summarise(count = n()) %>%
  ungroup()

hgtf_m_tax=hgtf_m_tax[,-20]
counts_df_m <- hgtf_m_tax %>%
  group_by(Category_code, geo_loc_name.x, `Depth (m).x`, DorR, Phylum, Class) %>%
  summarise(count = n()) %>%
  ungroup()


gene_label=c(donor='HGT donor', recipient="HGT recipient")

colnames(counts_df_l)=c("Category", "Env", "Depth", "DorR", "Phylum", "Class", "count")
colnames(counts_df_m)=c("Category", "Env", "Depth", "DorR", "Phylum", "Class", "count")

counts_df_l <- counts_df_l[counts_df_l$Category != "" & !is.na(counts_df_l$Category), ]
counts_df_m <- counts_df_m[counts_df_m$Category != "" & !is.na(counts_df_m$Category), ]



category_order <- c(
  "J", "A", "K", "L", # INFORMATION STORAGE AND PROCESSING
  "D", "V", "T", "M", "N", "Z", "W", "U", "O", # CELLULAR PROCESSES AND SIGNALING
  "C", "G", "E", "F", "H", "I", "P", "Q", # METABOLISM
  "R", "S" # POORLY CHARACTERIZED
)

counts_df_l$Category <- fct_relevel(counts_df_l$Category, category_order)
counts_df_m$Category <- fct_relevel(counts_df_m$Category, category_order)

size_legend_hgt=range(c(counts_df_l$count, counts_df_m$count))

counts_df_l$Class <- factor(counts_df_l$Class, levels = all_classes)
counts_df_m$Class <- factor(counts_df_m$Class, levels = all_classes)

counts_df_l$Class <- fct_relevel(counts_df_l$Class, "unclassified", after = Inf)
counts_df_m$Class <- fct_relevel(counts_df_m$Class, "unclassified", after = Inf)




hgtf_depth_l=ggplot(counts_df_l, aes(x = Depth, y=Category , size = count, fill=Class, colour = Class)) +
  geom_point(data = transform(counts_df_l, DorR = NULL), colour = "#cccccc") +
  geom_point(alpha=0.8, shape=21, stroke=0.5, position =position_dodge(.5), show.legend = T)+
  theme_linedraw(base_size = 12) +
  coord_flip()+
  scale_x_reverse(breaks = int_breaks)+
  scale_y_discrete(labels = label_wrap_gen(10))+
  facet_nested_wrap(DorR~Env, scales = "free_y", nrow = 2, ncol = 3, labeller=labeller(DorR=gene_label),
                    nest_line = element_line(colour="white", linetype = "dotted"))+
  labs(x = "Depth (m)", title = "Freshwater environments", tag = expression(bold(a))) +
  theme(legend.position = "none",
        plot.title = element_text(face="bold"),
        legend.spacing = unit(0.01, 'cm'),
        legend.key.size = unit(0.5, 'cm'),
        axis.title.x = element_blank(),
        strip.text = element_text(face="bold", size=12, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(1,1,1,12, unit = "pt"))+
  scale_size(range = c(1,10), breaks = c(1, 5, 10), limits = size_legend_hgt)+
  scale_fill_manual(values = all_class_colors, drop=F) +  
  scale_colour_manual(values = all_class_colors, drop=F) +
  guides(size=guide_legend(title="# HGT-mediated\ngenes"), colour='none',
         fill=guide_legend(title='Class', override.aes = list(size=4), ncol=3))

hgtf_depth_m=ggplot(counts_df_m, aes(x = Depth, y=Category , size = count, fill=Class, colour=Class)) +
  geom_point(data = transform(counts_df_m, DorR = NULL), colour = "#cccccc") +
  geom_point(alpha=0.8, shape=21, stroke=0.5, position =position_dodge(.5), show.legend = T)+
  theme_linedraw(base_size = 12) +
  coord_flip()+
  facet_nested_wrap(DorR~Env, scales = "free_y", nrow = 2, ncol = 2, labeller = labeller(DorR=gene_label),
                    nest_line = element_line(colour="white", linetype = "dotted"))+
  scale_x_reverse(breaks = int_breaks)+
  scale_y_discrete(labels = label_wrap_gen(10))+
  labs(x = "Depth (m)", title = "Marine environments", tag = expression(bold(b))) +
  theme(legend.position = "right",
        plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        legend.key.spacing = unit(0, 'cm'),
        legend.spacing = unit(0.01, 'cm'),
        legend.key.size = unit(0.5, 'cm'),
        strip.text = element_text(face="bold", size=12, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='#f2f2f2'))+
  scale_size(range = c(1,10), breaks = c(1, 5, 10), limits = size_legend_hgt)+
  scale_fill_manual(values = all_class_colors, drop=F) + 
  scale_colour_manual(values = all_class_colors, drop=F) +
  guides(size=guide_legend(title="# HGT-mediated\ngenes", order = 1), colour='none',
         fill=guide_legend(title='Class', override.aes = list(size=3), ncol=4))


blankplot=ggplot() + theme_void()
upper_plot=ggarrange(hgtf_depth_l, blankplot, common.legend = F, ncol = 2, widths = c(0.8, 0.2), align = "v")
bottom_plot=ggarrange(hgtf_depth_m)
ggarrange(upper_plot, bottom_plot, ncol = 1, common.legend = F, heights = c(0.495, 0.505), legend = "right")

ggsave("HGT_function_Depth.pdf", dpi=600, width = 15, height = 14, bg='white', device = cairo_pdf)

write.csv(rbind(counts_df_l, counts_df_m),"/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/Figure2.csv")



## Stats and summary
# Top 3 donor+recipient groups (collapsed at Phylum/Class level)
top_taxa_L <- counts_df_l %>%
  group_by(Phylum, Class) %>%
  summarise(total_HGT_genes = sum(count), .groups="drop") %>%
  arrange(desc(total_HGT_genes)) %>%
  slice_head(n = 3)

top_taxa_L

top_taxa_M <- counts_df_m %>%
  group_by(Phylum, Class) %>%
  summarise(total_HGT_genes = sum(count), .groups="drop") %>%
  arrange(desc(total_HGT_genes)) %>%
  slice_head(n = 3)

top_taxa_M


#Differences between donor vs recipient among taxonomic groups
agg_L <- counts_df_l %>%
  group_by(Class, DorR) %>%
  summarise(total_HGT_genes = sum(count), .groups = "drop")

table_fisher_L <- agg_L %>%
  pivot_wider(names_from = DorR, values_from = total_HGT_genes, values_fill = 0)

table_matrix_L <- as.matrix(table_fisher_L[,-1])

fisher_result_L <- fisher.test(table_matrix_L, simulate.p.value = T, B=999)
fisher_result_L

agg_L %>%
  group_by(DorR) %>%
  arrange(desc(total_HGT_genes)) %>%
  slice_head(n = 3)

agg_M <- counts_df_m %>%
  group_by(Class, DorR) %>%
  summarise(total_HGT_genes = sum(count), .groups = "drop")

table_fisher_M <- agg_M %>%
  pivot_wider(names_from = DorR, values_from = total_HGT_genes, values_fill = 0)

table_matrix_M <- as.matrix(table_fisher_M[,-1])

fisher_result_M <- fisher.test(table_matrix_M, simulate.p.value = T, B=999)
fisher_result_M

agg_M %>%
  group_by(DorR) %>%
  arrange(desc(total_HGT_genes)) %>%
  slice_head(n = 3)


#Distribution of HGT across COG categories
summarise_enrichment <- function(df) {
  # top 3 by raw counts
  top_cats <- df %>%
    group_by(Category) %>%
    summarise(total_HGT_genes = sum(count), .groups="drop") %>%
    arrange(desc(total_HGT_genes)) %>%
    slice_head(n = 3)
  
  # chi-sq test for uniform distribution
  cat_counts <- table(df$Category)
  chisq_res <- chisq.test(cat_counts)
  
  # residuals -> enrichment status
  enrichment <- data.frame(Category = names(chisq_res$residuals),
                           Residual = as.numeric(chisq_res$residuals)) %>%
    mutate(Status = case_when(
      Residual >  2 ~ "Enriched",
      Residual < -2 ~ "Depleted",
      TRUE          ~ "NS"
    ))
  
  # join counts + residuals
  summary <- df %>%
    group_by(Category) %>%
    summarise(total_HGT_genes = sum(count), .groups="drop") %>%
    left_join(enrichment, by = "Category") %>%
    arrange(desc(total_HGT_genes))
  
  list(
    top_categories = top_cats,
    chisq_summary  = summary,
    chisq_stat     = chisq_res$statistic,
    chisq_df       = chisq_res$parameter,
    chisq_p        = chisq_res$p.value
  )
}


res_L <- summarise_enrichment(counts_df_l)
res_M <- summarise_enrichment(counts_df_m)


res_L$top_categories
print(res_L$chisq_summary, n=24)
res_L$chisq_stat
res_L$chisq_df
res_L$chisq_p

res_M$top_categories
res_M$chisq_summary
res_M$chisq_stat
res_M$chisq_df
res_M$chisq_p


res_L$chisq_summary$Environment = "Freshwater"
res_M$chisq_summary$Environment = "Marine"
write.csv(rbind(res_L$chisq_summary, res_M$chisq_summary),"/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/COG_HGT_enrichment.csv")


######## Relationship between sampling depths and HGT activity #######
calculate_median_distance <- function(df, column) {
  if (nrow(df) > 1) {
    distances <- dist(df[[column]], method = "euclidean")
    return(median(as.numeric(distances)))
  } else {
    return(NA)  # No pairwise distances if only one point
  }
}

median_distance_L <- hgtf_l %>%
  group_by(group.y) %>%
  summarise(median_distance = calculate_median_distance(cur_data(), "geographic.location..depth."))

print(median_distance_L)

median_distance_M <- hgtf_m %>%
  group_by(group.y) %>%
  summarise(median_distance = calculate_median_distance(cur_data(), "Depth (m)"))

print(median_distance_M)

median_distance_L_c=merge(median_distance_L, hgt_table, by.x = "group.y", by.y = "Sample")
median_distance_M_c=merge(median_distance_M, hgt_tableM, by.x = "group.y", by.y = "Sample")


combined_median_distances <- bind_rows(
  mutate(median_distance_L_c, Source = "L"),
  mutate(median_distance_M_c, Source = "M")
)

combined_median_distances$group.y=factor(combined_median_distances$group.y, levels = legend_stress)

hgt_depth_link=ggplot(combined_median_distances, aes(x = median_distance, y=HGT_activity , colour=group.y, fill=group.y, shape=Source)) +
  geom_point(alpha=0.6, stroke=0.5, show.legend = T, size=5)+
  theme_linedraw(base_size = 12) +
  #scale_x_reverse(breaks = int_breaks)+
  #scale_y_discrete(labels = label_wrap_gen(35))+
  labs(x = "Median distance (m)", y="HGT activity (%)") +
  theme(legend.justification.inside = c(0.9,0.9),
        legend.location = "plot",
        legend.spacing = unit(0.01, 'cm'),
        legend.key.size = unit(0.5, 'cm'),
        panel.grid = element_line(color='#f2f2f2'))+
  scale_colour_viridis_d(option = 'turbo', end = 0.9, drop=F)+
  scale_fill_viridis_d(option = 'turbo', end = 0.9, drop=F)+
  scale_shape_manual(labels=c("L"="Freshwater", "M"="Marine"), values = c("L"=19, "M"=15))+
  guides(colour="none", fill="none", shape=guide_legend(title = "", position = "inside"))

hgt_depth_link
quartz.save("Link_HGT_Depth.pdf", type="pdf", dpi=600, width = 5, height = 5, bg='white')

write.csv(combined_median_distances,"/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/Figure_linkHGT_Depth.csv")


##### Resistome across depth #####
### Lakes ###
#RGI
rgi_lakes_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/summary_RGI_allbins_ARG.csv")
rgi_lakes_all <- rgi_lakes_all %>%
  mutate(group = sub("_.*", "", ORF_ID))

rgi_lakes_all=merge(rgi_lakes_all, meta, by.x="group", by.y='sample')  
rgi_lakes_all$Method="CARD"
rgi_lakes_all <- rgi_lakes_all %>%
  mutate(tax_bin = sub("_[^_]+$", "", ORF_ID))

#ResFinderFG2.0
#rf_lakes_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/summary_BLAST_allbins_RF.csv")
#RF_code=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/ResFinderFG_Code.csv")

#rf_lakes_all=merge(rf_lakes_all, RF_code, by.x='AB', by.y='code')
#write.csv(rf_lakes_all, "/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/BLAST_allbins_RF_ann.csv")
rf_lakes_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/BLAST_allbins_RF_ann.csv")

rf_lakes_all <- rf_lakes_all %>%
  mutate(group = sub("_.*", "", query))

rf_lakes_all=merge(rf_lakes_all, meta, by.x="group", by.y='sample')  
rf_lakes_all$Method="ResFinderFG2.0"
rf_lakes_all <- rf_lakes_all %>%
  mutate(tax_bin = sub("_[^_]+$", "", query))

#LatentARGs
#latent_lakes_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/summary_BLAST_allbins_latent.csv")
#latent_cluster=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Latent_cluster_ARG.csv")
#latent_seq=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/LatentARGs_seq.csv")
#latent_db=merge(latent_cluster, latent_seq, by.x='cluster', by.y='Cluster')

#latent_lakes_all=merge(latent_lakes_all, latent_db, by.x='target', by.y='Variant')
#write.csv(latent_lakes_all, "/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/BLAST_allbins_latent_ann.csv")
latent_lakes_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/BLAST_allbins_latent_ann.csv")

latent_lakes_all <- latent_lakes_all %>%
  mutate(group = sub("_.*", "", query))

latent_lakes_all=merge(latent_lakes_all, meta, by.x="group", by.y='sample')  
latent_lakes_all$Method="Latent ARGs"
latent_lakes_all <- latent_lakes_all %>%
  mutate(tax_bin = sub("_[^_]+$", "", query))

rgi_lakes_all <- rgi_lakes_all %>%
  separate_rows(Drug.Class, sep = "; ")
rf_lakes_all <- rf_lakes_all %>%
  separate_rows(Class, sep = ", ")
latent_lakes_all <- latent_lakes_all %>%
  separate_rows(Class.1, sep = ", ")


rgi_lakes_all=rgi_lakes_all[,c("ORF_ID", "tax_bin", "Drug.Class", "Method", "ID")]
rgi_lakes_all= merge(rgi_lakes_all, lakes_tax_meta_depth, by.x="tax_bin", by.y="bin")
rf_lakes_all= merge(rf_lakes_all, lakes_tax_meta_depth, by.x="tax_bin", by.y="bin")
latent_lakes_all= merge(latent_lakes_all, lakes_tax_meta_depth, by.x="tax_bin", by.y="bin")


rgi_lakes_all=rgi_lakes_all[,c("ORF_ID", "Drug.Class", "environment..feature.", "geographic.location..depth.", 
                               "Method", "Phylum", "Class", "ID", "Depth")]

rf_lakes_all=rf_lakes_all[,c("query", "Class.x", "environment..feature..x", "geographic.location..depth..x", 
                               "Method", "Phylum", "Class.y", "ID", "Depth")]

latent_lakes_all=latent_lakes_all[,c("query", "Class.1", "environment..feature..x", "geographic.location..depth..x", 
                             "Method", "Phylum", "Class.y", "target", "Depth")]


colnames(rgi_lakes_all)=c("ORF_ID","Drug_Class", "Env", "Depth", "Method","Phylum", "Class", "ARG", "Bin_depth")
colnames(rf_lakes_all)=c("ORF_ID","Drug_Class", "Env", "Depth", "Method","Phylum", "Class", "ARG", "Bin_depth")
colnames(latent_lakes_all)=c("ORF_ID","Drug_Class", "Env", "Depth", "Method","Phylum", "Class", "ARG", "Bin_depth")

Lakes_ARG=rbind(rgi_lakes_all, rf_lakes_all, latent_lakes_all)


### Marine ###
#RGI
rgi_m_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/summary_RGI_allbins_ARG.csv")
rgi_m_all <- rgi_m_all %>%
  mutate(group = sub("_.*", "", ORF_ID))

rgi_m_all=merge(rgi_m_all, metaM, by.x="group", by.y='Sample_name')  
rgi_m_all$Method="CARD"
rgi_m_all <- rgi_m_all %>%
  mutate(tax_bin = sub("_[^_]+$", "", ORF_ID))

#ResFinderFG2.0
#rf_m_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/summary_BLAST_allbins_RF.csv")
#rf_m_all=merge(rf_m_all, RF_code, by.x='AB', by.y='code')
#write.csv(rf_m_all, "/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/BLAST_allbins_RF_ann.csv")
rf_m_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/BLAST_allbins_RF_ann.csv")

rf_m_all <- rf_m_all %>%
  mutate(group = sub("_.*", "", query))

rf_m_all=merge(rf_m_all, metaM, by.x="group", by.y='Sample_name')
rf_m_all$Method="ResFinderFG2.0"
rf_m_all <- rf_m_all %>%
  mutate(tax_bin = sub("_[^_]+$", "", query))

#Latent ARGs
#latent_m_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/summary_BLAST_allbins_latent.csv")
#latent_m_all=merge(latent_m_all, latent_db, by.x='target', by.y='Variant')
#write.csv(latent_m_all, "/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/BLAST_allbins_latent_ann.csv")
latent_m_all=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/BLAST_allbins_latent_ann.csv")

latent_m_all <- latent_m_all %>%
  mutate(group = sub("_.*", "", query))

latent_m_all=merge(latent_m_all, metaM, by.x="group", by.y='Sample_name')  
latent_m_all$Method="Latent ARGs"
latent_m_all <- latent_m_all %>%
  mutate(tax_bin = sub("_[^_]+$", "", query))


rgi_m_all <- rgi_m_all %>%
  separate_rows(Drug.Class, sep = "; ")
rf_m_all <- rf_m_all %>%
  separate_rows(Class, sep = ", ")
latent_m_all <- latent_m_all %>%
  separate_rows(Class.1, sep = ", ")

rgi_m_all= merge(rgi_m_all, marine_tax_meta_depth, by.x="tax_bin", by.y="bin")
rf_m_all= merge(rf_m_all, marine_tax_meta_depth, by.x="tax_bin", by.y="bin")
latent_m_all= merge(latent_m_all, marine_tax_meta_depth, by.x="tax_bin", by.y="bin")


rgi_m_all=rgi_m_all[,c("ORF_ID", "Drug.Class", "geo_loc_name.x", "Depth (m).x",
                               "Method", "Phylum", "Class", "ID", "Depth")]

rf_m_all=rf_m_all[,c("query","Class.x", "geo_loc_name.x", "Depth (m).x", 
                             "Method", "Phylum", "Class.y", "ID", "Depth")]

latent_m_all=latent_m_all[,c("query","Class.1", "geo_loc_name.x", "Depth (m).x", 
                                     "Method", "Phylum", "Class.y", "target", "Depth")]



colnames(rgi_m_all)=c("ORF_ID","Drug_Class", "Env", "Depth", "Method","Phylum", "Class", "ARG", "Bin_depth")
colnames(rf_m_all)=c("ORF_ID","Drug_Class", "Env", "Depth", "Method","Phylum", "Class", "ARG", "Bin_depth")
colnames(latent_m_all)=c("ORF_ID","Drug_Class", "Env", "Depth", "Method","Phylum", "Class", "ARG", "Bin_depth")


Marine_ARG=rbind(rgi_m_all, rf_m_all, latent_m_all)


Lakes_ARG$Drug_Class[Lakes_ARG$Drug_Class=="disinfecting agents and antiseptics"] <- "disinfect. & antiseptics"
Marine_ARG$Drug_Class[Marine_ARG$Drug_Class=="disinfecting agents and antiseptics"] <- "disinfect. & antiseptics"



Lakes_ARG$Class <- factor(Lakes_ARG$Class, levels = all_classes)
Marine_ARG$Class <- factor(Marine_ARG$Class, levels = all_classes)

Marine_ARG$Class <- fct_relevel(Marine_ARG$Class, "unclassified", after = Inf)
Lakes_ARG$Class <- fct_relevel(Lakes_ARG$Class, "unclassified", after = Inf)


size_legend=range(c(Lakes_ARG$Bin_depth, Marine_ARG$Bin_depth))

Lakes_ARG$Method=factor(Lakes_ARG$Method, levels = c("CARD", "ResFinderFG2.0", "Latent ARGs"))
Marine_ARG$Method=factor(Marine_ARG$Method, levels = c("CARD", "ResFinderFG2.0", "Latent ARGs"))


set.seed(123)

Lakes_ARG_plot=ggplot(Lakes_ARG, aes(x = Depth, y=Drug_Class, size = Bin_depth, fill=Class, colour = Class)) +
  geom_point(data = transform(Lakes_ARG, Method = NULL), colour = "#cccccc", position =position_dodge(.5)) +
  geom_point(alpha=0.6, shape=21, stroke=0.5, position = position_dodge(.5), show.legend = T)+
  theme_linedraw(base_size = 10) +
  coord_flip()+
  scale_x_reverse(breaks = int_breaks)+
  scale_y_discrete(labels = label_wrap_gen(25))+
  facet_nested_wrap(Method~Env, scales = "free_y", nrow = 3, ncol = 3,
                    nest_line = element_line(colour="white", linetype = "dotted", linewidth = 0.5))+
  labs(x = "Depth (m)", title = "Freshwater environments", tag = expression(bold(a))) +
  theme(legend.position = "none",
        plot.title = element_text(face='bold'),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        strip.text = element_text(face="bold", size=10, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(1,1,1,10, unit = "pt"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = all_class_colors, drop=F) +  
  scale_colour_manual(values = all_class_colors, drop=F) +  
  scale_size_continuous(range=c(1,10), breaks=c(1,10,100,1000, 2000), limits=size_legend)+
  guides(size=guide_legend(title="Genome bin depth"), colour='none',
         fill=guide_legend(override.aes = list(size=4), ncol=4, title='Class'))


Marine_ARG_plot=ggplot(Marine_ARG, aes(x = Depth, y=Drug_Class , size = Bin_depth, fill=Class, colour = Class)) +
  geom_point(data = transform(Marine_ARG, Method = NULL), colour = "#cccccc", position =position_dodge(.5)) +
  geom_point(alpha=0.6, shape=21, stroke=0.5, position =position_dodge(.5), show.legend = T)+
  theme_linedraw(base_size = 10) +
  coord_flip()+
  scale_x_reverse(breaks = int_breaks)+
  scale_y_discrete(labels = label_wrap_gen(25))+
  facet_nested_wrap(Method~Env, scales = "free_y", nrow = 3, ncol = 2,
                    nest_line = element_line(colour="white", linetype = "dotted", linewidth = 0.5))+
  labs(x = "Depth (m)", title = "Marine environments", tag = expression(bold(b))) +
  theme(legend.position = "right",
        legend.justification = "right",
        legend.key.spacing = unit(0, 'cm'),
        plot.title = element_text(face='bold'),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing = unit(0.01, 'cm'),
        strip.text = element_text(face="bold", size=10, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='#f2f2f2'),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = all_class_colors, drop=F) +  
  scale_colour_manual(values = all_class_colors, drop=F) +
  scale_size_continuous(range=c(1,10), breaks=c(1,10,100,1000, 2000), limits=size_legend)+
  guides(size=guide_legend(title="Genome bin depth"), colour='none',
         fill=guide_legend(title='Class', override.aes = list(size=3), ncol=4))


blankplot=ggplot() + theme_void()
upper_plot=ggarrange(Lakes_ARG_plot, blankplot, common.legend = F, ncol = 2, widths = c(0.8, 0.2), align = "v")
ggarrange(upper_plot, Marine_ARG_plot, ncol = 1, common.legend = F, heights = c(0.495, 0.505), legend = "right")

ggsave("allARGs_Depth.pdf", dpi=600, width = 12, height = 13, bg='white', device = cairo_pdf)


### This is just to append the list of ARGs with the proper gene names
### Load the input files again before running this script below:
Marine_ARG_c <- Marine_ARG %>% mutate(ARG = as.character(ARG))
Lakes_ARG_c <- Lakes_ARG %>% mutate(ARG = as.character(ARG))

map_rgi <- bind_rows(
  rgi_m_all %>% select(ID, ARG_ID = Best_Hit_ARO),
  rgi_lakes_all %>% select(ID, ARG_ID = Best_Hit_ARO)
) %>% filter(!is.na(ID)) %>% mutate(ID = as.character(ID)) %>%
  distinct(ID, .keep_all = TRUE)

map_rf <- bind_rows(
  rf_m_all %>% select(ID, ARG_ID = ARG),
  rf_lakes_all %>% select(ID, ARG_ID = ARG)
) %>% filter(!is.na(ID)) %>% mutate(ID = as.character(ID)) %>%
  distinct(ID, .keep_all = TRUE)

map_latent <- bind_rows(
  latent_m_all %>% select(ID = target, ARG_ID = ARG),
  latent_lakes_all %>% select(ID = target, ARG_ID = ARG)
) %>% filter(!is.na(ID)) %>% mutate(ID = as.character(ID)) %>%
  distinct(ID, .keep_all = TRUE)

lookup <- bind_rows(map_rgi, map_rf, map_latent)

dups <- lookup %>% group_by(ID) %>% filter(n() > 1) %>% distinct(ID) %>% pull(ID)
if(length(dups) > 0) {
  stop("Found ID(s) present in multiple reference sources: ",
       paste(head(dups, 20), collapse = ", "),
       ". Each ARG should exist in only one reference df. Inspect these IDs.")
}

# safe one-to-one join
full_ARG_M <- Marine_ARG_c %>% left_join(lookup, by = c("ARG" = "ID"))
full_ARG_L <- Lakes_ARG_c %>% left_join(lookup, by = c("ARG" = "ID"))

write.csv(rbind(full_ARG_M, full_ARG_L),"/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/ARG_full.csv")




# Summary table data (unique ARGs per dataset):
length(unique(rgi_lakes_all$ARG))
length(unique(rgi_m_all$ARG))
length(unique(rf_lakes_all$ARG))
length(unique(rf_m_all$ARG))
length(unique(latent_lakes_all$ARG))
length(unique(latent_m_all$ARG))

#Checking whether the same ORF was mapped to genes across databases
arg_multiple_methods_M <- Marine_ARG %>%
  group_by(ORF_ID) %>%
  summarise(n_methods = n_distinct(Method)) %>%
  filter(n_methods > 1)

arg_multiple_d_methods_M <- Marine_ARG %>%
  semi_join(arg_multiple_methods_M, by = "ORF_ID")

write.csv(arg_multiple_d_methods_M, "/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/Marine_ARG_overlapmethods.csv")


arg_multiple_methods_L <- Lakes_ARG %>%
  group_by(ORF_ID) %>%
  summarise(n_methods = n_distinct(Method)) %>%
  filter(n_methods > 1)

arg_multiple_d_methods_L <- Lakes_ARG %>% #empty df = no overlapping ARGs in Freshwater data
  semi_join(arg_multiple_methods_L, by = "ORF_ID")


## Stats and summary
# Unique ARGs per environment
unique_args_L <- Lakes_ARG %>%
  group_by(Env) %>%
  summarise(n_unique_ARGs = n_distinct(ARG),
            n_unique_ORFs = n_distinct(ORF_ID),
            .groups = "drop")

unique_args_L

unique_args_M <- Marine_ARG %>%
  group_by(Env) %>%
  summarise(n_unique_ARGs = n_distinct(ARG),
            n_unique_ORFs = n_distinct(ORF_ID),
            .groups = "drop")

unique_args_M

#most common drug classes
drug_distribution_L <- Lakes_ARG %>%
  group_by(Env, Drug_Class) %>%
  summarise(n_unique_ARGs = n_distinct(ARG), 
            .groups = "drop") %>%
  arrange(Env, desc(n_unique_ARGs))

drug_distribution_L

drug_distribution_M <- Marine_ARG %>%
  group_by(Env, Drug_Class) %>%
  summarise(n_unique_ARGs = n_distinct(ARG), 
            .groups = "drop") %>%
  arrange(Env, desc(n_unique_ARGs))

drug_distribution_M

#taxonomic groups harboring ARGs
taxonomic_dist_L <- Lakes_ARG %>%
  group_by(Env, Phylum, Class) %>%
  summarise(n_unique_ARGs = n_distinct(ARG), .groups = "drop") %>%
  arrange(Env, desc(n_unique_ARGs))

taxonomic_dist_L

taxonomic_dist_M <- Marine_ARG %>%
  group_by(Env, Phylum, Class) %>%
  summarise(n_unique_ARGs = n_distinct(ARG), .groups = "drop") %>%
  arrange(Env, desc(n_unique_ARGs))


args_by_class_L <- Lakes_ARG %>%
  group_by(ARG, Drug_Class) %>%
  summarise(n_classes = n_distinct(Class),
            n_phyla = n_distinct(Phylum),
            .groups = "drop") %>%
  filter(n_classes > 1 | n_phyla > 1) %>%
  arrange(desc(n_classes))

args_by_class_M <- Marine_ARG %>%
  group_by(ARG, Drug_Class) %>%
  summarise(n_classes = n_distinct(Class),
            n_phyla = n_distinct(Phylum),
            .groups = "drop") %>%
  filter(n_classes > 1 | n_phyla > 1) %>%
  arrange(desc(n_classes))


###### Virus distribution and taxonomy ad lifestyle #######
vir_M=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/virus_contig_summary.csv", header = T)
vir_L=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/virus_contig_summary.csv", header = T)


vir_M <- vir_M %>%
  separate(taxonomy, into = c("Unranked", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";")
vir_M$Phylum[is.na(vir_M$Phylum) | vir_M$Phylum == ""] <- "unclassified virus"
vir_M$Class[is.na(vir_M$Class) | vir_M$Class == ""] <- "unclassified virus"


vir_L <- vir_L %>%
  separate(taxonomy, into = c("Unranked", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";")
vir_L$Phylum[is.na(vir_L$Phylum) | vir_L$Phylum == ""] <- "unclassified virus"
vir_L$Class[is.na(vir_L$Class) | vir_L$Class == ""] <- "unclassified virus"

#remove too short contigs (<5kb)
vir_M <- vir_M %>% filter(length >= 3000)
vir_L <- vir_L %>% filter(length >= 3000)

vir_M_meta=merge(vir_M, metaM, by.x="Sample", by.y='Sample_name')  
vir_L_meta=merge(vir_L, meta, by.x="Sample", by.y='sample_id')  

counts_vir_M <- vir_M_meta %>%
  group_by(group, geo_loc_name, `Depth (m)`, Phylum, Class) %>%
  summarise(count = n()) %>%
  ungroup()

counts_vir_L <- vir_L_meta %>%
  group_by(group, environment..feature., geographic.location..depth., Phylum, Class) %>%
  summarise(count = n()) %>%
  ungroup()


colnames(counts_vir_M)=c("Group", "Env", "Depth", "Phylum", "Class", "count")
colnames(counts_vir_L)=c("Group", "Env", "Depth", "Phylum", "Class", "count")


legend_tax_vir=sort(union(counts_vir_L$Class, counts_vir_M$Class))
size_legend_vir=range(c(counts_vir_L$count, counts_vir_M$count))

counts_vir_L$Class=factor(counts_vir_L$Class, levels = legend_tax_vir)
counts_vir_M$Class=factor(counts_vir_M$Class, levels = legend_tax_vir)


counts_vir_M$Phylum <- fct_relevel(counts_vir_M$Phylum, "unclassified virus", after = Inf)
counts_vir_M$Class <- fct_relevel(counts_vir_M$Class, "unclassified virus", after = Inf)

counts_vir_L$Phylum <- fct_relevel(counts_vir_L$Phylum, "unclassified virus", after = Inf)
counts_vir_L$Class <- fct_relevel(counts_vir_L$Class, "unclassified virus", after = Inf)


Marine_vir_plot=ggplot(counts_vir_M, aes(x = Depth, y=Phylum, size = count, fill=Class, colour = Class)) +
  geom_point(alpha=0.8, shape=21, stroke=0.5, position =position_dodge(.5), show.legend = T)+
  theme_linedraw(base_size = 10) +
  coord_flip()+
  scale_x_reverse(breaks = int_breaks)+
  scale_y_discrete(labels = label_wrap_gen(18))+
  facet_nested_wrap(~Env, scales = "free_y", nrow = 3, ncol = 2, labeller = as_labeller(marine_label))+
  labs(x = "Depth (m)", title = "Marine environments", tag = expression(bold(b))) +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing = unit(0.01, 'cm'),
        strip.text = element_text(face="bold", size=10, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='#f2f2f2'),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
  scale_size_continuous(range = c(1,10), breaks=c(1,10,100, 1000, 10000), limits = size_legend_vir)+
  scale_fill_viridis_d(option = 'turbo', end = 0.9, drop=F)+
  scale_colour_viridis_d(option = 'turbo', end = 0.9, drop=F)+
  guides(size=guide_legend(title="Viral contigs"), colour='none',
         fill=guide_legend(title='Class', override.aes = list(size=3), ncol=2))


Lake_vir_plot=ggplot(counts_vir_L, aes(x = Depth, y=Phylum, size = count, fill=Class, colour = Class)) +
  geom_point(alpha=0.8, shape=21, stroke=0.5, position =position_dodge(.5), show.legend = T)+
  theme_linedraw(base_size = 10) +
  coord_flip()+
  scale_x_reverse(breaks = int_breaks)+
  scale_y_discrete(labels = label_wrap_gen(18))+
  facet_nested_wrap(~Env, scales = "free_y", nrow = 1, ncol = 3, labeller = as_labeller(freshw_label))+
  labs(x = "Depth (m)", title = "Freshwater environments", tag = expression(bold(a))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing = unit(0.01, 'cm'),
        strip.text = element_text(face="bold", size=10, margin = margin(3,0,3,0)),
        panel.grid = element_line(color='#f2f2f2'),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        plot.margin = margin(1,1,1,10, unit = "pt")) +
  scale_size_continuous(range = c(1,10), breaks=c(1,10,100, 1000, 10000), limits = size_legend_vir)+
  scale_fill_viridis_d(option = 'turbo', end = 0.9, drop=F)+
  scale_colour_viridis_d(option = 'turbo', end = 0.9, drop=F)+
  guides(size=guide_legend(title="# viral contigs"), colour='none',
         fill=guide_legend(title='Class', override.aes = list(size=3), ncol=2))


blankplot=ggplot() + theme_void()
upper_plot=ggarrange(Lake_vir_plot, blankplot, common.legend = F, ncol = 2, widths = c(0.97, 0.03), align = "hv")
ggarrange(upper_plot, Marine_vir_plot, ncol = 1, common.legend = F, align = "hv")

ggsave("Virus_Depth.pdf", dpi=600,  width = 9.5, height = 10, bg='white', device = cairo_pdf)

#Size distribution
vir_L_meta$Phylum <- fct_relevel(vir_L_meta$Phylum, "unclassified virus", after = Inf)
vir_M_meta$Phylum <- fct_relevel(vir_M_meta$Phylum, "unclassified virus", after = Inf)

vir_L_meta_f <- vir_L_meta %>%
  group_by(Phylum) %>%
  ungroup() %>%
  droplevels()

vir_M_meta_f <- vir_M_meta %>%
  group_by(Phylum) %>%
  ungroup() %>%
  droplevels()

vir_L_meta_f$group=factor(vir_L_meta_f$group, levels = legend_stress)
vir_M_meta_f$group=factor(vir_M_meta_f$group, levels = legend_stress)


vir_size_L=ggplot(vir_L_meta_f, aes(x=fct_rev(Phylum), y=length/1000, color=group))+
  geom_point(position = position_dodge(0.8), size=1, stroke=0.5, alpha=0.8, show.legend = TRUE)+
  scale_color_viridis_d(drop=FALSE, option = "turbo", breaks=legend_stress, "Sampling location")+
  theme_linedraw(12)+
  geom_abline(intercept = 5, slope = 0, lty=2)+
  coord_flip()+
  theme(panel.grid = element_line(color='#f2f2f2'),
        legend.key.size = unit(0.5, 'cm'),
        plot.margin = margin(9,1,5.5,1, unit = "pt"))+
  xlab('Phylum')+ylab("Length of viral contigs (kb)")+
  labs(tag=expression(bold(a)), subtitle = expression(bold("Freshwater environments")))+
  guides(color=guide_legend(ncol=3, override.aes = list(size=3)))


vir_size_M=ggplot(vir_M_meta_f, aes(x=fct_rev(Phylum), y=length/1000, color=group))+
  geom_point(position = position_dodge(0.8), size=1, stroke=0.5, alpha=0.8, show.legend = TRUE)+
  scale_color_viridis_d(drop=FALSE, option = "turbo", breaks=legend_stress, "Sampling location")+
  theme_linedraw(12)+
  geom_abline(intercept = 5, slope = 0, lty=2)+
  coord_flip()+
  theme(panel.grid = element_line(color='#f2f2f2'),
        legend.key.size = unit(0.5, 'cm'),
        plot.margin = margin(9,1,5.5,1, unit = "pt"))+
  xlab('')+ylab("Length of viral contigs (kb)")+
  labs(tag=expression(bold(b)), subtitle = expression(bold("Marine environments")))+
  guides(color=guide_legend(ncol=3, override.aes = list(size=3)))

ggarrange(vir_size_L, vir_size_M, legend = "right", ncol = 2, common.legend = T)

ggsave("virus_size.png", dpi=200, width = 11, height = 5, bg='white', device = ragg::agg_png)


###### ARG-carrying MGEs #####
# Oceans
plasmid_contig_M=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/plasmid_contig_summary.csv")

blastp_M=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/blastp_MGE_ARG_results.csv")
genomad_M=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/genomad_summary.csv")
genomad_M=genomad_M[,-2]

genomad_M[genomad_M$gene == "k141_433317_2", ]


blastp_M <- blastp_M %>%
  mutate(Sample_target = sub("_.*", "", target))

blastp_M <- blastp_M %>%
  left_join(metaM %>% select(Sample_name, `Depth (m)`), 
            by = c("Sample" = "Sample_name")) %>%
  rename(Depth_query = `Depth (m)`)

blastp_M <- blastp_M %>%
  left_join(metaM %>% select(Sample_name, `Depth (m)`), 
            by = c("Sample_target" = "Sample_name")) %>%
  rename(Depth_target = `Depth (m)`)

blastp_M <- blastp_M %>%
  left_join(genomad_M %>% select(gene, Sample, Source, Contig), 
            by = c("query" = "gene", "Sample" = "Sample"))

#Lakes
plasmid_contig_L=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/plasmid_contig_summary.csv")

blastp_L=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/blastp_MGE_ARG_results.csv")
genomad_L=read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/genomad_summary.csv")
genomad_L=genomad_L[,-2]


blastp_L <- blastp_L %>%
  mutate(Sample_target = sub("_.*", "", target))

blastp_L <- blastp_L %>%
  left_join(meta %>% select(sample_id, geographic.location..depth.), 
            by = c("Sample" = "sample_id")) %>%
  rename(Depth_query = geographic.location..depth.)

blastp_L <- blastp_L %>%
  left_join(meta %>% select(sample_id, geographic.location..depth.), 
            by = c("Sample_target" = "sample_id")) %>%
  rename(Depth_target = geographic.location..depth.)


blastp_L <- blastp_L %>%
  left_join(meta %>% select(sample_id, environment..feature.), 
            by = c("Sample" = "sample_id")) %>%
  rename(Type = environment..feature.)


blastp_L <- blastp_L %>%
  left_join(genomad_L %>% select(gene, Sample, Source, Contig), 
            by = c("query" = "gene", "Sample" = "Sample"))


### Load the ARG data again (until before merging them with bin depth data!) to have their original headers
compiled_m <- bind_rows(
  rf_m_all %>% select(query, Class, ID), 
  rgi_m_all %>% rename(query = ORF_ID) %>% select(query, Class = Drug.Class, ID), 
  latent_m_all %>% select(query, Class = Class.1, ID = target)
)

compiled_l <- bind_rows(
  rf_lakes_all %>% select(query, Class, ID), 
  rgi_lakes_all %>% rename(query = ORF_ID) %>% select(query, Class = Drug.Class, ID), 
  latent_lakes_all %>% select(query, Class = Class.1, ID = target) 
)


blastp_M_ARG <- blastp_M %>%
  left_join(compiled_m, by = c("target" = "query"))

blastp_L_ARG <- blastp_L %>%
  left_join(compiled_l, by = c("target" = "query"))


blastp_M_ARG=blastp_M_ARG[-12,] #removed one detected gene by ResFinderFG because the very same gene was also identified by CARD


# Filter viruses to match the length-filtered version in vir_M
blastp_M_ARG <- blastp_M_ARG %>%
  filter((Source == "virus" & Contig %in% vir_M_meta$seq_name) |
           (Source != "virus"))

#  Filter viruses to match the length-filtered version in vir_L
blastp_L_ARG <- blastp_L_ARG %>%
  filter((Source == "virus" & Contig %in% vir_L_meta$seq_name) |
           (Source != "virus"))


# Plotting the results
legend_ab=union(sort(blastp_M_ARG$Group),sort(blastp_L_ARG$Group))

blastp_M_ARG$Group=factor(blastp_M_ARG$Group, levels = legend_stress)
blastp_L_ARG$Group=factor(blastp_L_ARG$Group, levels = legend_stress)

legend_group=union(sort(blastp_M_ARG$Group), sort(blastp_L_ARG$Group))


blastp_L_ARG_ss = blastp_L_ARG %>% distinct(ID, Group, .keep_all = TRUE)
blastp_M_ARG_ss = blastp_M_ARG %>% distinct(ID, Group, .keep_all = TRUE)

legend_class <- sort(union(
  union(blastp_M_ARG$Class, blastp_L_ARG$Class),
  union(blastp_M_ARG_ss$Class, blastp_L_ARG_ss$Class)))

mge_M=ggplot(blastp_M_ARG, aes(x=Depth_target, y=Depth_query, color = Group, shape=Class))+
  geom_abline (slope=1, linetype = "dashed", color="black")+
  geom_point(size=3, stroke=1.5, show.legend = T)+
  ylab("MGE-ARG : depth (m)")+xlab("MAG-ARG : depth (m)")+
  scale_color_viridis_d(drop=F, option = "turbo", breaks=legend_group)+
  scale_shape_manual(drop=FALSE, breaks = legend_class, values=0:8, limits=legend_class)+
  theme_linedraw(12)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        strip.text = element_text(face="bold", size=12, margin = margin(3,0,3,0)),
        legend.position = "none")+
  facet_wrap(~Source)+
  guides(shape=guide_legend(title='Drug class'),
         color=guide_legend(title='Sample location', ncol = 2))+
  labs(title = expression(bold(f)))

mge_L=ggplot(blastp_L_ARG, aes(x=Depth_target, y=Depth_query, color = Group, shape=Class))+
  geom_abline (slope=1, linetype = "dashed", color="black")+
  geom_point(size=3, stroke=1.5, show.legend = T)+
  ylab("MGE-ARG : depth (m)")+xlab("MAG-ARG : depth (m)")+
  scale_color_viridis_d(drop=FALSE, option = "turbo", breaks=legend_group)+
  scale_shape_manual(drop=FALSE, breaks = legend_class, values=0:8, limits=legend_class)+
  theme_linedraw(12)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        strip.text = element_text(face="bold", size=12, margin = margin(3,0,3,0)))+
  facet_wrap(~Source)+
  guides(shape=guide_legend(title='Drug class'),
         color=guide_legend(title='Sample location', ncol = 2))+
  labs(title = expression(bold(c)))

blankplot=ggplot() + theme_void()
onetooneplot <- ggarrange(blankplot, mge_L, blankplot, mge_M,  ncol = 1, common.legend = TRUE, legend = "right", 
                          align = "hv", heights = c(0.15, 0.85, 0.15, 0.85))

onetooneplot


#ARG/site
blastp_M_ARG_ss$Class=factor(blastp_M_ARG_ss$Class, levels = legend_class)
blastp_L_ARG_ss$Class=factor(blastp_L_ARG_ss$Class, levels = legend_class)

summary_marg_L <- blastp_L_ARG %>%
  group_by(Group, Source, ID, Class) %>%  
  summarise(n_rows = n(), .groups = "drop")

summary_marg_M <- blastp_M_ARG %>%
  group_by(Group, Source, ID, Class) %>%  
  summarise(n_rows = n(), .groups = "drop")

mge_sum_L=ggplot(summary_marg_L, aes(x=Source, shape=Class, color=Group))+
  geom_point(stat='count', position = position_dodge(0.6), size=4, stroke=1.5, show.legend = TRUE)+
  scale_shape_manual(drop=FALSE, breaks = legend_class, values=0:8, limits=legend_class)+
  scale_color_viridis_d(drop=FALSE, option = "turbo", breaks=legend_group)+
  scale_y_continuous(breaks=(seq(0, 3, 1)), limits = c(1,3))+
  theme_linedraw(12)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(9,3,5.5,1, unit = "pt"))+
  xlab('Source of MGE')+ylab("# mobilised\nunique ARG per site")+
  labs(title=expression(bold(b)))


mge_sum_M=ggplot(summary_marg_M, aes(x=Source, shape=Class, color = Group))+
  geom_count(stat='count', position = position_dodge(0.6), size=4, stroke=1.5, show.legend = TRUE)+
  scale_shape_manual(drop=FALSE, breaks = legend_class, values=0:8, limits=legend_class)+
  scale_y_continuous(breaks=(seq(0, 3, 1)), limits = c(1,3))+
  scale_color_viridis_d(drop=FALSE, option = "turbo", breaks=legend_group)+
  theme_linedraw(12)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(9,3,5.5,1, unit = "pt"))+
  xlab('Source of MGE')+ylab("# mobilised\nunique ARG per site")+
  labs(title=expression(bold(e)))

#sourcedata
summary_marg_M$Sample_location = "Marine"
summary_marg_L$Sample_location = "Freshwater"
write.csv(rbind(summary_marg_M, summary_marg_L), "/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/MGE_ARG.csv")

#Proportion of mobilizable ARGs in relation to MAG-ARG richness (unique instances!)
length(unique(blastp_M_ARG$ID)) # 5 mobile ARGs
length(unique(blastp_L_ARG$ID)) # 9 mobile ARGs

length(unique(compiled_m$ID)) # 39 ARGs
length(unique(compiled_l$ID)) # 144 ARGs
write.csv(data.frame(unique(compiled_l)), "/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/compiledARGs_L.csv")
write.csv(data.frame(unique(compiled_m)), "/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/compiledARGs_M.csv")

length(unique(blastp_M_ARG$ID))/length(unique(compiled_m$ID))*100 #12.82%
length(unique(blastp_L_ARG$ID))/length(unique(compiled_l$ID))*100 #6.25%


marine_pie <- data.frame(
  category=c("Mobilised ARGs", "Non-mobilised ARGs"),
  count=c(5, 34)
)

total <- 39   # all unique ARGs
marine_pie$fraction <- marine_pie$count / total


marine_pie$fraction <- c(marine_pie$count[1] / marine_pie$count[2], 1 - (marine_pie$count[1] / marine_pie$count[2]))
marine_pie$ymax <- cumsum(marine_pie$fraction)
marine_pie$ymin <- c(0, head(marine_pie$ymax, n=-1))
marine_pie$labelPosition <- (marine_pie$ymax + marine_pie$ymin) / 2
marine_pie$label <- paste0(marine_pie$category, "\n", sprintf("%.2f", marine_pie$fraction * 100), "% (n=", marine_pie$count, ")")

marine_pie_plot <- ggplot(marine_pie, aes(ymax=ymax, ymin=ymin, xmax=5, xmin=3, fill=category)) +
  geom_rect() +
  geom_label(x=4, aes(y=labelPosition, label=label), size=3.5) +
  scale_fill_brewer(palette="Blues") +
  coord_polar(theta="y") +
  xlim(c(2, 5)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(30,10,10,10, unit = "pt"))+
  labs(title=expression(bold(d)), tag = expression(bold("Marine\nenvironments")))


lake_pie <- data.frame(
  category=c("Mobilised ARGs", "Non-mobilised ARGs"),
  count=c(9, 135)
)

total <- 144   # all unique ARGs
lake_pie$fraction <- lake_pie$count / total
lake_pie$ymax <- cumsum(lake_pie$fraction)
lake_pie$ymin <- c(0, head(lake_pie$ymax, n=-1))
lake_pie$labelPosition <- (lake_pie$ymax + lake_pie$ymin) / 2
lake_pie$label <- paste0(lake_pie$category, "\n", sprintf("%.2f", lake_pie$fraction * 100), "% (n=", lake_pie$count, ")")

lake_pie_plot <- ggplot(lake_pie, aes(ymax=ymax, ymin=ymin, xmax=5, xmin=3, fill=category)) +
  geom_rect() +
  geom_label(x=4, aes(y=labelPosition, label=label), size=3.5) +
  scale_fill_brewer(palette="Blues") +
  coord_polar(theta="y") +
  xlim(c(2, 5)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(30,10,10,10, unit = "pt"))+
  labs(title=expression(bold(a)), tag = expression(bold("Freshwater\nenvironments")))


mge_summ_L=ggarrange(lake_pie_plot, mge_sum_L, legend = "none", ncol = 1, heights = c(0.6, 0.4), align = "v")
mge_summ_M=ggarrange(marine_pie_plot, mge_sum_M, legend = "none", ncol = 1, heights = c(0.6, 0.4), align = "v")
mge_summ=ggarrange(mge_summ_L, mge_summ_M, legend = 'none', ncol = 1)
ggarrange(mge_summ, onetooneplot, common.legend = T, legend = "right", widths = c(0.25, 0.75))

ggsave("MGE_ARG_combined.pdf", dpi=600, width = 15.5, height = 11, bg='white', device = cairo_pdf)

#source data
marine_pie$Sample_source="Marine"
lake_pie$Sample_source="Freshwater"

write.csv(rbind(marine_pie, lake_pie), "/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/pie.csv")




vir_L_meta <- vir_L_meta %>%
  mutate(seq_name = ifelse(grepl("provirus", seq_name), 
                           sub("provirus.*", "provirus", seq_name), 
                           seq_name))

vir_M_meta <- vir_M_meta %>%
  mutate(seq_name = ifelse(grepl("provirus", seq_name), 
                           sub("provirus.*", "provirus", seq_name), 
                           seq_name))


vir_mge_L=merge(blastp_L_ARG, vir_L_meta, by.x=c("Contig", "Sample"), by.y=c("seq_name", "Sample"), all.x = T)
vir_mge_L= vir_mge_L[vir_mge_L$Source != "plasmid", ]

vir_mge_M=merge(blastp_M_ARG, vir_M_meta, by.x=c("Contig", "Sample"), by.y=c("seq_name", "Sample"), all.x = T)
vir_mge_M= vir_mge_M[vir_mge_M$Source != "plasmid", ]

vir_mge_L$Group=factor(vir_mge_L$Group, levels = legend_stress)
vir_mge_M$Group=factor(vir_mge_M$Group, levels = legend_stress)

vir_mge_L$Class.x=factor(vir_mge_L$Class.x, levels = legend_class)
vir_mge_M$Class.x=factor(vir_mge_M$Class.x, levels = legend_class)


vir_mge_L <- vir_mge_L %>%
  mutate(Type = ifelse(grepl("Provirus", topology), "Prophages", "Phages"))
vir_mge_M <- vir_mge_M %>%
  mutate(Type = ifelse(grepl("Provirus", topology), "Prophages", "Phages"))

#Loading life cycle typing results
phatyp_L=read.table("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/phatyp_results/final_prediction/phatyp_prediction.tsv", header = T)
phatyp_M=read.table("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/phatyp_results/final_prediction/phatyp_prediction.tsv", header = T)

phatyp_L <- phatyp_L %>%
  rename(Phatype = TYPE)
phatyp_M <- phatyp_M %>%
  rename(Phatype = TYPE)

phatyp_L <- phatyp_L %>%
  mutate(Accession = gsub("-", "_", Accession))
phatyp_M <- phatyp_M %>%
  mutate(Accession = gsub("-", "_", Accession))

vir_mge_L <- vir_mge_L %>%
  mutate(Accession_c = paste(Contig, gsub("-", "_", coordinates), gsub("-", "_", Sample), sep = "_"))
vir_mge_L <- vir_mge_L %>%
  mutate(Accession_c = gsub("__", "_", Accession_c))

vir_mge_M <- vir_mge_M %>%
  mutate(Accession_c = paste(Contig, gsub("-", "_", coordinates), gsub("-", "_", Sample), sep = "_"))
vir_mge_M <- vir_mge_M %>%
  mutate(Accession_c = gsub("__", "_", Accession_c))

vir_mge_L_typ <- vir_mge_L %>%
  left_join(phatyp_L %>% select(Accession, Phatype, PhaTYPScore), 
            by = c("Accession_c" = "Accession"))
vir_mge_M_typ <- vir_mge_M %>%
  left_join(phatyp_M %>% select(Accession, Phatype, PhaTYPScore), 
            by = c("Accession_c" = "Accession"))




vir_mge_L_typ <- vir_mge_L_typ %>%
  mutate(Phatype = case_when(
    Phatype == "-" ~ "unknown", 
    TRUE ~ Phatype
  )) %>%
  filter(!is.na(Phatype),          # remove NA
         Phatype != "filtered",    # remove "filtered"
         PhaTYPScore >= 0.75)      # keep only >= 0.75

vir_mge_M_typ <- vir_mge_M_typ %>%
  mutate(Phatype = case_when(
    Phatype == "-" ~ "unknown", 
    TRUE ~ Phatype
  )) %>%
  filter(!is.na(Phatype),          # remove NA
         Phatype != "filtered",    # remove "filtered"
         PhaTYPScore >= 0.75)      # keep only >= 0.75

#Note: No marine instances after filtration!

vir_mge_L_typ$Phatype=factor(vir_mge_L_typ$Phatype, levels = c("temperate", "virulent"))

vir_mge_L_typ$Type <- tolower(vir_mge_L_typ$Type)


virmge_L=ggplot(vir_mge_L_typ, aes(x=Phylum, y=length/1000, color=Group, shape=Class.x))+
  geom_point(position = position_dodge(0.5), size=4, stroke=1.5,  show.legend = TRUE)+
  scale_color_viridis_d(drop=FALSE, option = "turbo", breaks=legend_group, name="Sample location")+
  scale_shape_manual(drop=FALSE, breaks = legend_class, values=0:8, limits=legend_class, "Drug class")+
  facet_grid2(vars(Type), vars(Phatype), render_empty = FALSE)+
  scale_x_discrete(labels = label_wrap_gen(10))+
  theme_linedraw(12)+
  coord_flip()+
  ylim(0,50)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.margin = margin(9,1,5.5,1, unit = "pt"),
        strip.text = element_text(face="bold"))+
  xlab('Phylum')+ylab("Length of viral contigs (kb)")+
  labs(subtitle = expression(bold("Freshwater environments")))

virmge_L



ggsave("MGE_virus_tax.pdf", dpi=600, width = 8, height = 6, bg='white', device = cairo_pdf)

write.csv(vir_mge_L_typ, "/Users/matev/Documents/Research/Chalmers/BlendARGs/Manuscript/Figures/source_data/phage_MGE_ARG.csv")


###### Plasmid distributions ######
plasmid_M_meta=merge(plasmid_contig_M, metaM, by.x="Sample", by.y='Sample_name')  
plasmid_L_meta=merge(plasmid_contig_L, meta, by.x="Sample", by.y='sample_id')  

classify_conjugation <- function(x) {
  case_when(
    is.na(x) | x == ""        ~ "No conjugation gene",
    grepl(";", x)             ~ "Multiple conjugation genes",
    TRUE                      ~ "Single conjugation gene"
  )
}


plasmid_M_meta <- plasmid_M_meta %>%
  mutate(conjugation_status = factor(classify_conjugation(conjugation_genes),
                                     levels = c("No conjugation gene", 
                                                "Single conjugation gene", 
                                                "Multiple conjugation genes")))

plasmid_L_meta <- plasmid_L_meta %>%
  mutate(conjugation_status = factor(classify_conjugation(conjugation_genes),
                                     levels = c("No conjugation gene", 
                                                "Single conjugation gene", 
                                                "Multiple conjugation genes")))

#Plasmid sizes
median(plasmid_M_meta$length)
median(plasmid_L_meta$length)
min(plasmid_M_meta$length)
min(plasmid_L_meta$length)
max(plasmid_M_meta$length)
max(plasmid_L_meta$length)

#Fraction of plasmids with conjugation genes
percentages_conj <- plasmid_M_meta %>%
  count(conjugation_status) %>%
  mutate(percentage = n / sum(n) * 100)

percentages_conj

percentages_conj <- plasmid_L_meta %>%
  count(conjugation_status) %>%
  mutate(percentage = n / sum(n) * 100)

percentages_conj

#Plot size distribution
plasmid_L_meta$group=factor(plasmid_L_meta$group, levels = legend_stress)
plasmid_M_meta$group=factor(plasmid_M_meta$group, levels = legend_stress)

plasmid_size_L=ggplot(plasmid_L_meta, aes(x=conjugation_status, y=log10(length/1000), color=group))+
  geom_point(position = position_dodge(0.8), size=1, stroke=0.5, alpha=0.6, show.legend = TRUE)+
  scale_x_discrete(labels = label_wrap_gen(15))+
  scale_color_viridis_d(drop=FALSE, option = "turbo", breaks=legend_stress, "Sampling location")+
  theme_linedraw(12)+
  coord_flip()+
  theme(panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(5,8,1,5, unit = "pt"),
        axis.title = element_text(size=12))+
  xlab('')+ylab(expression("Length of plasmid contigs log"[10] * "(kb)"))+
  labs(tag=expression(bold(a)), subtitle = expression(bold("Freshwater environments")))+
  guides(color=guide_legend(ncol=4, override.aes = list(size=3)))


plasmid_size_M=ggplot(plasmid_M_meta, aes(x=conjugation_status, y=log10(length/1000), color=group))+
  geom_point(position = position_dodge(0.8), size=1, stroke=0.5, alpha=0.6, show.legend = TRUE)+
  scale_x_discrete(labels = label_wrap_gen(15))+
  scale_color_viridis_d(drop=FALSE, option = "turbo", breaks=legend_stress, "Sampling location")+
  theme_linedraw(12)+
  coord_flip()+
  theme(panel.grid = element_line(color='#f2f2f2'),
        axis.text.y = element_blank(),
        plot.margin = margin(5,8,1,5, unit = "pt"),
        axis.title = element_text(size=12))+
  xlab('')+ylab(expression("Length of plasmid contigs log"[10] * "(kb)"))+
  labs(tag=expression(bold(b)), subtitle = expression(bold("Marine environments")))+
  guides(color=guide_legend(ncol=4, override.aes = list(size=3)))

ggarrange(plasmid_size_L, plasmid_size_M, legend = "right", ncol = 2, common.legend = T, align = "v")

ggsave("plasmid_size.pdf", dpi=600, width = 12, height = 5, bg='white', device = cairo_pdf)


###### FastANI ######
ocean_ani = read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Ocean/final_analysis/fastani_results/fastani_summary.csv", header = T)
ocean_ani$pair = paste(ocean_ani$source, ocean_ani$target, sep = "-")

lake_ani = read.csv("/Users/matev/Documents/Research/Chalmers/BlendARGs/Lakes/final_analysis/fastani_summary.csv", header = T)
lake_ani$pair = paste(lake_ani$source, lake_ani$target, sep = "-")

#calculate and add depth distances
depth_M=metaM[,c("Sample_name", "Depth (m)")]
distance_matrix_M <-as.matrix(dist(depth_M, method = "euclidean", upper = T, diag = T))

depth_L=meta[,c("sample", "geographic.location..depth.")]
distance_matrix_L <-as.matrix(dist(depth_L, method = "euclidean", upper = T, diag = T))

rownames(distance_matrix_M) <- depth_M$Sample_name
colnames(distance_matrix_M) <- depth_M$Sample_name
rownames(distance_matrix_L) <- depth_L$sample
colnames(distance_matrix_L) <- depth_L$sample

depth_M_dist <- melt(distance_matrix_M)
depth_M_dist$pair = paste(depth_M_dist$Var1, depth_M_dist$Var2, sep = "-")
depth_L_dist <- melt(distance_matrix_L)
depth_L_dist$pair = paste(depth_L_dist$Var1, depth_L_dist$Var2, sep = "-")


ocean_ani_d = merge(ocean_ani, depth_M_dist, by="pair", all.x = T)
ocean_ani_d$value[is.na(ocean_ani_d$value)] <- 0
ocean_ani_d <- ocean_ani_d[ocean_ani_d$bin_q != ocean_ani_d$bin_t, ]

lake_ani_d = merge(lake_ani, depth_L_dist, by="pair", all.x = T)
lake_ani_d$value[is.na(lake_ani_d$value)] <- 0
lake_ani_d <- lake_ani_d[lake_ani_d$bin_q != lake_ani_d$bin_t, ]

ocean_ani_d <- ocean_ani_d %>%
  mutate(PairID = ifelse(bin_q < bin_t, 
                         paste(bin_q, bin_t, sep = "-"), 
                         paste(bin_t, bin_q, sep = "-")))
ocean_ani_d <- ocean_ani_d[!duplicated(ocean_ani_d$PairID), ]
ocean_ani_d$PairID <- NULL

lake_ani_d <- lake_ani_d %>%
  mutate(PairID = ifelse(bin_q < bin_t, 
                         paste(bin_q, bin_t, sep = "-"), 
                         paste(bin_t, bin_q, sep = "-")))
lake_ani_d <- lake_ani_d[!duplicated(lake_ani_d$PairID), ]
lake_ani_d$PairID <- NULL


#Plotting
ocean_ani_d$Sample=factor(ocean_ani_d$Sample, levels = legend_stress)
lake_ani_d$Sample=factor(lake_ani_d$Sample, levels = legend_stress)

ani_M=ggplot(ocean_ani_d, aes(x=value, y=100-ANI, color = Sample, fill=Sample))+
  geom_point(size=2, stroke=1, shape=21, alpha=0.8, show.legend = TRUE)+
  theme_linedraw(12)+
  coord_flip()+
  scale_x_reverse()+
  ylim(0,26)+
  scale_color_viridis_d(drop=FALSE, option = "turbo")+
  scale_fill_viridis_d(drop=FALSE, option = "turbo")+
  #facet_wrap(~Sample)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(9,1,5.5,1, unit = "pt"),
        axis.title.y = element_blank())+
  xlab('Distance between water layers (m)')+ylab("MAG divergence (%)")+
  labs(tag=expression(bold(b)), subtitle = expression(bold("Marine environments")))+
  guides(fill=guide_legend(ncol=3, title="Sample location", by.row=T, override.aes = list(size=2)), color='none')

ani_L=ggplot(lake_ani_d, aes(x=value, y=100-ANI, color = Sample, fill=Sample))+
  geom_point(size=2, stroke=1, shape=21, alpha=0.8, show.legend = TRUE)+
  theme_linedraw(12)+
  coord_flip()+
  scale_x_reverse()+
  ylim(0,26)+
  scale_color_viridis_d(drop=FALSE, option = "turbo")+
  scale_fill_viridis_d(drop=FALSE, option = "turbo")+
  #facet_wrap(~Sample)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(9,1,5.5,1, unit = "pt"))+
  xlab('Distance between water layers (m)')+ylab("MAG divergence (%)")+
  labs(tag=expression(bold(a)), subtitle = expression(bold("Freshwater environments")))+
  guides(fill=guide_legend(ncol=3, title="Sample location", by.row=T, override.aes = list(size=2)), color='none')

ggarrange(ani_L, ani_M, ncol = 2, common.legend = T, widths = c(0.5, 0.5), legend = "right", align = "hv")
quartz.save("ANI_Depth.pdf", type="pdf", dpi=800, width = 12, height = 6, bg='white')




# Genetic incompatibility test by MAG divergence
high_sim_genomes_L <- unique(c(lake_ani_d$bin_q[lake_ani_d$ANI > 80],
                             lake_ani_d$bin_t[lake_ani_d$ANI > 80]))

high_sim_genomes_M <- unique(c(ocean_ani_d$bin_q[ocean_ani_d$ANI > 80],
                               ocean_ani_d$bin_t[ocean_ani_d$ANI > 80]))

all_genomes_L <- unique(tax_lake$user_genome)
all_genomes_M <- unique(tax_m$user_genome)

clean_name <- function(x) sub("\\.fa$", "", basename(x))
all_bins_L <- clean_name(all_genomes_L)
high_sim_bins_L <- clean_name(high_sim_genomes_L)

prop_high_sim_L <- length(intersect(all_bins_L, high_sim_bins_L)) / length(all_bins_L) * 100
prop_high_sim_L

all_bins_M <- clean_name(all_genomes_M)
high_sim_bins_M <- clean_name(high_sim_genomes_M)

prop_high_sim_M <- length(intersect(all_bins_M, high_sim_bins_M)) / length(all_bins_M) * 100
prop_high_sim_M


#Checking if the presence or lack of HGT is influenced by MAG divergence
hgt_ani_L=merge(hgts, lake_ani_d, by = "Sample", all.y = T)
hgt_ani_M=merge(hgtsM, ocean_ani_d, by = "Sample", all.y = T)

hgt_ani_L$HGT <- ifelse(hgt_ani_L$HGT_activity == 0, "Absence of HGT\n(N=12,193)", "Presence of HGT\n(N=13,593)")
hgt_ani_M$HGT <- ifelse(hgt_ani_M$HGT_activity == 0, "Absence of HGT\n(N=276)", "Presence of HGT\n(N=1,186)")

sum(hgt_ani_L$HGT_activity == 0)
length(hgt_ani_L$HGT_activity)-sum(hgt_ani_L$HGT_activity == 0)

sum(hgt_ani_M$HGT_activity == 0)
length(hgt_ani_M$HGT_activity)-sum(hgt_ani_M$HGT_activity == 0)

hgt_ani_M$Sample=factor(hgt_ani_M$Sample, levels = legend_stress)
hgt_ani_L$Sample=factor(hgt_ani_L$Sample, levels = legend_stress)

hgt_ani_L_plot=ggplot(hgt_ani_L, aes(y=100-ANI, x=HGT, color = Sample))+
  geom_point(show.legend = F, alpha=0.2, position = position_dodge2(0.6, padding = 0.2)) +
  stat_summary(fun.data = "mean_cl_boot", aes(color = Sample, shape="mean"), linewidth = 1, size = 1.5, 
               position = position_dodge(0.5), show.legend = T)+
  theme_linedraw(12)+
  ylim(0,26)+
  scale_color_viridis_d(drop=FALSE, option = "turbo")+
  scale_shape_manual(values = c("mean"=18), "Mean value")+
  #facet_wrap(~Sample)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(9,1,5.5,1, unit = "pt"),
        axis.title.x = element_blank())+
  ylab("MAG divergence (%)")+
  labs(tag=expression(bold(a)), subtitle = expression(bold("Freshwater environments")))+
  guides(color=guide_legend(ncol=3, title="Sample location", by.row=T, override.aes = list(size=1)),
         shape='none')


hgt_ani_M_plot=ggplot(hgt_ani_M, aes(y=100-ANI, x=HGT, color = Sample))+
  geom_point(show.legend = F, alpha=0.2, position = position_dodge2(0.6, padding=0.2)) +
  stat_summary(fun.data = "mean_cl_boot", aes(color = Sample, shape="mean"), linewidth = 1, size = 1.5, 
               position = position_dodge(0.5), show.legend = T)+
  theme_linedraw(12)+
  ylim(0,26)+
  scale_color_viridis_d(drop=FALSE, option = "turbo")+
  scale_shape_manual(values = c("mean"=18), "Mean value")+
  #facet_wrap(~Sample)+
  theme(panel.grid = element_line(color='#f2f2f2'),
        plot.margin = margin(9,1,5.5,1, unit = "pt"), 
        axis.title = element_blank())+
  ylab("MAG divergence (%)")+
  labs(tag=expression(bold(b)), subtitle = expression(bold("Marine environments")))+
  guides(color=guide_legend(ncol=3, title="Sample location", by.row=T, override.aes = list(size=1)),
         shape="none")

ggarrange(hgt_ani_L_plot, hgt_ani_M_plot, ncol = 2, common.legend = T, widths = c(0.5, 0.5), legend = "right", align = "hv")
quartz.save("ANI_HGT.pdf", type="pdf", dpi=800, width = 11, height = 6, bg='white')
