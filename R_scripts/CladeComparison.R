library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

rm(list=ls())
setwd("~/work/CRBE_STAGE2A/ecology")

parse_18S_taxonomy <- function(tax_string, list_clade_ITS) {
  
  split_tax <- unlist(strsplit(tax_string, "\\|"))
  
  if (split_tax[7] == "Dikarya" && !is.na(split_tax[7])) {
    i <- 1  
  } else {
    i <- 0
  }
  
  phylum_candidate <- split_tax[8 + 0]
  if (phylum_candidate %in% list_clade_ITS[[clade]]) {
    phylum <- phylum_candidate  
  } else {
    phylum <- split_tax[8 + i] 
  }
  
  ranks <- data.frame(
    Kingdom = split_tax[6], 
    Phylum = gsub("mycotina", "mycota", phylum), 
    Class = split_tax[10 + i],   
    Order = split_tax[11 + i],  
    Family = split_tax[12 + i],  
    Genus = split_tax[13 + i],   
    Species = split_tax[14 + i]  
  )
  
  ranks[] <- lapply(ranks, function(x) ifelse(grepl("incertae_sedis", x, ignore.case = TRUE), NA, x))
  ranks[] <- lapply(ranks, function(x) ifelse(grepl("phyta", x, ignore.case = TRUE), "Plantae", x))
  ranks[] <- lapply(ranks, function(x) ifelse(grepl("mycota", x, ignore.case = TRUE), paste0("Fungi_", x), x))
  return(ranks)
}


parse_ITS_taxonomy <- function(tax_string) {
  split_tax <- unlist(strsplit(tax_string, "\\|"))
  
  # Créer un vecteur vide par défaut
  get_tax_val <- function(prefix) {
    val <- grep(prefix, split_tax, value = TRUE)
    if (length(val) == 0) return(NA)
    return(gsub(".*__", "", val))
  }
  
  df <- data.frame(
    Kingdom = get_tax_val("k__"),
    Phylum = get_tax_val("p__"),
    Class = get_tax_val("c__"),
    Order = get_tax_val("o__"),
    Family = get_tax_val("f__"),
    Genus = get_tax_val("g__"),
    Species = get_tax_val("s__")
  )
  
  df[] <- lapply(df, function(x) ifelse(grepl("incertae_sedis", x, ignore.case = TRUE), NA, x))
  df[] <- lapply(df, function(x) ifelse(grepl("phyta", x, ignore.case = TRUE), "Plantae", x))
  df[] <- lapply(df, function(x) ifelse(grepl("mycota", x, ignore.case = TRUE), paste0("Fungi_", x), x))
  
  return(df)
}


table_18S <- read.table("OTU_table_OTU97_200filtered_quebec18S.txt", sep = "\t", header = TRUE)
table_ITS <- read.table("OTU_table_OTU97_filtered_quebecITS.txt", sep = "\t", header = TRUE)

clade <- "Phylum"

tax_ITS_raw <- do.call(rbind, lapply(table_ITS$taxonomy, function(tax_string) {
  split_tax <- unlist(strsplit(tax_string, "\\|"))
  df <- data.frame(
    Kingdom = grep("k__", split_tax, value = TRUE),
    Phylum = grep("p__", split_tax, value = TRUE),
    Class = grep("c__", split_tax, value = TRUE),
    Order = grep("o__", split_tax, value = TRUE),
    Family = grep("f__", split_tax, value = TRUE),
    Genus = grep("g__", split_tax, value = TRUE),
    Species = grep("s__", split_tax, value = TRUE)
  )
  df[] <- lapply(df, function(x) gsub(".*__", "", x))
  return(df)
}))

list_clade_ITS <- tax_ITS_raw %>%
  group_by(!!sym(clade))  %>%  
  summarise() 

tax_ITS <- do.call(rbind, lapply(table_ITS$taxonomy, parse_ITS_taxonomy)) %>% 
  mutate(Source = "ITS", OTU = table_ITS$OTU) %>%
  mutate(Source = "ITS", OTU = table_ITS$OTU, Abundance = table_ITS[, 2])


tax_18S <- do.call(rbind, lapply(table_18S$taxonomy, function(tax_string) {
  parse_18S_taxonomy(tax_string, list_clade_ITS)
})) %>% 
  mutate(Source = "18S", OTU = table_18S$OTU) %>%
  mutate(Source = "18S", OTU = table_18S$OTU, Abundance = table_18S[, 2])

num_otus_18S <- nrow(table_18S)
num_otus_ITS <- nrow(table_ITS)

total_abundance_18S <- sum(tax_18S$Abundance, na.rm = TRUE)
total_abundance_ITS <- sum(tax_ITS$Abundance, na.rm = TRUE)

clade_18S <- tax_18S %>%
  group_by(!!sym(clade)) %>%
  summarise(
    Count = n(),
    Total_Abundance = sum(Abundance, na.rm = TRUE)
  ) %>%
  filter(Count > 5) %>%
  mutate(
    ProportionOTUs = (Count / sum(Count)) * 100,
    ProportionAbundance = (Total_Abundance / sum(Total_Abundance)) * 100
  )

clade_ITS <- tax_ITS %>%
  group_by(!!sym(clade)) %>%
  summarise(
    Count = n(),
    Total_Abundance = sum(Abundance, na.rm = TRUE)
  ) %>%
  filter(Count > 5) %>%
  mutate(
    ProportionOTUs = (Count / sum(Count)) * 100,
    ProportionAbundance = (Total_Abundance / sum(Total_Abundance)) * 100
  )

all_clades <- unique(c(clade_18S[[clade]], clade_ITS[[clade]]))

clade_colors <- scales::hue_pal()(length(all_clades))
names(clade_colors) <- all_clades

#number of OTUs
plot_18S <- ggplot(clade_18S, aes(x = "", y = ProportionOTUs, fill = !!sym(clade))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = clade_colors) + 
  labs(title = paste("Proportions des", clade, "(18S)\nNombre d'OTUs :", num_otus_18S),
       fill = clade) +
  theme_void() +
  theme(legend.position = "right")

plot_ITS <- ggplot(clade_ITS, aes(x = "", y = ProportionOTUs, fill = !!sym(clade))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = clade_colors) + 
  labs(title = paste("Proportions des", clade, "(ITS)\nNombre d'OTUs :", num_otus_ITS),
       fill = clade) +
  theme_void() +
  theme(legend.position = "right")

grid.arrange(plot_18S, plot_ITS, ncol = 2)

#Abudances
plot_18S <- ggplot(clade_18S, aes(x = "", y = ProportionAbundance, fill = !!sym(clade))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = clade_colors) + 
  labs(title = paste("Proportions des", clade, "(18S)\nNombre d'OTUs :", num_otus_18S),
       fill = clade) +
  theme_void() +
  theme(legend.position = "right")

plot_ITS <- ggplot(clade_ITS, aes(x = "", y = ProportionAbundance, fill = !!sym(clade))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = clade_colors) + 
  labs(title = paste("Proportions des", clade, "(ITS)\nNombre d'OTUs :", num_otus_ITS),
       fill = clade) +
  theme_void() +
  theme(legend.position = "right")

grid.arrange(plot_18S, plot_ITS, ncol = 2)

