library(ggplot2)

rm(list=ls())
setwd("~/work/CRBE_STAGE2A/ecology")

### Charger la table FungalTraits
database <- read.table("13225_2020_466_MOESM4_ESM.csv", sep = ",", header = TRUE)
database$Genus <- database$GENUS
database$trophic_mode_fg <- database$primary_lifestyle

# Regrouper tous les types de saprotrophes en un seul groupe
database$trophic_mode_fg[grep("saprotroph", database$trophic_mode_fg, ignore.case = TRUE)] <- "saprotroph"

### Charger la table OTU
OTU_table <- read.table("OTU_table_OTU97_filtered_ITS.txt", sep = "\t", header = TRUE)

# Extraire les genres fongiques
split_tax <- unlist(strsplit(OTU_table$taxonomy, "\\|"))
OTU_table$Genus <- grep("g__", split_tax, value = TRUE)
OTU_table$Genus <- lapply(OTU_table$Genus, function(x) gsub(".*__", "", x))
OTU_table$Genus <- unlist(OTU_table$Genus)

# Ajouter la colonne de mode trophique par correspondance avec le genre
OTU_table$Trophic_mode <- NA
for (i in 1:nrow(OTU_table)) {
  if (OTU_table$Genus[i] %in% database$Genus) {
    trophic_modes <- sort(unique(database$trophic_mode_fg[which(database$Genus == OTU_table$Genus[i])]))
    if (length(trophic_modes) > 1 && any(is.na(trophic_modes))) {
      trophic_modes <- trophic_modes[!is.na(trophic_modes)]
    }
    OTU_table$Trophic_mode[i] <- paste0(trophic_modes, collapse = '-')
  }
}

# OTU NUMBER
trophic_counts <- as.data.frame(table(OTU_table$Trophic_mode))
colnames(trophic_counts) <- c("Trophic_mode", "Count")

ggplot(trophic_counts, aes(x = "", y = Count, fill = Trophic_mode)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Proportion des modes trophiques des OTUs") +
  theme(legend.position = "right")

#ABUNDANCE
abundance_by_mode <- aggregate(abundance ~ Trophic_mode, data = OTU_table, sum, na.rm = TRUE)

# Renommer les colonnes pour le graphique
colnames(abundance_by_mode) <- c("Trophic_mode", "Abundance")

# Graphique en secteurs basé sur l'abondance
ggplot(abundance_by_mode, aes(x = "", y = Abundance, fill = Trophic_mode)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Proportion des abondances par mode trophique") +
  theme(legend.position = "right")
