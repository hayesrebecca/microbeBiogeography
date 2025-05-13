
social_obligate_list <- c("Lactobacillaceae",
                          "Bifidobacteriaceae",
                          "Neisseriaceae",
                          "Orbaceae",
                          "Bartonella",
                          "Acetobacteraceae",
                          "Anaplasmataceae",
                          "Erwiniaceae",
                          "Hafniaceae"
)

solitary_obligate_list <- c("Acetobacteraceae",
                            "Bacillaceae",
                            "Burkholderiaceae",
                            "Clostridiaceae",
                            "Comamonadaceae",
                            "Enterobacteriaceae",
                            "Lachnospiraceae",
                            "Lactobacillaceae",
                            "Methylobacteriaceae",
                            "Moraxellaceae",
                            "Sphingomonadaceae", 
                            "Oxalobacteraceae",
                            "Anaplasmataceae",
                            "Erwiniaceae",
                            "Hafniaceae"
)


these_taxa_present <- function(meta_list, taxa){
  these_taxa <- meta_list %>%
    filter(grepl(taxa, bacteria))  # <-- use grepl here
  these_taxa$bacteria
}

## social hosts

bombus_social <- these_taxa_present(bombus_meta, paste(social_obligate_list, collapse="|"))
bombus_solitary <- these_taxa_present(bombus_meta, paste(solitary_obligate_list, collapse="|"))

apis_social <- these_taxa_present(apis_meta, paste(social_obligate_list, collapse="|"))
apis_solitary <- these_taxa_present(apis_meta, paste(solitary_obligate_list, collapse="|"))

## solitary hosts
anthophora_social <- these_taxa_present(anthophora_meta, paste(social_obligate_list, collapse="|"))
anthophora_solitary <- these_taxa_present(anthophora_meta, paste(solitary_obligate_list, collapse="|"))

melissodes_social <- these_taxa_present(melissodes_meta, paste(social_obligate_list, collapse="|"))
melissodes_solitary <- these_taxa_present(melissodes_meta, paste(solitary_obligate_list, collapse="|"))

andrena_social <- these_taxa_present(andrena_meta, paste(social_obligate_list, collapse="|"))
andrena_solitary <- these_taxa_present(andrena_meta, paste(solitary_obligate_list, collapse="|"))

megachile_social <- these_taxa_present(megachile_meta, paste(social_obligate_list, collapse="|"))
megachile_solitary <- these_taxa_present(megachile_meta, paste(solitary_obligate_list, collapse="|"))

hosts <- c("bombus", "apis", "anthophora", "melissodes", "andrena", "megachile")

host_microbes <- c()
for(i in hosts){
   soc <- paste0(i, "_social")
   host_microbes <- base::append(host_microbes, soc)
   sol <- paste0(i, "_solitary")
   host_microbes <- base::append(host_microbes, sol)
}


############################
## making summary table for each genus

apis_meta$genus <- rep("Apis", length(apis_meta$bacteria))
bombus_meta$genus <- rep("Bombus", length(bombus_meta$bacteria))
melissodes_meta$genus <- rep("Melissodes", length(melissodes_meta$bacteria))
anthophora_meta$genus <- rep("Anthophora", length(anthophora_meta$bacteria))
andrena_meta$genus <- rep("Andrena", length(andrena_meta$bacteria))
megachile_meta$genus <- rep("Megachile", length(megachile_meta$bacteria))

all_meta <- rbind(apis_meta, 
                  bombus_meta,
                  melissodes_meta,
                  anthophora_meta,
                  andrena_meta,
                  megachile_meta)

library(tidyverse)

library(knitr)

library(dplyr)
library(tidyr)
library(stringr)

combined_table <- all_meta %>%
  filter(percent_sites == 100) %>%
  # First, create a flag for pathogens
  mutate(PathogenObligate = ifelse(grepl("Anaplasmataceae|Erwiniaceae|Hafniaceae", bacteria), 1, 0)) %>%
  # Then create the Weak column based on all the Obligate columns
  mutate(Weak = ifelse(rowSums(across(c(SolitaryObligate, SocialObligate, BothObligate, PathogenObligate))) > 0, 0, 1)) %>%
  # Extract taxonomy levels
  mutate(
    ASV_Order = str_extract(bacteria, "o__[^;]+") %>% str_remove("o__"),
    ASV_Family = str_extract(bacteria, "f__[^;]+") %>% str_remove("f__"),
    ASV_Genus = str_extract(bacteria, "g__[^;]+") %>% str_remove("g__"),
    ASV_Species = str_extract(bacteria, "s__[^;]+") %>% str_remove("s__")
  ) %>%
  pivot_longer(
    cols = c(SolitaryObligate, SocialObligate, BothObligate, PathogenObligate, Weak),
    names_to = "Associate",
    values_to = "AssociateVal"
  ) %>%
  filter(AssociateVal == 1) %>%
  select(-AssociateVal, -obligate_color, -n_sites_present, -percent_sites, -bacteria, -Associate,
         genus, ASV_Order, ASV_Family, ASV_Genus, ASV_Species) %>%
  mutate(ASV_Species = gsub("(?<!_)_(?!_)", " ", ASV_Species, perl = TRUE)) %>%
  replace_na(list(ASV_Species="", ASV_Family="", ASV_Order="", ASV_Genus=""))

kable(combined_table, "latex")
