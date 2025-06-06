# Module: tree_functions

# This module includes functions to aid in plotting microbial phylogenetic trees for each bee host genus.

# Function: match_shared_ID
# This function takes two data frames and returns a filtered version of the first data frame.
# It keeps only the rows where the `UniqueID` column matches any `UniqueID` values 
# found in the second data frame.

# Arguments:
# - first_df: A data frame that contains a column named `UniqueID`. 
#   This is the data frame to be filtered.
# - second_df: A data frame that contains a column named `UniqueID`. 
#   This is the data frame used to identify matching `UniqueID` values.

# Returns:
# A filtered version of `first_df` containing only the rows with `UniqueID` values 
# that are also present in `second_df`.
match_shared_ID <- function(first_df, second_df){
  shared <- first_df$UniqueID[first_df$UniqueID %in% second_df$UniqueID]
  matched_df <- first_df %>%
    filter(UniqueID %in% shared)
  matched_df
}

# Function: match_shared_tiplabels
# This function matches the tip labels of a phylogenetic tree to the column names 
# of a presence-absence table and returns an updated table with only the matched columns.

# Arguments:
# - tree: A phylogenetic tree object that contains a vector of tip labels (`tip.label`).
# - pres_abs_table: A presence-absence table, where columns represent features, 
#   and there is a column named `UniqueID` to identify rows.

# Returns:
# A filtered version of the presence-absence table that retains only the columns 
# corresponding to the tip labels of the tree, along with the `UniqueID` column.
match_shared_tiplabels <- function(tree, pres_abs_table){ #input a phyloseq tree and a presence absense table
  tree_tips <- tree$tip.label #create object to store tree tip labels (strains)

  all_cols <- pres_abs_table %>% #filter the pres/abs table to remove unique ID
    select(-UniqueID) 
 
  match_cols <- all_cols[colnames(all_cols) %in% tree_tips] #match the features to the tree tips
  
  match_cols$UniqueID <- pres_abs_table$UniqueID # add back in the unique IDs
  match_cols #return the updated pres abs table with the tiplabels matched to the tree
}

# Function: collapse_identical_tips
# This function collapses identical tips in a phylogenetic tree by removing redundant tips with the same label, ensuring the tree retains a single representative branch for identical labels.

# Arguments:
# - phy: A phylogenetic tree object. The tree should contain tip labels stored in `phy$tip.label`.
# - tip_label: A character string representing the label of the tips in the tree that should be collapsed into a single branch.

# Returns:
# A new phylogenetic tree where all tips with the specified `tip_label` are collapsed into a single branch.

# Example:
# Input:
# phy <- read.tree(text = "((A,A),(B,B),(C,C));") # Example tree
# collapsed_phy <- collapse_identical_tips(phy, "A")
# print(collapsed_phy)
# Output:
# The tree will have a single tip labeled "A", with redundant tips removed.

# Details:
# - The function identifies all tips in the tree with the label matching `tip_label`.
# - It determines which tips can be collapsed based on their most recent common ancestor (MRCA).
# - Only one tip for each group of identical labels is retained, and others are removed.
#
collapse_identical_tips <- function(phy,tip_label){
  #matching_tips is initialized with the indices of tips in the phylogenetic tree (phy) whose labels match the provided tip_label. The function identifies all tips with the same label.
  matching_tips <- which(phy$tip.label==tip_label)
  nt <- length(phy$tip.label) # number of tips in tree
  nm <- length(matching_tips) # Number of tips matching the label
  keep <- numeric(nm) #keep is initialized as a numeric vector of length nm. It is used to keep track of which tips should be retained (1) and which tips should be dropped (0) in the new tree.
  
  #The while loop iterates through the indices of matching_tips to determine which tips to keep and which to drop.
  cur_tip <- 1
  #Inside the loop, the variable cur_tip is the current tip being considered, and next_tip is the tip immediately after cur_tip.
  while(cur_tip<=nm){
    if(cur_tip == nm){
      keep[cur_tip] <- 1
      break
    }
    next_tip <- cur_tip + 1
    #mrca_ (most recent common ancestor) is calculated using the getMRCA function for the tips identified by matching_tips[cur_tip] and matching_tips[next_tip]. This helps find the common ancestor of the current tip and the next tip.
    mrca_ <- getMRCA(phy,c(matching_tips[cur_tip],matching_tips[next_tip]))
    #descendants contains the indices of all descendants of the common ancestor, which includes both tips and internal nodes.
    descendants <- getDescendants(phy, mrca_)
    #descendant_tips is calculated to include only those indices from descendants that correspond to actual tips (i.e., indices less than or equal to nt).
    descendant_tips <- descendants[descendants<=nt]
    #The function checks if all descendant_tips are in the list of matching_tips. If they are, it means all these tips can be collapsed into a single branch, and they are marked to be kept.
    #The variable keep is updated accordingly, and cur_tip is advanced to skip the tips that have been collapsed.
    if(all(descendant_tips %in% matching_tips)){
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + length(descendant_tips)
    }else{ #If not all descendant_tips are in the list of matching_tips, it means that not all tips can be collapsed, so the current tip is marked to be kept, and cur_tip is incremented by 1.
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + 1
    }
  }
  #After the loop completes, the to_drop variable contains the indices of the tips that need to be dropped to collapse identical labels.
  to_drop <- matching_tips[!keep]
  #create a new phylogenetic tree (new_phy) by using the drop.tip function to remove the tips identified in to_drop.
  new_phy <- drop.tip(phy,to_drop)
  
  return(new_phy)
}

# Function: phylotree_heatmap_byGenus
# This function creates a phylogenetic tree with an associated heatmap that visualizes 
# bacterial presence/absence across different sites. The heatmap includes information 
# on bacterial families and the number of sites where a given taxon is present.

# Arguments:
# - tree.object: A phylogenetic tree object.
# - meta: Metadata containing site and genus/species information.
# - genus.or.spp: A string specifying whether to focus on 'Species' or 'Genus'.
# - this.species: The specific species or genus of interest.
# - presAbsTable: A table with bacterial presence/absence data.
# - site.order: A vector specifying the desired order of sites.
# - all_levels: A logical value indicating whether to include all taxonomic levels 
#   (default TRUE).
# - levels_to_drop: Taxonomic levels to exclude from the tree.
# - clade_names: Optional vector of clade names for custom groupings.
# - do_collapse: A logical value indicating whether to collapse identical tips in 
#   the tree (default FALSE).

# Returns:
# A list containing:
# 1. A `ggtree` object with the phylogenetic tree and heatmap.
# 2. A dataframe with features and site metadata.
# 3. The ordered tip labels of the tree.
#
phylotree_heatmap_byGenus <- function(tree.object,
                                      meta,
                                      genus.or.spp,
                                      this.species,
                                      presAbsTable,
                                      site.order,
                                      all_levels=TRUE,
                                      levels_to_drop,
                                      clade_names=NULL,
                                      do_collapse=FALSE,
                                      add_tip_labs=FALSE){
  #filter to include just the unique IDs in the specified genus
  if (genus.or.spp=='Species'){
    sp_ids <- meta %>%
      filter(GenusSpecies==this.species) %>%
      select(UniqueID)
    
    #pull out all sites that include the specified genus
    my_sites <- unique(meta$Site[meta$GenusSpecies==this.species])
  }
  if (genus.or.spp=='Genus'){
    sp_ids <- meta %>%
      filter(Genus==this.species) %>%
      select(UniqueID)
    
    #pull out all sites that include the specified genus
    my_sites <- unique(meta$Site[meta$Genus==this.species])
  }
  
  #remove tips from the tree that are not in the list of unique IDs in the specified genus
  trimmed_tree <- prune_samples(rownames(tree.object@sam_data) %in% sp_ids$UniqueID, tree.object)
  
  #remove taxa from the tree where there were less than zero observations
  pruned_tree <- prune_taxa(taxa_sums(trimmed_tree) > 0, trimmed_tree)
  
  #make labels from the taxonomic info
  feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')
  
  # ## convert to a phylo class which is more useful downstream
  gentree <- phy_tree(pruned_tree, errorIfNULL=TRUE)
  
  ## match the tip labs to the table with feature ID and Taxon
  gentree$tip.label  <-  feature.2.tax.16s$Taxon[match(gentree$tip.label,
                                                       feature.2.tax.16s$Feature.ID)]
  
  # Identify tips with labels exactly matching '16s:D_0__Bacteria'
  matching_tips <- grep('^16s:D_0__Bacteria$', gentree$tip.label)
  
  # Drop the matching tips
  gentree <- drop.tip(gentree, matching_tips)
  orig_len <- length(gentree$tip.label)
  print(paste('original # tips', orig_len))
  
  if (do_collapse == TRUE){
    #pull out unique tip labels
    groups <- unique(gentree$tip.label)
    
    #collapse branches that have the same label
    for (this.group in groups){
      gentree <- collapse_identical_tips(gentree, this.group)
    }
  }
  
  
  print(length(gentree$tip.label))
  if (do_collapse == TRUE){
    #pull out unique tip labels
    original_nodes <- gentree %>% as_tibble()
    groups <- unique(gentree$tip.label)
    
    
    # collapse branches that have the same label
    for (this.group in groups){
      gentree <- collapse_identical_tips(gentree, this.group)
    }
    new_nodes <- gentree %>% as_tibble()
  }
  print(length(gentree$tip.label))
  
  matched_presabs <- match_shared_tiplabels(gentree, presAbsTable)

  matched_pres_meta <- match_shared_ID(matched_presabs, meta)
  
  matched_id <- matched_pres_meta$UniqueID
  row.names(matched_pres_meta) <- matched_id
  if (genus.or.spp=='Species'){
    meta_match_sites <- match_shared_ID(meta, matched_pres_meta) %>%
      select(UniqueID, Site, GenusSpecies) %>%
      mutate(Site = factor(Site)) %>%
      filter(GenusSpecies==this.species) %>%
      select(!GenusSpecies) %>%
      group_by(UniqueID, Site) %>%
      mutate(n= n()) %>%
      pivot_wider(names_from=Site,
                  values_from = n,
                  names_expand = TRUE,
                  id_expand=TRUE) %>%
      pivot_longer(cols=2:length(colnames(.)),
                   names_to='Site',
                   values_to='Site_present') %>%
      filter(Site_present > 0) 
  }
  
  if (genus.or.spp == 'Genus') {
    
    # Get UniqueID and Site where the genus is present
    meta_match_sites <- match_shared_ID(meta, matched_pres_meta) %>%
      select(UniqueID, Site, Genus) %>%
      filter(Genus == this.species) %>%
      select(-Genus) %>%
      distinct(UniqueID, Site) %>%
      mutate(Site = factor(Site))
  }
  
  # Now get the total number of unique sites for this genus
  total_sites_genus <- meta_match_sites %>%
    pull(Site) %>%
    n_distinct()
  
  # Then match metadata and calculate percent presence
  features_site_metadata <- match_shared_ID(matched_pres_meta, meta_match_sites) %>%
    right_join(meta_match_sites, by = 'UniqueID') %>%
    pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
    filter(bact_pres == 1) %>%
    group_by(bacteria, Site) %>%
    summarise(present = 1, .groups = 'drop') %>%   # mark presence at each site
    group_by(bacteria) %>%
    summarise(
      n_sites_present = n(),                     # how many sites this bacteria is present at
      percent_sites = (n_sites_present / total_sites_genus) * 100
    ) %>%
    mutate(
      SocialObligate = as.numeric(str_detect(bacteria, "Bifidobacteriaceae|Neisseriaceae|Orbaceae|Bartonella")),
      SolitaryObligate = as.numeric(str_detect(bacteria, "Bacillaceae|Burkholderiaceae|Clostridiaceae|Comamonadaceae|Enterobacteriaceae|Methylobacteriaceae|Moraxellaceae|Sphingomonadaceae|Oxalobacteraceae")),
      BothObligate = as.numeric(str_detect(bacteria, "Lactobacillaceae|Acetobacteraceae"))
    )
  
  
  
  features_site_metadata <- features_site_metadata %>%
    mutate(obligate_color = case_when(
      # BothObligate
      grepl("Acetobacteraceae", bacteria) ~ "#e6194B",
      grepl("Lactobacillaceae", bacteria) ~ "#3cb44b",
      
      # SocialObligate
      grepl("Bartonella", bacteria) ~ "#ffe119",
      grepl("Bifidobacteriaceae", bacteria) ~ "#4363d8",
      grepl("Neisseriaceae", bacteria) ~ "#f58231",
      grepl("Orbaceae", bacteria) ~ "#911eb4",
      
      # SolitaryObligate
      grepl("Bacillaceae", bacteria) ~ "#42d4f4",
      grepl("Burkholderiaceae", bacteria) ~ "#f032e6",
      grepl("Clostridiaceae", bacteria) ~ "#bfef45",
      grepl("Comamonadaceae", bacteria) ~ "#fabed4",
      grepl("Enterobacteriaceae", bacteria) ~ "#469990",
      grepl("Methylobacteriaceae", bacteria) ~ "#dcbeff",
      grepl("Moraxellaceae", bacteria) ~ "#9A6324",
      grepl("Oxalobacteraceae", bacteria) ~ "#800000",
      grepl("Sphingomonadaceae", bacteria) ~ "#aaffc3",
      
      # Pathogens
      grepl("Wolbachia", bacteria) ~ "#808000",
      grepl("Erwinia", bacteria) ~ "#ffd8b1",
      grepl("Hafnia", bacteria) ~ "#000075",
      
      TRUE ~ NA_character_
    ))
  
  
  
  
  tree_tip_labs <- gentree$tip.label
  
  #dropping branches that weren't in the presence abs table
  final_drop <- gentree$tip.label[!(gentree$tip.label %in% features_site_metadata$bacteria)]
  
  
  if (length(final_drop) > 0){
    gentree <- drop.tip(gentree, final_drop)
  }
  
  ## save out order of tips
  is_tip <- gentree$edge[,2] <= length(gentree$tip.label)
  
  ordered_tips <- gentree$edge[is_tip, 2]
  
  tip.order <- gentree$tip.label[ordered_tips]
  
  if(add_tip_labs==FALSE){
  p <- ggtree(gentree, layout='rectangular') 
  p
  }else{
    p <- ggtree(gentree, layout='rectangular') + geom_tiplab()
    p
  }
  # Add logical columns to p$data for each family
  p$data$Acetobacteraceae_match <- grepl("Acetobacteraceae", p$data$label, fixed = TRUE)
  p$data$Lactobacillaceae_match <- grepl("Lactobacillaceae", p$data$label, fixed = TRUE)
  
  p$data$Bartonella_match <- grepl("Bartonella", p$data$label, fixed = TRUE)
  p$data$Bifidobacteriaceae_match <- grepl("Bifidobacteriaceae", p$data$label, fixed = TRUE)
  p$data$Neisseriaceae_match <- grepl("Neisseriaceae", p$data$label, fixed = TRUE)
  p$data$Orbaceae_match <- grepl("Orbaceae", p$data$label, fixed = TRUE)
  
  p$data$Bacillaceae_match <- grepl("Bacillaceae", p$data$label, fixed = TRUE)
  p$data$Burkholderiaceae_match <- grepl("Burkholderiaceae", p$data$label, fixed = TRUE)
  p$data$Clostridiaceae_match <- grepl("Clostridiaceae", p$data$label, fixed = TRUE)
  p$data$Comamonadaceae_match <- grepl("Comamonadaceae", p$data$label, fixed = TRUE)
  p$data$Enterobacteriaceae_match <- grepl("Enterobacteriaceae", p$data$label, fixed = TRUE)
  p$data$Methylobacteriaceae_match <- grepl("Methylobacteriaceae", p$data$label, fixed = TRUE)
  p$data$Moraxellaceae_match <- grepl("Moraxellaceae", p$data$label, fixed = TRUE)
  p$data$Oxalobacteraceae_match <- grepl("Oxalobacteraceae", p$data$label, fixed = TRUE)
  p$data$Sphingomonadaceae_match <- grepl("Sphingomonadaceae", p$data$label, fixed = TRUE)
  
  # Pathogens
  p$data$Wolbachia_match <- grepl("Wolbachia", p$data$label, fixed = TRUE)
  p$data$Erwinia_match <- grepl("Erwinia", p$data$label, fixed = TRUE)
  p$data$Hafnia_match <- grepl("Hafnia", p$data$label, fixed = TRUE)
  
  
  # Define your color vector
  taxa_colors <- c(
    Acetobacteraceae_match = '#e6194B',
    Lactobacillaceae_match = '#3cb44b',
    
    Bartonella_match = '#ffe119',
    Bifidobacteriaceae_match = '#4363d8',
    Neisseriaceae_match = '#f58231',
    Orbaceae_match = '#911eb4',
    
    Bacillaceae_match = '#42d4f4',
    Burkholderiaceae_match = '#f032e6',
    Clostridiaceae_match = '#bfef45',
    Comamonadaceae_match = '#fabed4',
    Enterobacteriaceae_match = '#469990',
    Methylobacteriaceae_match = '#dcbeff',
    Moraxellaceae_match = '#9A6324',
    Oxalobacteraceae_match = '#800000',
    Sphingomonadaceae_match = '#aaffc3',
    
    Wolbachia_match = '#808000',
    Erwinia_match = '#ffd8b1',
    Hafnia_match = '#000075'
  )
  
  
  # Plot
  p2 <- p +
    new_scale_fill() +
    coord_cartesian(clip="off") +
    geom_fruit(
      data=features_site_metadata,
      geom=geom_tile,
      pwidth=0.05,
      offset=0.1,
      mapping=aes(y=bacteria,
                  fill=percent_sites, width=0.1),
      show.legend=FALSE) +
    scale_fill_gradient(high = "black", low ="lightgrey") +
    new_scale_fill() +
    geom_fruit(
      data=features_site_metadata,
      geom=geom_tile,
      pwidth=0.05,
      offset=0.085,
      mapping=aes(y=bacteria,
                  fill=obligate_color, width=0.1),
      show.legend=FALSE) +
    scale_fill_identity() +
    new_scale_fill() +
    
    # BothObligate taxa
    geom_tippoint(aes(subset = Acetobacteraceae_match),
                  pch = 21, fill = taxa_colors["Acetobacteraceae_match"], size = 4) +
    geom_tippoint(aes(subset = Lactobacillaceae_match),
                  pch = 21, fill = taxa_colors["Lactobacillaceae_match"], size = 4) +
    
    # SocialObligate taxa
    geom_tippoint(aes(subset = Bartonella_match),
                  pch = 22, fill = taxa_colors["Bartonella_match"], size = 4) +
    geom_tippoint(aes(subset = Bifidobacteriaceae_match),
                  pch = 22, fill = taxa_colors["Bifidobacteriaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Neisseriaceae_match),
                  pch = 22, fill = taxa_colors["Neisseriaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Orbaceae_match),
                  pch = 22, fill = taxa_colors["Orbaceae_match"], size = 4) +
    
    # SolitaryObligate taxa
    geom_tippoint(aes(subset = Bacillaceae_match),
                  pch = 24, fill = taxa_colors["Bacillaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Burkholderiaceae_match),
                  pch = 24, fill = taxa_colors["Burkholderiaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Clostridiaceae_match),
                  pch = 24, fill = taxa_colors["Clostridiaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Comamonadaceae_match),
                  pch = 24, fill = taxa_colors["Comamonadaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Enterobacteriaceae_match),
                  pch = 24, fill = taxa_colors["Enterobacteriaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Methylobacteriaceae_match),
                  pch = 24, fill = taxa_colors["Methylobacteriaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Moraxellaceae_match),
                  pch = 24, fill = taxa_colors["Moraxellaceae_match"], size = 4) +
    geom_tippoint(aes(subset = Oxalobacteraceae_match),
                  pch = 24, fill = taxa_colors["Oxalobacteraceae_match"], size = 4) +
    geom_tippoint(aes(subset = Sphingomonadaceae_match),
                  pch = 24, fill = taxa_colors["Sphingomonadaceae_match"], size = 4) +
    
    # Pathogen taxa
    geom_tippoint(aes(subset = Wolbachia_match),
                  pch = 23, fill = taxa_colors["Wolbachia_match"], size = 4) +
    geom_tippoint(aes(subset = Erwinia_match),
                  pch = 23, fill = taxa_colors["Erwinia_match"], size = 4) +
    geom_tippoint(aes(subset = Hafnia_match),
                  pch = 23, fill = taxa_colors["Hafnia_match"], size = 4)
  
  
  
  ## list [[1]] is tree, [[2]] is metadata, [[3]] is tip.order

  
  return(list(p2, features_site_metadata, tip.order))
  
  
}

