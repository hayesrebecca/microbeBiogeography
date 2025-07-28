data <- spec.net

# what we need
# site to site comparison of host turnover
# columns should be Species
# rows should be sites
# values should be abundance
abund <- data %>%
    filter(WeightsMicrobe == 1) %>%
    select(Site, GenusSpecies) %>%
    group_by(Site, GenusSpecies) %>%
    arrange(Site) %>%
    filter(GenusSpecies != "") %>%
    summarize(HostAbund = n()) %>%
    pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
    mutate_all(~replace(., is.na(.), 0))
sitenames <- abund$Site
abund$Site <- NULL
rownames(abund) <- sitenames

#distance matrix of sites
geo <- data %>%
    filter(WeightsMicrobe == 1) %>%
    group_by(Site) %>%
    select(Site, Long, Lat) %>%
    distinct() %>%
    arrange(Site)
rownames(geo) <- sitenames
geo$Site <- NULL

  

#abundance data frame - bray curtis dissimilarity
dist.abund <- vegdist(abund, method = "bray")
dist.sim <- 1 - dist.abund

#geographic data frame - haversine distance in m (takes a df with lat and long and calculates dist)
d.geo <- distm(geo, fun = distHaversine)
dist.geo <- as.dist(d.geo)/1000

#abundance vs geographic mantel test
abund_geo  = mantel(dist.sim, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
print(abund_geo)

dist_decay_model <- betapart::decay.model(dist.sim,
                                          dist.geo,
                                          y.type="similarities",
                                          model.type = "exponential",
                                          perm=999)
dist_decay_model

bee.decay <- plot_decay_ggplot_single(dist_decay_model,
                                     xlab="Geographic Distance (km)",
                                     ylab="Bee Community Similarity",
                                     mod1color='navy',
                                     col = "black",
                                     lty = "solid",
                                     lwd = 1.5,
                                     cex = 1)
bee.decay
