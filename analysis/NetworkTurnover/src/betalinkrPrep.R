
## this script prepares networks for analysis using the betalinkr function from betapartite

# CH <- only_obligate_network$CH
# CH <- CH[,colnames(CH)!=""]
# HM <- only_obligate_network$HM
# JC <- only_obligate_network$JC
# MM <- only_obligate_network$MM 
# MM <- MM[,colnames(MM)!=""]
# PL <- only_obligate_network$PL
# PL <- PL[,colnames(PL)!=""]
# RP <- only_obligate_network$RP
# SC <- only_obligate_network$SC 
# SM <- only_obligate_network$SM

CH <- only_obligate_network_BM$CH
CH <- CH[,colnames(CH)!=""]
HM <- only_obligate_network_BM$HM
JC <- only_obligate_network_BM$JC
MM <- only_obligate_network_BM$MM 
MM <- MM[,colnames(MM)!=""]
PL <- only_obligate_network_BM$PL
PL <- PL[,colnames(PL)!=""]
SC <- only_obligate_network_BM$SC 
SM <- only_obligate_network_BM$SM


lower.order <- "Microbes"
higher.order <- "Pollinators"


obligate_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC),
                                          partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)

#View(obligate_poll_betalink)

colnames(obligate_poll_betalink) <- c("Site1",
                                      "Site2",
                                      "DissimilaritySpeciesComposition",
                                      "OnlySharedLinks",
                                      "WholeNetworkLinks",
                                      "SpeciesTurnoverLinks",
                                      paste("TurnoverAbsence",lower.order,sep=""),
                                      paste("TurnoverAbsence",higher.order,sep=""),
                                      "TurnoverAbsenceBoth")


geo <- unique(spec.net[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site

## add column for geographic distance between sites
obligate_poll_betalink$GeoDist <- apply(obligate_poll_betalink, 1, function(x){
  geo.dist[x["Site1"],  x["Site2"]]
})

dir.create("figures", showWarnings = FALSE)
dir.create("figures/final", showWarnings = FALSE)


#### species level networks
# CH1 <- only_transient_network$CH
# CH1 <- CH1[,colnames(CH1)!=""]
# HM1 <- only_transient_network$HM
# JC1 <- only_transient_network$JC
# MM1 <- only_transient_network$MM 
# MM1 <- MM1[,colnames(MM1)!=""]
# PL1 <- only_transient_network$PL
# PL1 <- PL1[,colnames(PL1)!=""]
# RP1 <- only_transient_network$RP
# SC1 <- only_transient_network$SC 
# SM1 <- only_transient_network$SM

CH1 <- only_transient_network_BM$CH
CH1 <- CH1[,colnames(CH1)!=""]
HM1 <- only_transient_network_BM$HM
JC1 <- only_transient_network_BM$JC
MM1 <- only_transient_network_BM$MM 
MM1 <- MM1[,colnames(MM1)!=""]
PL1 <- only_transient_network_BM$PL
PL1 <- PL1[,colnames(PL1)!=""]
#RP1 <- only_transient_network_BM$RP
SC1 <- only_transient_network_BM$SC 
SM1 <- only_transient_network_BM$SM


lower.order <- "Microbes"
higher.order <- "Pollinators"


transient_poll_betalink <- betalinkr_multi(webarray = webs2array(CH1, HM1, JC1, MM1, PL1, SM1, SC1),
                                           partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)

#View(transient_poll_betalink)

colnames(transient_poll_betalink) <- c("Site1",
                                       "Site2",
                                       "DissimilaritySpeciesComposition",
                                       "OnlySharedLinks",
                                       "WholeNetworkLinks",
                                       "SpeciesTurnoverLinks",
                                       paste("TurnoverAbsence",lower.order,sep=""),
                                       paste("TurnoverAbsence",higher.order,sep=""),
                                       "TurnoverAbsenceBoth")

## change Site1 and Site2 names to match obligate turnover betalinkr
transient_poll_betalink$Site1 <- gsub("1", "", transient_poll_betalink$Site1)
transient_poll_betalink$Site2 <- gsub("1", "", transient_poll_betalink$Site2)

## add column for geographic distance between sites
transient_poll_betalink$GeoDist <- apply(transient_poll_betalink, 1, function(x){
  geo.dist[x["Site1"],  x["Site2"]]
})

