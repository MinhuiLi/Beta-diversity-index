setwd("D:/000Document/200lmhworkingdirectory/20240220urban")
##vegan package vegdist() calculates community similarity/distance, details ?vegdist
# Below similarity (S) and distance metrics (D) conversion formulas, are all using S = 1-D
library(vegan)
library(pheatmap)
library(RColorBrewer) # Heatmap color sets 
library(phyloseq)

#loading species data
site <- read.csv('park_demo.csv', header = T, row.names = 1, sep = ',')
site <- t(site) 

###### Qualitative (binary) asymmetric similarity index: Jaccard similarity/difference, Sorensen similarity/difference, Simpson similarity/difference

#jaccard difference
jacc_dis <- vegdist(site, method = 'jaccard', binary = T)
jacc_dis
jacc <- as.matrix(jacc_dis) # Result is of dis type, if you want to output locally, you need to first convert it to matrix type
write.table(jacc, 'jaccard_distance.txt', col.names = NA, sep = '\t', quote = FALSE)
pheatmap(jacc, color = colorRampPalette(brewer.pal(7,'RdYlBu'))(100), main = "Jaccard Distance")

#jaccard similarity
jacc_sim <- 1 - jacc_dis
jacc_sim
jaccs <- as.matrix(jacc_sim)
write.table(jaccs, 'jaccard_similarity.txt', col.names = NA, sep = '\t', quote = FALSE)
pheatmap(jaccs, color = colorRampPalette(brewer.pal(7,'RdYlBu'))(100), main = "Jaccard Similarity")

###### Quantitative symmetric distance measurement: Euclidean distance, chord distance, Hellinger distance, chi-square distance

#Euclidean distance
eucl_dis <- vegdist(site, method = 'euclidean')
eucl_dis
eucl <- as.matrix(eucl_dis)
write.table(eucl, 'euclidean_distance.txt', col.names = NA, sep = '\t', quote = FALSE)
pheatmap(eucl, color = colorRampPalette(brewer.pal(7,'RdYlBu'))(100), main = "Euclidean Distance")

#Euclidean similarity
#This distance when converted to similarity needs to be normalized first, for example
eucl_dis_norm <- eucl_dis/max(eucl_dis)
eucl_sim <- 1 - eucl_dis_norm

#Chord distance
site_chord <- decostand(site, method = 'normalize')
site_chord
chord_dis <- vegdist(site_chord, method = 'euclidean')
chord_dis
chord <- as.matrix(chord_dis)
write.table(chord, 'chord_distance.txt', col.names = NA, sep = '\t', quote = FALSE)
pheatmap(chord, color = colorRampPalette(brewer.pal(7,'RdYlBu'))(100), main = "Chord Distance")

#Hellinger distance
site_hell <- decostand(site, method = 'hellinger')
hell_dis <- vegdist(site_hell, method = 'euclidean')
hell <- as.matrix(hell_dis) 
write.table(hell, 'hellinger_distance.txt', col.names = NA, sep = '\t', quote = FALSE)
pheatmap(hell, color = colorRampPalette(brewer.pal(7,'RdYlBu'))(100), main = "Hellinger Distance")

#Chi-square distance
#First do chi-square standardization, then calculate Euclidean distance
site_chi <- decostand(site, method = 'chi.square')
chi_dis <- vegdist(site_chi, 'euclidean')
chi_dis
chi <- as.matrix(chi_dis)
write.table(chi, 'chi-square_distance.txt', col.names = NA, sep = '\t', quote = FALSE)
pheatmap(chi, color = colorRampPalette(brewer.pal(7,'RdYlBu'))(100), main = "Chi-square Distance")

###### Quantitative asymmetric distance measurement: Bray-Curtis distance, Unifrac distance

#Bray-Curtis distance
bray_dis <- vegdist(site, method = 'bray')
bray_dis
bray <- as.matrix(bray_dis)
write.table(bray, 'bray-curtis_distance.txt', col.names = NA, sep = '\t', quote = FALSE)
pheatmap(bray, color = colorRampPalette(brewer.pal(7,'RdYlBu'))(100), main = "Bray-Curtis Distance")

#Bray-Curtis similarity
bray_sim <- 1 - bray_dis
bray_sim
brays <- as.matrix(bray_sim) 
write.table(brays, 'bray-curtis_similarity.txt', col.names = NA, sep = '\t', quote = FALSE)
pheatmap(brays, color = colorRampPalette(brewer.pal(7,'RdYlBu'))(100), main = "Bray-Curtis Similarity")
