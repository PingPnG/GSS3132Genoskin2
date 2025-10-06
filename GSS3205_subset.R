#############################################
#Need to add in the statistical analysis
#Need to add the ranking
#Need to add statistical figure into pptx
#############################################
rm(list=ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(cowplot)


args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]

rm(args)

file <-"GSS3205_CT.txt.ddct"
RefTime=48
TimePoints=c(48,72)

A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)
df <- A
#rownames(df) <-A[,1]
housekeeping <- c("ACTB", "GAPDH", "PPIA", "B2M")
inflammation <- c("DDIT4","dnajb9","DSG3","HBEGF","PLAUR","pmaip1","SESN2","slc30a1","stc2","TRIB3")
TissueDamage <- c("CLDN1","cxcl8","mmp10","MMP3","mt1g", "ddit3", "il13ra2", "MMP1")
Redox <- c("AKR1C1","akr1c2","HMOX1","NQO1", "osgin1", "slc7a11", "SQSTM1")
USS<- c( "DDIT4","dnajb9","DSG3","HBEGF","PLAUR","pmaip1","SESN2","slc30a1",
  "stc2","TRIB3","akr1c2","HMOX1","NQO1","osgin1","slc7a11","SQSTM1","ddit3","il13ra2",
  "MMP1","ARRDC3","asns","BCL3","btg1","cbs","chac1","fth1","GCLC","gclm",
  "gpt2","HERPUD1","KRT16","mt1x","NDRG1","pck2","pnrc1","psat1","psph","slc1a4",
  "smim14","srxn1","trim16","txnrd1","UPP1","VEGFA"
)
All <- c(
  "abca12",  "AKR1C1", "angptl4", "ARNTL2", "ARRDC3", "asns", "ass1",
  "atf3",  "BCL3", "btg1", "CCL20", "cebpb", "cers3", "chac1", "CLDN1",
  "csta", "CXCL1", "cxcl14", "cxcl8", "ddit3", "DDIT4", "dnajb9", "DSG3", "dst",
  "dusp1", "dusp10", "eppk1", "FAT4", "fgfr2", "FRAT1", "fth1", "gadd45a",
   "GCLC", "gclm", "gdf15", "gjb2", "gpt2", "HBEGF", "HERPUD1",
  "hist1h2ac", "HMOX1", "hspa1a", "hspa1b", "IL1a", "il1rl1", "il23a", "il6r",
  "irak2", "IRF7", "KRT16", "KRT6B", "KRTAP2-3", "lamp3", "MMP1", "mmp10",
  "MMP3", "mt1g", "mt1x", "NCF2", "NDRG1", "NQO1", "osgin1", "pck2", "PLAU",
  "PLAUR", "pmaip1", "pnrc1",  "PPP1R15A", "psat1", "psph", "PTGS2",
  "scnn1a", "serpine1", "serpine2", "slc1a4", "slc30a1", "slc7a11", "smim14",
  "sod2", "sprr1b", "SQSTM1", "srxn1", "stc2", "TGM1", "TIMP3", "tnc", "TNIP2",
  "TPBG", "TRIB3", "trim16", "UPP1", "VEGFA"
)
# Combine all other gene sets
other_genes <- unique(c(inflammation, TissueDamage, Redox, USS))

# Find genes in All but not in any other list
Rest <- setdiff(All, other_genes)
gene_sets <- list(Inflammation = inflammation, TissueDamage = TissueDamage, Redox = Redox, USS = USS, All=All, Rest=Rest)

for (nm in names(gene_sets)) {
  ddct_matrix <- df[(rownames(df) %in% gene_sets[[nm]]), ]
  print(nm)
  print(nrow(df))           # Before
  print(nrow(ddct_matrix))  # After

  orgID=colnames(ddct_matrix)
  #splitname<-strsplit(orgID, "[._]")
  #orgID_fixed <- gsub("2.8", "2d8", orgID)
  #orgID_fixed <- gsub("X", "", orgID_fixed)
  splitname<-strsplit(orgID, "[._]")
  Clen=length(splitname)
  trt=rep("NA", Clen)
  dup=rep("NA", Clen)
  Time=rep("NA", Clen)
  for(mm in  1:Clen ){
   trt[mm]=splitname[[mm]][1]
    dup[mm]=splitname[[mm]][3]
    Time[mm]=splitname[[mm]][2]
  }

  SID=paste0(trt, ".", Time,  ".", dup)
  trtTime=paste0(trt, ".", Time)
  colnames(ddct_matrix)=SID
  meta<- data.frame(orgID,Time,  dup, trt, SID, trtTime)

  library(officer)
  library(rvg)
  # Create PowerPoint
  doc <- read_pptx()
####################Now Let's Work on the PCA and Ranking of the material#######################
###############################################
  pca_result <- prcomp(t(ddct_matrix), center = TRUE, scale. = TRUE)
  loadings <- pca_result$scale

  write.csv(loadings, file = paste0(file,".", nm, ".ddct_PCAGeneImportance.csv") , row.names = TRUE)
  # Create a data frame with the PCA results
  pca_df <- as.data.frame(pca_result$x)
  pca_df$SID <- rownames(pca_df)
  combined_df <- merge(pca_df, meta, by = "SID", all.x = TRUE)
  write.csv(combined_df, file = paste0(file,".", nm, ".ddct_PCA.csv") , row.names = TRUE)


  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

  pc1_label <- paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)")
  pc2_label <- paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)")

  library(colorRamps)

  unique_trtTime <- unique(combined_df$trtTime)
  color_count <- length(unique_trtTime)
  custom_colors <- setNames(
  #colorRampPalette(brewer.pal(8, "Set3"))(color_count), # Generate as many colors as needed
  #paletteer_c("grDevices::rainbow", color_count) 
    primary.colors(color_count),
    unique_trtTime
  )

# Print the custom colors
#print(custom_colors)


avePCAData <- combined_df %>%
  group_by(trt, Time, trtTime) %>%
  summarise(
    avePC1 = mean(PC1, na.rm = TRUE),  # Calculate the average of PC1
    avePC2 = mean(PC2, na.rm = TRUE),   # Calculate the average of PC2
    se_PC1 = sd(PC1, na.rm = TRUE) / sqrt(n()),
    se_PC2 = sd(PC2, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup() 
avePCAData$Time <- as.numeric(avePCAData$Time)

#############################################################################

plot5 <- ggplot(avePCAData, aes(x = Time, y = avePC1, color = trt, group = trt)) +
  geom_line(linewidth = 1) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = avePC1 - se_PC1, ymax = avePC1 + se_PC1), width = 0.2) +
  labs(x = "Time", y = "Mean PCA PC1", title = paste(nm, " PC1:", pc1_label)) +
  theme(
    legend.position = "bottom",
    panel.background = element_blank(),  # Remove panel background
    plot.background = element_blank(),  # Remove plot background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
    
  ) 


jpeg(paste0(file, ".", nm, '.ddct.PC1.jpg'), width=1800, height=1800, res=300)
print(plot5)
dev.off()

avePCAData2 <- avePCAData %>%
  group_by(trt) %>%
  mutate(PC1_shift = avePC1 - avePC1[Time == RefTime]) %>%  # Shift all ddt values based on Time = 24
  ungroup()

plot52 <- ggplot(avePCAData2, aes(x = Time, y = PC1_shift, color = trt, group = trt)) +
  geom_line(linewidth=1) +  # Draw line for each trt
  geom_point(size=6) +  # Draw points for each data point
  geom_errorbar(aes(ymin = PC1_shift - se_PC1, ymax = PC1_shift + se_PC1), width = 0.2) +  # Error bars
  labs(x = "Time", y = "Mean PCA PC1 shift ", title = paste(nm, " PC1_shift:", pc1_label)) +
  theme(    legend.position = "bottom",
            panel.background = element_blank(),  # Remove panel background
            plot.background = element_blank(),  # Remove plot background
            panel.grid.major = element_blank(),  # Remove major grid lines
            panel.grid.minor = element_blank()   # Remove minor grid lines
  ) 
jpeg(paste0(file,".", nm,  '.ddct.PC1shift.jpg'), width=1800, height=1800, res=300)
print(plot52)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT PCA PC1 Ranking", location = ph_location_type(type = "title")) %>%
  ph_with(dml(ggobj = plot5), location = ph_location(left = 0.5, top = 1.8, width = 4.5, height = 4.5)) %>% 
  ph_with(dml(ggobj = plot52), location = ph_location(left = 0.5, top = 1.8, width = 4.5, height = 4.5))


plot6 <- ggplot(avePCAData, aes(x = Time, y = avePC2, color = trt, group = trt)) +
  geom_line() +  # Draw line for each trt
  geom_point() +  # Draw points for each data point
  geom_errorbar(aes(ymin = avePC2 - se_PC2, ymax = avePC2 + se_PC2), width = 0.2) +  # Error bars
  labs(x = "Time", y = "Mean PCA PC2 ", title = paste(nm, " PC2:", pc2_label)) +
  theme(    legend.position = "bottom",
            panel.background = element_blank(),  # Remove panel background
            plot.background = element_blank(),  # Remove plot background
            #panel.grid.major = element_blank(),  # Remove major grid lines
            panel.grid.minor = element_blank()   # Remove minor grid lines
  ) 
jpeg(paste0(file, ".", nm, '.ddct.PC2.jpg'), width=1800, height=1800, res=300)
print(plot6)
dev.off()

avePCAData3 <- avePCAData %>%
  group_by(trt) %>%
  mutate(PC2_shift = avePC2 - avePC2[Time == RefTime]) %>%  # Shift all ddt values based on Time = 24
  ungroup()

plot62 <- ggplot(avePCAData3, aes(x = Time, y = PC2_shift, color = trt, group = trt)) +
  geom_line() +  # Draw line for each trt
  geom_point() +  # Draw points for each data point
  geom_errorbar(aes(ymin = PC2_shift - se_PC2, ymax = PC2_shift + se_PC2), width = 0.2) +  # Error bars
  labs(x = "Time", y = "Mean PCA PC1 shift ", title = paste(nm, " PC1_shift:", pc1_label)) +
  theme(    legend.position = "bottom",
            panel.background = element_blank(),  # Remove panel background
            plot.background = element_blank(),  # Remove plot background
            #panel.grid.major = element_blank(),  # Remove major grid lines
            panel.grid.minor = element_blank()   # Remove minor grid lines
  )  
jpeg(paste0(file,".", nm, '.ddct.PC2shift.jpg'), width=1800, height=1800, res=300)
print(plot62)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT PCA PC1 Ranking", location = ph_location_type(type = "title")) %>%
  ph_with(dml(ggobj = plot6), location = ph_location(left = 0.5, top = 1.8, width = 4.5, height = 4.5)) %>% 
  ph_with(dml(ggobj = plot62), location = ph_location(left = 0.5, top = 1.8, width = 4.5, height = 4.5))

#############################################################################
#############################################################
avePCAData$Time <- factor(avePCAData$Time, levels = TimePoints, ordered = TRUE)

nudge_value_x <- (max(avePCAData$avePC1) - min(avePCAData$avePC1)) / 25
nudge_value_y <- (max(avePCAData$avePC2) - min(avePCAData$avePC2)) / 25 
p<-ggplot(avePCAData, aes(x = avePC1, y = avePC2, color = trt, shape = trt)) +
  geom_point(size = 6) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  xlab(pc1_label) +
  ylab(pc2_label) +
  coord_fixed() +
  theme_bw()  + 
  theme(legend.position = "top") +  # Position legend at the top
  geom_line(aes(group = trt), linetype = "solid") +  # Add a line based on 'trt'
  geom_text(data = avePCAData, aes(label = Time, x = avePC1, y = avePC2), size=2,
            nudge_x = nudge_value_x, nudge_y = - nudge_value_y ) +  # Further nudge the labels
  ggtitle(paste(nm, " ", file, " ddct PCA"))

img_path <- paste0(file, ".", nm, ".ddct.avePCA.jpg")

jpeg(img_path, width = 1800, height = 1800, res = 300)
print(p)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT Average PCA",
          location = ph_location_type(type = "title")) %>%
  ph_with(
    external_img(img_path, width = 8, height = 4.5),  # size in inches
    location = ph_location(left = 0.5, top = 1.8)
  )

############################
combined_df$Time <- factor(combined_df$Time, levels = TimePoints, ordered = TRUE)
nudge_value_x <- (max(pca_df$PC1) - min(pca_df$PC1)) / 25
nudge_value_y <- (max(pca_df$PC2) - min(pca_df$PC2)) / 25 
p1<-ggplot(combined_df, aes(x = PC1, y = PC2, color = trt, shape=trt)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  xlab(pc1_label) +
  ylab(pc2_label) +
  coord_fixed() +
  theme_bw()  +
  theme(legend.position = "bottom") +  # Position legend at the bottom
  stat_chull(aes(color = trt, fill = trt), alpha = 0.05, geom = "polygon") +
  geom_text(data = combined_df, aes(label = Time, x = PC1, y = PC2), size=2,
            nudge_x = nudge_value_x, nudge_y = - nudge_value_y ) + ggtitle(paste(nm, " ", file, " ddct PCA"))

img_path <- paste0(file, ".", nm, ".PCA2.jpg")

jpeg(img_path, width = 1800, height = 1800, res = 300)
print(p)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT PCA2",
          location = ph_location_type(type = "title")) %>%
  ph_with(
    external_img(img_path, width = 8, height = 4.5),  # size in inches
    location = ph_location(left = 0.5, top = 1.8)
  )


##########################################################################

p2<-ggplot(combined_df, aes(x = PC1, y = PC2, color = trt, shape=trt)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  xlab(pc1_label) +
  ylab(pc2_label) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "top") +  # Position legend at the bottom
  stat_chull(aes(color = trt, fill = trt), alpha = 0.05, geom = "polygon") +
  ggtitle(paste(nm, file, "PCA"))

img_path=paste0(file, ".", nm, '.PCA1.jpg')
jpeg(img_path, width = 1800, height = 1800, res = 300)
print(p)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT  PCA1",
          location = ph_location_type(type = "title")) %>%
  ph_with(
    external_img(img_path, width = 8, height = 4.5),  # size in inches
    location = ph_location(left = 0.5, top = 1.8)
  )
############################################################
flatten_dist_mat <- function(distance_mat) {
  myDist <- unlist(distance_mat)
  # Find minimum and maximum distances
  min_dist <- min(myDist)
  max_dist <- max(myDist)
  # Calculate similarity scores based on distances
  similarityScore <-  1 + 2 * (min_dist - myDist) / (max_dist - min_dist)
  return(similarityScore)
}
combined_df$trtTime=paste0(combined_df$trt,".", combined_df$Time)

for (myTrtTime in unique(combined_df$trtTime)){
  #myTrtTimeGSS = "Nil.24.GSS3132"
  # Calculate the reference point
  reference_point <- combined_df %>%
    filter(trtTime == myTrtTime) %>%  # Filter for control group
    summarise(
      avePC1 = mean(PC1),  # Calculate average of PC1
      avePC2 = mean(PC2)   # Calculate average of PC2
    )
  for (i in 1:dim(avePCAData)[1]) {
    point1<-avePCAData[i,c("avePC1","avePC2")]
    myDist <- sqrt(sum((point1 - reference_point)^2)) ###Eucledian distance
    avePCAData$EucDistance[i] <-myDist
  }
  avePCAData$similarityScore <- flatten_dist_mat(avePCAData$EucDistance)
  avePCAData <- avePCAData[order(avePCAData$similarityScore, decreasing = TRUE), ]
  
  for (i in 1:dim(pca_df)[1]) {
    point1<-pca_df[i,c("PC1","PC2")]
    myDist <- sqrt(sum((point1 - reference_point)^2)) ###Eucledian distance
    pca_df$EucDistance[i] <-myDist
  }
  pca_df$similarityScore <- flatten_dist_mat(pca_df$EucDistance)
  
  
  p <- ggplot(avePCAData, aes(y = reorder(trtTime, similarityScore), x = similarityScore)) +
    geom_bar(stat = "identity", aes(fill = trtTime)#, color = "black"
    ) +
    scale_fill_manual(values = custom_colors) +  # Apply the custom colors
    labs(y = "", x = "SimilarityScore", title = paste0("Rank to ", myTrtTime) )+
    theme(
      legend.position = "none",  # Remove the legend if not needed
      panel.background = element_rect(fill = "white", color = NA),  # White background
      plot.background = element_rect(fill = "white", color = NA),  # White plot area background
      panel.grid = element_blank(),  # Remove grid lines
      axis.line = element_line(color = "black")  # Add black axis lines
    )+
    ggtitle(paste(nm, file, "\n", myTrtTime, " Rank"))
  # 
  # Save the plot to a file
  file_name <- paste0(file,".",nm, ".", myTrtTime, ".jpg")
  #ggsave(file_name, plot = p, width = 1000, height = 800, res=300)
  jpeg(file_name, width=800, height=800, res=300)
  print(p)
  dev.off()
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = paste(file, myTrtTime, "Rank"), location = ph_location_type(type = "title")) %>%
    ph_with(dml(ggobj = p), location = ph_location(left = 0.5, top = 1.8, width = 8, height = 4.5))
  names(avePCAData)[names(avePCAData) == "similarityScore"] <- paste0("SimilarityScore.", myTrtTime)
  names(avePCAData)[names(avePCAData) == "EucDistance"] <- paste0("EucDistance.", myTrtTime)
  names(pca_df)[names(pca_df) == "similarityScore"] <- paste0("SimilarityScore.", myTrtTime)
  names(pca_df)[names(pca_df) == "EucDistance"] <- paste0("EucDistance.", myTrtTime)
}
library(openxlsx)
write.xlsx(avePCAData,paste0(file,".", nm, ".rank.xlsx"))
write.xlsx(pca_df,paste0(file,".", nm, ".individual.rank.xlsx"))
print(doc, target = paste0(file,".", nm, ".pptx"))

}




