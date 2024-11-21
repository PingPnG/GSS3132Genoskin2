################################################################
# create PCA plot from matrix
#Ping Hu
#11-21-2024
##################################################################
rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(cowplot)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]
#outname <- args[2]
rm(args)
#file <- "combinedDCT.newID.ddct.noUnt"
outname = file

A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)
data_matrix<- A[, 2:d[2]]
row.names(data_matrix) <- A[, 1]

################################################################################
splitname<-strsplit(colnames(data_matrix), "[._]")
Clen=length(splitname)
trt=rep("NA", Clen)
GSS=rep("NA", Clen)
dup=rep("NA", Clen)
Time=rep("NA", Clen)

for(mm in  1:Clen ){
  trt[mm]=splitname[[mm]][1]
  GSS[mm]=splitname[[mm]][3]
  dup[mm]=splitname[[mm]][4]
  Time[mm]=splitname[[mm]][2]
}

SID=paste0(trt, ".", Time, ".",GSS, ".", dup)
trtTimeGSS=paste0(trt, ".", Time, ".",GSS)
trtTime=paste0(trt, ".", Time)
trtGSS=paste0(trt, ".", GSS)
meta_data<- data.frame(Time, dup, GSS, trt, SID, trtTimeGSS, trtTime, trtGSS)

rm(trt, Time, A, Clen, d,dup, GSS, mm, SID, splitname, trtTime,trtTimeGSS, trtGSS, file)
list=ls()
list
######################################
pca_result <- prcomp(t(data_matrix), center = TRUE, scale. = TRUE)
loadings <- pca_result$scale
write.csv(loadings, file = paste0(outname,".PCAGeneImportance.csv") , row.names = TRUE)
# Create a data frame with the PCA results
pca_df <- as.data.frame(pca_result$x)
pca_df$SID <- rownames(pca_df)
combined_df <- merge(pca_df, meta_data, by = "SID", all.x = TRUE)
write.csv(combined_df, file = paste0(outname,".PCAData.csv") , row.names = TRUE)

variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

pc1_label <- paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)")
pc2_label <- paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)")


library(colorRamps)
 unique_trtTimeGSS <- unique(combined_df$trtGSS)
 color_count <- length(unique_trtTimeGSS)
 custom_colors <- setNames(
   primary.colors(color_count),
   unique_trtTimeGSS
 )
print(custom_colors)
 

 avePCAData <- combined_df %>%
   group_by(trt, Time, GSS) %>%
   summarise(
     avePC1 = mean(PC1, na.rm = TRUE),  # Calculate the average of PC1
     avePC2 = mean(PC2, na.rm = TRUE),   # Calculate the average of PC2
     se_PC1 = sd(PC1, na.rm = TRUE) / sqrt(n()),
     se_PC2 = sd(PC2, na.rm = TRUE) / sqrt(n())
   ) %>%
   ungroup() %>%
   mutate(trtTimeGSS = paste0(trt, Time, GSS)) 
 avePCAData$Time <- as.numeric(avePCAData$Time)
 avePCAData$trtGSS=paste0(avePCAData$trt, ".", avePCAData$GSS)
#############################################################################
 plot1 <- ggplot(avePCAData, aes(x = Time, y = avePC1, color = trt, group = trt)) +
   geom_line() +  # Draw line for each trt
   geom_point() +  # Draw points for each data point
   geom_errorbar(aes(ymin = avePC1 - se_PC1, ymax = avePC1 + se_PC1), width = 0.2) +  # Error bars
   facet_wrap(~ GSS) +  # Separate panels for each GSS
   labs(x = "Time", y = "Mean PCA PC1 ", title = paste(outname, "PC1:", pc1_label)) +
   theme_minimal() 
 jpeg(paste0(outname, 'PC1.jpg'), width=1800, height=1800, res=300)
 print(plot1)
 dev.off()
 
 avePCAData2 <- avePCAData %>%
   group_by(trt, GSS) %>%
   mutate(PC1_shift = avePC1 - avePC1[Time == 24]) %>%  # Shift all ddt values based on Time = 24
   ungroup()
 
 plot2 <- ggplot(avePCAData2, aes(x = Time, y = PC1_shift, color = trt, group = trt)) +
   geom_line() +  # Draw line for each trt
   geom_point() +  # Draw points for each data point
   geom_errorbar(aes(ymin = PC1_shift - se_PC1, ymax = PC1_shift + se_PC1), width = 0.2) +  # Error bars
   facet_wrap(~ GSS) +  # Separate panels for each GSS
   labs(x = "Time", y = "Mean PCA PC1 shift ", title = paste(outname, "PC1_shift24hr:", pc1_label)) +
   theme_minimal() 
 jpeg(paste0(outname, '.PC1shift24hr.jpg'), width=1800, height=1800, res=300)
 print(plot2)
 dev.off()
 
 
 plot3 <- ggplot(avePCAData, aes(x = Time, y = avePC2, color = trt, group = trt)) +
   geom_line() +  # Draw line for each trt
   geom_point() +  # Draw points for each data point
   geom_errorbar(aes(ymin = avePC2 - se_PC2, ymax = avePC2 + se_PC2), width = 0.2) +  # Error bars
   facet_wrap(~ GSS) +  # Separate panels for each GSS
   labs(x = "Time", y = "Mean PCA PC2 ", title = paste(outname, "PC2:", pc2_label)) +
   theme_minimal() 
 jpeg(paste0(outname, '.PC2.jpg'), width=1800, height=1800, res=300)
 print(plot3)
 dev.off()
 
 avePCAData3 <- avePCAData %>%
   group_by(trt, GSS) %>%
   mutate(PC2_shift = avePC2 - avePC2[Time == 24]) %>%  # Shift all ddt values based on Time = 24
   ungroup()
 
 plot4 <- ggplot(avePCAData3, aes(x = Time, y = PC2_shift, color = trt, group = trt)) +
   geom_line() +  # Draw line for each trt
   geom_point() +  # Draw points for each data point
   geom_errorbar(aes(ymin = PC2_shift - se_PC2, ymax = PC2_shift + se_PC2), width = 0.2) +  # Error bars
   facet_wrap(~ GSS) +  # Separate panels for each GSS
   labs(x = "Time", y = "Mean PCA PC2 shift ", title = paste(outname, "PC2_shift24hr:", pc2_label)) +
   theme_minimal() 
 jpeg(paste0(outname, '.PC2shift24hr.jpg'), width=1800, height=1800, res=300)
 print(plot4)
 dev.off()
 
 #############################################################
 avePCAData$Time <- factor(avePCAData$Time, levels = c(24, 48, 72, 120), ordered = TRUE)
 jpeg(paste0(outname, '.avePCA.jpg'), width=1800, height=1800, res=300)
 nudge_value_x <- (max(avePCAData$avePC1) - min(avePCAData$avePC1)) / 25
 nudge_value_y <- (max(avePCAData$avePC2) - min(avePCAData$avePC2)) / 25 
 ggplot(avePCAData, aes(x = avePC1, y = avePC2, color = trtGSS, shape = trtGSS)) +
   geom_point(size = 3) +
   scale_color_manual(values = custom_colors) +  # Apply custom colors
   xlab(pc1_label) +
   ylab(pc2_label) +
   coord_fixed() +
   theme_bw() +facet_wrap(~ GSS) + 
   theme(legend.position = "top") +  # Position legend at the top
   #stat_chull(aes(color = trtGSS), alpha = 0.05, geom = "polygon") +  # Convex hull for groups
   geom_line(aes(group = trtGSS), linetype = "solid") +  # Add a line based on 'trt'
   geom_text(data = avePCAData, aes(label = Time, x = avePC1, y = avePC2), size=2,
             nudge_x = nudge_value_x, nudge_y = - nudge_value_y ) +  # Further nudge the labels
   ggtitle(paste(outname, "average PCA"))
 dev.off()
 ############################
 combined_df$Time <- factor(combined_df$Time, levels = c(24, 48, 72, 120), ordered = TRUE)
 
jpeg(paste0(outname, '.PCA_GSS.jpg'), width=1800, height=1800, res=300)
 nudge_value_x <- (max(pca_df$PC1) - min(pca_df$PC1)) / 25
 nudge_value_y <- (max(pca_df$PC2) - min(pca_df$PC2)) / 25 
ggplot(combined_df, aes(x = PC1, y = PC2, color = trtGSS, shape=trtGSS)) +
  geom_point(size = 1) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  xlab(pc1_label) +
  ylab(pc2_label) +
  coord_fixed() +
  theme_bw() +facet_wrap(~ GSS) +
  theme(legend.position = "bottom") +  # Position legend at the bottom
  stat_chull(aes(color = trtGSS, fill = trtGSS), alpha = 0.05, geom = "polygon") +
  geom_text(data = combined_df, aes(label = Time, x = PC1, y = PC2), size=2,
            nudge_x = nudge_value_x, nudge_y = - nudge_value_y ) + ggtitle(paste(outname, "PCA"))
dev.off()

##########################################################################
jpeg(paste0(outname, '.PCA.jpg'), width=1800, height=1800, res=300)
ggplot(combined_df, aes(x = PC1, y = PC2, color = trtGSS, shape=trtGSS)) +
  geom_point(size = 1) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  xlab(pc1_label) +
  ylab(pc2_label) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "top") +  # Position legend at the bottom
  stat_chull(aes(color = trtGSS, fill = trtGSS), alpha = 0.05, geom = "polygon") +
  ggtitle(paste(outname, "PCA"))
dev.off()

