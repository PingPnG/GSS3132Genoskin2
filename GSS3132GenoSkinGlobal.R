rm(list=ls())
#install.packages('bigmemory')
#BiocManager::install("NMF", dependencies = TRUE)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(cowplot)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]
outname <- args[2]
rm(args)
file <- "combinedDCT.newID.txt"
outname <- "GSS3132DCT"

A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)
heatmap_data <- A[, 2:d[2]]
row_names <- A[, 1]
row.names(heatmap_data) <- A[, 1]


ZZ=as.numeric(min(heatmap_data))
if(ZZ<=0){
  data_matrix <- as.matrix(heatmap_data) -ZZ+1
}else{
  data_matrix <- as.matrix(heatmap_data)
}

orgID=colnames(data_matrix)
splitname<-strsplit(orgID, "[._]")
Clen=length(splitname)
trt=rep("NA", Clen)
GSS=rep("NA", Clen)
dup=rep("NA", Clen)
Time=rep("NA", Clen)
donor=rep("NA", Clen)
for(mm in  1:Clen ){
  trt[mm]=splitname[[mm]][2]
  GSS[mm]=splitname[[mm]][1]
  dup[mm]=splitname[[mm]][5]
  Time[mm]=splitname[[mm]][3]
  donor[mm]=splitname[[mm]][4]
}

SID=paste0(trt, ".", Time, ".",GSS, ".", dup)
trtTime=paste0(trt, ".", Time)
trtTimeGSS=paste0(trt, ".", Time, ".", GSS)

meta<- data.frame(orgID,Time, donor, dup, GSS, trt, SID, trtTime, trtTimeGSS)
data_matrix_noRNA <-data_matrix[,Time != "6"]
meta_noRNA <-meta[Time!="6", ]
save(data_matrix_noRNA, file = "data_matrix_noRNA.RData")
save(meta_noRNA, file = "meta_noRNA.RData")
ddct_matrix <- matrix(NA, 
                       nrow = nrow(heatmap_data), 
                       ncol = ncol(heatmap_data)-6)
rownames(ddct_matrix) <- rownames(heatmap_data)
#colnames(ddct_matrix) <- colnames(heatmap_data)


#####This is the original dct Data######Should I do a bar plot with error bar
for (i in  1:d[1]){
  genename=row_names[i]
  testdata1<-data.frame(Time, donor, dup, GSS, trt, SID, trtTime, trtTimeGSS, orgID)
  testdata1$Time <- as.numeric(testdata1$Time)
  testdata1$dct=as.numeric(heatmap_data[i,])
  Time6 <- testdata1 %>%
    filter(Time == 6)
  testdata <-testdata1 %>%
    filter(Time != 6)
 
 
  unt_data <- testdata %>%
    filter(trt == "Unt") %>%
    select(Time, dup, dct, GSS) %>%
    rename(unt_dct = dct)
  
  # Join the original dataset with the "Unt" reference data to match on Time and donor
  testdata2 <- testdata %>%
    left_join(unt_data, by = c("Time", "dup", "GSS")) %>%
    mutate(ddct = dct - unt_dct) %>%  # Calculate ddct as the difference
    select(-unt_dct)  # Remove the temporary column
  #####enter into ddct matrix
  
  if(genename == rownames(ddct_matrix)[i]){
        ddct_matrix[i,] <-testdata2$ddct
        if(i==1){
          colnames(ddct_matrix) <-testdata2$SID
        }else{
          all_match <- all(colnames(ddct_matrix) %in% testdata2$SID)
          if (all_match) {
            #print("All column names of ddct_matrix match orgID in testdata2.")
          } else {
            print(paste(genename," Not all column names of ddct_matrix match orgID in testdata2."))
          }
        }
  }
  
  testdata2$Time <- as.numeric(testdata2$Time)
  summary_testdata2 <- testdata2 %>%
    group_by(Time, trt, GSS) %>%
    summarize(
      mean_ddct = mean(ddct, na.rm = TRUE),
      se_ddct = sd(ddct, na.rm = TRUE) / sqrt(n())
    )
  plot5 <- ggplot(summary_testdata2, aes(x = Time, y = mean_ddct, color = trt, group = trt)) +
    geom_line() +  # Draw line for each trt
    geom_point() +  # Draw points for each data point
    geom_errorbar(aes(ymin = mean_ddct - se_ddct, ymax = mean_ddct + se_ddct), width = 0.2) +  # Error bars
    facet_wrap(~ GSS) +  # Separate panels for each GSS
    labs(x = "Time", y = "Mean DDCT=Fold", title = paste("Fold", genename)) +
    theme_minimal()
  summary_testdata2 <- summary_testdata2 %>%
    group_by(trt, GSS) %>%
    mutate(ddct_shift = mean_ddct - mean_ddct[Time == 24]) %>%  # Shift all ddt values based on Time = 24
    ungroup()
  plot6 <- ggplot(summary_testdata2, aes(x = Time, y = ddct_shift, color = trt, group = trt)) +
    geom_line() +  # Draw line for each trt
    geom_point() +  # Draw points for each data point
    geom_errorbar(aes(ymin = ddct_shift - se_ddct, ymax = ddct_shift + se_ddct), width = 0.2) +  # Error bars
    facet_wrap(~ GSS) +  # Separate panels for each GSS
    labs(x = "Time", y = "Shifted Mean DDCT =FOLD", title = paste("Fold", genename)) +  # Set title to genename
    theme_minimal()
  
  summary_testdata1 <- testdata1 %>%
    group_by(trtTimeGSS) %>%
    summarize(
      mean_dct = mean(dct, na.rm = TRUE),
      se_dct = sd(dct, na.rm = TRUE) / sqrt(n())
    )
  
  plot1 <- ggplot(summary_testdata1, aes(x = trtTimeGSS, y = mean_dct)) +
    geom_bar(stat = "identity", fill = "skyblue") +   # Bar plot
    geom_errorbar(aes(ymin = mean_dct - se_dct, ymax = mean_dct + se_dct), width = 0.2) +
    labs(x = "trtTimeGSS", y = "Mean dct=Signal", title = paste("Signal", genename)) +
    coord_flip() +  # Make the plot horizontal
    theme_minimal()
  Time6 <- Time6 %>%
    mutate(
      category = ifelse(trt %in% c("Geraniol", "TNFa", "SLS"), "irritant", "non-irritant")
    )
  # Perform t-test
  t_test_result <- t.test(dct ~ category, data = Time6)
  
  # Perform Wilcoxon test
  wilcox_test_result <- wilcox.test(dct ~ category, data = Time6)
  
  # Print the p-values
  #cat("T-test p-value:", t_test_result$p.value, "\n")
  #cat("Wilcoxon test p-value:", wilcox_test_result$p.value, "\n")
  
  
  plot4 <- ggplot(Time6, aes(x = trt, y = dct, fill = category)) +
    geom_bar(stat = "identity") +  # Bar plot
    scale_fill_manual(values = c("irritant" = "red", "non-irritant" = "blue")) +  # Set colors for each category
    labs(x = "RNASeq_Treat", y = "dct=Signal", title = paste("Signal",genename)) +
    coord_flip() +  # Make the plot horizontal
    theme_minimal() +
    annotate("text", x = Inf, y = Inf, label = paste("T-test p-value:", round(t_test_result$p.value, 4)), 
             hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
    annotate("text", x = Inf, y = Inf, label = paste("Wilcoxon p-value:", round(wilcox_test_result$p.value, 4)), 
             hjust = 1.1, vjust = 3, size = 4, color = "black")
  
  testdata$Time <- as.numeric(testdata$Time)
  
  summary_testdata <- testdata %>%
    group_by(Time, trt, GSS) %>%
    summarize(
      mean_dct = mean(dct, na.rm = TRUE),
      se_dct = sd(dct, na.rm = TRUE) / sqrt(n())
    )
  
  plot2 <- ggplot(summary_testdata, aes(x = Time, y = mean_dct, color = trt, group = trt)) +
    geom_line() +  # Draw line for each trt
    geom_point() +  # Draw points for each data point
    geom_errorbar(aes(ymin = mean_dct - se_dct, ymax = mean_dct + se_dct), width = 0.2) +  # Error bars
    facet_wrap(~ GSS) +  # Separate panels for each GSS
    labs(x = "Time", y = "Mean DCT=signal", title = paste("Signal", genename)) +
    theme_minimal()
  
  # Shift ddt values so that the value at Time = 24 becomes 0 for each trt and GSS group
  
  
  summary_testdata <- summary_testdata %>%
    group_by(trt, GSS) %>%
    mutate(dct_shift = mean_dct - mean_dct[Time == 24]) %>%  # Shift all ddt values based on Time = 24
    ungroup()
  
 
  # Create the line plot with error bars, pivoted on GSS
  plot3 <- ggplot(summary_testdata, aes(x = Time, y = dct_shift, color = trt, group = trt)) +
    geom_line() +  # Draw line for each trt
    geom_point() +  # Draw points for each data point
    geom_errorbar(aes(ymin = dct_shift - se_dct, ymax = dct_shift + se_dct), width = 0.2) +  # Error bars
    facet_wrap(~ GSS) +  # Separate panels for each GSS
    labs(x = "Time", y = "Shifted Mean DCT=Signal", title = paste("Signal", genename)) +  # Set title to genename
    theme_minimal()
  
  # Display the plot

  
#  output_filename <- paste0(genename, ".jpg")
#  combined_plot <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 1) 
#  ggsave(output_filename, combined_plot, width = 8, height = 24, units = "in", dpi = 300, bg = "white")
  
}

# Save ddct_matrix as a tab-delimited text file
write.table(ddct_matrix, file = paste0(file, ".ddct"), sep = "\t", row.names = TRUE, col.names = TRUE)
save(ddct_matrix, file = "ddct_matrix.RData")
################################################################################
splitname<-strsplit(colnames(ddct_matrix), "[._]")
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
meta_ddct<- data.frame(Time, dup, GSS, trt, SID)
ddct_matrix_noUnt <- ddct_matrix[, trt != "Unt", drop = FALSE]
meta_ddct_noUnt <-meta_ddct[trt!="Unt", ]
# Print the filtered matrix
#print(ddct_matrix_filtered)
save(ddct_matrix_noUnt, file = "ddct_matrix_noUnt.RData")
save(meta_ddct_noUnt, file = "meta_ddct_noUnt.RData")
rm(trt, Time, A, Clen, d, donor, dup, GSS, mm, orgID, row_names, SID, splitname, ZZ, trtTime,trtTimeGSS, all_match,i,genename, 
   plot1, plot2, plot3, plot4, plot5, plot6, summary_testdata, summary_testdata1, summary_testdata2, t_test_result, 
   testdata, testdata1, testdata2, Time6, unt_data, wilcox_test_result, meta_ddct, ddct_matrix, heatmap_data)
list=ls()
list
###########################Now need to analysis the meta_ddct_noUnt


###############################################
pca_result <- prcomp(t(ddct_matrix_noUnt), center = TRUE, scale. = TRUE)
loadings <- pca_result$scale
#write.xlsx(loadings, paste0(outname,".PCAGeneImportance.xlsx"))
write.csv(loadings, file = paste0(outname,".ddct_noUnt_PCAGeneImportance.csv") , row.names = TRUE)
# Create a data frame with the PCA results
pca_df <- as.data.frame(pca_result$x)
pca_df$SID <- rownames(pca_df)
combined_df <- merge(pca_df, meta_ddct_noUnt, by = "SID", all.x = TRUE)
write.csv(combined_df, file = paste0(outname,".ddct_noUnt_PCA.csv") , row.names = TRUE)

#combined_df$Time <- as.numeric(combined_df$Time)

variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

pc1_label <- paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)")
pc2_label <- paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)")


 # custom_colors <- c(#"Untreated" = "black",  
 #                    "Nil" = "cyan", 
 #                    "PS" = "green","WD" = "red" #, 
 #                    #"Geraniol" ="hotpink", "Petrolatum" = "purple", 
 #                    #"SheaButterAAK"="purple",  "SheabutterRita"="purple", 
 #                    #"SLS" ="hotpink",           "TNFa" ="hotpink"
 # )
 #library(RColorBrewer) 
combined_df$trtGSS=paste0(combined_df$trt, ".", combined_df$GSS)

library(colorRamps)
 # Extract unique values from avePCAData$trtTimeGSS
 unique_trtTimeGSS <- unique(combined_df$trtGSS)
 
 # Generate a custom color palette with more than 21 colors
 # Using colorRampPalette to create a large set of distinct colors
 color_count <- length(unique_trtTimeGSS)
 custom_colors <- setNames(
   #colorRampPalette(brewer.pal(8, "Set3"))(color_count), # Generate as many colors as needed
   #paletteer_c("grDevices::rainbow", color_count) 
   primary.colors(color_count),
   unique_trtTimeGSS
 )
 
 # Print the custom colors
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
 plot5 <- ggplot(avePCAData, aes(x = Time, y = avePC1, color = trt, group = trt)) +
   geom_line() +  # Draw line for each trt
   geom_point() +  # Draw points for each data point
   geom_errorbar(aes(ymin = avePC1 - se_PC1, ymax = avePC1 + se_PC1), width = 0.2) +  # Error bars
   facet_wrap(~ GSS) +  # Separate panels for each GSS
   labs(x = "Time", y = "Mean PCA PC1 ", title = paste("PC1:", pc1_label)) +
   theme_minimal() 
 jpeg(paste0(outname, '.ddct_noUnt.PC1.jpg'), width=1800, height=1800, res=300)
 print(plot5)
 dev.off()
 
 avePCAData2 <- avePCAData %>%
   group_by(trt, GSS) %>%
   mutate(PC1_shift = avePC1 - avePC1[Time == 24]) %>%  # Shift all ddt values based on Time = 24
   ungroup()
 
 plot52 <- ggplot(avePCAData2, aes(x = Time, y = PC1_shift, color = trt, group = trt)) +
   geom_line() +  # Draw line for each trt
   geom_point() +  # Draw points for each data point
   geom_errorbar(aes(ymin = PC1_shift - se_PC1, ymax = PC1_shift + se_PC1), width = 0.2) +  # Error bars
   facet_wrap(~ GSS) +  # Separate panels for each GSS
   labs(x = "Time", y = "Mean PCA PC1 shift ", title = paste("PC1_shift:", pc1_label)) +
   theme_minimal() 
 jpeg(paste0(outname, '.ddct_noUnt.PC1shift.jpg'), width=1800, height=1800, res=300)
 print(plot52)
 dev.off()
 
 
 plot6 <- ggplot(avePCAData, aes(x = Time, y = avePC2, color = trt, group = trt)) +
   geom_line() +  # Draw line for each trt
   geom_point() +  # Draw points for each data point
   geom_errorbar(aes(ymin = avePC2 - se_PC2, ymax = avePC2 + se_PC2), width = 0.2) +  # Error bars
   facet_wrap(~ GSS) +  # Separate panels for each GSS
   labs(x = "Time", y = "Mean PCA PC2 ", title = paste("PC2:", pc2_label)) +
   theme_minimal() 
 jpeg(paste0(outname, '.ddct_noUnt.PC2.jpg'), width=1800, height=1800, res=300)
 print(plot6)
 dev.off()
 
 avePCAData3 <- avePCAData %>%
   group_by(trt, GSS) %>%
   mutate(PC2_shift = avePC2 - avePC2[Time == 24]) %>%  # Shift all ddt values based on Time = 24
   ungroup()
 
 plot62 <- ggplot(avePCAData3, aes(x = Time, y = PC2_shift, color = trt, group = trt)) +
   geom_line() +  # Draw line for each trt
   geom_point() +  # Draw points for each data point
   geom_errorbar(aes(ymin = PC2_shift - se_PC2, ymax = PC2_shift + se_PC2), width = 0.2) +  # Error bars
   facet_wrap(~ GSS) +  # Separate panels for each GSS
   labs(x = "Time", y = "Mean PCA PC1 shift ", title = paste("PC1_shift:", pc1_label)) +
   theme_minimal() 
 jpeg(paste0(outname, '.ddct_noUnt.PC2shift.jpg'), width=1800, height=1800, res=300)
 print(plot62)
 dev.off()
 
 #############################################################################
 #############################################################
 avePCAData$Time <- factor(avePCAData$Time, levels = c(24, 48, 72, 120), ordered = TRUE)
 jpeg(paste0(outname, '.ddct_noUnt.avePCA.jpg'), width=1800, height=1800, res=300)
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
   ggtitle(paste(outname, "PCA"))
 dev.off()
 ############################
 combined_df$Time <- factor(combined_df$Time, levels = c(24, 48, 72, 120), ordered = TRUE)
 
jpeg(paste0(outname, '.PCA2.jpg'), width=1800, height=1800, res=300)
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
jpeg(paste0(outname, '.PCA1.jpg'), width=1800, height=1800, res=300)
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
combined_df$trtTimeGSS=paste0(combined_df$trtTime,".", combined_df$GSS)
for (myTrtTimeGSS in unique(combined_df$trtTimeGSS)){
    #myTrtTimeGSS = "Nil.24.GSS3132"
    # Calculate the reference point
    reference_point <- combined_df %>%
      filter(trtTimeGSS == myTrtTimeGSS) %>%  # Filter for control group
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

      
      p <- ggplot(avePCAData, aes(y = reorder(trtTimeGSS, similarityScore), x = similarityScore)) +
        geom_bar(stat = "identity", aes(fill = trtTimeGSS)#, color = "black"
        ) +
        scale_fill_manual(values = custom_colors) +  # Apply the custom colors
        labs(y = "", x = "SimilarityScore", title = paste0("Rank to ", myTrtTimeGSS) )+
        theme(
          legend.position = "none",  # Remove the legend if not needed
          panel.background = element_rect(fill = "white", color = NA),  # White background
          plot.background = element_rect(fill = "white", color = NA),  # White plot area background
          panel.grid = element_blank(),  # Remove grid lines
          axis.line = element_line(color = "black")  # Add black axis lines
        )+
        ggtitle(paste(outname, "\n", myTrtTimeGSS, " Rank"))
      # 
      # Save the plot to a file
      file_name <- paste0(outname,".", myTrtTimeGSS, ".jpg")
      #ggsave(file_name, plot = p, width = 1000, height = 800, res=300)
      jpeg(file_name, width=800, height=800, res=300)
      print(p)
      dev.off()
      names(avePCAData)[names(avePCAData) == "similarityScore"] <- paste0("SimilarityScore.", myTrtTimeGSS)
      names(avePCAData)[names(avePCAData) == "EucDistance"] <- paste0("EucDistance.", myTrtTimeGSS)
      names(pca_df)[names(pca_df) == "similarityScore"] <- paste0("SimilarityScore.", myTrtTimeGSS)
      names(pca_df)[names(pca_df) == "EucDistance"] <- paste0("EucDistance.", myTrtTimeGSS)
}
library(openxlsx)
write.xlsx(avePCAData,paste0(outname,".rank.xlsx"))
write.xlsx(pca_df,paste0(outname,".individual.rank.xlsx"))
