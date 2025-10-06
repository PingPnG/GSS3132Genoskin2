#############################################
#Need to add in the statistical analysis
#Need to add the ranking
#Need to add statistical figure into pptx
#############################################
rm(list=ls())

#BiocManager::install("NMF", dependencies = TRUE)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(cowplot)

truefc<-function(VVV){
  #print(VVV)
  if (is.finite(VVV )){
    XXX=VVV
    if(VVV==0){
      XXX=NA
    }else if(VVV<1){
      XXX=-1/VVV
    }
    return(XXX)
  }else{
    return("NA")
  }
}

args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]
#outname <- args[2]
rm(args)

#outname <- "GSS3196"
file <-"GSS3205_CT.txt"

A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)
CT <- A[,2:d[2]]
rownames(CT) <-A[,1]
#housekeeping <- c("ACTB", "GAPDH", "PPIA", "B2M")
housekeeping <- c("GAPDH", "PPIA", "B2M")
CT_housekeeping <- CT[rownames(CT) %in% housekeeping, , drop = FALSE]
housekeeping_col_means <- colMeans(CT_housekeeping, na.rm = TRUE)
dct <- sweep(CT, 2, housekeeping_col_means, FUN = "-")

####change the datato positive
ZZ=as.numeric(min(dct))
if(ZZ<=0){
  data_matrix <- as.matrix(dct) -ZZ+1
}else{
  data_matrix <- as.matrix(dct)
}

orgID=colnames(data_matrix)
orgID_fixed <- gsub("\\._", "_", orgID)
splitname<-strsplit(orgID_fixed, "[._]")
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
colnames(data_matrix)=SID
write.table(data_matrix, file = paste0(file, ".all.dct"), sep = "\t", row.names = TRUE, col.names = TRUE)
meta<- data.frame(orgID,Time,  dup, trt, SID, trtTime)

#save(data_matrix, file = "data_matrix.RData")
#save(meta, file = "meta.RData")
ddct_matrix <- matrix(NA, 
                       nrow = nrow(dct), 
                       ncol = ncol(dct))
rownames(ddct_matrix) <- rownames(dct)
####Use this to hold all the statistical results
big_results_list <- list()

trt_groups=unique(trtTime)
trt_pairs <- combn(trt_groups, 2, simplify = FALSE)
trt_pairs_filtered <- Filter(function(pair) {
  is_untreated_pair <- any(grepl("Untreated", pair))
  time1 <- sub(".*\\.", "", pair[1])
  time2 <- sub(".*\\.", "", pair[2])
  
  base1 <- sub("\\..*", "", pair[1])
  base2 <- sub("\\..*", "", pair[2])
  (is_untreated_pair && (time1 == time2))||(base1==base2)
  
}, trt_pairs)
trt_pairs_filtered <- lapply(trt_pairs_filtered, function(pair) {
  if (grepl("Untreated", pair[1]) && !grepl("Untreated", pair[2])) {
    return(rev(pair))  # swap positions
  } else {
    return(pair)
  }
})


library(officer)
library(rvg)
# Create PowerPoint
doc <- read_pptx()
#####This is the original dct Data######Should I do a bar plot with error bar
for (i in  1:d[1]){
  genename=rownames(dct)[i]
  testdata<-data.frame(Time,  dup, trt, SID, trtTime, orgID)
  testdata$Time <- as.numeric(testdata$Time)
  testdata$dct=as.numeric(dct[i,])
  ####As instance not match, so remove dup in here
  # Step 1: Extract and average DCT for Untreated group by Time
  unt_data <- testdata %>%
    filter(trt == "Untreated") %>%
    select(Time, dup, dct) %>%
    rename(unt_dct = dct)
  
  avg_dct_by_time <- unt_data %>%
    group_by(Time) %>%
    summarise(avg_dct = mean(unt_dct, na.rm = TRUE), .groups = "drop")
  
  # Step 2: Join average Untreated DCT back to the full dataset by Time
  testdata2 <- testdata %>%
    left_join(avg_dct_by_time, by = "Time") %>%
    mutate(ddct = dct - avg_dct)
  ###Enter into the ddct matrix
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
    group_by(trtTime, trt, Time) %>%
    summarize(
      mean_ddct = mean(ddct, na.rm = TRUE),
      se_ddct = sd(ddct, na.rm = TRUE) / sqrt(n()),
      mean_dct = mean(dct, na.rm = TRUE),
      se_dct = sd(dct, na.rm = TRUE) / sqrt(n()),
      count=n()
    )
  
  results_list <- list()   
  KP_dct = kruskal.test(testdata2$dct ~ testdata2$trt)$p.value
  KP_ddct = kruskal.test(testdata2$ddct ~ testdata2$trt)$p.value
  results_list[["genename"]] <- genename   
 
  for ( test in trt_groups) { 
    results_list[[paste0("mean_ddct.", test)]] <- summary_testdata2$mean_ddct[summary_testdata2$trtTime == test]
    results_list[[paste0("se_ddct.", test)]] <- summary_testdata2$se_ddct[summary_testdata2$trtTime == test]
    results_list[[paste0("mean_dct.", test)]] <- summary_testdata2$mean_dct[summary_testdata2$trtTime == test]
    results_list[[paste0("se_dct.", test)]] <- summary_testdata2$se_dct[summary_testdata2$trtTime == test]
    results_list[[paste0("count.", test)]] <- summary_testdata2$count[summary_testdata2$trtTime == test]
  }
  results_list[["KP_dct"]] <- KP_dct
  results_list[["KP_ddct"]] <- KP_ddct
  
  for (pair in trt_pairs_filtered)  {
    trt1 <- pair[1]
    trt2 <- pair[2]
    # Extract the gene values for each treatment group
    ddct1 <- testdata2[testdata2$trtTime == trt1, "ddct"]
    ddct2 <- testdata2[testdata2$trtTime == trt2, "ddct"]
    dct1 <- testdata2[testdata2$trtTime == trt1, "dct"]
    dct2 <- testdata2[testdata2$trtTime == trt2, "dct"]
    if (length(ddct1) > 1 && length(ddct2) > 1) {
      mDDCT1 <- mean(ddct1, na.rm = TRUE)
      mDDCT2 <- mean(ddct2, na.rm = TRUE)
      mDCT1 <- mean(dct1, na.rm = TRUE)
      mDCT2 <- mean(dct2, na.rm = TRUE)     
      fc_DDCT <- mDDCT1 / mDDCT2
      TFC_DDCT <- truefc(fc_DDCT)
      fc_DCT <- mDCT1 / mDCT2
      TFC_DCT <- truefc(fc_DCT)
      if (length(unique(ddct1)) > 1 && length(unique(ddct2)) > 1) {
        wilcox_ddct <- wilcox.test(ddct1, ddct2, paired = FALSE)$p.value
      } else {
        wilcox_ddct <- 1
      }
      if (length(unique(dct1)) > 1 && length(unique(dct2)) > 1) {
        wilcox_dct <- wilcox.test(dct1, dct2, paired = FALSE)$p.value
      } else {
        wilcox_dct <- 1
      }
      if (length(unique(ddct1)) > 1 && length(unique(ddct2)) > 1) {
        t_ddct <- t.test(ddct1, ddct2, paired = FALSE)$p.value
      } else {
        t_ddct <- 1
      }
      if (length(unique(dct1)) > 1 && length(unique(dct2)) > 1) {
        t_dct <- t.test(dct1, dct2, paired = FALSE)$p.value
      } else {
        t_dct <- 1
      }
      
      # Store results with named keys
      label <- paste0(trt1, "_vs_", trt2)
      results_list[[paste0("wilcox_ddct.", label)]] <- wilcox_ddct
      results_list[[paste0("tp_ddct.", label)]] <- t_ddct
      results_list[[paste0("wilcox_dct.", label)]] <- wilcox_dct
      results_list[[paste0("tp_dct.", label)]] <- t_dct
      results_list[[paste0("fc_ddct.", label)]] <- fc_DDCT
      results_list[[paste0("truefold_ddct.", label)]] <- TFC_DDCT
      results_list[[paste0("fc_dct.", label)]] <- fc_DCT
      results_list[[paste0("truefold_dct.", label)]] <- TFC_DCT
    }
  }
  big_results_list[[genename]] <- results_list
  
  
  plot5 <- ggplot(summary_testdata2, aes(x = Time, y = mean_ddct, color = trt, group = trt)) +
    geom_line(linewidth=3) +  # Draw line for each trt
    geom_point(size=6) +  # Draw points for each data point
    geom_errorbar(aes(ymin = mean_ddct - se_ddct, ymax = mean_ddct + se_ddct), width = 0.2) +  # Error bars
    labs(x = "Time", y = "Mean DDCT=Fold", title = paste("Fold", genename)) +
    theme_minimal()
  
  summary_testdata2 <- summary_testdata2 %>%
    group_by(trt) %>%
    mutate(ddct_shift = mean_ddct - mean_ddct[Time == 48]) %>%  # Shift all ddt values based on Time = 48
    ungroup()
  plot6 <- ggplot(summary_testdata2, aes(x = Time, y = ddct_shift, color = trt, group = trt)) +
    geom_line(linewidth=1) +  # Draw line for each trt
    geom_point(size=6) +  # Draw points for each data point
    geom_errorbar(aes(ymin = ddct_shift - se_ddct, ymax = ddct_shift + se_ddct), width = 0.2) +  # Error bars
    labs(x = "Time", y = "Shifted Mean DDCT =FOLD", title = paste("Fold", genename)) +  # Set title to genename
    theme_minimal()
  
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = paste0(genename, " DDCT"), location = ph_location_type(type = "title")) %>%
    ph_with(dml(ggobj = plot5), location = ph_location(left = 0.5, top = 1.8, width = 4.5, height = 4.5)) %>%
    # Add right plot (bxp2)
    ph_with(dml(ggobj = plot6), location = ph_location(left = 5.2, top = 1.8, width = 4.5, height = 4.5)) #%>%

  plot2 <- ggplot(summary_testdata2, aes(x = Time, y = mean_dct, color = trt, group = trt)) +
    geom_line(linewidth=1) +  # Draw line for each trt
    geom_point(size=6) +  # Draw points for each data point
    geom_errorbar(aes(ymin = mean_dct - se_dct, ymax = mean_dct + se_dct), width = 0.2) +  # Error bars
    labs(x = "Time", y = "Mean DCT=signal", title = paste("Signal", genename)) +
    theme_minimal()
  
  # Shift ddt values so that the value at Time = 48 becomes 0 for each trt and GSS group
  summary_testdata2 <- summary_testdata2 %>%
    group_by(trt) %>%
    mutate(dct_shift = mean_dct - mean_dct[Time == 48]) %>%  # Shift all ddt values based on Time = 48
    ungroup()
  
 
  # Create the line plot with error bars, pivoted on GSS
  plot3 <- ggplot(summary_testdata2, aes(x = Time, y = dct_shift, color = trt, group = trt)) +
    geom_line(linewidth=1) +  # Draw line for each trt
    geom_point(size=6) +  # Draw points for each data point
    geom_errorbar(aes(ymin = dct_shift - se_dct, ymax = dct_shift + se_dct), width = 0.2) +  # Error bars
    labs(x = "Time", y = "Shifted Mean DCT=Signal", title = paste("Signal", genename)) +  # Set title to genename
    theme_minimal()
 
  # Display the plot
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = paste0(genename, " DCT"), location = ph_location_type(type = "title")) %>%
    ph_with(dml(ggobj = plot2), location = ph_location(left = 0.5, top = 1.8, width = 4.5, height = 4.5)) %>%
    # Add right plot (bxp2)
    ph_with(dml(ggobj = plot3), location = ph_location(left = 5.2, top = 1.8, width = 4.5, height = 4.5)) #%>%
  
  # output_filename <- paste0(genename, ".jpg")
  # combined_plot <- plot_grid(plot1, plot2, plot3, plot5, plot6, ncol = 1) 
  # ggsave(output_filename, combined_plot, width = 8, height = 24, units = "in", dpi = 300, bg = "white")

}

big_results_list <- lapply(big_results_list, function(df) {
  for (col in names(df)) {
    # Convert character to numeric if needed
    if (is.character(df[[col]])) {
      df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    }
    
    # Replace NA/NaN with 0 for columns starting with "truefold"
    if (grepl("^truefold", col)) {
      df[[col]][is.na(df[[col]]) | is.nan(df[[col]])] <- 0
    }
    
    # Replace NA/NaN with 1 for columns starting with "KP", "wilcox", "tp" or "fdr"
    if (grepl("^(KP|wilcox|tp|fdr)", col)) {
      df[[col]][is.na(df[[col]]) | is.nan(df[[col]])] <- 1
    }
  }
  df
})

# Save ddct_matrix as a tab-delimited text file
write.table(ddct_matrix, file = paste0(file, ".ddct"), sep = "\t", row.names = TRUE, col.names = TRUE)
save(ddct_matrix, file = "ddct_matrix.RData")

# Convert big results list to a single data frame

big_results_table <- bind_rows(big_results_list, .id = "genename")


pData <- big_results_table %>% select(starts_with("Kruskal"), starts_with("wilcox"), starts_with("tp"))
pData[is.na(pData)] <- 1
d=dim(pData)
df<-data.frame(matrix(ncol=d[2], nrow=d[1]))
colnames(df)=paste0("fdr.",colnames(pData))

for (i in 1:d[2]) {
  p <- as.numeric(pData[[i]])       # Use [[i]] or as.numeric() to avoid list issues
  df[, i] <- p.adjust(p, method = "fdr")
}
final_result_table <- cbind(big_results_table, df)
write.table(final_result_table, file = paste0(file, ".stat"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)     

####################Now Let's Work on the PCA and Ranking of the material#######################
###############################################
pca_result <- prcomp(t(ddct_matrix), center = TRUE, scale. = TRUE)
loadings <- pca_result$scale

write.csv(loadings, file = paste0(file,".ddct_PCAGeneImportance.csv") , row.names = TRUE)
# Create a data frame with the PCA results
pca_df <- as.data.frame(pca_result$x)
pca_df$SID <- rownames(pca_df)
combined_df <- merge(pca_df, meta, by = "SID", all.x = TRUE)
write.csv(combined_df, file = paste0(file,".ddct_PCA.csv") , row.names = TRUE)

#combined_df$Time <- as.numeric(combined_df$Time)

variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

pc1_label <- paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)")
pc2_label <- paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)")

library(colorRamps)
# Extract unique values from avePCAData$trtTimeGSS
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
  geom_line(linewidth=1) +  # Draw line for each trt
  geom_point(size=6) +  # Draw points for each data point
  geom_errorbar(aes(ymin = avePC1 - se_PC1, ymax = avePC1 + se_PC1), width = 0.2) +  # Error bars
  labs(x = "Time", y = "Mean PCA PC1 ", title = paste("PC1:", pc1_label)) +
  theme_minimal() 


jpeg(paste0(file, '.ddct.PC1.jpg'), width=1800, height=1800, res=300)
print(plot5)
dev.off()

avePCAData2 <- avePCAData %>%
  group_by(trt) %>%
  mutate(PC1_shift = avePC1 - avePC1[Time == 48]) %>%  # Shift all ddt values based on Time = 24
  ungroup()

plot52 <- ggplot(avePCAData2, aes(x = Time, y = PC1_shift, color = trt, group = trt)) +
  geom_line(linewidth=1) +  # Draw line for each trt
  geom_point(size=6) +  # Draw points for each data point
  geom_errorbar(aes(ymin = PC1_shift - se_PC1, ymax = PC1_shift + se_PC1), width = 0.2) +  # Error bars
  labs(x = "Time", y = "Mean PCA PC1 shift ", title = paste("PC1_shift:", pc1_label)) +
  theme_minimal() 
jpeg(paste0(file, '.ddct.PC1shift.jpg'), width=1800, height=1800, res=300)
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
  labs(x = "Time", y = "Mean PCA PC2 ", title = paste("PC2:", pc2_label)) +
  theme_minimal() 
jpeg(paste0(file, '.ddct.PC2.jpg'), width=1800, height=1800, res=300)
print(plot6)
dev.off()

avePCAData3 <- avePCAData %>%
  group_by(trt) %>%
  mutate(PC2_shift = avePC2 - avePC2[Time == 48]) %>%  # Shift all ddt values based on Time = 24
  ungroup()

plot62 <- ggplot(avePCAData3, aes(x = Time, y = PC2_shift, color = trt, group = trt)) +
  geom_line() +  # Draw line for each trt
  geom_point() +  # Draw points for each data point
  geom_errorbar(aes(ymin = PC2_shift - se_PC2, ymax = PC2_shift + se_PC2), width = 0.2) +  # Error bars
  labs(x = "Time", y = "Mean PCA PC1 shift ", title = paste("PC1_shift:", pc1_label)) +
  theme_minimal() 
jpeg(paste0(file, '.ddct.PC2shift.jpg'), width=1800, height=1800, res=300)
print(plot62)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT PCA PC1 Ranking", location = ph_location_type(type = "title")) %>%
  ph_with(dml(ggobj = plot6), location = ph_location(left = 0.5, top = 1.8, width = 4.5, height = 4.5)) %>% 
  ph_with(dml(ggobj = plot62), location = ph_location(left = 0.5, top = 1.8, width = 4.5, height = 4.5))

#############################################################################
#############################################################
avePCAData$Time <- factor(avePCAData$Time, levels = c(48, 72), ordered = TRUE)
jpeg(paste0(file, '.ddct.avePCA.jpg'), width=1800, height=1800, res=300)
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
  ggtitle(paste(file, " ddct PCA"))
print(p)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT Average PCA", location = ph_location_type(type = "title")) %>%
  ph_with(dml(ggobj = p), location = ph_location(left = 0.5, top = 1.8, width = 8, height = 4.5))
############################
combined_df$Time <- factor(combined_df$Time, levels = c(48, 72), ordered = TRUE)

jpeg(paste0(file, '.PCA2.jpg'), width=1800, height=1800, res=300)
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
            nudge_x = nudge_value_x, nudge_y = - nudge_value_y ) + ggtitle(paste(file, " ddct PCA"))
print(p1)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT PCA", location = ph_location_type(type = "title")) %>%
  ph_with(dml(ggobj = p1), location = ph_location(left = 0.5, top = 1.8, width = 8, height = 4.5))
##########################################################################
jpeg(paste0(file, '.PCA1.jpg'), width=1800, height=1800, res=300)
p2<-ggplot(combined_df, aes(x = PC1, y = PC2, color = trt, shape=trt)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  xlab(pc1_label) +
  ylab(pc2_label) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "top") +  # Position legend at the bottom
  stat_chull(aes(color = trt, fill = trt), alpha = 0.05, geom = "polygon") +
  ggtitle(paste(file, "PCA"))
print(p2)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "All Gene DDCT PCA", location = ph_location_type(type = "title")) %>%
  ph_with(dml(ggobj = p2), location = ph_location(left = 0.5, top = 1.8, width = 8, height = 4.5))

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
  #myTrtTimeGSS = "Nil.48.GSS3132"
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
    ggtitle(paste(file, "\n", myTrtTime, " Rank"))
  # 
  # Save the plot to a file
  file_name <- paste0(file,".", myTrtTime, ".jpg")
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
write.xlsx(avePCAData,paste0(file,".rank.xlsx"))
write.xlsx(pca_df,paste0(file,".individual.rank.xlsx"))
print(doc, target = paste0(file, ".pptx"))






