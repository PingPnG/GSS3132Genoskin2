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
colnames(data_matrix)=SID
write.table(data_matrix, file = paste0(file, ".all.dct"), sep = "\t", row.names = TRUE, col.names = TRUE)
meta<- data.frame(orgID,Time, donor, dup, GSS, trt, SID, trtTime, trtTimeGSS)
data_matrix_noRNA <-data_matrix[,Time != "6"]
write.table(data_matrix_noRNA, file = paste0(file, ".dct"), sep = "\t", row.names = TRUE, col.names = TRUE)
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
  
      
      KP_dct = kruskal.test(testdata2$dct ~ testdata2$trt)$p.value
      KP_ddct = kruskal.test(testdata2$ddct ~ testdata2$trt)$p.value
      WP_WDvPS_dct = wilcox.test(testdata2$dct[testdata2$trt=="WD"], testdata2$dct[testdata2$trt=="PS"], paired=TRUE)$p.value
      TP_WDvPS_dct = t.test(testdata2$dct[testdata2$trt=="WD"], testdata2$dct[testdata2$trt=="PS"], paired=TRUE)$p.value
      WP_WDvPS_ddct = wilcox.test(testdata2$ddct[testdata2$trt=="WD"], testdata2$ddct[testdata2$trt=="PS"], paired=TRUE)$p.value
      TP_WDvPS_ddct = t.test(testdata2$ddct[testdata2$trt=="WD"], testdata2$ddct[testdata2$trt=="PS"], paired=TRUE)$p.value
      WD_PS=mean(testdata2$ddct[testdata2$trt=="WD"]) - mean(testdata2$ddct[testdata2$trt=="PS"])
      result <- paste(genename, "All", "All", WD_PS, KP_dct, WP_WDvPS_dct, TP_WDvPS_dct, KP_ddct, WP_WDvPS_ddct, TP_WDvPS_ddct,sep = ",")
      title <- paste("genename", "All", "All", "WD_PS", "KP_dct", "WP_WDvPS_dct", "TP_WDvPS_dct", "KP_ddct", "WP_WDvPS_ddct", "TP_WDvPS_ddct",sep = ",")
      for (j in  unique(testdata2$Time)){
        testdata3 <-testdata2[testdata2$Time == j,]
        d=dim(testdata3)
        if(d[1] >0){
          if(length(unique(testdata3$trt)) >1){
            KP_dct = kruskal.test(testdata3$dct ~ testdata3$trt)$p.value
            KP_ddct = kruskal.test(testdata3$ddct ~ testdata3$trt)$p.value
            WP_WDvPS_dct = wilcox.test(testdata3$dct[testdata3$trt=="WD"], testdata3$dct[testdata3$trt=="PS"], paired=TRUE)$p.value
            TP_WDvPS_dct = t.test(testdata3$dct[testdata3$trt=="WD"], testdata3$dct[testdata3$trt=="PS"], paired=TRUE)$p.value
            WP_WDvPS_ddct = wilcox.test(testdata3$ddct[testdata3$trt=="WD"], testdata3$ddct[testdata3$trt=="PS"], paired=TRUE)$p.value
            TP_WDvPS_ddct = t.test(testdata3$ddct[testdata3$trt=="WD"], testdata3$ddct[testdata3$trt=="PS"], paired=TRUE)$p.value
            WD_PS=mean(testdata3$ddct[testdata3$trt=="WD"]) - mean(testdata3$ddct[testdata3$trt=="PS"])
            result <- paste(result, j, "All", WD_PS,KP_dct, WP_WDvPS_dct, TP_WDvPS_dct, KP_ddct, WP_WDvPS_ddct, TP_WDvPS_ddct,sep = ",")
            title <- paste(title, j, "All", "WD_PS", "KP_dct", "WP_WDvPS_dct", "TP_WDvPS_dct", "KP_ddct", "WP_WDvPS_ddct", "TP_WDvPS_ddct",sep = ",")
            
          }
        }
        for (k in  unique(testdata2$GSS)){
          testdata3 <-testdata2[testdata2$Time == j & testdata2$GSS==k,]
          d=dim(testdata3)
          if(d[1] >0){
            if(length(unique(testdata3$trt)) >1){
              KP_dct = kruskal.test(testdata3$dct ~ testdata3$trt)$p.value
              KP_ddct = kruskal.test(testdata3$ddct ~ testdata3$trt)$p.value
              WP_WDvPS_dct = wilcox.test(testdata3$dct[testdata3$trt=="WD"], testdata3$dct[testdata3$trt=="PS"], paired=TRUE)$p.value
              TP_WDvPS_dct = t.test(testdata3$dct[testdata3$trt=="WD"], testdata3$dct[testdata3$trt=="PS"], paired=TRUE)$p.value
              WP_WDvPS_ddct = wilcox.test(testdata3$ddct[testdata3$trt=="WD"], testdata3$ddct[testdata3$trt=="PS"], paired=TRUE)$p.value
              TP_WDvPS_ddct = t.test(testdata3$ddct[testdata3$trt=="WD"], testdata3$ddct[testdata3$trt=="PS"], paired=TRUE)$p.value
              WD_PS=mean(testdata3$ddct[testdata3$trt=="WD"]) - mean(testdata3$ddct[testdata3$trt=="PS"])
              result <- paste(result, j, k,WD_PS, KP_dct, WP_WDvPS_dct, TP_WDvPS_dct, KP_ddct, WP_WDvPS_ddct, TP_WDvPS_ddct,sep = ",")
              title <- paste(title, j, k, "WD_PS", "KP_dct", "WP_WDvPS_dct", "TP_WDvPS_dct", "KP_ddct", "WP_WDvPS_ddct", "TP_WDvPS_ddct",sep = ",")
            }
          }
        }
      }
  
  
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
  Irritant_No=mean(Time6$dct[Time6$category == "irritant"]) - mean(Time6$dct[Time6$category == "non-irritant"])
  # Print the p-values
  #cat("T-test p-value:", t_test_result$p.value, "\n")
  #cat("Wilcoxon test p-value:", wilcox_test_result$p.value, "\n")
  
  result<-paste(result, "6hr", "RNASeq", Irritant_No, t_test_result$p.value,wilcox_test_result$p.value,sep = ",")
  title<-paste(title, "6hr", "RNASeq", "Irritant_No","t_test_pvalue","wilcox_test_pvalue",sep = ",")
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
  print(result)
  if(i==1){
    print(title)
  }
}

# Save ddct_matrix as a tab-delimited text file
#write.table(ddct_matrix, file = paste0(file, ".ddct"), sep = "\t", row.names = TRUE, col.names = TRUE)
#save(ddct_matrix, file = "ddct_matrix.RData")
