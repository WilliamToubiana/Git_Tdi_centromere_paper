library(ggplot2)
library(cowplot)
library(hash)
library(stringr)
library(modeest)
library(gtools)

##### CODE
### merges df by the first column
merge_df_by_col <- function(match_vec,match_rep_name){
  
  ### so merge requires that the columns NOT to match on have different names or it tries to match on those too
  ### and even if I try to specify by.y and by.x it still doesn't work
  ### so I will bang a prefix on all but the first (matching) column
  ### add first df and replace names
  
  A1 = match_vec[1]
  A1 = eval(parse(text=A1))
  match_coln <- colnames(A1)[1]
  
  new_col_n <- c()
  A1_col <- colnames(A1)
  for(el in A1_col){
    match_pre = match_rep_name[1]
    new_names <- paste(match_pre,el,sep = "")
    new_col_n = c(new_col_n,new_names )
  }
  
  colnames(A1) <- new_col_n
  colnames(A1)[1] <- match_coln
  df_merged <- A1
  
  for(i in 1:(length(match_vec))){
    curr_df_name <- match_vec[i]
    curr_df <- eval(parse(text=curr_df_name))
    match_coln <- colnames(curr_df)[1]
    curr_df_col <- colnames(curr_df)
    new_col_n <- c()
    
    for(el in curr_df_col){
      match_pre = match_rep_name[i]
      new_names <- paste(match_pre,el,sep = "")
      new_col_n = c(new_col_n,new_names )
      #print(new_names)
    }
    colnames(curr_df) <- new_col_n
    colnames(curr_df)[1] <- match_coln
    df_merged = merge(df_merged,curr_df)
  }
  return(df_merged)
}


plot_cov_along_chr <- function(df, chr_want){
  df_chr <- subset(df, df$scaf_name == chr_want)
  print(head(df_chr))
  
  p1 <- ggplot(df_chr, aes(x=mid_pos , cov)) +
    theme_bw() +
    geom_line() +
    ggtitle(paste("scaf", chr_want)) +
    xlab("Position")
  
  return(p1)	
}

peakfinder <- function(d){
  dh <- hist(d,plot=FALSE, breaks=1000)
  ins <- dh[["counts"]]
  nbins <- length(ins)
  ss <- which(rank(ins)%in%seq(from=nbins,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}

## data Tdi 
### Notes
### good genome. X == Tdi_LRv5a_scf3
### all males look right


setwd("output")
Tdi_F_ReSeq_Di02         <- read.table("Tdi_F_ReSeq_Di02_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di02$mid_pos <-  (Tdi_F_ReSeq_Di02$end + Tdi_F_ReSeq_Di02$start) / 2
Tdi_F_ReSeq_Di04         <- read.table("Tdi_F_ReSeq_Di04_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di04$mid_pos <-  (Tdi_F_ReSeq_Di04$end + Tdi_F_ReSeq_Di04$start) / 2
Tdi_F_ReSeq_Di06         <- read.table("Tdi_F_ReSeq_Di06_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di06$mid_pos <-  (Tdi_F_ReSeq_Di06$end + Tdi_F_ReSeq_Di06$start) / 2
Tdi_F_ReSeq_Di08         <- read.table("Tdi_F_ReSeq_Di08_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di08$mid_pos <-  (Tdi_F_ReSeq_Di08$end + Tdi_F_ReSeq_Di08$start) / 2
Tdi_F_ReSeq_Di10         <- read.table("Tdi_F_ReSeq_Di10_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di10$mid_pos <-  (Tdi_F_ReSeq_Di10$end + Tdi_F_ReSeq_Di10$start) / 2

Tdi_M_18_3997         <- read.table("Tdi_M_18-3997_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_M_18_3997$mid_pos <-  (Tdi_M_18_3997$end + Tdi_M_18_3997$start) / 2
Tdi_M_18_3998         <- read.table("Tdi_M_18-3998_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_M_18_3998$mid_pos <-  (Tdi_M_18_3998$end + Tdi_M_18_3998$start) / 2

### add some info
Tdi_F_ReSeq_Di02$samp_name <- rep("Di02", length(Tdi_F_ReSeq_Di02[,1]))
Tdi_F_ReSeq_Di04$samp_name <- rep("Di04", length(Tdi_F_ReSeq_Di04[,1]))
Tdi_F_ReSeq_Di06$samp_name <- rep("Di06", length(Tdi_F_ReSeq_Di06[,1]))
Tdi_F_ReSeq_Di08$samp_name <- rep("Di08", length(Tdi_F_ReSeq_Di08[,1]))
Tdi_F_ReSeq_Di10$samp_name <- rep("Di10", length(Tdi_F_ReSeq_Di10[,1]))
Tdi_M_18_3997$samp_name    <- rep("18_3997", length(Tdi_M_18_3997[,1]))
Tdi_M_18_3998$samp_name    <- rep("18_3998", length(Tdi_M_18_3998[,1]))

Tdi_F_ReSeq_Di02$sex <- rep("F", length(Tdi_F_ReSeq_Di02[,1]))
Tdi_F_ReSeq_Di04$sex <- rep("F", length(Tdi_F_ReSeq_Di04[,1]))
Tdi_F_ReSeq_Di06$sex <- rep("F", length(Tdi_F_ReSeq_Di06[,1]))
Tdi_F_ReSeq_Di08$sex <- rep("F", length(Tdi_F_ReSeq_Di08[,1]))
Tdi_F_ReSeq_Di10$sex <- rep("F", length(Tdi_F_ReSeq_Di10[,1]))
Tdi_M_18_3997$sex    <- rep("M", length(Tdi_M_18_3997[,1]))
Tdi_M_18_3998$sex    <- rep("M", length(Tdi_M_18_3998[,1]))


#############################################################################################
### chr lens

Tdi_chr_lens <- aggregate(end~scaf_name, Tdi_F_ReSeq_Di02, FUN=max)
head(Tdi_chr_lens[order(Tdi_chr_lens$end,decreasing=T),], n = 20)
### big drop after 12


plot_all_chr_Tdi <- function(df){
  out_plota <- plot_grid(
    plot_cov_along_chr(df, "Tdi_LRv5a_scf1"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf2"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf3"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf4"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf5"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf6"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf7"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf8"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf9"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf10"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf11"),
    plot_cov_along_chr(df, "Tdi_LRv5a_scf12"), ncol = 1)
  
  return(out_plota)
  
}

png(filename = "Tdi_F_ReSeq_Di02_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di02)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_F_ReSeq_Di04_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di04)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_F_ReSeq_Di06_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di06)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_F_ReSeq_Di08_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di08)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_F_ReSeq_Di10_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di10)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_M_18_3997_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_M_18_3997 )
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_M_18_3998_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_M_18_3998 )
dev.off()
getwd() ## where has my plot gone....



##### median cov window

want_scafs_Tdi <- c(
  "Tdi_LRv5a_scf1",
  "Tdi_LRv5a_scf2",
  "Tdi_LRv5a_scf3",
  "Tdi_LRv5a_scf4",
  "Tdi_LRv5a_scf5",
  "Tdi_LRv5a_scf6",
  "Tdi_LRv5a_scf7",
  "Tdi_LRv5a_scf8",
  "Tdi_LRv5a_scf9",
  "Tdi_LRv5a_scf10",
  "Tdi_LRv5a_scf11",
  "Tdi_LRv5a_scf12")




Tdi_F_ReSeq_Di02_ws     <- Tdi_F_ReSeq_Di02[Tdi_F_ReSeq_Di02$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di02_ws_med <- aggregate(Tdi_F_ReSeq_Di02_ws$cov, list(Tdi_F_ReSeq_Di02_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di02_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di02_med_window_cov")
Tdi_F_ReSeq_Di04_ws     <- Tdi_F_ReSeq_Di04[Tdi_F_ReSeq_Di04$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di04_ws_med <- aggregate(Tdi_F_ReSeq_Di04_ws$cov, list(Tdi_F_ReSeq_Di04_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di04_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di04_med_window_cov")
Tdi_F_ReSeq_Di06_ws     <- Tdi_F_ReSeq_Di06[Tdi_F_ReSeq_Di06$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di06_ws_med <- aggregate(Tdi_F_ReSeq_Di06_ws$cov, list(Tdi_F_ReSeq_Di06_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di06_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di06_med_window_cov")
Tdi_F_ReSeq_Di08_ws     <- Tdi_F_ReSeq_Di08[Tdi_F_ReSeq_Di08$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di08_ws_med <- aggregate(Tdi_F_ReSeq_Di08_ws$cov, list(Tdi_F_ReSeq_Di08_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di08_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di08_med_window_cov")
Tdi_F_ReSeq_Di10_ws     <- Tdi_F_ReSeq_Di10[Tdi_F_ReSeq_Di10$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di10_ws_med <- aggregate(Tdi_F_ReSeq_Di10_ws$cov, list(Tdi_F_ReSeq_Di10_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di10_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di10_med_window_cov")
Tdi_M_18_3997_ws     <- Tdi_M_18_3997[Tdi_M_18_3997$scaf_name %in% want_scafs_Tdi,]
Tdi_M_18_3997_ws_med <- aggregate(Tdi_M_18_3997_ws$cov, list(Tdi_M_18_3997_ws$scaf_name), FUN=median)
colnames(Tdi_M_18_3997_ws_med) <- c("scaf", "Tdi_M_18_3997_med_window_cov")
Tdi_M_18_3998_ws     <- Tdi_M_18_3998[Tdi_M_18_3998$scaf_name %in% want_scafs_Tdi,]
Tdi_M_18_3998_ws_med <- aggregate(Tdi_M_18_3998_ws$cov, list(Tdi_M_18_3998_ws$scaf_name), FUN=median)
colnames(Tdi_M_18_3998_ws_med) <- c("scaf", "Tdi_M_18_3998_med_window_cov")

match_dfs <- c("Tdi_F_ReSeq_Di02_ws_med", "Tdi_F_ReSeq_Di04_ws_med", "Tdi_F_ReSeq_Di06_ws_med", "Tdi_F_ReSeq_Di08_ws_med", "Tdi_F_ReSeq_Di10_ws_med", "Tdi_M_18_3997_ws_med", "Tdi_M_18_3998_ws_med")
match_name <- c("","", "", "", "", "", "")
Tdi_LRv5a_med_window_cov <- merge_df_by_col(match_dfs,match_name)
write.csv(Tdi_LRv5a_med_window_cov, "Tdi_LRv5a_med_window_cov.csv")




###############################################################################
### male vs female cov

#############################################################################
### get norm cov

norm_cov <- function(df){
  modal_cov <- mlv(df$cov, method = "shorth")
  df$cov_n <- df$cov / modal_cov
  return(df)
}

Tdi_F_ReSeq_Di02 <- norm_cov(Tdi_F_ReSeq_Di02)
Tdi_F_ReSeq_Di04 <- norm_cov(Tdi_F_ReSeq_Di04)
Tdi_F_ReSeq_Di06 <- norm_cov(Tdi_F_ReSeq_Di06)
Tdi_F_ReSeq_Di08 <- norm_cov(Tdi_F_ReSeq_Di08)
Tdi_F_ReSeq_Di10 <- norm_cov(Tdi_F_ReSeq_Di10)
Tdi_M_18_3997    <- norm_cov(Tdi_M_18_3997)
Tdi_M_18_3998    <- norm_cov(Tdi_M_18_3998)
Tdi_all_norm     <- rbind(Tdi_F_ReSeq_Di02, Tdi_F_ReSeq_Di04, Tdi_F_ReSeq_Di06, Tdi_F_ReSeq_Di08, Tdi_F_ReSeq_Di10, Tdi_F_ReSeq_Di02, Tdi_M_18_3997, Tdi_M_18_3998)

head(Tdi_all_norm)
tail(Tdi_all_norm)



plot_MF_cov_all <- function(df, chr_want){
  df_chr <- subset(df, df$scaf_name == chr_want)
  print(head(df_chr))
  ggplot(df_chr, aes(x=mid_pos, cov_n, col = samp_name)) +
    theme_bw() +
    geom_line(alpha=0.5) +
    ggtitle(paste("scaf", chr_want)) + 
    #scale_colour_manual(values=c("royalblue2", "red3")) +
    xlab("Position") + ylim(0,3)
}

plot_MF_cov_all_smooth <- function(df, chr_want){
  df_chr <- subset(df, df$scaf_name == chr_want)
  print(head(df_chr))
  p1 <- ggplot(df_chr, aes(x=mid_pos, cov_n, col = samp_name)) +
    theme_bw() +
    geom_smooth(method="auto", se=FALSE, fullrange=FALSE, level=0.95) +
    ggtitle(paste("scaf", chr_want)) + 
    #scale_colour_manual(values=c("royalblue2", "red3")) +
    xlab("Position")
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
  
}

png(filename = "All_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 40, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf1"),
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf2"),
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf3"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf4"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf5"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf6"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf7"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf8"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf9"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf10"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf11"), 
  plot_MF_cov_all(Tdi_all_norm, "Tdi_LRv5a_scf12"), 
  ncol = 1)
dev.off()
getwd() ## where has my plot gone....

plot_MF_cov_all_Tdimfcolour <- function(df, chr_want){
  df_chr <- subset(df, df$scaf_name == chr_want)
  print(head(df_chr))
  ggplot(df_chr, aes(x=mid_pos, cov_n, col = samp_name)) +
    theme_bw() +
    geom_line(alpha=0.5) +
    ggtitle(paste("scaf", chr_want)) + 
    scale_colour_manual(values=c("royalblue1", "royalblue2", "red1", "red2", "red3", "red4", "pink1")) +
    xlab("Position") + ylim(0,3)
}


png(filename = "All_to_Tdi_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window_mfcolour.png", width = 20, height = 40, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf1"),
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf2"),
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf3"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf4"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf5"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf6"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf7"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf8"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf9"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf10"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf11"), 
  plot_MF_cov_all_Tdimfcolour(Tdi_all_norm, "Tdi_LRv5a_scf12"), 
  ncol = 1)
dev.off()
getwd() ## where has my plot gone....



### scaf cov window plot

cov_plot_hist_Tdi <- function(df, max_cov) {
  df1 <- df
  samp = deparse(substitute(df))
  print(samp)
  p1 <- ggplot(df1) +
    theme_bw() +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf1"),     aes(cov, color="Tdi_LRv5a_scf1"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf2"),     aes(cov, color="Tdi_LRv5a_scf2"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf3"),     aes(cov, color="Tdi_LRv5a_scf3"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf4"),     aes(cov, color="Tdi_LRv5a_scf4"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf5"),     aes(cov, color="Tdi_LRv5a_scf5"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf6"),     aes(cov, color="Tdi_LRv5a_scf6"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf7"),     aes(cov, color="Tdi_LRv5a_scf7"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf8"),     aes(cov, color="Tdi_LRv5a_scf8"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf9"),     aes(cov, color="Tdi_LRv5a_scf9"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf10"),    aes(cov, color="Tdi_LRv5a_scf10"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf11"),    aes(cov, color="Tdi_LRv5a_scf11"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tdi_LRv5a_scf12"),    aes(cov, color="Tdi_LRv5a_scf12"), size = 1, stat="density") +
    scale_color_manual(name = "Scaf", values = c("Tdi_LRv5a_scf1" = "#1B9E77",
                                                 "Tdi_LRv5a_scf2" = "#D95F02",
                                                 "Tdi_LRv5a_scf3" = "red3",
                                                 "Tdi_LRv5a_scf4" = "#A6761D",
                                                 "Tdi_LRv5a_scf5" = "yellow2",
                                                 "Tdi_LRv5a_scf6" = "#666666",
                                                 "Tdi_LRv5a_scf7" = "lightblue",
                                                 "Tdi_LRv5a_scf8" = "royalblue2",
                                                 "Tdi_LRv5a_scf9" = "darkorchid",
                                                 "Tdi_LRv5a_scf10" = "#7570B3",
                                                 "Tdi_LRv5a_scf11" = "black",
                                                 "Tdi_LRv5a_scf12" = "grey33")) + 
    
    xlim(c(0, max_cov)) + 
    ggtitle(samp) +
    xlab("Coverage") 
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}


Tdi_legend_sep <- cowplot::get_legend(cov_plot_hist_Tdi(Tdi_M_18_3997, 80)$p1)

Tdi_cov_plot <- plot_grid(
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di02, 30)$p2,
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di04, 30)$p2,
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di06, 30)$p2,
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di08, 30)$p2,
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di10, 30)$p2,
  cov_plot_hist_Tdi(Tdi_M_18_3997, 80)$p2,
  cov_plot_hist_Tdi(Tdi_M_18_3998, 80)$p2,
  ncol = 2)


png(filename = "Tdi_cov_plot.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tdi_cov_plot, Tdi_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?



################################################################################################################
#### class windows X A


Tdi_all_norm_c <- as.data.frame(cbind(Tdi_F_ReSeq_Di02$scaf_name, Tdi_F_ReSeq_Di02$start, Tdi_F_ReSeq_Di02$end, Tdi_F_ReSeq_Di02$mid_pos, 
                                      Tdi_F_ReSeq_Di02$cov_n, Tdi_F_ReSeq_Di04$cov_n, Tdi_F_ReSeq_Di06$cov_n, Tdi_F_ReSeq_Di08$cov_n, Tdi_F_ReSeq_Di10$cov_n, Tdi_M_18_3997$cov_n, Tdi_M_18_3998$cov_n))
colnames(Tdi_all_norm_c) <- c("scaf_name", "start", "end", "mid_pos", "Tdi_F_ReSeq_Di02_cov_n", "Tdi_F_ReSeq_Di04_cov_n", "Tdi_F_ReSeq_Di06_cov_n", "Tdi_F_ReSeq_Di08_cov_n", "Tdi_F_ReSeq_Di10_cov_n", "Tdi_M_18_3997_cov_n", "Tdi_M_18_3998_cov_n")

Tdi_all_norm_c$start <- as.numeric(Tdi_all_norm_c$start )
Tdi_all_norm_c$end   <- as.numeric(Tdi_all_norm_c$end )
Tdi_all_norm_c$mid_pos <- as.numeric(Tdi_all_norm_c$mid_pos)
Tdi_all_norm_c$Tdi_F_ReSeq_Di02_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di02_cov_n)
Tdi_all_norm_c$Tdi_F_ReSeq_Di04_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di04_cov_n)
Tdi_all_norm_c$Tdi_F_ReSeq_Di06_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di06_cov_n)
Tdi_all_norm_c$Tdi_F_ReSeq_Di08_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di08_cov_n)
Tdi_all_norm_c$Tdi_F_ReSeq_Di10_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di10_cov_n)
Tdi_all_norm_c$Tdi_M_18_3997_cov_n    <- as.numeric(Tdi_all_norm_c$Tdi_M_18_3997_cov_n)
Tdi_all_norm_c$Tdi_M_18_3998_cov_n    <- as.numeric(Tdi_all_norm_c$Tdi_M_18_3998_cov_n)

Tdi_all_norm_c$female_cov_n = 
  (Tdi_all_norm_c$Tdi_F_ReSeq_Di02_cov_n +
  Tdi_all_norm_c$Tdi_F_ReSeq_Di04_cov_n +
  Tdi_all_norm_c$Tdi_F_ReSeq_Di06_cov_n +
  Tdi_all_norm_c$Tdi_F_ReSeq_Di08_cov_n +
  Tdi_all_norm_c$Tdi_F_ReSeq_Di10_cov_n) / 5
  
Tdi_all_norm_c$male_cov_n = 
  (Tdi_all_norm_c$Tdi_M_18_3997_cov_n + Tdi_all_norm_c$Tdi_M_18_3998_cov_n) / 2

Tdi_all_norm_c$MF <- log2(Tdi_all_norm_c$male_cov_n / Tdi_all_norm_c$female_cov_n)

Tdi_all_norm_c_cut <- subset(Tdi_all_norm_c, Tdi_all_norm_c$MF < -0.5) 
Tdi_Sex_chr_peak <- peakfinder(Tdi_all_norm_c_cut$MF)
Tdi_Auto_peak    <- peakfinder(Tdi_all_norm_c$MF)

Tdi_Sex_chr_peak
Tdi_Auto_peak    


Tdi_all_norm_c$XA_s <- ifelse(Tdi_all_norm_c$MF < Tdi_Sex_chr_peak + 0.1 & Tdi_all_norm_c$MF > Tdi_Sex_chr_peak - 0.1, "X", "A")
Tdi_all_norm_c$XA_l <- ifelse(Tdi_all_norm_c$MF < Tdi_Auto_peak - 0.5, "X", "A")

write.csv(Tdi_all_norm_c, "Tdi_all_norm_c_windows_100000.csv", row.names=FALSE)

png(filename = "Tdi_all_norm_c_windows_100000_hist.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tdi_all_norm_c$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tdi_Sex_chr_peak, col='red', lwd=1)
abline(v=Tdi_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....

###########################################################################################################################################################
##### plot MF by scaf dotplot


### need to get offset first

Tdi_chr_lens_s <- Tdi_chr_lens[gtools::mixedorder(Tdi_chr_lens$scaf_name), ]
head(Tdi_chr_lens_s)

offset_v <- c(0)
c_total <- 0
for(i in seq(1, length(Tdi_chr_lens_s[,1]) -1)){
  c_len <- Tdi_chr_lens_s[i,]$end
  print(c_len)
  c_total <- c_total + c_len
  offset_v <- c(offset_v, c_total)
  
}

Tdi_chr_lens_s$offset <- offset_v

## add to dict

Tdi_offset_dict <- hash()
for(i in seq(1:length(Tdi_chr_lens_s[,1]))){
  chr_n <- Tdi_chr_lens_s$scaf_name[i]
  offset_n <- Tdi_chr_lens_s$offset[i]
  Tdi_offset_dict [[chr_n]] <-  offset_n
}

## add to cov data
Tdi_chr_offset <- c()
for(i in seq(1:length(Tdi_all_norm_c[,1]))){
  scaf_n <- Tdi_all_norm_c$scaf_name[i]
  off_n <- Tdi_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tdi_chr_offset <- c(Tdi_chr_offset, off_n)
}

Tdi_all_norm_c <- cbind(Tdi_all_norm_c, Tdi_chr_offset)
Tdi_all_norm_c$genome_pos <- Tdi_all_norm_c$mid_pos + Tdi_all_norm_c$Tdi_chr_offset

Tdi_all_norm_c$scaf_class_1 <- ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf1", Tdi_all_norm_c$scaf_name,
                                      ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf2", Tdi_all_norm_c$scaf_name,
                                             ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf3", Tdi_all_norm_c$scaf_name,
                                                    ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf4", Tdi_all_norm_c$scaf_name,
                                                           ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf5", Tdi_all_norm_c$scaf_name,
                                                                  ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf6", Tdi_all_norm_c$scaf_name,
                                                                         ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf7", Tdi_all_norm_c$scaf_name,
                                                                                ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf8", Tdi_all_norm_c$scaf_name,
                                                                                       ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf9", Tdi_all_norm_c$scaf_name,
                                                                                              ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf10", Tdi_all_norm_c$scaf_name,
                                                                                                     ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf11", Tdi_all_norm_c$scaf_name,
                                                                                                            ifelse(Tdi_all_norm_c$scaf_name == "Tdi_LRv5a_scf12", Tdi_all_norm_c$scaf_name,"other"))))))))))))

Tdi_all_norm_c$scaf_class_1o <- ordered(Tdi_all_norm_c$scaf_class_1, levels= c("Tdi_LRv5a_scf1", "Tdi_LRv5a_scf2", "Tdi_LRv5a_scf3", "Tdi_LRv5a_scf4", "Tdi_LRv5a_scf5", "Tdi_LRv5a_scf6", "Tdi_LRv5a_scf7", "Tdi_LRv5a_scf8", "Tdi_LRv5a_scf9", "Tdi_LRv5a_scf10","Tdi_LRv5a_scf11", "Tdi_LRv5a_scf12", "other"))


head(Tdi_all_norm_c)
tail(Tdi_all_norm_c)


## all
ggplot(Tdi_all_norm_c, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2))

## drop 'other'


Tdi_all_norm_c_LGs <- subset(Tdi_all_norm_c, Tdi_all_norm_c$scaf_class_1 != "other")



png(filename = "Tdi_cov_dotplot_LGs.png", width  = 12, height = 7, units = "in", bg = "white", res = 300)
plot_grid(Tdi_cov_plot, Tdi_legend_sep, rel_widths  = c(1, 0.5))
ggplot(Tdi_all_norm_c_LGs, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none")
dev.off()
getwd() ## where has my plot gone....?




