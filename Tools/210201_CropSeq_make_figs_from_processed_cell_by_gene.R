options(stringsAsFactors = FALSE)

Output_dir = "/Volumes/labs/Deshpande_Lab/_LAB_SHARE/Karina/HOX_Project/CROPSeq/ShendureLab/PlottingResults/"
sample="316"
params = "log10_violin_gene_subset"
pseudocount = .1

## upload cell by gene table
processed_cell_by_gene_name = paste0(Output_dir,"JB_",sample,"_processed_cell_by_gene_tbl")
processed_cell_by_gene = read.table(processed_cell_by_gene_name, sep = "\t", header = TRUE)

## upload general gRNA info
general_gRNA_file_name = paste0(Output_dir,"library_complex_gRNAs_reoriented.txt")
general_gRNA_file = read.table(general_gRNA_file_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
general_gRNA_file[which(general_gRNA_file$gene == "CCDC101"),]$gene = "SGF29"
sub_general_gRNA_file = general_gRNA_file[,c(1,3)]
colnames(sub_general_gRNA_file) = c("gRNA","Target_Gene")

## upload cell - gRNA assignment file
gRNA_file_name = paste0(Output_dir, "JB_", sample, "_cell_gRNA_assignments")
gRNA_cell_file = read.table(gRNA_file_name, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
gRNA_cell_file_wTargetGene = merge(gRNA_cell_file, sub_general_gRNA_file)

## count UMIs per cell
num_reads_per_cell = as.data.frame(rowSums(processed_cell_by_gene))
colnames(num_reads_per_cell) = "Read_count"
num_reads_per_cell$Read_count = as.numeric(num_reads_per_cell$Read_count)

## make UMI count boxplot
umi_plot = ggplot(num_reads_per_cell, aes(y=Read_count, x ="1")) + 
  geom_boxplot() + theme_classic() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("UMIs per Cell") +
  ggtitle(paste0("Sample: ", sample,"\nNum Cells: ", nrow(num_reads_per_cell), "\nMedian: ", median(num_reads_per_cell$Read_count))) +
  theme(plot.title = element_text(size = 10))

pdf(paste0(Output_dir,paste0("JB_",sample,"UMIs_per_cell.pdf"), sep = ""), width = 3, height = 5)
umi_plot
dev.off()

genes_per_cell = as.data.frame(apply(processed_cell_by_gene, 1, function(c)sum(c!=0)))
colnames(genes_per_cell) = "Gene_Count"

##make gene count boxplot
gene_count_plot = ggplot(genes_per_cell, aes(y=Gene_Count, x ="1")) + 
  geom_boxplot() + theme_classic() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Genes per Cell") +
  ggtitle(paste0("Sample: ", sample,"\nNum Cells: ", nrow(num_reads_per_cell), "\nMedian: ", median(genes_per_cell$Gene_Count))) +
  theme(plot.title = element_text(size = 10))

pdf(paste0(Output_dir,paste0("JB_",sample,"_genes_per_cell.pdf"), sep = ""), width = 3, height = 5)
gene_count_plot
dev.off()

##scale data to the median read count
median_count = median(num_reads_per_cell$Read_count)
processed_cell_by_gene_scaled = processed_cell_by_gene/num_reads_per_cell$Read_count*median_count
processed_cell_by_gene_scaled_SAVED = processed_cell_by_gene_scaled

##Add pseudocount & log
processed_cell_by_gene_scaled_log = log(processed_cell_by_gene_scaled+pseudocount,base = 10)

## Add gRNA infomation to cell-by-gene
processed_cell_by_gene_scaled_wgRNAs = processed_cell_by_gene_scaled_log[which(rownames(processed_cell_by_gene_scaled_log) %in% gRNA_cell_file$cell_bc),]
gRNA_cell_file_wTargetGene = gRNA_cell_file_wTargetGene[match(rownames(processed_cell_by_gene_scaled_wgRNAs), gRNA_cell_file_wTargetGene$cell_bc),]

processed_cell_by_gene_scaled_wgRNAs$gRNA = gRNA_cell_file_wTargetGene$gRNA
processed_cell_by_gene_scaled_wgRNAs$Target_Gene = gRNA_cell_file_wTargetGene$Target_Gene

barplot(table(processed_cell_by_gene_scaled_wgRNAs$Target_Gene), las=2)
barplot(table(processed_cell_by_gene_scaled_wgRNAs$gRNA), las=2)

Num_cells_per_target_gene = as.data.frame(table(processed_cell_by_gene_scaled_wgRNAs$Target_Gene))

genes_targeted = unique(sub_general_gRNA_file$Target_Gene)

#genes_w_likely_effects = c("GFP_gene", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10", "HOXA11", "HOXA13", "MEIS1")
genes_w_likely_effects = c("HOXA9", "HOXA10", "HOXA11", "HOXA13", "MEIS1")


## Plot Targetted Genes as boxplots (we can change these to violin plots as well if needed)
library(grid)
grid.newpage()

for(i in 1:length(genes_targeted)){
  gene_of_interest = genes_targeted[i]
  if(gene_of_interest %in% colnames(processed_cell_by_gene_scaled_wgRNAs) == TRUE){
    processed_cell_by_gene_scaled_wgRNAs$color = "not_targeted"
    processed_cell_by_gene_scaled_wgRNAs[which(processed_cell_by_gene_scaled_wgRNAs$Target_Gene == gene_of_interest),]$color = "targeted"
    ggplot1 = ggplot(processed_cell_by_gene_scaled_wgRNAs, aes_string(y=gene_of_interest, x = "Target_Gene", fill = "color")) + 
      geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab(paste0("Scaled Read Count")) + 
      ggtitle(gene_of_interest)
    if(i != length(genes_targeted)){
      ggplot1 = ggplot1 + theme(axis.title.x=element_blank())
    }
    if(i == 1){
      ggplot2 = ggplot1
    } else if(i == 2){
      test_grid = rbind(ggplotGrob(ggplot2), ggplotGrob(ggplot1), size = "last")
    } else {
      test_grid = rbind(test_grid, ggplotGrob(ggplot1), size = "last")
    }
  }
}

pdf(paste0(Output_dir,paste0("JB_",sample,"_Targeted_Genes_BoxPlot_",params,".pdf"), sep = ""), width = 7, height = 60)
grid.draw(test_grid)
dev.off()

### Plot genes of interest as violin plots

grid.newpage()

#processed_cell_by_gene_scaled_wgRNAs$Target_Gene = factor(processed_cell_by_gene_scaled_wgRNAs$Target_Gene, levels= c("AFF2","KMT2A","AFF4", "CSNK2B", "DOT1L", "BRD4", "KAT7", "CSNK2A1", "MLLT3", "SGF29", "MLLT1", "MLLT6", "MLLT10", "NTC"))
processed_cell_by_gene_scaled_wgRNAs$Target_Gene = factor(processed_cell_by_gene_scaled_wgRNAs$Target_Gene, levels= c("NTC", "DOT1L", "MLLT1", "MLLT3", "BRD4"))

for(i in 1:length(genes_w_likely_effects)){
  gene_of_interest = genes_w_likely_effects[i]
  if(gene_of_interest %in% colnames(processed_cell_by_gene_scaled_wgRNAs) == TRUE){
    #processed_cell_by_gene_scaled_wgRNAs$color = "not_targeted"
    #processed_cell_by_gene_scaled_wgRNAs[which(processed_cell_by_gene_scaled_wgRNAs$Target_Gene == gene_of_interest),]$color = "targeted"
    ggplot1 = ggplot(processed_cell_by_gene_scaled_wgRNAs, aes_string(y=gene_of_interest, x = "Target_Gene")) + 
      #geom_boxplot(fill = "#F8766D") 
      geom_violin(fill = "skyblue3") +
      geom_boxplot(width=0.1, color="black", alpha=0.6, size = .5) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab(paste0("log10(Scaled Read Count + .1)")) + 
      ggtitle(gene_of_interest)
    #geom_jitter()
    if(i != length(genes_targeted)){
      ggplot1 = ggplot1 + theme(axis.title.x=element_blank())
    }
    if(i == 1){
      ggplot2 = ggplot1
    } else if(i == 2){
      test_grid = rbind(ggplotGrob(ggplot2), ggplotGrob(ggplot1), size = "last")
    } else {
      test_grid = rbind(test_grid, ggplotGrob(ggplot1), size = "last")
    }
  }
}


pdf(paste0(Output_dir,paste0("JB_",sample,"_Genes_Of_Interest_subset2",params,".pdf"), sep = ""), width = 16, height = 40)
grid.draw(test_grid)
dev.off()

## Make gRNAs boxplots

gRNAs = general_gRNA_file$sgrna

for(i in 1:length(gRNAs)){
  gene_of_interest = gRNAs[i]
  if(gene_of_interest %in% colnames(processed_cell_by_gene_scaled_wgRNAs) == TRUE){
    #processed_cell_by_gene_scaled_wgRNAs$color = "not_targeted"
    #processed_cell_by_gene_scaled_wgRNAs[which(processed_cell_by_gene_scaled_wgRNAs$Target_Gene == gene_of_interest),]$color = "targeted"
    ggplot1 = ggplot(processed_cell_by_gene_scaled_wgRNAs, aes_string(y=gene_of_interest, x = "gRNA")) + 
      geom_boxplot(fill = "#00BFC4") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab(paste0("Scaled Read Count")) + 
      ggtitle(gene_of_interest)
    if(i != length(genes_targeted)){
      ggplot1 = ggplot1 + theme(axis.title.x=element_blank())
    }
    if(i == 1){
      ggplot2 = ggplot1
    } else if(i == 2){
      test_grid = rbind(ggplotGrob(ggplot2), ggplotGrob(ggplot1), size = "last")
    } else {
      test_grid = rbind(test_grid, ggplotGrob(ggplot1), size = "last")
    }
  }
}

pdf(paste0(Output_dir,paste0("JB_",sample,"_gRNAs_Plot_",params,".pdf"), sep = ""), width = 7, height = 100)
grid.draw(test_grid)
dev.off()


###   WORK IN PROGRESS BEYOND THIS POINT   ###

## Now let's generate a table of medians
library(dplyr)
summarized_medians = group_by(processed_cell_by_gene_scaled_wgRNAs, Target_Gene) %>%
  summarise(
    count = n(),
    median_log10 = median(MEIS1, na.rm = TRUE),
    median = 10^median(MEIS1, na.rm = TRUE)-.1
  )

table_of_medians = data.frame()
table_of_medians_horizontal = data.frame()

for(i in 1:length(genes_w_likely_effects)){
  gene_of_interest = genes_w_likely_effects[i]
  summarized_medians = group_by(processed_cell_by_gene_scaled_wgRNAs, Target_Gene) %>%
    summarise(
      count = n(),
      median_log10 = median(.data[[gene_of_interest]], na.rm = TRUE),
      median = 10^median(.data[[gene_of_interest]], na.rm = TRUE)-.1
    )
  summarized_medians = as.data.frame(summarized_medians)
  summarized_medians$gene_of_interest = gene_of_interest
  summarized_medians$NTC_median = summarized_medians[which(summarized_medians$Target_Gene == "NTC"),4]
  summarized_medians$NTC_median_log10 = log(summarized_medians$NTC_median +pseudocount, base = 10)
  table_of_medians = rbind(table_of_medians, summarized_medians)
  #table_of_medians$gene_of_int = gene_of_interest
  table_of_medians_horizontal = rbind(table_of_medians_horizontal, c(summarized_medians$median, gene_of_interest))
  colnames(table_of_medians_horizontal) = c(summarized_medians$Target_Gene, "Gene_of_Interest")
  
  #rownames(table_of_medians_horizontal) = c(rownames(table_of_medians_horizontal), gene_of_interest)
}

rownames(table_of_medians_horizontal) = table_of_medians_horizontal$Gene_of_Interest
table_of_medians_horizontal= table_of_medians_horizontal[,c(1:(ncol(table_of_medians_horizontal)-1))]
table_of_medians_vertical = t(table_of_medians_horizontal)
test = apply(table_of_medians_vertical,2, as.numeric)
rownames(test) = rownames(table_of_medians_vertical)
table_of_medians_vertical = as.data.frame(test)
#write.csv(table_of_medians,"/Volumes/labs/Deshpande_Lab/_LAB_SHARE/Karina/HOX_Project/CROPSeq/ShendureLab/PlottingResults/tableofmedians.csv",row.names = TRUE)
#table_of_medians_nonZero = table_of_medians[which(table_of_medians$NTC_median != 0),]
#table_of_medians_nonZero$fold_change = table_of_medians_nonZero$median/table_of_medians_nonZero$NTC_median
#table_of_medians_nonZero$log2_fold_change = log(table_of_medians_nonZero$fold_change, base = 2)

pseudocount = 1
## Ok, let's attempt to make a heatmap of medians ...
matrix_of_medians = as.matrix(table_of_medians_vertical)
matrix_of_medians_log10 = log(matrix_of_medians + pseudocount, base = 10)
test_dist = dist(matrix_of_medians_log10, method = "euclidean")
#test_dist = dist(matrix_of_medians, method = "euclidean")

hclust_obj = hclust(test_dist, method = "ward.D2")

plot(hclust_obj)

#Hmmm... we'll come back here ....

## Let's try to do this for all genes

all_gene_medians_tbl = data.frame()

for(i in 1:length(target_genes)){
  target_gene = target_genes[i]
  temp_medians = apply(processed_cell_by_gene_scaled_wgRNAs[which(processed_cell_by_gene_scaled_wgRNAs$Target_Gene == target_gene),c(1:(ncol(processed_cell_by_gene_scaled_wgRNAs)-32))],2,median)
  all_gene_medians_tbl = rbind(all_gene_medians_tbl,temp_medians)
}
rownames(all_gene_medians_tbl) = target_genes
colnames(all_gene_medians_tbl) = colnames(processed_cell_by_gene_scaled_wgRNAs[which(processed_cell_by_gene_scaled_wgRNAs$Target_Gene == target_gene),c(1:(ncol(processed_cell_by_gene_scaled_wgRNAs)-32))])

all_gene_medians_tbl_change_pseudocount = all_gene_medians_tbl
all_gene_medians_tbl_change_pseudocount = 10^all_gene_medians_tbl_change_pseudocount - .1
all_gene_medians_tbl_change_pseudocount = log(all_gene_medians_tbl_change_pseudocount + 1)

test_dist = dist(all_gene_medians_tbl_change_pseudocount, method = "euclidean")
hclust_obj = hclust(test_dist, method = "ward.D2")
plot(hclust_obj)





## Perform Wilcoxon Rank Sum Test (in progress)
gene_of_interest = "HOXA13"
group1_targ_gene = "KMT2A"
group2_targ_gene = "KAT7"

group1_log = processed_cell_by_gene_scaled_wgRNAs[which(processed_cell_by_gene_scaled_wgRNAs$Target_Gene == group1_targ_gene),which(colnames(processed_cell_by_gene_scaled_wgRNAs) %in% gene_of_interest)]
group2_log = processed_cell_by_gene_scaled_wgRNAs[which(processed_cell_by_gene_scaled_wgRNAs$Target_Gene == group2_targ_gene),which(colnames(processed_cell_by_gene_scaled_wgRNAs) %in% gene_of_interest)]

group1_non_log = (10^group1_log)+pseudocount
group2_non_log = (10^group2_log)+pseudocount

wilcox.test(group1_non_log, group2_non_log)$p.value
wilcox.test(group1_log, group2_log, conf.int = TRUE)$p.value

wilcox.test(group1_log, group2_log, conf.int = TRUE)

wilcox.test(group1_non_log, group2_non_log, conf.int = TRUE)


  
