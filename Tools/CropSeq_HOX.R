options(stringsAsFactors = FALSE)

library(ggplot2)
library(Seurat)

#install.packages("berryFunctions")
library(berryFunctions)

upload_rds_from_file = TRUE
seurat_table_file_name = "seur_tbl_200303.rds"

upload_DE_genes_from_file = TRUE
DE_gene_by_perturbation_file_name = "JB_BothRepsWgRNAs_DE_genes_by_perturbed_gene.txt"
DE_gene_by_gRNA_file_name = "JB_BothRepsWgRNAs_DE_genes_by_gRNA.txt"

Output_dir = "/modifyaccordingly/"
sample="BothRepsWgRNAs"

processed_cell_by_gene_wgRNA_REAL_genes = read.table(paste0(Output_dir,"JB_",sample,"_processed_cell_by_gene_wgRNA_REAL_genes.txt", sep = ""), sep = "\t", stringsAsFactors = FALSE)
gRNA_cell_file_wTargetGene = read.table(paste0(Output_dir,"JB_",sample,"_gRNA_cell_file_wTargetGene.txt", sep = ""), sep = "\t", stringsAsFactors = FALSE, header = TRUE)
gRNA_cell_file_wTargetGene = gRNA_cell_file_wTargetGene[match(rownames(processed_cell_by_gene_wgRNA_REAL_genes), gRNA_cell_file_wTargetGene$cell_bc),]

if(upload_rds_from_file == FALSE){
  seur_tbl = CreateSeuratObject(counts = t(processed_cell_by_gene_wgRNA_REAL_genes), min.cells = 3, min.features = 200, project = "CROPSeq_collab")
  VlnPlot(seur_tbl, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  seur_tbl <- NormalizeData(seur_tbl, normalization.method = "LogNormalize", scale.factor = 10000)
  
  seur_tbl@meta.data$cell_bc = gRNA_cell_file_wTargetGene$cell_bc
  seur_tbl@meta.data$gRNA = gRNA_cell_file_wTargetGene$gRNA
  seur_tbl@meta.data$Target_Gene = gRNA_cell_file_wTargetGene$Target_Gene
}else{
  seur_tbl = readRDS(file=paste0(Output_dir, seurat_table_file_name))
}

##### RIDGE PLOTS #####

my_levels <- rev(c("SGF29", "KAT7", "AFF2", "AFF4", "MLLT1", "MLLT3", "MLLT6", "MLLT10", "DOT1L", "KMT2A", "BRD4", "CSNK2B", "CSNK2A1", "NTC")) 
seur_tbl@meta.data$Target_Gene <- factor(x = seur_tbl@meta.data$Target_Gene, levels = my_levels)
Idents(seur_tbl) <- "Target_Gene"

features_to_plot = c("HOXA9", "HOXA13", "MEIS1", "BMI1", "LYZ")

ridge_to_save = RidgePlot(object = seur_tbl, features = features_to_plot, ncol = 1) + theme(legend.position = 'none')

pdf(paste0(Output_dir,paste0("JB_",sample,"_Updated_Ridge_Plots_210302.pdf"), sep = ""), width = 5, height = 20)
ridge_to_save
dev.off()

#####to calculate DE genes #####

# DE genes split by perturbed genes
Idents(seur_tbl) <- "Target_Gene"
perturbed_genes = unique(seur_tbl@meta.data$Target_Gene)
perturbed_genes = as.character(perturbed_genes[which(!(perturbed_genes %in% "NTC"))])

if(upload_DE_genes_from_file == FALSE){
  tbl_of_DE_genes = data.frame()
  
  for(i in 1:length(perturbed_genes)){
    temp_gene = perturbed_genes[i]
    print(temp_gene)
    if(is.error(FindMarkers(seur_tbl, ident.1 = temp_gene, ident.2 = "NTC", min.pct = 0.25)) == FALSE){ ##this deals w/ errors which appear when there are no DE genes between the perturbed gene sample and the NTC
      test.markers <- FindMarkers(seur_tbl, ident.1 = temp_gene, ident.2 = "NTC", min.pct = 0.25)
      temp_tbl = as.data.frame(cbind(rep(temp_gene, nrow(test.markers)), rep("NTC", nrow(test.markers)), rownames(test.markers)))
      colnames(temp_tbl) = c("Group_1", "Group_2", "DE_gene")
      
      combined_tbl = cbind(temp_tbl, test.markers)
      rownames(combined_tbl) = NULL
      tbl_of_DE_genes = as.data.frame(rbind(tbl_of_DE_genes, combined_tbl))
    }
  }
  tbl_of_DE_genes_by_perturbed_gene = tbl_of_DE_genes
} else {
  tbl_of_DE_genes_by_perturbed_gene = read.table(paste0(Output_dir, DE_gene_by_perturbation_file_name), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

#write.table(tbl_of_DE_genes_by_perturbed_gene, paste0(Output_dir,"JB_",sample,"_DE_genes_by_perturbed_gene.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

# DE genes split by gRNA
Idents(seur_tbl) <- "gRNA"
NTCs = c("CROPSeq_18", "CROPSeq_19")
gRNAs = as.character(unique(seur_tbl@meta.data$gRNA))
gRNAs = gRNAs[which(!(gRNAs %in% NTCs))]

if(upload_DE_genes_from_file == FALSE){
  tbl_of_DE_genes = data.frame()
  for(i in 1:length(gRNAs)){
    temp_gene = gRNAs[i]
    print(temp_gene)
    if(is.error(FindMarkers(seur_tbl, ident.1 = temp_gene, ident.2 = NTCs, min.pct = 0.25)) == FALSE){ ##this deals w/ errors which appear when there are no DE genes between the perturbed gene sample and the NTC
      test.markers <- FindMarkers(seur_tbl, ident.1 = temp_gene, ident.2 = NTCs, min.pct = 0.25)
      temp_tbl = as.data.frame(cbind(rep(temp_gene, nrow(test.markers)), rep("NTC", nrow(test.markers)), rownames(test.markers)))
      colnames(temp_tbl) = c("Group_1", "Group_2", "DE_gene")
      
      combined_tbl = cbind(temp_tbl, test.markers)
      rownames(combined_tbl) = NULL
      tbl_of_DE_genes = as.data.frame(rbind(tbl_of_DE_genes, combined_tbl))
    }
  }
  tbl_of_DE_genes_by_gRNA = tbl_of_DE_genes
} else {
  tbl_of_DE_genes_by_gRNA = read.table(paste0(Output_dir, DE_gene_by_gRNA_file_name), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

#write.table(tbl_of_DE_genes_by_gRNA, paste0(Output_dir,"JB_",sample,"_DE_genes_by_gRNA.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

##### Make UMAP ##### ---- ****
Idents(seur_tbl) <- "Target_Gene"

if(upload_rds_from_file == FALSE){
  all.genes <- rownames(seur_tbl)
  seur_tbl <- ScaleData(seur_tbl, features = all.genes)
  seur_tbl <- FindVariableFeatures(seur_tbl, selection.method = "vst", nfeatures = 2000)
  seur_tbl <- RunPCA(seur_tbl, features = VariableFeatures(object = seur_tbl))
  DimPlot(seur_tbl, reduction = "pca")
  
  ElbowPlot(seur_tbl) ## to determine number of PCs to use
  
  seur_tbl <- RunUMAP(seur_tbl, dims = 1:20)
  DimPlot(seur_tbl, reduction = "umap")
}

##Save by running:
#saveRDS(seur_tbl, file = paste0(Output_dir, "seur_tbl_200303.rds"))

### To reupload this seurat object, run this:
#seur_tbl = readRDS(file=paste0(Output_dir, "seur_tbl_200303.rds"))

#####dendrogram & heatmap #####

##scale the data so all cells have the same total read count:
num_reads_per_cell = as.data.frame(rowSums(processed_cell_by_gene_wgRNA_REAL_genes))
colnames(num_reads_per_cell) = "Read_count"
num_reads_per_cell$Read_count = as.numeric(num_reads_per_cell$Read_count)

median_count = median(num_reads_per_cell$Read_count)
processed_cell_by_gene_scaled = processed_cell_by_gene_wgRNA_REAL_genes/num_reads_per_cell$Read_count*median_count
head(rowSums(processed_cell_by_gene_scaled))

## Center & scale the data: (different normalization than done before, but yields similar results)
#subtract mean expression from each column
mean_exp = colMeans(processed_cell_by_gene_scaled)
tbl_centered = as.data.frame(t(t(processed_cell_by_gene_scaled) - mean_exp))

#and scale by standard deviation
col_stdev = apply(processed_cell_by_gene_scaled, 2, sd)
tbl_centered_sdScaled = as.data.frame(t(t(tbl_centered)/col_stdev))

tbl_centered_sdScaled$Target_Gene = gRNA_cell_file_wTargetGene$Target_Gene

##calculate the median experession for each gene in each group & make a dendrogram & heatmaps

library("pheatmap")

make_dendrogram = function(tbl_of_medians_w_Target_Gene_Col){
  all_gene_medians_tbl = data.frame()
  
  genes_targeted = as.character(unique(tbl_of_medians_w_Target_Gene_Col$Target_Gene))
  
  for(i in 1:length(genes_targeted)){
    target_gene = genes_targeted[i]
    temp_medians = apply(tbl_of_medians_w_Target_Gene_Col[which(tbl_of_medians_w_Target_Gene_Col$Target_Gene == target_gene),c(1:(ncol(tbl_of_medians_w_Target_Gene_Col)-1))],2,median)
    all_gene_medians_tbl = rbind(all_gene_medians_tbl,temp_medians)
  }
  
  rownames(all_gene_medians_tbl) = genes_targeted
  colnames(all_gene_medians_tbl) = colnames(tbl_of_medians_w_Target_Gene_Col[which(tbl_of_medians_w_Target_Gene_Col$Target_Gene == target_gene),c(1:(ncol(tbl_of_medians_w_Target_Gene_Col)-1))])
  
  print("making dendrogram")
  test_dist = dist(all_gene_medians_tbl, method = "euclidean")
  hclust_obj = hclust(test_dist, method = "ward.D2")
  
  # print("making heatmap")
  # col_sums = colSums(all_gene_medians_tbl)
  # col_sums_nonZero = col_sums[col_sums != 0]
  # all_gene_medians_tbl_nonZero = all_gene_medians_tbl[which(colnames(all_gene_medians_tbl) %in% names(col_sums_nonZero))]
  # pheatmap_obj = pheatmap(t(all_gene_medians_tbl_nonZero), fontsize = .5)
  return(list(hclust_obj, all_gene_medians_tbl))
}

##All genes
hclust_pheatmap_all_genes = make_dendrogram(tbl_centered_sdScaled)
hclust_obj_all_genes = hclust_pheatmap_all_genes[[1]]
plot(hclust_obj_all_genes)

pheatmap_obj_all_genes = hclust_pheatmap_all_genes[[2]]
col_var = apply(pheatmap_obj_all_genes, 2, var)
col_var_nonZero = col_var[col_var != 0]
pheatmap_obj_all_genes_nonZero = pheatmap_obj_all_genes[which(colnames(pheatmap_obj_all_genes) %in% names(col_var_nonZero))]
pheatmap_all_genes = pheatmap(t(pheatmap_obj_all_genes_nonZero), fontsize_row = .5, fontsize_col = 14, angle_col = 90, cutree_cols = 4)

pdf(paste0(Output_dir,paste0("JB_",sample,"_pheatmap_ALL_gene_test.pdf"), sep = ""), width = 10, height = 14)
pheatmap_all_genes
dev.off()

#plot(pheatmap_obj_all_genes)

##only include genes which are DE across samples
tbl_centered_sdScaled_justDEgenes = tbl_centered_sdScaled[,which(colnames(tbl_centered_sdScaled) %in% tbl_of_DE_genes_by_perturbed_gene$DE_gene)]
tbl_centered_sdScaled_justDEgenes$Target_Gene = gRNA_cell_file_wTargetGene$Target_Gene

hclust_pheatmap_DE_genes = make_dendrogram(tbl_centered_sdScaled_justDEgenes)
hclust_obj_DE_genes = hclust_pheatmap_DE_genes[[1]]
plot(hclust_obj_DE_genes)

medians_tbl_DE_genes = hclust_pheatmap_DE_genes[[2]]
pheatmap_DE_genes = pheatmap(t(medians_tbl_DE_genes), fontsize_row = .5, fontsize_col = 8, angle_col = 90, cutree_cols = 4)

pdf(paste0(Output_dir,paste0("JB_",sample,"pheatmap_DE_gene_test.pdf"), sep = ""), width = 4, height = 8)
pheatmap_DE_genes
dev.off()


