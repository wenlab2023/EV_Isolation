wd <- "~/R/"
setwd(wd)

wdcode <- paste0("code/") 
ver <- 'Figure4/'
wdresult <- paste0("result/",ver)
wddata <- paste0("data/",ver)
wdplot <- paste0("plot/",ver)
source(paste0(wdcode,"function.R"))
load_packages()
load_plot_packages() 
library(Seurat)
library(BPCells)
library(ggplot2)


metadata <- readRDS(paste0(wddata,'metadata.RDS'))
#a single type of protein   -------------------------------------------


combined_data1.1 <- read.csv(paste0(wddata,"F4DEPs.csv"),check.names=F)
combined_data1.1 <- combined_data1.1[,-1]
unique(combined_data1.1$cate)
combined_data1.1 %>% dplyr::filter(cate=='AF4') %>% 
  ggplot(aes(x=log2FC ,y=avg_percent,color=spec1))+
  geom_point(size=1.5,alpha=0.7)+
  scale_color_manual(values=c('#B5252E','#416997','#D4978F','#DAA83D','#ADC482'))+
  labs(y='log2(FC)',x = "Average percent",title = 'AF4 samples',colour = "")+
  theme2_noback+themelegendnone+
  geom_text_repel(aes(x=log2FC, y=avg_percent, label=new_column),
                  segment.size = 0.1,nudge_y = 0,
                  vjust=1, color="black", size=2, max.overlaps=20) +    
  scale_y_continuous(limits = c(0, 0.4), expand = c(0.05, 0),breaks = seq(0, 0.4, by = 0.05))+
  scale_x_continuous(limits = c(0, 2), expand = c(0.05, 0),breaks = seq(0, 2, by = 0.5))

# scale_x_continuous(limits = c(0, 0.4), expand = c(0.1, 0),breaks = seq(0, 0.4, by = 0.1), minor_breaks = seq(0, 0.4, by = 0.025)) +  
# scale_y_continuous(limits = c(0, 2), expand = c(0.1, 0),breaks = seq(0, 2, by = 0.5), minor_breaks = seq(0, 2, by = 0.25)) 

ggsave(paste0(wdplot,"a.AF4.pdf"), width = 5, height =12, units = "cm")


combined_data1.1 %>%dplyr::filter(cate=='Exo-CMDS') %>% 
  ggplot(aes(x=log2FC ,y=avg_percent,color=spec1))+
  geom_point(size=1.5,alpha=0.7)+
  scale_color_manual(values=c('#B5252E','#416997','#D4978F','#DAA83D','#ADC482'))+
  labs(y='log2(FC)',x = "average percent",title = 'EXO-CMDS',colour = "")+
  theme2_noback+themelegendnone+
  geom_text_repel(aes(x=log2FC, y=avg_percent, label=new_column),
                  segment.size = 0.1,nudge_y = 0,
                  vjust=1, color="black", size=2, max.overlaps=40) +    
  scale_y_continuous(limits = c(0, 0.02), expand = c(0.05, 0),breaks = seq(0, 0.02, by = 0.005))+
  scale_x_continuous(limits = c(0, 2), expand = c(0.05, 0),breaks = seq(0, 2, by = 0.5))

ggsave(paste0(wdplot, "a.CMDS.pdf"), width = 5, height =5, units = "cm")


combined_data1.1 %>%dplyr::filter(cate=='SEC') %>% 
  ggplot(aes(x=log2FC ,y=avg_percent,color=spec1))+
  geom_point(size=1.5,alpha=0.7)+
  scale_color_manual(values=c('#416997','#D4978F','#DAA83D','#ADC482'))+
  labs(y='log2(FC)',x = "average percent",title = 'SEC',colour = "")+
  theme2_noback+themelegendnone+
  geom_text_repel(aes(x=log2FC, y=avg_percent, label=new_column),
                  segment.size = 0.1,nudge_y = 0,
                  vjust=1, color="black", size=2, max.overlaps=20) +    
  scale_y_continuous(limits = c(0, 0.02), expand = c(0.05, 0),breaks = seq(0, 0.02, by = 0.005))+
  scale_x_continuous(limits = c(0, 1.5), expand = c(0.05, 0),breaks = seq(0, 1.5, by = 0.25))
ggsave(paste0(wdplot, "a.SEC.pdf"), width = 5, height =5, units = "cm")


combined_data1.1 %>%dplyr::filter(cate=='Plasma') %>% 
  ggplot(aes(x=log2FC,y=avg_percent,color=spec1))+
  geom_point(size=1.5,alpha=0.7)+
  scale_color_manual(values=c('#B5252E','#416997','#D4978F','#DAA83D','#ADC482'))+
  labs(y='log2(FC)',x = "average percent",title = 'Plasma-control',colour = "")+
  theme2_noback+themelegendnone+
  geom_text_repel(aes(x=log2FC, y=avg_percent, label=new_column),
                  segment.size = 0.1,nudge_y = 0,
                  vjust=1, color="black", size=2, max.overlaps=20) +    
  scale_x_continuous(limits = c(0, 1.8), expand = c(0.05, 0),breaks = seq(0, 1.8, by = 0.5))+
  scale_y_continuous(limits = c(0, 0.09), expand = c(0.05, 0),breaks = seq(0, 0.09, by = 0.015))  
ggsave(paste0(wdplot, "a.control.pdf"), width = 5, height =7, units = "cm")




# b single EV ------------------------


options(future.globals.maxSize = 1e9)
options(stringsAsFactors = F) 

load(paste0(wddata,'singleEV_data.Rdata'))
rm(inputdata,sparse_bind_matrix,sparse_inputdada)

obj <- CreateSeuratObject(counts = sparse_bind,project = "seurat", min.cells=0, min.features=0, names.delim = "_")
obj <- NormalizeData(obj,normalization.method = "LogNormalize",scale.factor = 1e+07)

obj <- FindVariableFeatures(obj)
obj <- SketchData( 
  object = obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
obj

DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:15) 
obj <- FindClusters(obj, resolution = 0.16, algorithm = 1)
obj2 <- RunUMAP(obj, dims = 1:15, return.model = T) 

color <- c('#98B9DA','#C97231','#A4B841','#264653','#D79E62','#9EC8D7','#E1C367','#475194','#DAA83D', #8
           '#BC5D50','#D4978F','#E2BAB4','#85A642',#12
           '#ADC482','#BCD3EA','#75B2C5','#DCB655','#E1C959','#B4D3E0',#18
           '#E9E475','#5B9CB8','#3E6B9B',#21
           '#4D83AA','#DADFE8','#D79E62','#9EC8D7'
)

DimPlot(obj2,label=T,label.size = 5,
        cols = color,
        repel=F,
        reduction = "umap" #"umap"'tsne'
) 

obj3 <- ProjectData(
  object = obj2,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:15,
  refdata = list(cluster_full = "seurat_clusters")
)


obj3$cluster_full <- factor(obj3$cluster_full, levels = c(0:15))
meta <- obj3@meta.data
meta <- left_join(meta,metadata, by=c('orig.ident'='seurat'))
rownames(meta) <- rownames(obj3@meta.data)
obj3@meta.data <- meta


#2 bind cluster 

DefaultAssay(obj3) <- "RNA"

metadata_obj <- obj3@meta.data

metadata_obj$cluster_full[metadata_obj$cluster_full == 13] <- 10
metadata_obj$cluster_full[metadata_obj$cluster_full == 14] <- 13
metadata_obj$cluster_full[metadata_obj$cluster_full == 15] <- 14
metadata_obj$cluster_full <- factor(metadata_obj$cluster_full, levels = c(0:14))

table(metadata_obj$cluster_full)

obj3@meta.data <- metadata_obj

saveRDS(obj3,file = paste0(wdresult,'obj3.RDS'))





##bumap all-------------------------


color <- c('#98B9DA','#C97231','#899E2F','#DAA83D',#3 DAA83D
           '#4D83AA','#E1C959','#475194','#9EC8D7', #7 #264653
           '#BC5D50','#264653','#D79E62','#E1C959',#11'#E2BAB4'
           '#5B9CB8','#cca7cb','#E2E9B8'
           
           ,'#DCB655','#E1C959','#B4D3E0',#18
           '#E9E475','#DADFE8','#D79E62',#21 DADFE8  D79E62
           '#4D83AA','#5B9CB8','#5B9CB8','#9EC8D7'
)

p1 <- DimPlot(obj3, label = T, label.size = 2, reduction = "full.umap", 
              group.by = "cluster_full", alpha = 0.4,pt.size = 1,
              cols = color) + 
  labs(x='',y='',title='')
#NoLegend()+
# themedim
p1
ggsave(paste0(wdplot,"b.pdf"),plot=p1 ,width = 18, height =13, units = "cm")


#c  umap methods-----------------

DimPlot(obj3, label = T, label.size = 2,cols = color, reduction = "full.umap",
        split.by = "method1",pt.size = 1,ncol=2,
        group.by = "cluster_full", alpha = 1.2,raster=T) + 
  labs(x='',y='',title='')+
  #NoLegend()+
  theme2_noback
ggsave(paste0(wdplot,"c.pdf") ,width = 18, height =14, units = "cm")



#d cell fraction -----------------
obj3 <- readRDS(wdresult,'obj3.RDS')

bar_input <- obj3@meta.data
sumall <- length((bar_input$cluster_full))

result <- bar_input %>%
  dplyr::group_by(method1) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::count(cluster_full) %>%
  dplyr::mutate(percentage = n / sum(n) ) %>%
  dplyr::ungroup()

head(result)
sum(result$percentage)

result <- dplyr::filter(result,cluster_full %in% c(5,7,6,1,9,11,2,8,12,13,14)) 
result <- left_join(result,metadata,by=c('orig.ident'='seurat'))
result$cluster_full <- factor(result$cluster_full, levels =  c(5,7,6,1,9,11,2,8,12,13,14) )

bar_input_af4 <- dplyr::filter(result,method1 == 'AF4')
bar_input_mf <- dplyr::filter(result,method1 == 'MF')
bar_input_sec <- dplyr::filter(result,method1 == 'SEC')
bar_input_control <- dplyr::filter(result,method1 == 'Control')


color <- c( '#E1C959','#9EC8D7',  '#475194','#C97231',#3 DAA83D
            '#264653','#E1C959', '#899E2F', '#BC5D50', '#5B9CB8', '#cca7cb', '#E2E9B8')


p1 <- ggplot(bar_input_af4, 
             aes(x = name_new1, y = percentage, fill = factor(cluster_full))) + 
  geom_bar(stat = "identity", position = "stack", width = 0.7)+
  RotatedAxis()+
  scale_fill_manual(values=color)+ 
  labs(x='',y='',fill='')+
  theme1_noback +   scale_y_continuous(limits = c(0, 0.9), expand = c(0,0), breaks = seq(0, 0.9, 0.1))+
  theme(legend.position = 'none') + labs(title = 'AF4')
#library(ggalt)
p1
ggsave(paste0(wdplot,"d.AF4.pdf"),plot=p1,width = 5.5, height = 4.5, units = "cm")


p2 <- ggplot(bar_input_mf, 
             aes(x = name_new1, y = percentage, fill = factor(cluster_full))) + 
  geom_bar(stat = "identity", position = "stack", width = 0.7)+
  RotatedAxis()+
  scale_fill_manual(values=color)+ 
  labs(x='',y='',fill='')+
  theme1_noback + scale_y_continuous(limits = c(0, 0.9), expand = c(0,0), breaks = seq(0, 0.9, 0.1))+
  theme(legend.position = 'none')+ labs(title = 'EXO-CMDS')
p2
ggsave(paste0(wdplot,"d.CMDS.pdf"),plot=p2,width = 5.5, height = 4.5, units = "cm")


p3 <- ggplot(bar_input_sec, 
             aes(x = name_new1, y = percentage, fill = factor(cluster_full))) + 
  geom_bar(stat = "identity", position = "stack", width = 0.7)+
  RotatedAxis()+
  scale_fill_manual(values=color)+ 
  labs(x='',y='',fill='')+
  theme1_noback + scale_y_continuous(limits = c(0, 0.9), expand = c(0,0), breaks = seq(0, 0.9, 0.1))+
  theme(legend.position = 'none') + labs(title = 'SEC')
p3
ggsave(paste0(wdplot,"d.sec.pdf"),plot=p3,width = 5.5, height = 4.5, units = "cm")


p4 <- ggplot(bar_input_control, 
             aes(x = name_new1, y = percentage, fill = factor(cluster_full))) + 
  geom_bar(stat = "identity", position = "stack", width = 0.7)+
  RotatedAxis()+
  scale_fill_manual(values=color)+ 
  labs(x='',y='',fill='')+
  theme1_noback + 
  scale_y_continuous(limits = c(0, 0.9), expand = c(0,0), breaks = seq(0, 0.9, 0.1))+
  theme(legend.position = 'none')+ labs(title = 'Plasma-control')
p4
ggsave(paste0(wdplot,"d.control.pdf"),plot=p4,width = 5.5, height =4.5, units = "cm")




#e findmarker-----------------

Idents(obj3) <- obj3$cluster_full
DefaultAssay(obj3) <- "RNA" 
obj3 <- ScaleData(obj3)

nk.markers <- FindAllMarkers(obj3,
                             logfc.threshold = 0,
                             test.use = "wilcox",
                             slot = "data",
                             min.pct = 0.2, 
                             only.pos = T)
# write.xlsx(nk.markers,file=paste0(wdresult,"DEPs.xlsx"),
#            sheetName="1", append=T, row.names=F,col.names = T,showNA = T)

gene <- unique(nk.markers$gene)
nk.markers2 <- nk.markers %>%
  group_by(cluster) %>%
  arrange(desc(pct.1)) %>%
slice_head(n = 3)

clusters <- c(0:14)


expression_data <- GetAssayData(obj3, layer = "data")[unique(nk.markers2$gene), ]
metadata_obj <- obj3@meta.data
metadata_obj <- rownames_to_column(metadata_obj)

expression_long <- as.data.frame(as.matrix(expression_data)) %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "cell", value = "expression", -gene)

expression_long <- left_join(expression_long,metadata_obj,by=c('cell'='rowname'))

summary_data <- expression_long %>%
  dplyr::group_by(cluster_full, gene) %>%
  dplyr::summarise(
    avg_exp = mean(expression),
    pct_exp = sum(expression > 0) / n(),
    .groups = "drop"  #
    
  ) 

summary_data$gene <-  factor(summary_data$gene, levels = unique(nk.markers$gene))
summary_data$pct_exp[summary_data$pct_exp==0] <- NA

p <- ggplot(summary_data, aes(x = cluster_full, y = gene)) +
  geom_point(aes(size = pct_exp, color = avg_exp)) +
  scale_size_continuous(
    range = c(0, 5),
    breaks = c(0,  0.25,  0.5,0.75,1),  
    labels = c("", "25%", "50%", "75%", "100%")
  ) +   
  scale_color_gradient2(low = "#F5F3B1", mid = "#C55548", high = "#1F1243", midpoint = 5) +
  theme_minimal() +
  labs(
    #title = "DotPlot of Gene Expression",
    x = "",y = "", color = "Average Expression", size = "Percentage Expression" ) +
  coord_flip() +
  scale_x_discrete(labels = c('PBA_c0','PBA_c1','PBA_c2','PBA_c3','PBA_c4','PBA_c5','PBA_c6','PBA_c7',
                              'PBA_c8','PBA_c9','PBA_c10','PBA_c11','PBA_c12','PBA_c13','PBA_c14'))+
  theme(axis.text.x = element_text( size = 9,angle =90,vjust =0.5,hjust=1),
        axis.text.y = element_text( size = 9),
        legend.position = 'none',
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.3, "cm"), 
        legend.text = element_text(size = 9),
        axis.line.y = element_line(colour="black", size=0.25),
        axis.line.x = element_line(colour="black", size=0.25),
        
        axis.ticks = element_line(colour = "black", size = 0.25),
        axis.ticks.length.x=unit(0.08, "cm"),
        axis.ticks.length.y=unit(0.08, "cm")
  ) 
p

count_evs <- as.data.frame(table(obj3@meta.data[["cluster_full"]]))
colnames(count_evs) <- c('cluster','counts')

p2 <- ggplot(count_evs, aes(x=cluster, y=counts)) +
  geom_bar(data = count_evs, aes(x=cluster, y=counts),
           stat = "identity", position = position_dodge2(width=1), width = 0.6, color='black', fill='#6888a5', linewidth=0.2) +
  theme1_noback +
  theme(axis.text.x = element_text( size = 9,angle =90,vjust = 0.5,hjust=1 ),
        axis.text.y = element_text( size = 9),
        axis.title.x = element_text( size = 9),
        axis.ticks = element_line(colour = "black", size = 0.25),
        axis.ticks.length.x=unit(0.08, "cm"),
        axis.ticks.length.y=unit(0, "cm")
  ) +
  scale_x_discrete(labels = c('','','','','','','','','','','','','','',''))+
  labs(x='', y='EV number')+
  rotate()+
  scale_y_continuous(limits = c(0e+00, 4e+05), expand = c(0,0)) 
p2

library(patchwork)

p %>%
  insert_right(p2,width=0.14)
ggsave(paste0(wdplot,"e.pdf") ,width = 17, height =8, units = "cm")


#s3------------------------------
p <- list()

for (i in unique(nk.markers$cluster)) {
  a <- dplyr::filter(nk.markers,cluster==i)
  gene <- a$gene
  for (j in gene) {
    p[[j]] <- FeaturePlot(obj3,features= j ,raster=T,
                          cols = c("lightgrey", "#d37b6d"), alpha = 3)+
      theme4_backline
  }
}
combined_plot <- do.call(grid.arrange, c(p, ncol = 5))

ggsave(paste0( wdplot,j,".png"),combined_plot,width = 24, height = 5, units = "cm")

#f marker bar cluster----------------------

for (i in c('CD9' ,'CD81')) {
  
  summary_data2 <- dplyr::filter(summary_data, gene == i)# %>%
  # arrange(desc(avg_exp))
  summary_data2$cluster <- (paste0('PBA_','c',summary_data2$cluster_full))
  summary_data2$cluster <- factor(summary_data2$cluster,levels = summary_data2$cluster)
  
  summary_data2$cluster
  custom_colors <- c(
    "PBA_c0" = "#98B9DA", "PBA_c1" = "#C97231", "PBA_c2" = "#899E2F", "PBA_c3" = "#DAA83D","PBA_c4" = "#4D83AA","PBA_c5" = "#E9E475", 
    "PBA_c6" = "#475194",  "PBA_c7" = "#9EC8D7","PBA_c8" = "#BC5D50", "PBA_c9" = "#264653",  "PBA_c10" = "#D79E62",  "PBA_c11" = "#E1C959",
    "PBA_c12" = "#5B9CB8", "PBA_c13" = "#cca7cb",  "PBA_c14" = "#E2E9B8"
  )
  
  summary_data2 %>% 
    ggplot( aes(x=cluster, y=avg_exp,fill=cluster)) + 
    geom_bar(
      stat = "identity", position = position_dodge2(width=1),width = 0.6,color='black',linewidth=0.2)+
    theme_bw()+
    scale_fill_manual(values = custom_colors)+
    theme2_noback+themetext30+themelegendnone+labs(x='')+
    labs(title =i,x='',y='Average Expression' )
  ggsave(paste0(wdplot,'f',i,".pdf") ,width = 7.5, height =4.5, units = "cm")
}



#g cluster bulk corr ---------------------------------

library(pheatmap)
library(ComplexHeatmap)


bind <- NULL
bind <- data.frame(gene=rownames(obj3))
normalized_data <- GetAssayData(obj3, layer = "counts") 
normalized_data <- obj3@assays$RNA$scale.data
for(i in 0:14){
  cellname <- rownames(dplyr::filter(obj3@meta.data,cluster_full == i))
  selected_cells_data <- normalized_data[, colnames(normalized_data) %in% cellname]
  # gene_mean_expression <- apply(selected_cells_data, 1, mean, trim = 0.01)
  gene_mean_expression <- apply(selected_cells_data, 1, mean)
  #gene_expression_df <- data.frame(gene = rownames(selected_cells_data), mean_expression = gene_mean_expression)
  gene_expression_df <- data.frame( mean_expression = gene_mean_expression)
  colnames(gene_expression_df) <- i
  bind <- cbind(bind,gene_expression_df)
  
}
bind <- bind[,-1]
colnames(bind) <- c( "c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10" ,"c11" ,"c12", "c13" ,"c14")
bind[is.na(bind)] <- 0
corr_cluster <- cor(bind,method='spearman')


annotation_col <- as.data.frame(rownames(corr_cluster))
colnames(annotation_col) <- 'sample'
annotation_col$Cluster <- c('Low technical specificity','SEC specific','Low technical specificity', #2
                            'Low technical specificity','Low technical specificity','AF4 specific', #5
                            'CMDS specific','AF4 specific','Low technical specificity', #8
                            'SEC specific','Low technical specificity','SEC specific', #11
                            'Low technical specificity','Low technical specificity','Low technical specificity') #14

annotation_col <- column_to_rownames(annotation_col,var='sample')

ann_colors = list(
  Cluster=c('AF4 specific'='#B5252E','CMDS specific'='#416997','SEC specific'='#DAA83D','Low technical specificity'='grey')
)    


bk <- c(seq(-0.4, -0.01, by = 0.02), seq(0, 0.8, by = 0.05))
p <- pheatmap::pheatmap(
  corr_cluster,
  color = c(colorRampPalette(colors = c("#1F1243", "white"))(length(bk) / 2),
            colorRampPalette(colors = c("white", "#B5252E"))(length(bk) / 2)),
  clustering_method = "ward.D2",
  show_colnames = TRUE,
  show_rownames = TRUE,
  annotation_col = annotation_col,
  #annotation_row = annotation_col,
  annotation_colors = ann_colors,
  breaks = bk,
  cellwidth = 12,
  cellheight = 12,
  fontsize = 13
)
ggsave(paste0(wdplot,"g.pdf"),plot= p,width = 20, height =12, units = "cm")



# S3 --------------------------------

bar_input <- obj3@meta.data

cluster_percent <- bar_input %>%
  dplyr::group_by(orig.ident, cluster_full) %>%
  dplyr::summarise(cell_count = n(), .groups = 'drop') %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(total_cells = sum(cell_count),
                percent = cell_count / total_cells ) %>%
  dplyr::ungroup() %>%
  dplyr::select(orig.ident, cluster_full, percent)

cluster_percent <- left_join(cluster_percent,metadata,by=c('orig.ident'='seurat'))


cluster_percent$method2 <- factor(cluster_percent$method2, levels = c('AF4','MF','SEC','Ori'))
comparisons <- list(c('AF4', 'MF'),c('AF4', 'Ori') ,c('AF4', 'SEC'),c('MF', 'Ori'),c('MF', 'SEC'),c('Ori', 'SEC'))
library(ggsignif)

p <- list()
for(i in c(0:14)){
  
  data <- dplyr::filter(cluster_percent,cluster_full == i)
  
  p_values <- sapply(comparisons, function(comp) {
    wilcox.test(data$percent[data$method2 == comp[1]], 
                data$percent[data$method2 == comp[2]])$p.value
  })
  
  
  fdr_values <- p.adjust(p_values, method = "fdr")
  
  significant_comparisons <- comparisons#[fdr_values < 0.05]
  cluster <- paste0('c',i)
  
  p[[cluster]] <- ggplot(data, aes(x=method2, y=percent, fill=method2)) +
    geom_jitter( width = 0.1, size = 1.7,shape=21,alpha=0.8)+
    geom_boxplot(width = 0.5,colour = "black", outlier.size=0,fill='white', outlier.color = 'white', position = position_dodge(width = 0.9),alpha=0.1)+
    scale_fill_manual(values = c('#B5252E','#416997','#ADC482','#DAA83D')
                      ,name="")+
    labs(y='EV fraction',x = "",color = "",title=paste0('PBA_c',i))+ 
    theme2_noback+ themelegendnone+ themetext30+
    geom_signif(
      comparisons = significant_comparisons,  
      step_increase = 0.07,
      map_signif_level = T,
      test = wilcox.test,
      tip_length = 0,
      vjust = 0.7,
      size=0.2,
      textsize = 5 
    )+  
    scale_x_discrete(labels = c('AF4', 'Exo-CMDS', 'SEC', 'Plasma-control')) 
  p[[cluster]]
}
library(gridExtra)
combined_plot <- do.call(grid.arrange, c(p, ncol = 5))
ggsave(paste0(wdplot,'s2',".pdf"),combined_plot, width = 21, height =24, units = "cm")



