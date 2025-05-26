wd <- "~/R/"

setwd(wd)

wdcode <- paste0("code/") 
ver <- 'Figure6/'
wdresult <- paste0("result/",ver)
wddata <- paste0("data/",ver)
wdplot <- paste0("plot/",ver)
source(paste0(wdcode,"function.R"))
load_packages()
load_plot_packages() 


log_data <- readRDS(paste0(wddata,"log_data.rds"))
metadata <- readRDS(file=paste0(wddata,'metadata.RDS'))


#a  -------------------------------
log_data[is.na(log_data)] <- 0

impute_data1 <- log_data
cluster <- 17

row_cluster <- hclust(dist(impute_data1),method="ward.D2")
cutree_result <- as.data.frame(cutree(row_cluster, k = cluster))
colnames(cutree_result) <- c('cluster')
cutree_result$gene <- rownames(cutree_result)
df <- dft(cutree_result$cluster)

clusterlist <- list()
for (i in unique(cutree_result$cluster)) {
  a <- dplyr::filter(cutree_result,cluster==i)
  i <- paste0('MS',as.character(i))
  clusterlist[[i]] <- a$gene
}

annotation_col <- as.data.frame(rownames(impute_data1))
colnames(annotation_col) <- 'gene'
annotation_col <- left_join(annotation_col,cutree_result,by='gene')
annotation_col <- data.frame(annotation_col[,2])
rownames(annotation_col) <- rownames(impute_data1)
colnames(annotation_col) <- 'Cluster'
annotation_col$Cluster <- as.character(annotation_col$Cluster)


ann_colors = list(
  Cluster=c('1'='#5E2B69','2'='#416997','3'='#1AA799','4'='#AECD5A','5'='#F0D64D',
            '6'='#B5252E','7'='#F7D499','8'='#EEC6D4','9'='#0093D0','10'='#F0AD82',
            '11' = '#85A642', '12'='#264653','13'='#BC5D50','14'='#706d94','15'='#C97231'  ,
            '16' = '#E1C959', '17'='#D79E62','18'='#4D83AA','19'='#DAA83D','20'='#98B9DA' ,
            '21' = '#85A642', '22'='#264653','23'='#BC5D50','24'='#706d94','25'='#C97231' )
  # Group=c("N"="blue","T"="red")
)

library(pheatmap)
pdf(paste0(wdplot,"a.pdf"), width = 12.5, height = 4.5)
bk <- c(seq(15,20,by=1),seq(20.1,33,by=1))
pheatmap::pheatmap(
  t(impute_data1),
  #scale = "column",#
  cluster_rows =T, 
  cluster_cols = T, 
  # #color = colorRampPalette(c('#FA7F6F','#E9C46A','#8983bf')),
  color = c(colorRampPalette(colors = c("#F5F3B1","#C55548"))(length(bk)/2),
            colorRampPalette(colors = c("#C55548","#1F1243"))(length(bk)/2)),#EBEBEB,BD6379 F9E855 48938B 3B0F4F
  clustering_method = "ward.D2",
  cutree_col = cluster,
  show_colnames = F,
  show_rownames = T,
  # annotation_row = annotation_row,
  #scale = 'column',
  #annotation_row = annotation_col,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_col = annotation_col$Cluster,
  breaks= bk
)
dev.off() 

saveRDS(clusterlist,file=paste0(wdresult,'MScluster.RDS'))


#b -----------
clusterlist <- readRDS(paste0(wddata,'MScluster.RDS'))
deplist <- readRDS(file=paste0(wddata,'MS_dep.RDS'))

names(clusterlist) <- c("MS_c1","MS_c2" , "MS_c3",  "MS_c4" , "MS_c5"  ,"MS_c6" , "MS_c7" , "MS_c8" , "MS_c9" ,
                        "MS_c10", "MS_c11" ,"MS_c12", "MS_c13" ,"MS_c14", "MS_c15" ,"MS_c16" ,"MS_c17")

datacbind <- list_list_overlap_jaccard_p(clusterlist,deplist,bg_num=1588,FDR=T)

inpudata <- as.matrix(datacbind[['jaccardindex']])
inpudata <- melt(inpudata)
result_p <- melt(as.matrix(datacbind[['p']]))
colnames(inpudata) <- c('Var1','Var2','value')
colnames(result_p) <- c('Var1','Var2','value')

fig <- list()
heat <- ggplot(inpudata, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = value), color = "white", size = 0.2) +
  # scale_fill_manual(values =c(  "white", "#F08080", "#CD0000"))+
  scale_fill_gradient2(name = 'Jaccard index',
                       low = "#F5F3B1",
                       mid = "#C55548",
                       high = "#1F1243",
                       midpoint = 0.5,
                       # space = 'Lab',
                       na.value = "#E6E6E6",
                       guide = "colorbar"
                       # aesthetics = "colour"
  ) + #fill colour
  geom_text(data = result_p, aes(x = Var2, y = Var1, label = value), size = 5,vjust = 0.8)  +   
  theme(
    panel.grid = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x = element_text(size = 6, angle = 30, hjust = 1),
    axis.ticks.x=element_blank(), 
    axis.title.y=element_blank(),
    legend.position = "none",
    legend.key.size = unit(0.3, "cm"), 
    legend.text = element_text(size = 6.5),
    axis.text.y = element_text(size = 6.5),
    # axis.text.y = element_blank(),
    axis.ticks.y=element_blank(),
    axis.line = element_blank()
  )+themelegendright
heat
ggsave(paste0(wdplot,"b.pdf"), width = 7, height =6, units = "cm")

#c ------------------------------


load(paste0(wddata,'F6data.R'))


pba_af <- PBA[,1:6]
pba_mf <- PBA[,7:12]
pba_sec <- PBA[,13:18]

log_data[is.na(log_data)] <- 0
data_af <- log_data[,1:6]
data_mf <- log_data[,7:12]
data_sec <- log_data[,13:18]
data_af <- data_af[rowSums(data_af > 0) >= 4, ]
data_mf <- data_mf[rowSums(data_mf > 0) >= 4, ]
data_sec <- data_sec[rowSums(data_sec > 0) >= 4, ]
data_af[data_af==0] <- NA
data_mf[data_mf==0] <- NA
data_sec[data_sec==0] <- NA

marker <- c('CD81','APOE') 

p1 <- list()
for (i in marker)  {
  
  ms1 <- melt(cbind(data_af[i,],data_mf[i,],data_sec[i,])  ) %>% 
    mutate(catems='MS',tech=c(rep("AF4", 6), rep("MF", 6), rep("SEC", 6)))
  pba1 <- melt(PBA[i,]) %>% mutate(catepbs='PBA')
  inputdata <- left_join(ms1,pba1,by='variable') 
  
  p1[[i]] <- ggplot(inputdata, aes(x = value.x, y = value.y) )+
    geom_point(aes(fill= tech),size = 1.5, alpha = 1,shape=21) +
    theme2_noback +theme(legend.position = 'right') +
    labs(title = i, y = "PBA - log2(NPX)", x = "MS - log2(Intensity)", legend = 'spearman corr') +
    scale_fill_manual(values=c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','PBA-control'))+
    theme2_noback+ themelegendnone+
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.5) + 
    stat_cor(method = "pearson", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
             # label.x = 20, label.y = 15, 
             size = 2.7) 
  
}

library(gridExtra)
combined_plot <- do.call(grid.arrange, c(p1, ncol = 2) )
ggsave(paste0(wdplot,'d',".pdf"),combined_plot,width =14, height =6, units = "cm")


# d ---------------------------

cluster_df <- bind_rows(lapply(names(clusterlist), function(cluster) {
  data.frame(Protein = clusterlist[[cluster]], Cluster =  sub("MS", "MS_c", cluster))
}))
group_df <- bind_rows(lapply(names(deplist), function(group) {
  data.frame(Protein = deplist[[group]], Group = group)
}))
plotdata <- inner_join(cluster_df, group_df, by = "Protein")

library(webshot)
library(networkD3)
sankey <- plotdata %>%
  dplyr::rename(source = Cluster, target = Group) %>%  
  dplyr::mutate(value = 1) %>%  
  dplyr::group_by(source, target) %>%
  dplyr::summarise(value = sum(value), .groups = "drop")  

nodes <- data.frame(name = unique(c(sankey$source, sankey$target)))

sankey$IDsource <- match(sankey$source, nodes$name) - 1
sankey$IDtarget <- match(sankey$target, nodes$name) - 1


library(networkD3)

custom_colors <- c(
  "AF4_detected only" = "#D1C1E1", 
  "AF4_elevated" = "#A778B4",
  "CMDS_detected only" = "#5568B8",
  "CMDS_elevated" = "#69B190",
  "SEC_detected only" = "#BEBC48",
  "SEC_elevated" = "#DF4828",
  "MS_c1" = '#5E2B69',       
  "MS_c2" = '#416997',
  "MS_c3" = '#1AA799', 
  "MS_c4" = '#AECD5A',
  "MS_c5" = '#F0D64D',
  "MS_c6" = '#B5252E',
  "MS_c7" = '#F7D499',
  "MS_c8" = '#EEC6D4',               
  "MS_c9" = '#0093D0',
  "MS_c10" = '#F0AD82',
  "MS_c11" = '#85A642',
  "MS_c12" = '#264653',
  "MS_c13" = '#BC5D50',
  "MS_c14" = '#706d94',
  "MS_c15" = '#C97231',              
  "MS_c16" = '#E1C959',
  "MS_c17" = '#D79E62'
)

color_scale <- paste0(
  'd3.scaleOrdinal()',
  '.domain(["', paste(names(custom_colors), collapse = '", "'), '"])',
  '.range(["', paste(custom_colors, collapse = '", "'), '"])'
)

nodes$group <- as.character(nodes$name) 


sankeyplot <- sankeyNetwork(
  Links = sankey, Nodes = nodes,
  Source = "IDtarget", Target = "IDsource",
  Value = "value", NodeID = "name",
  fontSize = 12, nodeWidth = 30,
  colourScale = color_scale, 
  NodeGroup = "group"       
)
webshot::install_phantomjs()
saveNetwork(sankeyplot,file=paste0(wdplot,"sankeyplot.html"))
webshot::webshot(paste0(wdplot,"sankeyplot.html"),file=paste0(wdplot,"/csankeyplot3.pdf"),vwidth=300,vheight=500)

library(dplyr) 
library(org.Hs.eg.db) 
library(clusterProfiler) 
library(DOSE)
library(ggplot2) 
library(RColorBrewer)
bindlist <- c(deplist,clusterlist)
#GO
go_enrich_function <- function( gene , backgroundgene){  
  OrgDb <- 'org.Hs.eg.db'
  go_up <- enrichGO(gene=gene,
                    OrgDb='org.Hs.eg.db',
                    universe = backgroundgene,
                    keyType = "SYMBOL", 
                    ont="CC", ## CC MF
                    pAdjustMethod = "BH",
                    qvalueCutoff=0.05,
                    pvalueCutoff=0.05,
  )
  r_go_up <- go_up@result
  r_go_up <- filter(r_go_up, pvalue < 0.05)
  return(r_go_up)
}

kegg_enrich_func <- function(gene){
  diff_entrez <- bitr(gene,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")
  KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                          organism = "hsa", 
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          minGSSize = 10,
                          maxGSSize = 500)
  
  KEGG_diff <- setReadable(KEGG_diff,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID")
  KEGG_diff <- filter(KEGG_diff, pvalue < 0.05)
  
}

p <- list()
p0 <- list()
resultbindgo <- NULL
for (i in names(bindlist)) {
  re1 <- go_enrich_function(bindlist[[i]],rownames(log_data))
  re1$cate <- i
  re1 <- re1 %>% 
    dplyr::filter(Count > 5 ) %>% 
    head(n = 3)
  resultbindgo <- rbind(resultbindgo,re1)
  
  
  custom_colors <- c(
    "AF4_detected only" =  "#D1C1E1FF",  "AF4_elevated" = "#A778B4FF","CMDS_detected only" = "#5568B8FF","CMDS_elevated"="#69B190FF",
    "SEC_detected only"="#BEBC48FF","SEC_elevated"= "#DF4828FF","MS_c1" ='#5E2B69',       
    "MS_c2"= '#416997',"MS_c3"='#1AA799' ,"MS_c4"='#AECD5A' ,"MS_c5"='#F0D64D' ,"MS_c6"='#B5252E' ,"MS_c7"= '#F7D499',"MS_c8"='#EEC6D4' ,              
    "MS_c9"='#0093D0' ,"MS_c10"= '#F0AD82',"MS_c11"= '#85A642',"MS_c12"='#264653' ,"MS_c13"='#BC5D50' ,"MS_c14"= '#706d94',"MS_c15" ='#C97231' ,             
    "MS_c16"= '#E1C959',"MS_c17"=  '#D79E62'  
  )
  
  p[[i]] <- ggplot(data = re1, 
                   aes(x = -log10(pvalue), y = Description, fill = cate)) +
    geom_bar(stat = "identity", width = 0.7, position = position_dodge(preserve = "single")) + 
    labs(x = "-log10(p-value)") +
    scale_y_discrete(position = "right") + 
    theme2_nobackplus + 
    scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
    scale_fill_manual(values = custom_colors,name='')+
    theme(
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(),   
      axis.title.y = element_blank(),  
      axis.ticks.y.right = element_blank(),           
      axis.ticks.y.left = element_line(),            
      axis.line.y.right = element_blank(),      
      axis.line.y.left = element_line() ,
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )+themelegendleft 
  
  
  p0[[i]] <- ggplot(data = re1, 
                    aes(x = -log10(pvalue), y = Description, fill = cate)) +
    geom_bar(stat = "identity", width = 0.7, position = position_dodge(preserve = "single")) + 
    labs(x = "-log10(p-value)") +
    scale_y_discrete(position = "right") + 
    theme2_nobackplus + 
    scale_x_continuous(limits = c(0, 0), expand = c(0, 0)) +
    scale_fill_manual(values = custom_colors,name='')+
    theme(
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(), 
      axis.text.y = element_text(size = 7, hjust = 0),  
      axis.title.y = element_blank(),  
      axis.ticks.y.right = element_blank(),                
      axis.ticks.y.left = element_line(),                  
      axis.line.y.right = element_blank(), 
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.line.y.left = element_line()                    
    )+themelegendleft 
  
  
}

for (i in names(p)) {
  p[[i]] <- p[[i]]+theme(legend.position = 'left', 
                         legend.title = element_text(size = 6),
                         legend.key.size = unit(legendsize, "cm"),  
                         legend.text = element_text(size = 6,angle=70),
                         panel.background = element_rect(fill=NA, colour=NA))
  p0[[i]] <- p0[[i]]+theme(legend.position = 'left', 
                           legend.title = element_text(size = 6),
                           legend.key.size = unit(legendsize, "cm"),  
                           legend.text = element_text(size = 6,angle=70),
                           panel.background = element_rect(fill=NA, colour=NA))
}
library(cowplot)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[6]],p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], 
          p[[12]], p[[13]], p[[14]], p[[16]],p[[17]], p[[18]], p[[19]], p[[21]],p[[22]],p[[23]],
          ncol = 1,  align = "v", axis = "tb")
ggsave(paste0(wdplot,"c1.pdf"), width = 7, height =20, units = "cm")

plot_grid(p0[[1]], p0[[2]], p0[[3]], p0[[4]], p0[[6]],p0[[7]], p0[[8]], p0[[9]], p0[[10]], p0[[11]], 
          p0[[12]], p0[[13]], p0[[14]], p0[[16]],p0[[17]], p0[[18]], p0[[19]], p0[[21]],p0[[22]],p0[[23]],
          ncol = 1,  align = "v",
          axis = "tb")
ggsave(paste0(wdplot,"c2.pdf"), width = 10, height =20, units = "cm")


#end -----------


# e f   -------------------------------------

clusterlist <- readRDS(paste0(wddata,'MScluster.RDS'))
deplist <- readRDS(file=paste0(wddata,'MS_dep.RDS'))

library(circlize)
library(viridis)
library(reshape2)

ac <- NULL
bc <- NULL 
source <- NULL
target <- NULL

#bindlist <- c(deplist,clusterlist)
names(clusterlist) <- c("MS_c1","MS_c2" , "MS_c3",  "MS_c4" , "MS_c5"  ,"MS_c6" , "MS_c7" , "MS_c8" , "MS_c9" ,
                        "MS_c10", "MS_c11" ,"MS_c12", "MS_c13" ,"MS_c14", "MS_c15" ,"MS_c16" ,"MS_c17")
tissuegene <- readRDS(paste0(wddata,'tissuegene.RDS'))
singlecellgene <- readRDS(paste0(wddata,'singelcellgene.RDS'))

## tissue spec -----------------------------
for (y in names(clusterlist)) {
  for (j in names(tissuegene)) {
    
    a <- clusterlist[[y]]
    b <- as.character(tissuegene[[j]])
    
    for (i in a) {
      source <- c(source,rep(y,length(b)))
      target <- c(target,rep(j,length(b)))
      ac <- c(ac,rep(i,length(b)))
      bc <- c(bc,b)
      
    }
  }}

source <- as.data.frame(source)
source <- cbind(source,as.data.frame(target))
source <- cbind(source,as.data.frame(ac))
source <- cbind(source,as.data.frame(bc))
source$score <- ifelse(source$ac == source$bc, 1, NA)
colnames(source) <- c('source', 'target', 'ligand', 'receptor' , 'value')
source <- as.data.frame(source)

source5 <- na.omit(source[,c(1,2,5)])
catenumb <- dft(source5$target) %>% 
  dplyr::filter(Freq > 1) 
source5 <- dplyr::filter(source5,target %in% catenumb$file) 
source5$target[source5$target == 'skin 1'] <- 'Skin'
source5$target <- sub("^(.)", "\\U\\1", source5$target, perl = TRUE)
length(unique(source5$source))+length(unique(source5$target)) 

source5 <- readRDS(paste0('result/F6/tissuecircle.RDS'))
gaplength <- c(rep(4,length(unique(source5$source))-1),14,rep(4,length(unique(source5$target))-1),14)


mycolor <- c('#BC5D50FF','#D4978FFF','#E2BAB4FF','#C97231FF','#D79E62FF','#E1C367FF','#E9E475FF','#85A642FF','#A4B841FF',#9
             '#ADC482FF','#75B2C5FF','#5B9CB8FF','#98B9DAFF','#3E6B9BFF',#14
             '#C27877','#e2a9a3','#EDC4BB',"#DDD8EFFF",'#D4C0AC','#B89E8A','#9E7A64','#7D576F','#957e95',"#C3A8D1FF",
             '#1E446C')
{
pdf(paste0(wdplot,"e.pdf"), width = 12, height = 12)

circos.clear()
circos.par(start.degree = 175, 
           gap.after = gaplength,
           #gap.degree = 4, 
           track.margin = c(-0.05, 0.1), 
           points.overflow.warning = FALSE,
           canvas.xlim = c(-2, 2),  
           canvas.ylim = c(-2, 2),
           points.overflow.warning = FALSE)


#par(cex = 2, mar = c(0, 0, 0, 0))
par(cex = 1, ps = 21, mar = c(0, 0, 0, 0))
chordDiagram(
  x = source5, 
  grid.col = mycolor, 
  transparency = 0.6,
  annotationTrack = "grid",  
  annotationTrackHeight = c(0.07, 0.01), 
  diffHeight = -0.05, 
  link.arr.type = "big.arrow",
  preAllocateTracks = list(track.height = 0.02)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 0.1,
      labels = sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 1 
    )
  },
  bg.border = NA
)
name <- c(namet(source5$source),namet(source5$target))
dev.off() 
}

## singlecell spec -------------------

ac <- NULL
bc <- NULL 
source <- NULL
target <- NULL

for (y in names(clusterlist)) {
  for (j in names(singlecellgene)) {
    
    a <- clusterlist[[y]]
    b <- as.character(singlecellgene[[j]])
    
    for (i in a) {
      source <- c(source,rep(y,length(b)))
      target <- c(target,rep(j,length(b)))
      ac <- c(ac,rep(i,length(b)))
      bc <- c(bc,b)
      
    }
  }}

source <- as.data.frame(source)
source <- cbind(source,as.data.frame(target))
source <- cbind(source,as.data.frame(ac))
source <- cbind(source,as.data.frame(bc))
source$score <- ifelse(source$ac == source$bc, 1, NA)
colnames(source) <- c('source', 'target', 'ligand', 'receptor' , 'value')
source <- as.data.frame(source)

source5 <- na.omit(source[,c(1,2,5)])
catenumb <- dft(source5$target) %>% 
  dplyr::filter(Freq > 1) 
source5 <- dplyr::filter(source5,target %in% catenumb$file) 
source5$target <- sub("^(.)", "\\U\\1", source5$target, perl = TRUE)
length(unique(source5$source))+length(unique(source5$target)) 
#saveRDS(source5,file=paste0(wdresult,'singlecellcircle.RDS'))

gaplength <- c(rep(4,length(unique(source5$source))-1),14,rep(4,length(unique(source5$target))-1),14)

mycolor <- c('#BC5D50FF','#D4978FFF','#E2BAB4FF','#C97231FF','#D79E62FF','#E1C367FF','#E9E475FF','#85A642FF','#A4B841FF',#9
             '#ADC482FF','#75B2C5FF','#5B9CB8FF','#98B9DAFF','#3E6B9BFF',#14
             '#C27877','#e2a9a3','#EDC4BB',"#DDD8EFFF",'#D4C0AC','#B89E8A','#9E7A64','#7D576F','#957e95',"#C3A8D1FF",
             '#1E446C','#2a577e',"#A778B4FF", #24
             '#BC5D50FF','#D4978FFF'#,'#E2BAB4FF','#C97231FF','#D79E62FF'#,'#E1C367FF','#E9E475FF'
)

{
  pdf(paste0(wdplot,"f.pdf"), width = 12, height = 12)
  
  circos.clear()
  circos.par(start.degree = 170, 
             gap.after = gaplength,
             #gap.degree = 4, 
             track.margin = c(-0.05, 0.1), 
             points.overflow.warning = FALSE,
             canvas.xlim = c(-2, 2),  
             canvas.ylim = c(-2, 2),
             points.overflow.warning = FALSE)
  
  
  #par(cex = 2, mar = c(0, 0, 0, 0))
  par(cex = 1, ps = 21, mar = c(0, 0, 0, 0))
  chordDiagram(
    x = source5, 
    grid.col = mycolor, 
    transparency = 0.6,
    annotationTrack = "grid",  
    annotationTrackHeight = c(0.07, 0.01), 
    diffHeight = -0.05, 
    link.arr.type = "big.arrow",
    preAllocateTracks = list(track.height = 0.02)
  )
  
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.index = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      circos.text(
        x = mean(xlim),
        y = ylim[1] + 0.1,
        labels = sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1 
      )
    },
    bg.border = NA
  )
  name <- c(namet(source5$source),namet(source5$target))
  dev.off() 
}
#end-------------
#S5 other protines ---------------------

load(paste0(wddata,'F6data.R'))

pro <- intersect(rownames(log_data),rownames(PBA))
p1 <- list()
for (i in pro)  {
  
  ms1 <- melt(cbind(data_af[i,],data_mf[i,],data_sec[i,])  ) %>% 
    mutate(catems='MS',tech=c(rep("AF4", 6), rep("MF", 6), rep("SEC", 6)))
  pba1 <- melt(PBA[i,]) %>% mutate(catepbs='PBA')
  inputdata <- left_join(ms1,pba1,by='variable') 
  
  p1[[i]] <- ggplot(inputdata, aes(x = value.x, y = value.y) )+
    geom_point(aes(fill= tech),size = 1.5, alpha = 1,shape=21) +
    theme2_noback +theme(legend.position = 'right') +
    labs(title = i, y = "PBA - log2(NPX)", x = "MS - log2(Intensity)", legend = 'spearman corr') +
    scale_fill_manual(values=c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','PBA-control'))+
    theme2_noback+ themelegendnone+
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.5) + 
    stat_cor(method = "pearson", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
             # label.x = 20, label.y = 15, 
             size = 2.7) 

}

library(gridExtra)
combined_plot <- do.call(grid.arrange, c(p1, ncol = 6) )

ggsave(paste0(wdplot,'S5',".pdf"),combined_plot,width =30, height =42, units = "cm")

