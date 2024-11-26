wd <- "~/R/"
setwd(wd)

wdcode <- paste0("code/")
ver <- 'Figure3/'
wdresult <- paste0("result/",ver)
wddata <- paste0("data/",ver)
wdplot <- paste0("plot/",ver)
source(paste0(wdcode,"function.R"))
load_packages()
load_plot_packages() 

# a PBA EV number ---------------------------

load(paste0(wddata,'F3adata.Rdata'))

data2$method1 <- factor(data2$method1, levels = c('AF4','MF','SEC','Control'))

meancount <- data2 %>%
  dplyr::group_by(method1) %>%
  dplyr::summarize(CV = mean(EV_count), .groups = 'drop')

ggplot(data2, aes(x=method1,y=as.numeric(EV_count),fill=storage))+
  geom_boxplot(width = 0.55,colour = "black", outlier.size=0, position = position_dodge(width = 0.9),size=0.3)+
  geom_jitter(aes(color = storage) ,position = position_dodge(width = 0.9), size = 0.5,color='black')+
  scale_fill_manual(values = c('#5568B8FF','#BEBC48'),   #'#B5252E', '#416997', '#DAA83D', '#ADC482'
                    name = "Storage Type",
                    labels = c("Freeze-thaw",'Non-freeze-thaw')) +
  labs(y='EV number',x = "",color = "")+
  scale_y_continuous(limits = c(1.3e+06, 5.2e+06), expand = c(0.1,0))+ 
  theme2_noback+ 
  scale_x_discrete(labels = c('AF4','Exo-CMDS','SEC','Plasma-control'))+
  themelegendnone

ggsave(paste0(wdplot,"a.pdf"), width = 10, height =6, units = "cm")

# b pca -------------------
#TMM标准化
load(paste0(wddata,'F3bdata+metadata.Rdata'))


library(edgeR)
data <- DGEList(counts = data2)
TMM <- calcNormFactors(data, method="TMM")
normalized_counts <- cpm(TMM, log = F)


wd_method <- 'PCA'

data_pca2 <- prcomp(t(as.matrix(normalized_counts)), scale=T)
percentage <- round(data_pca2$sdev^2 / sum(data_pca2$sdev^2) * 100,2)
percentage <- paste(colnames(as.data.frame(data_pca2[["x"]])),"(", paste(as.character(percentage), "%", ")", sep=""))

data_input <-as.data.frame(data_pca2[["x"]])
data_input$id <- rownames(data_input)
data_input <- left_join(data_input, metadata,by=c('id'='name'))

data_input$method1 <- factor(data_input$method1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))

data_input %>%
  ggplot(aes(x=PC1,y=PC2,fill=method1))+
  geom_point(size=4,shape = 21,alpha=0.8)+
  scale_fill_manual(values=c('#B5252E','#416997','#ADC482','#DAA83D'),
                    labels=c('AF4','EXO-CMDS','SEC','Plasma-control'))+
  labs(y=percentage[2],x = percentage[1],title = '',colour = "Method")+
  theme4_backline + themelegendnone+
  scale_x_continuous(limits = c(-20, 25), expand = c(0, 0),breaks = seq(-20, 25, by = 10), minor_breaks = seq(-20, 25, by = 5)) +  
  scale_y_continuous(limits = c(-20, 15), expand = c(0, 0),breaks = seq(-20, 20, by = 10), minor_breaks = seq(-20, 20, by = 5)) 

ggsave(paste0(wdplot,"b.pdf"), width = 7, height =6, units = "cm")



# c heatmap ---------------

corp <- cor(normalized_counts,method = 'spearman')#spearman  pearson

cor_matrix <- cor(normalized_counts,method = 'spearman')

annotation_col <- as.data.frame(rownames(cor_matrix))
colnames(annotation_col) <- 'sample'
annotation_col <- left_join(annotation_col,metadata,by=c('sample'='name'))
colnames(annotation_col) <- c('sample','Cluster','Frze')

annotation_col <- column_to_rownames(annotation_col,var='sample')

ann_colors = list(
  Cluster=c('AF4'='#B5252E','EXO-CMDS'='#416997','SEC'='#DAA83D','Plasma-control'='#ADC482'),
  Frze=c('Non-freeze-thaw'='#BEBC48','Freeze-thaw'='#5568B8FF') #'#3E6B9B','#DAA83D'
)    

library(pheatmap)
bk <- c(seq(0.6,0.7,by=0.01),seq(0.71,1,by=0.01))
pheatmap::pheatmap(
  corp,
  color = c(colorRampPalette(colors = c("#F5F3B1","#C55548"))(length(bk)/2),
            colorRampPalette(colors = c("#C55548","#1F1243"))(length(bk)/2)),#EBEBEB,BD6379 白 F9E855 48938B 3B0F4F
  clustering_method = "ward.D2",
  cutree_col = 4,
  cutree_row = 4,
  show_colnames = F,
  show_rownames = F,
  # annotation_row = annotation_row,
  #scale = 'column',
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  #  labels_col = annotation_col$Cluster,
  breaks= bk,
  cellwidth = 4, 
  cellheight = 4, 
  fontsize = 12, 
  border_color = NA,
  filename = paste0(wdplot,"c.pdf")
)
dev.off()

 # d vln --------------------

cor_matrix2 <- melt(cor_matrix)
colnames(cor_matrix2) <- c('Var1','Var2','value')
cor_matrix2 <- left_join(cor_matrix2,metadata,by=c('Var1'='name' ))
colnames(cor_matrix2) <- c("Var1","Var2","value","cate1", "storage1")
cor_matrix2 <- left_join(cor_matrix2,metadata,by=c('Var2'='name' ))
colnames(cor_matrix2) <- c("Var1","Var2","value","cate1", "storage1","cate2", "storage2")

cor_matrix2 <- dplyr::filter(cor_matrix2,cate1==cate2)

cv_values <- cor_matrix2 %>%
  dplyr::group_by(cate1) %>%
  dplyr::summarize(CV = cal_cv(value), .groups = 'drop')

mean_values <- cor_matrix2 %>%
  dplyr::group_by(cate1) %>%
  dplyr::summarize(MEAN = mean(value), .groups = 'drop')

cor_matrix2$cate1 <- factor(cor_matrix2$cate1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))

# r1 <- mutl_anova_test(cor_matrix2,'value',c('cate1'))
# r1

custom_labels <- data.frame(
  cate1 = c('AF4','EXO-CMDS','SEC','Plasma-control'),
  label = c('CV=0.0315','CV=0.0100','CV=0.0373','CV=0.0267')
)

ggplot(cor_matrix2, aes(x=cate1,y=as.numeric(value),fill=cate1))+
  geom_violin(trim = FALSE, color = "transparent") + 
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482')
                    ,name="")+
  scale_color_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482')
                     ,name="")+
  geom_boxplot(width = 0.1,colour = "white", size=0.5,outlier.size=0,  outlier.alpha = 0, position = position_dodge(width = 0.9),alpha=0.3)+
  labs(y='Spearman correlation coefficient',x = "",color = "")+  
  scale_y_continuous(limits = c(0.8, 1.05), expand = c(0,0), breaks = seq(0.8, 1, 0.1))+  
  theme2_noback+ themelegendnone+ themey90+
  geom_text(data = custom_labels, aes(x=cate1, y=1.02, label=label), color="black", size=2, vjust=0)


ggsave(paste0(wdplot,"d.pdf"), width = 9, height =6, units = "cm")


#e ev number------------------

result_count <- read_csv(paste0(wddata,'EV_nummber.csv'))
result_count <- result_count[,-c(1,2)]

plotdata <- NULL
for (i in colnames(result_count)) {
  result_count1 <- result_count[,i]
  sum1 <- result_count1[1,1]
  sum2 <- result_count1[2,1]
  sum3 <- result_count1[3,1]
  sum4 <- result_count1[4,1]
  sum5 <- sum(na.omit(result_count1[,i])) - sum1 - sum2 -sum3
  
  sum_bind <- rbind(sum1,sum2,sum3,sum5)
  sum_bind$cate <- c('1','2','3','>3')
  sum_bind$percentage <- (sum_bind[,i]/sum(na.omit(result_count1[,i])))*100
  sum_bind$id <- i
  colnames(sum_bind) <- c('value','cate','percentage','id')
  
  plotdata <- rbind(plotdata,sum_bind)
}

plotdata <- as.data.frame(plotdata)
plotdata$cate <- as.character(plotdata$cate)
colnames(plotdata) <- c('value', 'cate', 'percentage', 'id')

plotdata <- left_join(plotdata,metadata,by=c('id'='name'))
methods <- unique(plotdata$method1)

plotdata <- plotdata %>%
  dplyr::group_by(id) %>% 
  dplyr::mutate(total_value = sum(value, na.rm = TRUE), 
                percentage = (value / total_value) * 100,
                log10value = log10(value)) %>%
  ungroup()

r2 <- mutl_anova_test(plotdata,'value',c('method1','cate'))
r2 <- dplyr::filter(r2,X2==X4)

# View the updated plotdata
head(plotdata)

unique(plotdata$cate)
plotdata$method1 <- factor(plotdata$method1 , levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))
plotdata$cate <- factor(plotdata$cate , levels = c( "1" , "2"  ,"3" , ">3"))

mean_plotdata <- plotdata %>% 
  dplyr::group_by(method1,cate) %>%
  dplyr::summarize(mean = (mean(value)), .groups = 'drop')



ggplot(mean_plotdata, aes(x = cate, y = log10(mean), fill = method1)) +
  geom_bar(stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9) +
  labs(x = '', y = 'log10(mean(EV number))', title = '',fill='') +
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','EXO-CMDS','SEC','Plasma-control')) +
  scale_y_continuous(limits = c(0, 8), expand = c(0,0)) +
  theme1_noback + themelegendnone
#scale_x_discrete(labels = name1)

ggsave(paste0(wdplot,'e',".pdf"), width =6, height = 6, units = "cm")


# f evmarker ----------------------

marker <- c('CD9','CD63','CD81')

data2 <- normalized_counts[marker,]


data2.1 <- melt(as.matrix(log2(data2)))
colnames(data2.1) <- c('protein','name','value')
data2.1 <- left_join(data2.1,metadata,by=c('name'='name'))

unique(data2.1$method1)

value_median <- data2.1 %>%
  dplyr::group_by(method1, protein) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE), .groups = 'drop')


data2.1$protein <- factor(data2.1$protein, levels = marker)
value_median$protein <- factor(value_median$protein, levels = marker)
data2.1$method1 <- factor(data2.1$method1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))
value_median$method1 <- factor(value_median$method1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))

# anova <- mutl_anova_testplus2( data2.1,'value','method1')
# anova <- do.call(rbind, anova)
# anova <- dplyr::filter(anova,protein %in% marker)
# write.xlsx(anova,file=paste0(wdresult,'marker_evpercent_anova.xlsx'),
#            sheetName="logTMM-MARKERS", append=T, row.names=F,col.names = T,showNA = T)


ggplot(data2.1, aes(x=protein, y=value,fill=method1)) +
  geom_bar(data = value_median, aes(x=protein, y=median_value,fill=method1),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  #scale_fill_nejm()+
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','EXO-CMDS','SEC','Plasma-control'),name='')+
  geom_jitter(data=data2.1, aes(x=protein, y=value,fill=method1), position=position_dodge(0.8),size=0.5,alpha=0.9) +
  geom_errorbar(data=data2.1, aes(x=protein, y=value,fill=method1), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  theme1_noback+ themelegendnone+
  scale_y_continuous(limits = c(0, 20), expand = c(0,0),breaks = seq(0, 20, 5))+ 
  
  labs(x='',y='Mean log2(NPX)',fill = '')

ggsave(paste0(wdplot,"f.pdf"), width = 9, height =5.5, units = "cm")


# g cd9 cd81 cd63 marker percent --------------------


marker <- c("CD9", "CD63", "CD81")

markerdf <- readRDS(paste0(wddata,'protein_ev_count.RDS'))

markerdf$percent <- (markerdf$count/markerdf$evall)*100

markerdf <- dplyr::filter(markerdf,variable %in% marker)
markerdf <- left_join(markerdf,metadata,by='name')

percent_median <- markerdf %>%
  dplyr::group_by(method1, variable) %>%
  dplyr::summarise(median_percent = median(percent, na.rm = TRUE), .groups = 'drop')

markerdf$variable <- factor(markerdf$variable, levels = marker)
percent_median$variable <- factor(percent_median$variable, levels = marker)

colnames(markerdf)[1] <- 'protein'
# anova <- mutl_anova_testplus2( markerdf,'percent','method1')
# anova <- do.call(rbind, anova)
# write.xlsx(anova,file=paste0(wdresult,'marker_evpercent_anova.xlsx'),
#            sheetName="misev", append=T, row.names=F,col.names = T,showNA = T)


markerdf$method1 <- factor(markerdf$method1, levels =c('AF4','EXO-CMDS','SEC','Plasma-control'))
percent_median$method1 <- factor(percent_median$method1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))


ggplot(markerdf, aes(x=protein, y=percent, fill=method1)) +
  geom_bar(data = percent_median, aes(x=variable, y=median_percent,fill=method1),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','PBA-control'),name='')+
  geom_jitter(data=markerdf, aes(x=protein, y=percent), position=position_dodge(0.8),size=0.5,alpha=0.9) +
  geom_errorbar(data=markerdf, aes(x=protein, y=percent), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  scale_y_continuous(limits = c(0, 8), expand = c(0,0),breaks = seq(0, 8, 2))+ 
  theme1_noback+ themelegendnone+
  labs(x='',y='Percentage of detected EVs (%)',fill = '')

ggsave(paste0(wdplot,"g.pdf"), width = 18, height =6, units = "cm")


#  S1 a  ------------------


result_count <- read_csv(paste0(wddata,'EV_nummber.csv'))
result_count <- result_count[,-c(1,2)]

singleproteincount <- result_count[1,]
sum(singleproteincount)

plotdata <- NULL
for (i in colnames(result_count)) {
  result_count1 <- result_count[,i]
  sum1 <- result_count1[1,1]
  # sum2 <- result_count1[2,1]
  # sum3 <- result_count1[3,1]
  #sum4 <- result_count1[4,1]
  sum5 <- sum(na.omit(result_count1[,i])) - sum1 #- sum2 -sum3
  
  sum_bind <- rbind(sum1,sum5)
  sum_bind$cate <- c('Single protein expression','≥2 proteins expression')
  sum_bind$percentage <- (sum_bind[,i]/sum(na.omit(result_count1[,i])))*100
  sum_bind$id <- i
  colnames(sum_bind) <- c('value','cate','percentage','id')
  
  plotdata <- rbind(plotdata,sum_bind)
}
plotdata <- as.data.frame(plotdata)
plotdata$cate <- as.character(plotdata$cate)
colnames(plotdata) <- c('value', 'cate', 'percentage', 'id')

plotdata <- left_join(plotdata,metadata,by=c('id'='name'))
methods <- unique(plotdata$method1)

i <- 'AF4'
data2 <- dplyr::filter(plotdata, `method1` == i)
data2$num <- final_sequence <- c(rep(c(1, 1, 2, 2, 3, 3), times = 4), rep(c(4, 4, 5, 5, 6, 6), times = 2))
data2$name3 <- paste0(  data2$storage,'_',  data2$num )

p<- ggplot(data2, aes(x = id, y = as.numeric(percentage$AF_1a_1), fill = cate)) +
  geom_col(position = "stack", width = 0.8) +
  labs(x = '', y = 'Percentage(%)', title = i) +
  scale_fill_manual(values = c("#8f5362", "#e0a981", "#ecd09c", "#6888a5", "#6888a5")) +
 # scale_y_continuous(limits = c(0, 100.1), expand = c(0,0), breaks = seq(0, 100, 20)) +
  theme1_backline +  scale_x_discrete(labels = data2$name3)+themetext30+themelegendnone

p
ggsave(paste0(wdplot,'S1',i,".pdf"),plot =  p, width =5, height = 6, units = "cm")


i <- 'EXO-CMDS'
data2 <- dplyr::filter(plotdata, `method1` == i)
data2$num <- final_sequence <- c(rep(c(1, 1, 2, 2, 3, 3), times = 4), rep(c(4, 4, 5, 5, 6, 6), times = 2))
data2$name3 <- paste0(  data2$storage,'_',  data2$num )

# 生成图形
p<- ggplot(data2, aes(x = id, y = as.numeric(percentage$AF_1a_1), fill = cate)) +
  geom_col(position = "stack", width = 0.8) +
  labs(x = '', y = 'Percentage(%)', title = i) +
  scale_fill_manual(values = c("#8f5362", "#e0a981", "#ecd09c", "#6888a5", "#6888a5")) +
  scale_y_continuous(limits = c(0, 100.1), expand = c(0,0), breaks = seq(0, 100, 20)) +
  theme1_backline + 
  scale_x_discrete(labels = data2$name3)+themetext30+themelegendnone

p
ggsave(paste0(wdplot,'S1',i,".pdf"),plot =  p, width =5, height = 6, units = "cm")


i <- 'SEC'
data2 <- dplyr::filter(plotdata, `method1` == i)
data2$num <- final_sequence <- c(rep(c(1, 1, 2, 2, 3, 3), times = 4), rep(c(4, 4, 5, 5, 6, 6), times = 2))
data2$name3 <- paste0(  data2$storage,'_',  data2$num )

# 生成图形
p<- ggplot(data2, aes(x = id, y = as.numeric(percentage$AF_1a_1), fill = cate)) +
  geom_col(position = "stack", width = 0.8) +
  labs(x = '', y = 'Percentage(%)', title=i) +
  scale_fill_manual(values = c("#8f5362", "#e0a981", "#ecd09c", "#6888a5", "#6888a5")) +
  scale_y_continuous(limits = c(0, 100.1), expand = c(0,0), breaks = seq(0, 100, 20)) +
  theme1_backline + 
  scale_x_discrete(labels = data2$name3)+themetext30+themelegendnone

p
ggsave(paste0(wdplot,'S1',i,".pdf"),plot =  p, width =5, height = 6, units = "cm")



i <- 'Plasma-control'
data2 <- dplyr::filter(plotdata, `method1` == i)
data2$num <- final_sequence <- c(rep(c(1, 1, 2, 2, 3, 3), times = 4), rep(c(4, 4, 5, 5, 6, 6), times = 2))
data2$name3 <- paste0(  data2$storage,'_',  data2$num )

# 生成图形
p<- ggplot(data2, aes(x = id, y = as.numeric(percentage$AF_1a_1), fill = cate)) +
  geom_col(position = "stack", width = 0.8) +
  labs(x = '', y = 'Percentage(%)', title = i) +
  scale_fill_manual(values = c("#8f5362", "#e0a981", "#ecd09c", "#6888a5", "#6888a5")) +
  scale_y_continuous(limits = c(0, 100.1), expand = c(0,0), breaks = seq(0, 100, 20)) +
  theme1_backline + 
  scale_x_discrete(labels = data2$name3)+themetext30+themelegendnone

p
ggsave(paste0(wdplot,'S1',i,".pdf"),plot =  p, width =5, height = 6, units = "cm")

#  S1 b ------------------
marker <- c("CD63", "CD9", "CD81", "CD82", "CD47", "GNA12", "GNA13", "GNA14", "GNA11", "GNA15", 
            "TSAP6", "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGA8", 
            "ITGA9","ITGA10", "ITGA11", "ITGAE", "ITGAV", "ITGA2B", "ITGAD", "ITGAL", "ITGAX")

marker <- c("ITGAM","ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8", "TFR2", "LAMP1" ,
            "LAMP2", "SDC1", "SDC2", "SDC3", "SDC4", "BSG", "ADAM10", "GPC1", "NT5E", "CD59", 
            "TSG101", "CHMP1", "CHMP2", "CHMP3", "CHMP4", "CHMP5", "CHMP6", "CHMP7", "CHMP1A","CHMP1B", 
            "CHMP2A", "CHMP2B", "CHMP4A", "CHMP4B", "CHMP4C", "PDCD6IP", "VPS4A", "VPS4B", 
            "ARRDC1", "FLOT1", "FLOT2","CAV1", "CAV2", "CAV3", "SDCBP", "HSPA8", "HSP90AB1", 
            "ACTA1", "ACTA2", "ACTC1", "ACTB","ACTG1", "ACTG2", "TUBA1A", "TUBA1B", "TUBA1C",   "TUBA3C", "TUBA3D", "TUBA3E", "TUBA4A", 
            "TUBA8", "TUBB", "TUBB1", "TUBB2A", "TUBB2B", "TUBB3", "TUBB4A", "TUBB4B", "TUBB6", "TUBG1", "TUBG2", "GAPDH")

normalized_counts <- as.data.frame(normalized_counts)
data2 <- normalized_counts[marker,]
data2 <- na.omit(data2)


data2.1 <- melt(as.matrix(log2(data2)))
colnames(data2.1) <- c('protein','name','value')
data2.1 <- left_join(data2.1,metadata,by=c('name'))

unique(data2.1$method1)

value_median <- data2.1 %>%
  dplyr::group_by(method1, protein) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE), .groups = 'drop')

data2.1$protein <- factor(data2.1$protein, levels = marker)
value_median$protein <- factor(value_median$protein, levels = marker)
data2.1$method1 <- factor(data2.1$method1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))
value_median$method1 <- factor(value_median$method1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))

p1 <- ggplot(data2.1, aes(x=protein, y=value,fill=method1)) +
  geom_bar(data = value_median, aes(x=protein, y=median_value,fill=method1),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  #scale_fill_nejm()+
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','Plasma-control'),name='')+
  geom_jitter(data=data2.1, aes(x=protein, y=value,fill=method1), position=position_dodge(0.8),size=0.5,alpha=0.9) +
  geom_errorbar(data=data2.1, aes(x=protein, y=value,fill=method1), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  theme1_noback+ themelegendnone+
  scale_y_continuous(limits = c(0, 20), expand = c(0,0),breaks = seq(0, 20, 5))+  
  labs(x='',y='Mean log2(NPX)',fill = '')
p2 <- ggplot(data2.1, aes(x=protein, y=value,fill=method1)) +
  geom_bar(data = value_median, aes(x=protein, y=median_value,fill=method1),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  #scale_fill_nejm()+
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','Plasma-control'),name='')+
  geom_jitter(data=data2.1, aes(x=protein, y=value,fill=method1), position=position_dodge(0.8),size=0.5,alpha=0.9) +
  geom_errorbar(data=data2.1, aes(x=protein, y=value,fill=method1), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  theme1_noback+ themelegendnone+
  scale_y_continuous(limits = c(0, 20), expand = c(0,0),breaks = seq(0, 20, 5))+ 
  labs(x='',y='Mean log2(NPX)',fill = '')


markerdf <- readRDS(paste0(wddata,'protein_ev_count.RDS'))
markerdf$percent <- (markerdf$count/markerdf$evall)*100
markerdf <- dplyr::filter(markerdf,variable %in% marker)
markerdf <- left_join(markerdf,metadata,by='name')
percent_median <- markerdf %>%
  dplyr::group_by(method1, variable) %>%
  dplyr::summarise(median_percent = median(percent, na.rm = TRUE), .groups = 'drop')
markerdf$variable <- factor(markerdf$variable, levels = marker)
percent_median$variable <- factor(percent_median$variable, levels = marker)
colnames(markerdf)[1] <- 'protein'
markerdf$method1 <- factor(markerdf$method1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))
percent_median$method1 <- factor(percent_median$method1, levels = c('AF4','EXO-CMDS','SEC','Plasma-control'))

p3 <- ggplot(markerdf, aes(x=protein, y=percent, fill=method1)) +
  geom_bar(data = percent_median, aes(x=variable, y=median_percent,fill=method1),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','PBA-control'),name='')+
  geom_jitter(data=markerdf, aes(x=protein, y=percent), position=position_dodge(0.8),size=0.5,alpha=0.9) +
  geom_errorbar(data=markerdf, aes(x=protein, y=percent), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  scale_y_continuous(limits = c(0, 8), expand = c(0,0),breaks = seq(0, 8, 2))+ 
  theme1_noback+ themelegendnone+
  labs(x='',y='Percentage of detected EVs (%)',fill = '')

p4 <- ggplot(markerdf, aes(x=protein, y=percent, fill=method1)) +
  geom_bar(data = percent_median, aes(x=variable, y=median_percent,fill=method1),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','PBA-control'),name='')+
  geom_jitter(data=markerdf, aes(x=protein, y=percent), position=position_dodge(0.8),size=0.5,alpha=0.9) +
  geom_errorbar(data=markerdf, aes(x=protein, y=percent), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  scale_y_continuous(limits = c(0, 8), expand = c(0,0),breaks = seq(0, 8, 2))+ 
  theme1_noback+ themelegendnone+
  labs(x='',y='Percentage of detected EVs (%)',fill = '')

p1/p2/p3/p4+plot_layout(heights = c(1,1,1,1))

ggsave(paste0(wdplot,'S1b.pdf'), width = 20, height = 25, units = "cm")



