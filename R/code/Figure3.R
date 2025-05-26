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

# a --------------------------
load(paste0(wddata,'F3data1.R'))
metadata <- open_xlsx(paste0(wddata,'PBA metadata.xlsx'),1)

library(edgeR)
data <- DGEList(counts = oneisoed)
TMM <- calcNormFactors(data, method="TMM")
oneiso <- cpm(TMM, log = F)

wd_method <- 'PCA'

data_pca2 <- prcomp(t(as.matrix(oneiso)), scale=T)
percentage <- round(data_pca2$sdev^2 / sum(data_pca2$sdev^2) * 100,2)
percentage <- paste(colnames(as.data.frame(data_pca2[["x"]])),"(", paste(as.character(percentage), "%", ")", sep=""))

data_input <-as.data.frame(data_pca2[["x"]])
data_input$id <- rownames(data_input)
data_input <- left_join(data_input, metadata,by=c('id'='rawname'))

data_input$method1 <- factor(data_input$`Isolation method`, levels = c('AF4','Exo-CMDS','SEC','Plasma-control'))

data_input %>%
  ggplot(aes(x=PC1,y=PC2,fill=method1))+
  geom_point(size=3,shape = 21,alpha=0.8)+
  scale_fill_manual(values=c('#B5252E','#416997','#DAA83D','#ADC482'),
                    labels=c('AF4','Exo-CMDS','SEC','Plasma-control'))+
  labs(y=percentage[2],x = percentage[1],title = '',colour = "Method")+
  theme4_backline + themelegendnone+
  scale_x_continuous(limits = c(-20, 25), expand = c(0, 0),breaks = seq(-20, 25, by = 10), minor_breaks = seq(-20, 25, by = 5)) +  
  scale_y_continuous(limits = c(-20, 15), expand = c(0, 0),breaks = seq(-20, 20, by = 10), minor_breaks = seq(-20, 20, by = 5)) +
  themelegendright

ggsave(paste0(wdplot,"a.pdf"), width = 9, height =6, units = "cm")



# b --------------------------------------------

corp <- cor(oneiso,method = 'spearman')

cor_matrix <- cor(oneiso,method = 'spearman')

annotation_col <- as.data.frame(rownames(cor_matrix))
colnames(annotation_col) <- 'sample'
annotation_col <- left_join(annotation_col,metadata,by=c('sample'='rawname'))
annotation_col <- annotation_col[,c('sample','Isolation method','Freeze-thaw condition')]
colnames(annotation_col) <- c('sample','Cluster','Frze')

annotation_col <- column_to_rownames(annotation_col,var='sample')

ann_colors = list(
  Cluster=c('AF4'='#B5252E','Exo-CMDS'='#416997','SEC'='#DAA83D','Plasma-control'='#ADC482'),
  Frze=c('Non-freeze-thaw'='#BEBC48','Freeze-thaw'='#5568B8FF') 
)   
library(pheatmap)


bk <- c(seq(0.6,0.7,by=0.01),seq(0.71,1,by=0.01))
pheatmap::pheatmap(
  corp,
  color = c(colorRampPalette(colors = c("#F5F3B1","#C55548"))(length(bk)/2),
            colorRampPalette(colors = c("#C55548","#1F1243"))(length(bk)/2)),
  clustering_method = "ward.D2",
  cutree_col = 4,
  cutree_row = 4,
  show_colnames = F,
  show_rownames = F,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  breaks= bk,
  cellwidth = 4, 
  cellheight = 4, 
  fontsize = 12, 
  border_color = NA,
  filename = paste0(wdplot,"b.pdf")
)
dev.off()

corp <- cor(oneiso,method = 'spearman')
corp <- melt(corp)
colnames(corp) <- c('X1','X2','value')
corp <- left_join(corp,metadata,by=c('X1'='rawname'))
corp <- left_join(corp,metadata,by=c('X2'='rawname'))
corp2 <- corp %>% dplyr::filter(`Isolation method.x` == 'Plasma-control' & `Isolation method.y` == 'SEC') 

mean(corp2$value)# 0.841716

# c  --------------------

cor_matrix2 <- melt(cor_matrix)
colnames(cor_matrix2) <- c('Var1','Var2','value')
cor_matrix2 <- left_join(cor_matrix2,metadata,by=c('Var1'='rawname' ))
cor_matrix2 <- cor_matrix2[,c(1,2,3,5,4)]
colnames(cor_matrix2) <- c("Var1","Var2","value","cate1", "storage1")
cor_matrix2 <- left_join(cor_matrix2,metadata,by=c('Var2'='rawname' ))
cor_matrix2 <- cor_matrix2[,c(1,2,3,4,5,7,6)]
colnames(cor_matrix2) <- c("Var1","Var2","value","cate1", "storage1","cate2", "storage2")

cor_matrix2 <- dplyr::filter(cor_matrix2,cate1==cate2)

cv_values <- cor_matrix2 %>%
  dplyr::group_by(cate1) %>%
  dplyr::summarize(CV = cal_cv(value), .groups = 'drop')
cv_values$CV <- round(cv_values$CV,3)

mean_values <- cor_matrix2 %>%
  dplyr::group_by(cate1) %>%
  dplyr::summarize(MEAN = mean(value), .groups = 'drop')

cor_matrix2$cate1 <- factor(cor_matrix2$cate1, levels = c('AF4','Exo-CMDS','SEC','Plasma-control'))


ggplot(cor_matrix2, aes(x=cate1,y=as.numeric(value),fill=cate1))+
  geom_violin(trim = T, color = "transparent") + 
  
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482')
                    ,name="")+
  scale_color_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482')
                     ,name="")+
  geom_boxplot(width = 0.1,colour = "white", size=0.5,outlier.size=0,  outlier.alpha = 0, position = position_dodge(width = 0.9),alpha=0.3)+
  
  
  labs(y='Spearman correlation coefficient',x = "",color = "")+ 
  scale_y_continuous(limits = c(0.8, 1.05), expand = c(0,0), breaks = seq(0.8, 1, 0.1))+
  theme2_noback+ themelegendnone+ themey90+
  geom_text(data = cv_values, aes(x=cate1, y=1.02, label=paste0('CV=',CV)), color="black", size=2, vjust=0)


ggsave(paste0(wdplot,"c.pdf"), width = 6.5, height =5.7, units = "cm")

#d --------------------------------------

load(paste0(wddata,'F3data2.R'))


plotdata <- NULL
for (i in colnames(result_count)) {
  result_count1 <- result_count[,i]
  sum1 <- result_count1[1]
  sum2 <- result_count1[2]
  sum3 <- result_count1[3]
  sum5 <- sum(result_count1,na.rm = T) - sum1 - sum2 -sum3 
  
  sum_bind <- as.data.frame(rbind(sum1,sum2,sum3,sum5))
  colnames(sum_bind) <- 'number'
  sum_bind$procate <- c('Single protein expression','2 proteins expression','3 proteins expression','≥4 proteins expression')
  
  sum_bind$percentage <- (sum_bind[,1]/sum(na.omit(result_count1)))*100
  sum_bind$id <- i

  plotdata <- rbind(plotdata,sum_bind)
}
plotdata <- left_join(plotdata,metadata,by=c('id'='rawname'))
plotdata$procate <- factor(plotdata$procate, levels = c( "Single protein expression","2 proteins expression","3 proteins expression","≥4 proteins expression"   ))

value_median <- plotdata %>%
  dplyr::group_by(`Isolation method`, procate) %>%
  dplyr::summarise(median_value = median(percentage, na.rm = TRUE), .groups = 'drop')  

plotdata %>% dplyr::group_by(`Isolation method`,procate) %>% dplyr::summarise(meanv=mean(percentage))
plotdata %>% dplyr::group_by(procate) %>% dplyr::summarise(meanv=mean(percentage))

library(ggbreak)

ggplot(plotdata, aes(x=`Isolation method`, y=percentage,fill=procate)) +
  geom_bar(data = value_median, aes(x=`Isolation method`, y=median_value,fill=procate),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  scale_fill_manual(values = c('#475272','#68A39D','#D9B381','#ADC482'),name='')+
  geom_jitter(data=plotdata, aes(x=`Isolation method`, y=percentage,fill=procate), position=position_dodge(0.8),size=0.5,alpha=0.9) +
  geom_errorbar(data=plotdata, aes(x=`Isolation method`, y=percentage,fill=procate), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  theme1_noback+ themelegendnone+
  scale_y_break(breaks = c(15,75),space = 0.2,scales = 0.4,expand = c(0,0))+themelegendright+
  labs(x='',y='Percentage (%)',fill = '')+themetext30
ggsave(paste0(wdplot,"d.pdf"), width = 15, height =7, units = "cm")

# f-------------------

marker <- c('CD9','CD63','CD81','ITGA4','ITGB1','ADAM10','LAMP1','GPC1','NT5E','ITGA1')

data2 <- oneisoed[marker,]

data2.1 <- melt(as.matrix((data2)))
colnames(data2.1) <- c('protein','name','value')
data2.1 <- left_join(data2.1,metadata,by=c('name'='rawname'))

unique(data2.1$`Isolation method`)

value_median <- data2.1 %>%
  dplyr::group_by(`Isolation method`, protein) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE), .groups = 'drop')


data2.1$protein <- factor(data2.1$protein, levels = marker)
value_median$protein <- factor(value_median$protein, levels = marker)
data2.1$`Isolation method` <- factor(data2.1$`Isolation method`, levels = c('AF4','Exo-CMDS','SEC','Plasma-control'))
value_median$`Isolation method` <- factor(value_median$`Isolation method`, levels = c('AF4','Exo-CMDS','SEC','Plasma-control'))


data2.1.2 <- data2.1 %>% dplyr::group_by(`Isolation method`,protein) %>% dplyr::summarise(meanvalue=max(value))

ggplot(data2.1, aes(x=protein, y=value*100,fill=`Isolation method`)) +
  geom_bar(data = value_median, aes(x=protein, y=median_value*100,fill=`Isolation method`),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','EXO-CMDS','SEC','Plasma-control'),name='')+
  geom_jitter(data=data2.1, aes(x=protein, y=value*100,fill=`Isolation method`), position=position_dodge(0.8),size=0.1,alpha=0.9) +
  geom_errorbar(data=data2.1, aes(x=protein, y=value*100,fill=`Isolation method`), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  theme1_noback+ themelegendnone+ themetext30+
 scale_y_continuous(limits = c(0, 10), expand = c(0,0),breaks = seq(0, 10, 2))+ 
  
  labs(x='',y='Percentage of detected EVs (%)',fill = '')

ggsave(paste0(wdplot,"f.pdf"), width = 9, height =5.5, units = "cm")


#  S1  ------------------


marker <- c( "ITGA2", "ITGA3", "ITGA5", "ITGA6", "ITGA8", 
            "ITGA9", "ITGA11", "ITGAV", "ITGA2B", "ITGAL", "ITGAX")
marker <- c("ITGAM", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8",
            "LAMP2", "SDC1", "CAV1")
data2 <- oneisoed[marker,]


data2.1 <- melt(as.matrix((data2)))
colnames(data2.1) <- c('protein','name','value')
data2.1 <- left_join(data2.1,metadata,by=c('name'='rawname'))
data2.1 <- na.omit(data2.1)
unique(data2.1$`Isolation method`)

value_median <- data2.1 %>%
  dplyr::group_by(`Isolation method`, protein) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE), .groups = 'drop')


data2.1$protein <- factor(data2.1$protein, levels = marker)
value_median$protein <- factor(value_median$protein, levels = marker)
data2.1$`Isolation method` <- factor(data2.1$`Isolation method`, levels = c('AF4','Exo-CMDS','SEC','Plasma-control'))
value_median$`Isolation method` <- factor(value_median$`Isolation method`, levels = c('AF4','Exo-CMDS','SEC','Plasma-control'))

data2.1.2 <- data2.1 %>% dplyr::group_by(`Isolation method`,protein) %>% dplyr::summarise(meanvalue=max(value))

ggplot(data2.1, aes(x=protein, y=value*100,fill=`Isolation method`)) +
  geom_bar(data = value_median, aes(x=protein, y=median_value*100,fill=`Isolation method`),
           stat = "identity", position = position_dodge2(width=1),width = 0.8,color='black',linewidth=0.3,alpha=0.9)+
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','EXO-CMDS','SEC','Plasma-control'),name='')+
  geom_jitter(data=data2.1, aes(x=protein, y=value*100,fill=`Isolation method`), position=position_dodge(0.8),size=0.1,alpha=0.9) +
  geom_errorbar(data=data2.1, aes(x=protein, y=value*100,fill=`Isolation method`), position=position_dodge(0.8),stat="boxplot", width=0.65, coef=0) +
  theme1_noback+ themelegendnone+ themetext30+
  scale_y_continuous(limits = c(0, 2), expand = c(0,0),breaks = seq(0, 2, 0.5))+ 
  
  labs(x='',y='Percentage of detected EVs (%)',fill = '')

ggsave(paste0(wdplot,"s1.pdf"), width = 10, height =7, units = "cm")
