wd <- "~/R/"
setwd(wd)

wdcode <- paste0("code/") 
ver <- 'Figure5/'
wdresult <- paste0("result/",ver)
wddata <- paste0("data/",ver)
wdplot <- paste0("plot/",ver)
source(paste0(wdcode,"function.R"))
load_packages()
load_plot_packages() 

log_data <- readRDS(paste0(wddata,"log_data.rds"))
metadata <- readRDS(paste0(wddata,"metadata.RDS"))
#saveRDS(metadata,file=paste0(wddata,'metadata.RDS'))
# b log2intensity ----------------------

log_data[is.na(log_data)] <- 0
data_af <- log_data[,1:6]
data_mf <- log_data[,7:12]
data_sec <- log_data[,13:18]

mean <- as.data.frame(rowMeans(log_data,na.rm = T))
mean$protein <- rownames(mean)
colnames(mean) <- c('Mean','protein')
mean <- mean[order(-mean$Mean), ]
mean$rank <- seq(1,length(mean$Mean),1)

data_af <- melt(as.matrix(data_af))
data_mf <- melt(as.matrix(data_mf))
data_sec <- melt(as.matrix(data_sec))
colnames(data_af) <- c("Var1","Var2","value")
colnames(data_mf) <- c("Var1","Var2","value")
colnames(data_sec) <- c("Var1","Var2","value")

af <- data_af %>% 
  dplyr::group_by(Var1) %>%
  dplyr::summarise(
    mean_value = mean(value),
    lower_ci = mean(value) - qt(0.975, df=length(value)-1)*sd(value)/sqrt(length(value)),
    upper_ci = mean(value) + qt(0.975, df=length(value)-1)*sd(value)/sqrt(length(value))
  )%>%
  mutate(cate ='AF4')
mf <- data_mf %>% 
  dplyr::group_by(Var1) %>%
  dplyr::summarise(
    mean_value = mean(value),
    lower_ci = mean(value) - qt(0.975, df=length(value)-1)*sd(value)/sqrt(length(value)),
    upper_ci = mean(value) + qt(0.975, df=length(value)-1)*sd(value)/sqrt(length(value))
  )%>%
  mutate(cate ='EXO-CMDS')
sec <- data_sec %>% 
  dplyr::group_by(Var1) %>%
  dplyr::summarise(
    mean_value = mean(value),
    lower_ci = mean(value) - qt(0.975, df=length(value)-1)*sd(value)/sqrt(length(value)),
    upper_ci = mean(value) + qt(0.975, df=length(value)-1)*sd(value)/sqrt(length(value))
  ) %>%
  mutate(cate ='SEC')


bind <- rbind(af,mf,sec)
bind <- left_join(bind,mean,by=c('Var1'='protein'))
bind <- bind %>% 
  mutate(note = ifelse(Var1 %in%c('APOE','CD9','CD81','CD63','TSG101','PDCD6IP',"FLOT1"), Var1, NA))



ggplot(bind, aes(x=rank, y=mean_value, color=cate)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = cate), linetype = 1, size = 0.02, alpha = 0.2) +
  geom_smooth(aes(color = cate), alpha = 1.2, method = "loess", se = FALSE, span = 0.25, size = 0.35) +
  scale_color_manual(values = c("#d37b6d", "#6888a5", "#DAA83D", "#706d94")) +
  scale_fill_manual(values = c("#d37b6d", "#6888a5", "#DAA83D", "#706d94")) +
  theme_bw() +
  labs(x = 'Rank', y = 'log2(Intensity)', color = '') +
  scale_x_continuous(limits = c(0, 1550), expand = c(0, 0)) +   
  scale_y_continuous(limits = c(10, 38), expand = c(0, 0)) +   
  theme2_noback + themelegendnone +
  geom_vline(data = bind %>% filter(!is.na(note)), aes(xintercept = rank), linetype = "dashed", color = "grey", size = 0.1) +
  geom_text(data = bind %>% filter(!is.na(note)), aes(x = rank, y = 35, label = note), color = "black",
            alpha=0.3, angle = 0, hjust = 0.5, size = 1.7)
ggsave(paste0(wdplot,"b intensity.pdf"),  width = 9, height = 6, units = "cm")

# c freeze thaw ---------------------------------

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

d1 <- melt(as.matrix(data_af)) 
d2 <- melt(as.matrix(data_mf)) 
d3 <- melt(as.matrix(data_sec)) 
colnames(d1) <- c("Var1","Var2","value")
colnames(d2) <- c("Var1","Var2","value")
colnames(d3) <- c("Var1","Var2","value")
d <- left_join(rbind(d1,d2,d3),metadata,by=c('Var2'='new_name1' ))


ggplot(d, aes(x=method1,y=as.numeric(value),fill=storage))+
  geom_violin(trim = FALSE, color = "transparent",alpha=0.9) +  
  geom_boxplot(width = 0.2,colour = "white", outlier.size=0, position = position_dodge(width = 0.9),size=0.3,alpha=0.2)+
  scale_fill_manual(values = c('#5568B8','#BEBC48'),   
                    name = "",
                    labels = c("Freeze-thaw",'Non-freeze-thaw')) +
  labs(y='MS - log2(Intensity)',x = "",color = "")+  
  theme2_noback+ themelegendnone+
  scale_x_discrete(labels = c('AF4','Exo-CMDS','SEC','Ori-plasma'))
ggsave(paste0(wdplot,"c.pdf"), width = 8, height =6, units = "cm")





# d,e,f evmarker ----------------------

evmarker <- c('CD9','CD63','CD81')
marker1 <- c('APOA1','APOA2','APOB')
marker2 <- c('TGFBI','HSPA13','LDHA','LDHB')

msmarker <- function(marker1)    {
    log_data[is.na(log_data)] <- 0
    data_af <- log_data[,1:6]
    data_mf <- log_data[,7:12]
    data_sec <- log_data[,13:18]
    
    data_af2 <- rownames_to_column(data_af)
    data_af2 <- filter(data_af2,rowname %in% marker1)%>% 
      column_to_rownames(var = "rowname") %>% 
      as.matrix() %>% 
      melt()
    colnames(data_af2) <- c('gene','id','value')
    
    data_mf2 <- rownames_to_column(data_mf)
    data_mf2 <- filter(data_mf2,rowname %in% marker1) %>% 
      column_to_rownames(var = "rowname") %>% 
      as.matrix() %>% 
      melt()
    colnames(data_mf2) <- c('gene','id','value')
    
    data_sec2 <- rownames_to_column(data_sec)
    data_sec2 <- filter(data_sec2,rowname %in% marker1)%>% 
      column_to_rownames(var = "rowname") %>% 
      as.matrix() %>% 
      melt()
    colnames(data_sec2) <- c('gene','id','value')
    
    
    bind_data <- rbind(data_af2,data_mf2)
    bind_data <- rbind(bind_data,data_sec2)
    bind_data <- left_join(bind_data,metadata,by=c('id'='new_name1'))
    
    median_value1 <- data_af2 %>%
      dplyr::group_by(gene) %>%
      dplyr::summarize(
        value = ifelse(sum(value == 0) >= 3, 0, median(value)),
        method1 = 'AF4'
      )
    
    median_value2 <- data_mf2 %>%
      dplyr::group_by(gene) %>%
      dplyr::summarize(
        value = ifelse(sum(value == 0) >= 3, 0, median(value)),
        method1 = 'MF'
      )
    
    median_value3 <- data_sec2 %>%
      dplyr::group_by(gene) %>%
      dplyr::summarize(
        value = ifelse(sum(value == 0) >= 3, 0, median(value)),
        method1 = 'SEC'
      )
    
    median_value <- rbind(median_value1,median_value2)
    median_value <- rbind(median_value,median_value3)
    
    bind_data$gene <-  factor(bind_data$gene, levels = marker1)
    median_value$gene <- factor(median_value$gene , levels = marker1)
    
    
    e <- ggplot(bind_data, aes(x=gene, y=value,fill=method1)) +
      geom_bar(data = median_value, aes(x=gene, y=value,fill=method1),
               stat = "identity", position = position_dodge2(width=1),width = 0.65,color='black',linewidth=0.3)+
      scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','PBA-control'),name='')+
      geom_errorbar(data=bind_data, aes(x=gene, y=value,fill=method1), position=position_dodge(0.65),
                    stat="boxplot", width=0.6, coef=0,  inherit.aes = F) +
      geom_jitter(data=bind_data, aes(x=gene, y=value,fill=method1), position=position_dodge(0.65), color = "black",size=1,shape=20) +
      scale_y_continuous(limits = c(0, 40), expand = c(0,0))+  
      theme1_noback+ themelegendnone+
      labs(x='',y='log2(Intensity)',title = '')
    return(e)
}
d <- msmarker(evmarker)
e <- msmarker(marker1)
f <- msmarker(marker2)

pb <- d+e+f+ plot_layout(widths = c(3,3,4))

ggsave(paste0(wdplot,"def.pdf"), pb,width = 15, height =6, units = "cm")





#g pca --------------------------------------
log_data_nona <- log_data_nona <- na.omit(log_data)
data_pca2 <- prcomp(t(as.matrix(log_data_nona)), scale=T)
percentage <- round(data_pca2$sdev^2 / sum(data_pca2$sdev^2) * 100,2)
percentage <- paste(colnames(as.data.frame(data_pca2[["x"]])),"(", paste(as.character(percentage), "%", ")", sep=""))

data_input <-as.data.frame(data_pca2[["x"]])
data_input$id <- rownames(data_input)
data_input <- left_join(data_input, metadata,by=c('id'='new_name1'))

data_input %>%
  ggplot(aes(x=PC1,y=PC2,fill=method1))+
  geom_point(size=5,shape = 21,alpha=0.8)+
  scale_fill_manual(values=c('#B5252E','#416997','#DAA83D','#ADC482'))+
  labs(y=percentage[2],x = percentage[1],title = '',fill = "Method")+
  theme4_backline+ themelegendnone +
  scale_x_continuous(limits = c(-60, 35), expand = c(0, 0),breaks = seq(-60, 35, by = 20), minor_breaks = seq(-60, 35, by = 10)) +
  scale_y_continuous(limits = c(-30, 30), expand = c(0, 0),breaks = seq(-30, 30, by = 20), minor_breaks = seq(-30, 30, by = 10))
ggsave(paste0(wdplot,"g.pdf"), width = 6, height =6, units = "cm")



#h spearman correlation ----------------------------

log_data[log_data==0] <- NA

p <- na.omit(log_data)
corp <- cor(p,method = 'spearman')

annotation_col <- as.data.frame(rownames(corp))
colnames(annotation_col) <- 'sample'
annotation_col$Cluster <- c(rep('AF4',6),rep('MF',6),rep('SEC',6))
annotation_col$Freeze <- c(rep('Non-freeze-thaw',2),rep('Freeze-thaw',2),rep('Freeze-thaw',2),rep('Non-freeze-thaw',2),rep('Freeze-thaw',2),
                           rep('Freeze-thaw',2),rep('Non-freeze-thaw',2),rep('Freeze-thaw',2),rep('Freeze-thaw',2))
annotation_col <- column_to_rownames(annotation_col,var='sample')

ann_colors = list( 
  Cluster=c('AF4'='#B5252E','MF'='#416997','SEC'='#DAA83D'),
  Freeze=c('Non-freeze-thaw'='#BEBC48','Freeze-thaw'='#5568B8')
)    

bk <- c(seq(0.2,0.59,by=0.1),seq(0.6,1,by=0.1))
heatmap <-  pheatmap::pheatmap(
  corp,
  color = c(colorRampPalette(colors = c("#F5F3B1","#C55548"))(length(bk)/2),
            colorRampPalette(colors = c("#C55548","#1F1243"))(length(bk)/2)),
  clustering_method = "ward.D2",
  cutree_col = 3,
  cutree_row = 3,
  show_colnames = F,
  show_rownames = F,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  breaks= bk,
  cellwidth = 12, 
  cellheight = 12, 
  fontsize = 12
)
ggsave(paste0(wdplot,'h.pdf'),plot = heatmap, width = 20, height =16, units = "cm")


#i CV corr -----------------------
data2 <- melt(as.matrix(corp))
colnames(data2) <- c('Var1','Var2','value')
data2$cate1 <- sub("_.*", "", data2$Var1)
data2$cate1_storage <- sub("^[^_]*_([^_]*)_.*", "\\1", data2$Var1)
data2$cate2 <- sub("_.*", "", data2$Var2)
data2$cate2_storage <- sub("^[^_]*_([^_]*)_.*", "\\1", data2$Var2)

data2 <- dplyr::filter(data2,cate1==cate2)


data2 <- data2 %>%
  mutate(cate1 = case_when(
    cate1 == "AF4" ~ "AF4",
    cate1 == "MF" ~ "Exo-CMDS",
    cate1 == "SEC" ~ "SEC",
    TRUE ~ cate1
  ))
unique(data2$cate1)
data2$cate1 <- factor(data2$cate1, levels = c('AF4','Exo-CMDS','SEC'))

mean_values <- data2 %>%
  dplyr::group_by(cate1) %>%
  dplyr::summarize(CV = mean(value), .groups = 'drop')

cv_values <- data2 %>%
  dplyr::group_by(cate1) %>%
  dplyr::summarize(CV = cal_cv(value), .groups = 'drop')

# r1 <- mutl_anova_test(data2,'value',c('cate1'))
# r1

custom_labels <- data.frame(
  cate1 = c('AF4','Exo-CMDS','SEC'),
  label = c('CV=0.0197','CV=0.00756','CV=0.0211')
)


ggplot(data2, aes(x=cate1,y=as.numeric(value),fill=cate1))+
  geom_violin(trim = FALSE, color = "transparent") + 
  
  scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482')
                    ,name="")+
  scale_color_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482')
                     ,name="")+
  geom_boxplot(width = 0.2,colour = "white", size=0.3,outlier.size=0,  outlier.alpha = 0, position = position_dodge(width = 0.9),alpha=0.3)+
  labs(y='Spearman correlation coefficient',x = "",color = "")+ 
  scale_y_continuous(limits = c(0.9, 1.06), expand = c(0,0), breaks = seq(0.9, 1, 0.05))+ 
  theme2_noback+ themelegendnone+ themey90+
  geom_text(data = custom_labels, aes(x=cate1, y=1.04, label=label), color="black", size=2, vjust=0.3,angle=-20)
ggsave(paste0(wdplot,"i.pdf"), width = 7, height =6, units = "cm")



### DEPs -----------------------------------

dda_gene_cate <-  readRDS(paste0(wddata, 'method_detected_only.RDS'))
setdiff(rownames(log_data),c(dda_gene_cate[['af_spec']],dda_gene_cate[['mf_spec']], dda_gene_cate[['sec_spec']]))
other_gene <- setdiff(rownames(log_data),dda_gene_cate[['af_spec']])
other_gene <- setdiff(other_gene,dda_gene_cate[['mf_spec']])
other_gene <- setdiff(other_gene,dda_gene_cate[['sec_spec']])#641
impute_data1 <- log_data[other_gene,]

impute_data1 <- as.matrix(impute_data1)
impute_data1[impute_data1==0] <- NA

library(limma)
impute_data2 <- cbind(impute_data1,impute_data1[,1:6])
colnames(impute_data2) <- 1:24
sample_info <- data.frame(
  condition = factor(c(rep('control',18),rep('treat',6))),
  row.names = colnames(impute_data2)
)
design <- model.matrix(~ 0 + sample_info$condition)  
colnames(design) <- levels(sample_info$condition) 
design
fit <- lmFit(impute_data2, design)
contrast_matrix <- makeContrasts(treated_vs_control = treat - control, levels = design)
contrast_matrix
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_contrasts <- eBayes(fit_contrasts)
results <- topTable(fit_contrasts, number=Inf, adjust="BH", sort.by="P")
head(results)
results1 <- as.data.frame(results) %>% 
  mutate(cate = case_when(
    P.Value < 0.05 & logFC > 0 ~ "up",     
    P.Value < 0.05 & logFC < 0 ~ "down",   
    TRUE ~ "NA"                                
  ))

p1 <- ggplot(results1, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = cate), size = 1.5,alpha=0.7) +
  scale_color_manual(values = c("#1F1243", "grey","#B5252E", "red")) +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)", title = "AF4") +
  theme2_noback+themelegendnone
p1




impute_data2 <- cbind(impute_data1,impute_data1[,7:12])
colnames(impute_data2) <- 1:24
sample_info <- data.frame(
  condition = factor(c(rep('control',18),rep('treat',6))),
  row.names = colnames(impute_data2)
)
design <- model.matrix(~ 0 + sample_info$condition)  
colnames(design) <- levels(sample_info$condition) 
design
fit <- lmFit(impute_data2, design)
contrast_matrix <- makeContrasts(treated_vs_control = treat - control, levels = design)
contrast_matrix
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_contrasts <- eBayes(fit_contrasts)
results <- topTable(fit_contrasts, number=Inf, adjust="BH", sort.by="P")
head(results)
results2 <- as.data.frame(results) %>% 
  mutate(cate = case_when(
    P.Value < 0.05 & logFC > 0 ~ "up",    
    P.Value < 0.05 & logFC < 0 ~ "down",  
    TRUE ~ "NA"                               
  ))
p2 <- ggplot(results2, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = cate), size = 1.5,alpha=0.7) +
  scale_color_manual(values = c("#1F1243", "grey","#B5252E", "red")) +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)", title = "Exo-CMDS") +
  theme2_noback+themelegendnone
p2


impute_data2 <- cbind(impute_data1,impute_data1[,13:18])
colnames(impute_data2) <- 1:24
sample_info <- data.frame(
  condition = factor(c(rep('control',18),rep('treat',6))),
  row.names = colnames(impute_data2)
)
design <- model.matrix(~ 0 + sample_info$condition)
colnames(design) <- levels(sample_info$condition) 
design
fit <- lmFit(impute_data2, design)
contrast_matrix <- makeContrasts(treated_vs_control = treat - control, levels = design)
contrast_matrix
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_contrasts <- eBayes(fit_contrasts)
results <- topTable(fit_contrasts, number=Inf, adjust="BH", sort.by="P")
head(results)
results3 <- as.data.frame(results) %>% 
  mutate(cate = case_when(
    P.Value < 0.05 & logFC > 0 ~ "up",    
    P.Value < 0.05 & logFC < 0 ~ "down",  
    TRUE ~ "NA"                            
  ))
p3 <- ggplot(results3, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = cate), size = 1.5,alpha=0.7) +
  scale_color_manual(values = c("#1F1243", "grey","#B5252E", "red")) +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)", title = "SEC") +
  theme2_noback+themelegendright
p3

p1+p2+p3

write.xlsx(results1,file=paste0(wdresult,"DEPs.xlsx"), sheetName="AF4", append=T, row.names=F,col.names = T,showNA = T)
write.xlsx(results2,file=paste0(wdresult,"DEPs.xlsx"), sheetName="EXO-CMDS", append=T, row.names=F,col.names = T,showNA = T)
write.xlsx(results3,file=paste0(wdresult,"DEPs.xlsx"), sheetName="SEC", append=T, row.names=F,col.names = T,showNA = T)

deplist <- list()
deplist[['AF4_detected only']] <- dda_gene_cate[['af_spec']]
deplist[['AF4_elevated']] <- results1 %>% dplyr::filter(cate == 'up') %>% row.names
deplist[['CMDS_detected only']] <- dda_gene_cate[['mf_spec']]
deplist[['CMDS_elevated']] <- results2 %>% dplyr::filter(cate == 'up') %>% row.names
deplist[['SEC_detected only']] <- dda_gene_cate[['sec_spec']]
deplist[['SEC_elevated']] <-results3 %>% dplyr::filter(cate == 'up') %>% row.names


saveRDS(deplist,file=paste0(wdresult,'MS_dep.RDS'))


#f 三角图 ----------------------------------
deplist <- readRDS(file=paste0(wdresult,'MS_dep.RDS'))
library(ggtern)

daaf <- log_data[dda_gene_cate[['data_afgene']],1:6] %>% 
  mutate(AF4 = rowMeans(across(1:6),na.rm=T) )
damf <- log_data[dda_gene_cate[['data_mfgene']],7:12]%>% 
  mutate(`Exo-CMDS` = rowMeans(across(1:6),na.rm=T) )
dasec <- log_data[dda_gene_cate[['data_secgene']],13:18]%>% 
  mutate(SEC = rowMeans(across(1:6),na.rm=T) )

names(dda_gene_cate) <- c("data_afgene" , "data_mfgene" , "data_secgene" ,
                          "overlapgene","AF4_detected_only","CMDS_detected_only","SEC_detected_only" )


log_data[is.na(log_data)] <- 0
data_af <- log_data[,1:6]
data_mf <- log_data[,7:12]
data_sec <- log_data[,13:18]
data_af <- data_af[rowSums(data_af > 0) >= 4, ]%>% 
  mutate(AF4 = rowMeans(across(1:6),na.rm=T) )
data_mf <- data_mf[rowSums(data_mf > 0) >= 4, ]%>% 
  mutate(`Exo-CMDS` = rowMeans(across(1:6),na.rm=T) )
data_sec <- data_sec[rowSums(data_sec > 0) >= 4, ]%>% 
  mutate(SEC = rowMeans(across(1:6),na.rm=T) )


data2 <- log_data
data2$cate <- 'cate'
data2$RowNames <- rownames(data2)
data_af$RowNames <- rownames(data_af)
data_mf$RowNames <- rownames(data_mf)
data_sec$RowNames <- rownames(data_sec)

data2 <- data2[,19:20]
data2 <- left_join(data2,data_af[,7:8],by='RowNames')
data2 <- left_join(data2,data_mf[,7:8],by='RowNames')
data2 <- left_join(data2,data_sec[,7:8],by='RowNames')

data2[is.na(data2)] <- 0
data2 <- data2 %>% 
  mutate(exclusive = case_when(
    RowNames %in% deplist[['AF4_detected only']]  ~ "AF4-detected only",   
    RowNames %in% deplist[['CMDS_detected only']]  ~ "CMDS-detected only", 
    RowNames %in% deplist[['SEC_detected only']]  ~ "SEC-detected only", 
    RowNames %in% deplist[['AF4_elevated']]  ~ "AF4-elevated",   
    RowNames %in% deplist[['CMDS_elevated']]  ~ "CMDS-elevated", 
    RowNames %in% deplist[['SEC_elevated']]  ~ "SEC-elevated", 
    TRUE ~ NA                                
  ))

colnames(data2)
unique(data2$exclusive)
#data3 <- dplyr::filter(data2,exclusive %in% c('AF4-elevated','CMDS-elevated','SEC-elevated',NA ) )


p <- ggtern(data = data2, aes(AF4, `Exo-CMDS`, SEC)) + 
  geom_mask()+
  geom_point(aes(color = exclusive), size = 2, alpha = 0.7) +
  theme_bw() + 
  scale_color_manual(values = c("#9B62A7FF", "#5568B8FF",   
                                "#59A5A9FF", "#BEBC48FF", "#E78C35FF", "#721E17FF", 'grey')) +
  theme(axis.title.x = element_text(size = 7, face = 'bold'),
        axis.text.x = element_text(size = 7, angle = anglenum, vjust = vj, hjust = hj),
        axis.title.y = element_text(size = 7, face = 'bold'),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 7),          
        plot.subtitle = element_text(size = 7),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank( )
  ) +
  theme_arrowdefault()+themelegendtop

p
ggsave(paste0(wdplot,"j.pdf"), width = 12, height =9, units = "cm")




# S4 ev marker ----------------------

marker1 <- c("CD63", "CD9", "CD81", "CD82", "CD47", "GNA12", "GNA13", "GNA14", "GNA11", "GNA15", 
            "TSAP6", "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGA8", "ITGA9","ITGA10", "ITGA11", "ITGAE", "ITGAV", "ITGA2B", "ITGAD", "ITGAL", "ITGAX", "ITGAM", 
            "ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8", "TFR2", "LAMP1", "LAMP2", "SDC1", "SDC2", "SDC3", "SDC4", "BSG")

marker2 <- c(  "ADAM10", "GPC1", "NT5E", "CD59", 
              "TSG101", "CHMP1", "CHMP2", "CHMP3", "CHMP4", "CHMP5", "CHMP6", "CHMP7", "CHMP1A","CHMP1B", "CHMP2A", "CHMP2B", "CHMP4A", "CHMP4B", "CHMP4C", "PDCD6IP", "VPS4A", "VPS4B", 
              "ARRDC1", "FLOT1", "FLOT2")

marker3 <- c("CAV1", "CAV2", "CAV3", "SDCBP", "HSPA8", "HSP90AB1", 
            "ACTA1", "ACTA2", "ACTC1", "ACTB","ACTG1", "ACTG2", "TUBA1A", "TUBA1B", "TUBA1C",   "TUBA3C", "TUBA3D", "TUBA3E", "TUBA4A", 
            "TUBA8", "TUBB", "TUBB1", "TUBB2A", "TUBB2B", "TUBB3", "TUBB4A", "TUBB4B", "TUBB6", "TUBG1", "TUBG2", "GAPDH")


msmarker <- function(marker1)    {
  log_data[is.na(log_data)] <- 0
  data_af <- log_data[,1:6]
  data_mf <- log_data[,7:12]
  data_sec <- log_data[,13:18]
  
  data_af2 <- rownames_to_column(data_af)
  data_af2 <- filter(data_af2,rowname %in% marker1)%>% 
    column_to_rownames(var = "rowname") %>% 
    as.matrix() %>% 
    melt()
  colnames(data_af2) <- c('gene','id','value')
  
  data_mf2 <- rownames_to_column(data_mf)
  data_mf2 <- filter(data_mf2,rowname %in% marker1) %>% 
    column_to_rownames(var = "rowname") %>% 
    as.matrix() %>% 
    melt()
  colnames(data_mf2) <- c('gene','id','value')
  
  data_sec2 <- rownames_to_column(data_sec)
  data_sec2 <- filter(data_sec2,rowname %in% marker1)%>% 
    column_to_rownames(var = "rowname") %>% 
    as.matrix() %>% 
    melt()
  colnames(data_sec2) <- c('gene','id','value')
  
  
  bind_data <- rbind(data_af2,data_mf2)
  bind_data <- rbind(bind_data,data_sec2)
  bind_data <- left_join(bind_data,metadata,by=c('id'='new_name1'))
  
  median_value1 <- data_af2 %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(
      value = ifelse(sum(value == 0) >= 3, 0, median(value)),
      method1 = 'AF4'
    )
  
  median_value2 <- data_mf2 %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(
      value = ifelse(sum(value == 0) >= 3, 0, median(value)),
      method1 = 'MF'
    )
  
  median_value3 <- data_sec2 %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(
      value = ifelse(sum(value == 0) >= 3, 0, median(value)),
      method1 = 'SEC'
    )
  
  median_value <- rbind(median_value1,median_value2)
  median_value <- rbind(median_value,median_value3)
  
  bind_data$gene <-  factor(bind_data$gene, levels = marker1)
  median_value$gene <- factor(median_value$gene , levels = marker1)
  
  
  e <- ggplot(bind_data, aes(x=gene, y=value,fill=method1)) +
    geom_bar(data = median_value, aes(x=gene, y=value,fill=method1),
             stat = "identity", position = position_dodge2(width=1),width = 0.65,color='black',linewidth=0.3)+
    scale_fill_manual(values = c('#B5252E','#416997','#DAA83D','#ADC482'),labels=c('AF4','Exo-CMDS','SEC','PBA-control'),name='')+
    geom_errorbar(data=bind_data, aes(x=gene, y=value,fill=method1), position=position_dodge(0.65),
                  stat="boxplot", width=0.6, coef=0,  inherit.aes = F) +
    geom_jitter(data=bind_data, aes(x=gene, y=value,fill=method1), position=position_dodge(0.65), color = "black",size=1,shape=20) +
    scale_y_continuous(limits = c(0, 40), expand = c(0,0))+
    theme1_noback+ themelegendnone+ themetext30+
    labs(x='',y='log2(Intensity)',title = '')
  return(e)
}

d <- msmarker(marker1)
e <- msmarker(marker2)
f <- msmarker(marker3)

pb <- d+e+f+ plot_layout(heights = c(1,1,1))

ggsave(paste0(wdplot,"S4.pdf"), pb,width = 18, height =25, units = "cm")












