
load_packages <- function() {
  library(readxl)
  library(xlsx)
  #library(showtext)
  library(dplyr)
  library(tidyverse)
  # library('DESeq2')
  #library('stringr')
  library(plyr)
  #library(lubridate)
  library(zoo)
  #library(cowplot)
  library(writexl)
  # library(R.utils)
  library(reshape)
  #library(colourpicker)
  
}

# library(Seurat)
# library(ggplot2)
# library(clustree)#
# library(cowplot)
# library(dplyr)
# library(data.table)
# library(paletteer) #
# library(vcd)#
# library(gplots)
# library(limma)

load_plot_packages <- function() {
  
  library(ggpubr)
  library(ggplot2)
  library(ggsignif)
  library(ggrepel)
  library(ggplot2)
  library(patchwork)
  library(ggalt)
  library(umap)
  library(aplot)
  library(ggsci)
  
  
  
}
library(ggplot2)

namet <- function(file){ #
  
  a <- names(table(file))
  
  return(a)
}


dft <- function(file){ #
  
  a <- as.data.frame(table(file))
  
  return(a)
}


open_xlsx <- function(file,i){
  
  a <- as.data.frame(read_excel(file,sheet=i))
  
  return(a)
}
cal_cv <- function(x){
  y = na.omit(x)
  return(sd(y)/mean(y))
}

#theme new  -----------------------

xytitle <- 8
textsize <- 8
anglenum <- 0 # 45  90
vj <- 0 # 1.1  0 
hj <- 0.5 # 1.1  0.5
legendsize <- 0.4
legendposit <- '' #right left none
titlesize <- 9

theme4_backline <- theme(axis.title.x = element_text(size=xytitle,face='bold'),
                         axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj ),
                         axis.title.y = element_text(size=xytitle,face='bold'),
                         axis.text.y = element_text( size = textsize),
                         plot.title = element_text(size = titlesize),          
                         plot.subtitle = element_text(size = textsize),       
                         panel.background = element_rect(fill=NA, colour=NA),
                         legend.background = element_rect(fill=NA, colour=NA),
                         
                         panel.grid.major = element_line(color = "#EDEDED", size = 0.5), 
                         panel.grid.minor = element_line(color = "#EDEDED", size = 0.25),
                         panel.border = element_rect(colour = "black", fill = NA, size = 0.2),
                         #panel.background = element_rect(fill = NA, color = "black"),
                         strip.background = element_blank(),
                         axis.line = element_blank( ),
                         
                         axis.ticks = element_line(colour = "black", size = 0.25),
                         axis.ticks.length.x=unit(0.08, "cm"),
                         axis.ticks.length.y=unit(0.08, "cm")
)


theme2_noback <- theme(axis.title.x = element_text(size=xytitle,face='bold'),
                       axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj ),
                       axis.title.y = element_text(size=xytitle,face='bold'),
                       axis.text.y = element_text( size = textsize),
                       plot.title = element_text(size = titlesize),          
                       plot.subtitle = element_text(size = textsize),       
                       panel.grid.major = element_blank( ),
                       panel.grid.minor = element_blank( ),
                       panel.border = element_blank( ),
                       panel.background = element_rect(fill=NA, colour=NA),
                       legend.background = element_rect(fill=NA, colour=NA),
                       
                       axis.line = element_line(colour = "black", size = 0.15, 
                                                arrow = arrow(type = "closed", #open
                                                              angle = 12,
                                                              length = unit(0.2, "cm"))
                       ), #箭头                  
                       axis.ticks = element_line(colour = "black", size = 0.25),
                       axis.ticks.length.x=unit(0.08, "cm"),
                       axis.ticks.length.y=unit(0.08, "cm"),
                       
)



theme2_noback <- theme(axis.title.x = element_text(size=xytitle,face='bold'),
                       axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj ),
                       axis.title.y = element_text(size=xytitle,face='bold'),
                       axis.text.y = element_text( size = textsize),
                       plot.title = element_text(size = titlesize),          
                       plot.subtitle = element_text(size = textsize),       
                       panel.grid.major = element_blank( ),
                       panel.grid.minor = element_blank( ),
                       panel.border = element_blank( ),
                       panel.background = element_rect(fill=NA, colour=NA),
                       legend.background = element_rect(fill=NA, colour=NA),
                       
                       axis.line = element_line(colour = "black", size = 0.15, 
                                                arrow = arrow(type = "closed", #open
                                                              angle = 12,
                                                              length = unit(0.2, "cm"))
                       ), #箭头                  
                       axis.ticks = element_line(colour = "black", size = 0.25),
                       axis.ticks.length.x=unit(0.08, "cm"),
                       axis.ticks.length.y=unit(0.08, "cm"),
                       
)


theme2_nobackplus <- theme(
  axis.title.x = element_text(size = xytitle, face = 'bold'),
  axis.text.x = element_text(size = textsize, angle = anglenum, vjust = vj, hjust = hj),
  axis.title.y = element_text(size = xytitle, face = 'bold'),
  axis.text.y = element_text(size = textsize),
  plot.title = element_text(size = titlesize),
  plot.subtitle = element_text(size = textsize),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_rect(fill = NA, colour = NA),
  
  legend.background = element_rect(fill = NA, colour = NA),
  
  axis.line.x = element_line(
    colour = "darkgrey",
    size = 0.15,
    linetype = "dashed"
  ),
  axis.line.y = element_line(
    colour = "black",
    size = 0.15
  ),

  axis.ticks = element_line(colour = "black", size = 0.25),
  axis.ticks.length.x = unit(0.08, "cm"),
  axis.ticks.length.y = unit(0.08, "cm")
)


theme1_noback <- theme(axis.title.x = element_text(size=xytitle,face='bold'),
                       axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj ),
                       axis.title.y = element_text(size=xytitle,face='bold'),
                       axis.text.y = element_text( size = textsize),
                       plot.title = element_text(size = titlesize),          
                       plot.subtitle = element_text(size = textsize),       
                       panel.grid.major = element_blank( ),
                       panel.grid.minor = element_blank( ),
                       panel.border = element_blank( ),
                       
                       panel.background = element_rect(fill=NA, colour=NA),
                       legend.background = element_rect(fill=NA, colour=NA),
                       
                       axis.line.x = element_blank(),  
                       axis.line.y = element_line(colour="black", size=0.25),
                       axis.ticks = element_line(colour = "black", size = 0.25),
                       axis.ticks.length.x=unit(0.08, "cm"),
                       axis.ticks.length.y=unit(0.08, "cm")
)


theme1_backline <- theme(axis.title.x = element_text(size=xytitle,face='bold'),
                         axis.text.x = element_text( size = textsize, angle = anglenum,vjust = vj,hjust= hj ),
                         axis.title.y = element_text(size=xytitle,face='bold'),
                         axis.text.y = element_text( size = textsize),
                         plot.title = element_text(size = titlesize),          
                         plot.subtitle = element_text(size = textsize),       
                         # panel.grid.major = element_blank( ),
                         # panel.grid.minor = element_blank( ),
                         panel.border = element_blank( ),
                         
                         panel.background = element_rect(fill=NA, colour=NA),
                         legend.background = element_rect(fill=NA, colour=NA),
                         
                         panel.grid.major = element_line(color = "#EDEDED", size = 0.5), 
                         panel.grid.minor = element_line(color = "#EDEDED", size = 0.25),
                         # panel.border = element_rect(colour = "black", fill = NA, size = 0.2),
                         
                         axis.line.x = element_blank(),  
                         axis.line.y = element_line(colour="black", size=0.25),
                         
                         axis.ticks = element_line(colour = "black", size = 0.25),
                         axis.ticks.length.x=unit(0.08, "cm"),
                         axis.ticks.length.y=unit(0.08, "cm"),
                         
                         legend.title = element_text(size = 7), 
                         legend.text = element_text(size = 6)
)


themelegendright <- theme(legend.position = 'right', 
                          legend.title = element_text(size = textsize),
                          legend.key.size = unit(legendsize, "cm"),  
                          legend.text = element_text(size = textsize),
                          panel.background = element_rect(fill=NA, colour=NA))
themelegendleft <- theme(legend.position = 'left', 
                         legend.title = element_text(size = textsize),
                         legend.key.size = unit(legendsize, "cm"),  
                         legend.text = element_text(size = textsize,angle=90),
                         panel.background = element_rect(fill=NA, colour=NA))
themelegendtop <- theme(legend.position = 'top', 
                        legend.title = element_text(size = textsize),
                        legend.key.size = unit(legendsize, "cm"),  
                        legend.text = element_text(size = textsize),
                        panel.background = element_rect(fill=NA, colour=NA))
themelegendnone <- theme(legend.position = 'none', 
                         legend.title = element_text(size = textsize),
                         legend.key.size = unit(legendsize, "cm"),  
                         legend.text = element_text(size = textsize),
                         panel.background = element_rect(fill=NA, colour=NA))

themetext30 <- theme(axis.text.x = element_text( size = textsize, angle = 30,
                                                 vjust = 0.7,hjust= 0.7 )
)

themetext90 <- theme(axis.text.x = element_text( size = textsize, angle = 90,hjust = 1, vjust = 0.5 )
)

themey90 <- theme( axis.text.y = element_text( size = 7, angle = 90,vjust = 0,hjust= 0.5))



#-------------
mutl_anova_test <- function(df_ancova,dependent,independent){      
  
  fit <- aov(as.formula(paste(dependent, "~", paste(independent, collapse = "*"))),
             data=df_ancova)
  
  result2 <- TukeyHSD(fit)
  result2 <- do.call(rbind, TukeyHSD(fit))
  result2 <- rownames_to_column(as.data.frame(result2))
  
  result <- data.frame(do.call('rbind', strsplit(as.character(result2$rowname), "[:-]")))
  result <- cbind(result,result2)
  
  #result <- subset(result, X2 == X4)
  result$mark <- ifelse(result$`p adj` < 0.05, "*", NA)
  
  return(result)
  
}

mutl_anova_testplus <- function(inputdata,dependent,independent){      
  columns <- c('source','protein', 'Sum','Mean','F','pvalue')
  scores <- data.frame(matrix(ncol = length(columns), nrow = 0))
  result2 <- list()
  for (i in unique(inputdata$protein)) {
    df_ancova <- inputdata[inputdata$protein == i, ]
    fit <- aov(as.formula(paste(dependent, "~", paste(independent, collapse = "*"))),
               data=df_ancova)
    
    result <- summary(fit)
    result <- as.data.frame(result[[1]])
    result <- rownames_to_column(result)
    result$Df <- i
    colnames(result) <- columns
    result <- as.data.frame(apply(result, 2, function(x) gsub(" ", "", x)))
    result <- filter(result,source != 'Residuals')
    scores <- rbind(scores,result)
    # TukeyHSD(fit)
    result2[[i]] <- do.call(rbind, TukeyHSD(fit))
    result2[[i]] <-as.data.frame(result2[[i]] )
    result2[[i]]$protein <- i
    result2[[i]]$cate <- rownames(result2[[i]])
    
    
  }
  scores$pvalue <- as.numeric(scores$pvalue)
  
  set.seed(123)
  adj <- p.adjust(scores$`pvalue`, method = 'fdr', n = length(scores$`pvalue`))
  scores$qvalue <- adj
  # scores$rejected <- adj < fdrancova_result
  #return(result2)
  return(scores)
  
}



mutl_anova_testplus2 <- function(inputdata,dependent,independent){      
  columns <- c('source','protein', 'Sum','Mean','F','pvalue')
  scores <- data.frame(matrix(ncol = length(columns), nrow = 0))
  result2 <- list()
  for (i in unique(inputdata$protein)) {
    df_ancova <- inputdata[inputdata$protein == i, ]
    unique_levels <- unique(df_ancova$method1)
    
    if (length(unique_levels) < 2) {
      print(i)
    }else{
      
      fit <- aov(as.formula(paste(dependent, "~", paste(independent, collapse = "*"))),
                 data=df_ancova)
      
      result <- summary(fit)
      result <- as.data.frame(result[[1]])
      result <- rownames_to_column(result)
      result$Df <- i
      colnames(result) <- columns
      result <- as.data.frame(apply(result, 2, function(x) gsub(" ", "", x)))
      result <- filter(result,source != 'Residuals')
      scores <- rbind(scores,result)
      # TukeyHSD(fit)#
      result2[[i]] <- do.call(rbind, TukeyHSD(fit))
      result2[[i]] <-as.data.frame(result2[[i]] )
      result2[[i]]$protein <- i
      result2[[i]]$cate <- rownames(result2[[i]])
    }
    
  }
  scores$pvalue <- as.numeric(scores$pvalue)
  
  set.seed(123)
  adj <- p.adjust(scores$`pvalue`, method = 'fdr', n = length(scores$`pvalue`))
  scores$qvalue <- adj
  # scores$rejected <- adj < fdrancova_result
  return(result2)
  #return(scores)
  
}


# end--------------


## genelist overlap-------------
two_gene_list_overlap_pvalue<-function(gene_list_1,gene_list_2,bg_num){
  over_gene<-intersect(gene_list_1,gene_list_2)
  p_hyper=1-phyper(length(over_gene)-1,length(gene_list_1),bg_num-length(gene_list_1),length(gene_list_2))
  return(p_hyper)
}
jaccard_index <- function(x, y) {
  intersection <- length(intersect(x, y))
  union <- length(base::union(x, y))
  if (union == 0) {
    return(NA)
  } else {
    return(intersection / union)
  }
}


list_overlap <- function(af4_cluster,mf_cluster,bg_num){
  
  af_mf <- matrix(nrow = length(af4_cluster), 
                  ncol = length(mf_cluster), 
                  dimnames = list(names(af4_cluster),
                                  names(mf_cluster)))
  for (i in 1:length(af4_cluster)) {
    for (j in 1:length(mf_cluster)) {
      af_mf[i, j] <- two_gene_list_overlap_pvalue(af4_cluster[[i]], mf_cluster[[j]], bg_num)  
    }
  }
  
  return(af_mf)
  
}

list_overlap_jaccard <- function(af4_cluster,mf_cluster,bg_num){
  
  af_mf <- matrix(nrow = length(af4_cluster), 
                  ncol = length(mf_cluster), 
                  dimnames = list(names(af4_cluster),
                                  names(mf_cluster)))
  for (i in 1:length(af4_cluster)) {
    for (j in 1:length(mf_cluster)) {
      af_mf[i, j] <- jaccard_index(af4_cluster[[i]], mf_cluster[[j]])  
    }
  }
  return(af_mf)
}



two_list_overlap_num <- function(wgcna_genes,marker_list){
  
  result_matrix <- matrix(0, nrow = length(marker_list), ncol = length(wgcna_genes))
  rownames(result_matrix) <- names(marker_list) 
  colnames(result_matrix) <- names(wgcna_genes)
  
  for (i in names(marker_list)) {
    for (j in names(wgcna_genes)) {
      result_matrix[i, j] <- length(intersect(marker_list[[i]], wgcna_genes[[j]]))
    }
  }
  
  return(result_matrix)
} 


two_list_overlap_gene <- function(wgcna_genes,marker_list){
  
  result_matrix <- matrix(0, nrow = length(marker_list), ncol = length(wgcna_genes))
  rownames(result_matrix) <- names(marker_list) 
  colnames(result_matrix) <- names(wgcna_genes)
  
  for (i in names(marker_list)) {
    for (j in names(wgcna_genes)) {
      if(length(intersect(marker_list[[i]], wgcna_genes[[j]]))>0){
        result_matrix[i, j] <- paste(intersect(marker_list[[i]], wgcna_genes[[j]]), collapse = ", ")
        
      }else{
        result_matrix[i, j] <- ''
        
      }
      
      
    }
  }
  
  return(result_matrix)
} 


