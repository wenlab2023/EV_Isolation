wd <- "~/Desktop/EV_project/"
setwd(wd)

wdcode <- paste0("code/") 
ver <- 'F2/'
wdresult <- paste0("result/",ver)
wddata <- paste0("data/",ver)
wdplot <- paste0("plot/",ver)
source(paste0(wdcode,"function.R"))
load_packages()
load_plot_packages() 


# b -------------------

load(paste0(wddata,'F2bdata.Rdata'))

data1$storage <- factor(data1$storage, levels = c('none','froz'))

data1 %>%
  ggplot(aes(x = method, y = value, fill = storage)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.5, 
               size = 0.5, colour = 'black', 
               outlier.colour = 'white', outlier.size = 0.1,alpha=0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),alpha=0.8,
             shape = 21, colour = 'black', size = 1.3) + 
  scale_fill_manual(values = c("#B2182B","#264653","#CDAD00")) +  
  labs(y = 'Total particles/ml', x = "", fill = '') +  
  scale_y_continuous(limits = c(4e+10, 2.8e+11)) +
  theme2_noback +
  geom_signif(
    y_position = c(2.4e+11,2.55e+11), xmin = c(1,2), xmax = c(2.95,2.95),
    annotation = c("*", "*"), tip_length = 0,size = 0.2,textsize=3
  ) +themelegendnone

ggsave(paste0(wdplot,'b',"total particles.pdf"), width = 13, height =6, units = "cm")




# c   ------------------------

mean_plot <-  as.data.frame(read_excel(paste0(wddata,"NTA_inputdata.xlsx"),sheet=1))

i <- 'AF4'
mean_plot_cate <- dplyr::filter(mean_plot,method==i)
mean_plot_cate <- mean_plot_cate %>% arrange(X1) %>%
  mutate(cum_cells = cumsum(as.numeric(value)), total_cells = sum(as.numeric(value)), cum_freq = cum_cells / total_cells)
lower_bound <- min(mean_plot_cate$X1[mean_plot_cate$cum_freq > 0.025])
upper_bound <- max(mean_plot_cate$X1[mean_plot_cate$cum_freq < 0.975])

ggplot(mean_plot_cate, aes(x = as.numeric(X1), y = as.numeric(value), color = as.character(name))) +
  geom_line(size = 0.5, alpha = 0.7) +
  geom_vline(xintercept = lower_bound, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = upper_bound, color = "darkgrey", linetype = "dashed") +
  scale_color_manual(values = c("#B2182B", "#B2182B", "#264653", "#264653", "#264653", "#264653"), name = "") +
  labs(title = paste0(i), x = "Size(nm)", y = "Particles/ml")+
  scale_x_continuous(limits=c(0,400),
                     breaks=seq(0,600,100))+
  scale_y_continuous(labels = function(x) sprintf("%.1e", x))+
  theme2_noback+themelegendnone
ggsave(paste0(wdplot,'c',i,".pdf"), width = 7, height = 5.5, units = "cm")


i <- 'MF'
mean_plot_cate <- dplyr::filter(mean_plot,method==i)
mean_plot_cate <- mean_plot_cate %>% arrange(X1) %>%
  mutate(cum_cells = cumsum(as.numeric(value)), total_cells = sum(as.numeric(value)), cum_freq = cum_cells / total_cells)
lower_bound <- min(mean_plot_cate$X1[mean_plot_cate$cum_freq > 0.025])
upper_bound <- max(mean_plot_cate$X1[mean_plot_cate$cum_freq < 0.975])

ggplot(mean_plot_cate,aes(x=as.numeric(X1),y=(as.numeric(value)),color=as.character(name)))+
  geom_line(size=0.5,alpha=0.7)+
  geom_vline(xintercept = lower_bound, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = upper_bound, color = "darkgrey", linetype = "dashed") +
  scale_color_manual(values =c("#B2182B", "#B2182B","#264653","#264653","#264653","#264653"),name="")+
  #theme_bw()+
  labs(title = 'EXO-CMDS',x="Size(nm)",y="Particles/ml")+
  scale_x_continuous(limits=c(0,400), 
                     breaks=seq(0,600,100))+ 
  scale_y_continuous(labels = function(x) sprintf("%.1e", x))+
  # scale_y_continuous(limits=c(0,6e+09), 
  #                    breaks=seq(0,6e+09,1e+09))+
  # geom_vline(xintercept = 30, linetype = "dashed",size=0.3) +
  # geom_vline(xintercept = 150, linetype = "dashed",size=0.3)+
  theme2_noback+themelegendnone
ggsave(paste0(wdplot,'c',i,".pdf"), width = 7, height = 5.5, units = "cm")



i <- 'SEC'
mean_plot_cate <- dplyr::filter(mean_plot,method==i)    
mean_plot_cate <- mean_plot_cate %>% arrange(X1) %>%
  mutate(cum_cells = cumsum(as.numeric(value)), total_cells = sum(as.numeric(value)), cum_freq = cum_cells / total_cells)
lower_bound <- min(mean_plot_cate$X1[mean_plot_cate$cum_freq > 0.025])
upper_bound <- max(mean_plot_cate$X1[mean_plot_cate$cum_freq < 0.975])

ggplot(mean_plot_cate,aes(x=as.numeric(X1),y=(as.numeric(value)),color=as.character(name)))+
  geom_line(size=0.5,alpha=0.7)+
  geom_vline(xintercept = lower_bound, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = upper_bound, color = "darkgrey", linetype = "dashed") +
  scale_color_manual(values =c("#B2182B", "#B2182B","#264653","#264653","#264653","#264653"),name="")+
  #theme_bw()+
  labs(title = 'SEC',x="Size(nm)",y="Particles/ml")+
  scale_x_continuous(limits=c(0,400), 
                     breaks=seq(0,600,100))+ 
  scale_y_continuous(labels = function(x) sprintf("%.1e", x))+
  # scale_y_continuous(limits=c(0,6e+09),
  #                    breaks=seq(0,6e+09,1e+09))+ 
  # geom_vline(xintercept = 30, linetype = "dashed",size=0.3) +
  # geom_vline(xintercept = 150, linetype = "dashed",size=0.3)+
  theme2_noback+themelegendnone
ggsave(paste0(wdplot,'c',i,".pdf"), width = 7, height = 5.5, units = "cm")

ggsave(paste0(wdplot,'1_NTA/',i,"_concentration.pdf"), width = 7, height = 5.5, units = "cm")



