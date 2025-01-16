#plot nice map and table for kdr sample counts
#Joel Odero 05.06.24
#######
install.packages('reshape2')
install.packages("githubinstall")
library("githubinstall")
library('reshape2')
githubinstall("starmie")
library(tidyverse)
library(sf)
library(RColorBrewer)
library("ggrepel")
library("ggspatial")
library("ggsci")
library("data.table")
library(starmie)
library(conStruct)
library('scales')


#Admixture 

setwd("/Users/2576093o/OneDrive - University of Glasgow/Tz_funestus_analysis")
df_samples <- read.csv('df_samples.csv')
df_samples_filtered <- df_samples[!df_samples$admin1_name == "Kigoma", ]
dir = '/Users/2576093o/OneDrive - University of Glasgow/Tz_funestus_analysis/4857.TZ//'
fhand = '2RL.50000.0.05.0.0.TZ.'
admixoutput = list()
i = 0
for (a in seq(1,9)) {
  i<-i+1
  pf = paste0(dir, fhand, a, '.P')
  qf = paste0(dir, fhand, a, '.Q')
  lf = paste0(dir,'log',a,'.out')
  my_admix <- loadAdmixture(qf, pf, lf)
  admixoutput[[i]] <-  my_admix
}
#turn into admixture object multi
#plot cv diagnostics
admixlist_all = starmie::admixList(admixoutput)
bestK(admixlist_all)

#plot all K for the first six or so

starmie::plotMultiK(admixlist_all[1:9], populations = df_samples_filtered[,c('V1', 'admin1_name')], plot = TRUE)

#plot all K for K2-6
starmie::plotMultiK(admixlist_all[2:6], populations = df_samples[,c('V1','admin1_name')], plot = TRUE)

#plot all K for K2,K3 and K4
starmie::plotMultiK(admixlist_all[2:4], populations = df_samples[,c('V1','admin1_name')], plot = TRUE)

#plot admix only for good ks
bar_k2 = starmie::plotBar(admixoutput[[2]],plot = FALSE)
bar_k3 = starmie::plotBar(admixoutput[[3]],plot = FALSE)
bar_k4 = starmie::plotBar(admixoutput[[4]],plot = FALSE)

#bind df
bar_k2 = cbind(bar_k2,df_samples[,c('partner_sample_id','cohort_admin1_year', 'admin1_name')])
bar_k3 = cbind(bar_k3,df_samples[,c('partner_sample_id','cohort_admin1_year', 'admin1_name')])
bar_k4 = cbind(bar_k4,df_samples[,c('partner_sample_id','cohort_admin1_year', 'admin1_name')])

#reorder factors for plotting
bar_k2$admin1_name <- factor(bar_k2$admin1_name, levels = c("Pwani", "Tanga", "Morogoro", "Ruvuma", "Lindi", "Mtwara", "Dodoma","Kagera", 'Katavi',"Mwanza","Kigoma"))
bar_k3$admin1_name <- factor(bar_k3$admin1_name, levels = c("Pwani", "Tanga", "Morogoro", "Ruvuma", "Lindi", "Mtwara", "Dodoma","Kagera", 'Katavi',"Mwanza","Kigoma"))
bar_k4$admin1_name <- factor(bar_k4$admin1_name, levels = c("Pwani", "Tanga", "Morogoro", "Ruvuma", "Lindi", "Mtwara", "Dodoma","Kagera", 'Katavi',"Mwanza","Kigoma"))

#plotting
k2plot<-
  ggplot(bar_k2, aes(factor(partner_sample_id), value, fill = factor(Cluster))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~admin1_name, switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  labs(title = "K=2", y = "Ancestry", x = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.5)) +
  scale_fill_brewer(palette = "Spectral") +
  theme(
    panel.spacing.x = unit(0.05, "lines"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 18), 
    legend.position = "none",
    strip.text.x = element_text(angle = 90, face = 'bold', size ='22'),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 12, r = 80, b = 12, l = 12),  # Increased right margin
    axis.title.x = element_text(face = "bold", size = '25'),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = '25'),
    plot.title = element_text(face = "bold", size = 20)# Bold plot title
  )
k2plot 
p1<-k2plot+theme(strip.text.x = element_blank())
p1
ggsave("plot_k2.tiff", plot = k2plot, width = 20, height = 5, dpi = 300, bg = 'white')

k3plot<-
  ggplot(bar_k3, aes(factor(partner_sample_id), value, fill = factor(Cluster))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~admin1_name, switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  labs(title = "K=3", y = "Ancestry", x = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.5)) +
  scale_fill_brewer(palette = "Spectral") +
  theme(
    panel.spacing.x = unit(0.05, "lines"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 18), 
    legend.position = "none",
    strip.text.x = element_text(angle = 90, face = 'bold', size ='22'),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 12, r = 80, b = 12, l = 12),  # Increased right margin
    axis.title.x = element_text(face = "bold", size = '25'),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = '25'),# Bold y-axis title
    plot.title = element_text(face = "bold", size = 20)# Bold plot title
  )
k3plot 
p2<- k3plot+theme(strip.text.x = element_blank())
p2
ggsave("plot_k3.tiff", plot = k3plot, width = 20, height = 5, dpi = 300, bg = 'white')
ggsave("plot_k3.png", plot = k3plot, width = 20, height = 5, dpi = 300, bg = 'white')

k4plot<-
  ggplot(bar_k4, aes(factor(partner_sample_id), value, fill = factor(Cluster))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~admin1_name, switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  labs(title = "K=4", y = "Ancestry", x = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1.5)) +
  scale_fill_brewer(palette = "Spectral") +
  theme(
    panel.spacing.x = unit(0.05, "lines"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 18), 
    legend.position = "none",
    strip.text.x = element_text(angle = 90, face = 'bold', size ='22'),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 12, r = 80, b = 12, l = 12),  # Increased right margin
    axis.title.x = element_text(face = "bold", size = '25'),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = '25'),   # Bold y-axis title
    plot.title = element_text(face = "bold", size = 20)# Bold plot title
  )
k4plot
ggsave("plot_k4.tiff", plot = k4plot, width = 20, height = 5, dpi = 300, bg = 'white')

# Combine plots
install.packages('patchwork')
library(patchwork)

# Combine plots
combined_plot <- p1 / p2 / k4plot + plot_layout(ncol = 1)
combined_plot

# Save combined plot
ggsave("ADMIXTURE_combined_plot.png", plot = combined_plot, width = 20, height = 15, dpi = 300, bg = 'white')
ggsave("ADMIXTURE_combined_plot.tiff", plot = combined_plot, width = 20, height = 15, dpi = 300, bg = 'white')

