library(tidyverse)
library(gridExtra)
library('ggplot2')

setwd("/Users/2576093o/OneDrive - University of Glasgow/Tz_funestus_analysis")
df_samples <- read.csv('df_samples.csv')

#load roh data for each chrom
df_3rl <- read.csv('ROH_3RL.csv')
df_2rl <- read.csv('ROH_2RL.csv')
df_x <- read.csv('ROH_X.csv')

#roh for chrom 3rl

summary_df_3rl <- df_3rl %>% group_by(sample_id) %>% 
  summarise(mean_seglen = mean(roh_length), 
            sum_seglen = sum(roh_length), 
            f_roh = sum(roh_length)/210990000, 
            var_rohlen = var(roh_length), 
            n_segs = length(roh_length))

summary_df_3rl <- left_join(summary_df_3rl, df_samples)
view(summary_df_3rl)

summary_df_3rl$admin1_name <- factor(summary_df_3rl$admin1_name, 
                                     levels = c("Kagera", "Mwanza", "Katavi",
                                                'Kigoma','Dodoma','Morogoro', 
                                                'Pwani', 'Ruvuma', 'Tanga', 'Lindi', 
                                                'Mtwara'))

#do this for froh and nroh
ggplot(summary_df_3rl, aes(x=admin1_name, y=n_segs, colour=admin1_name))+
  xlab(bquote('3RL')) +
  geom_boxplot()

#joint distribution chrom 3rl
plot_3rl <- ggplot(summary_df_3rl, aes(x = f_roh, y = n_segs, colour = admin1_name)) +
  ylab(bquote("Count ROH")) +
  xlab("Frequency ROH") +
  ggtitle("3RL") +
  geom_point() +
  theme_classic()
plot_3rl
#**********************************************************************************************

#roh for chrom 2rl
summary_df_2rl <- df_2rl %>% group_by(sample_id) %>% 
  summarise(mean_seglen = mean(roh_length), 
            sum_seglen = sum(roh_length), 
            f_roh = sum(roh_length)/210990000, 
            var_rohlen = var(roh_length), 
            n_segs = length(roh_length))

summary_df_2rl <- left_join(summary_df_2rl, df_samples)
view(summary_df_2rl)

summary_df_2rl$admin1_name <- factor(summary_df_2rl$admin1_name, 
                                     levels = c("Kagera", "Mwanza", 
                                                "Katavi",'Kigoma','Dodoma',
                                                'Morogoro', 'Pwani', 'Ruvuma', 
                                                'Tanga', 'Lindi', 'Mtwara'))

#do this for froh and nroh
ggplot(summary_df_2rl, aes(x=admin1_name, y=n_segs, colour=admin1_name))+
  xlab(bquote('2RL')) +
  geom_boxplot()

#joint distribution chrom 2rl
plot_2rl <- ggplot(summary_df_2rl, aes(x = f_roh, y = n_segs, colour = admin1_name)) +
  ylab(bquote("Count ROH")) +
  xlab("Frequency ROH") +
  ggtitle("2RL") +
  geom_point() +
  theme_classic()
plot_2rl
#***************************************************************************************************
#roh for chrom x
summary_df_x <- df_x %>% group_by(sample_id) %>% 
  summarise(mean_seglen = mean(roh_length), 
            sum_seglen = sum(roh_length), 
            f_roh = sum(roh_length)/210990000, 
            var_rohlen = var(roh_length), 
            n_segs = length(roh_length))

summary_df_x <- left_join(summary_df_x, df_samples)
view(summary_df_x)

summary_df_x$admin1_name <- factor(summary_df_x$admin1_name, 
                                     levels = c("Kagera", "Mwanza", 
                                                "Katavi",'Kigoma','Dodoma',
                                                'Morogoro', 'Pwani', 'Ruvuma', 
                                                'Tanga', 'Lindi', 'Mtwara'))
#do this for froh and nroh
ggplot(summary_df_x, aes(x=admin1_name, y=n_segs, colour=admin1_name))+
  xlab(bquote('X')) +
  geom_boxplot()

#joint distribution chrom x
plot_x <- ggplot(summary_df_x, aes(x = f_roh, y = n_segs, colour = admin1_name)) +
  ylab(bquote("Count ROH")) +
  xlab("Frequency ROH") +
  ggtitle("X") +
  geom_point() +
  theme_classic()
plot_x

#***********************************************************************************************

#combine roh plots for 3 individual chroms
combined_plot <- grid.arrange(plot_2rl, plot_3rl, plot_x, nrow = 1)


#combine roh for all chroms
ROH_comb <- bind_rows(df_3rl, df_2rl, df_x)
view(ROH_comb)

summary_df_comb <- ROH_comb %>% group_by(sample_id) %>% 
  summarise(mean_seglen = mean(roh_length), 
            sum_seglen = sum(roh_length), 
            f_roh = sum(roh_length)/210990000, 
            var_rohlen = var(roh_length), 
            n_segs = length(roh_length))


summary_df_comb <- left_join(summary_df_comb, df_samples)
view(summary_df_comb)

summary_df_comb$admin1_name <- factor(summary_df_comb$admin1_name, 
                                   levels = c("Kagera", "Mwanza", 
                                              "Katavi",'Kigoma','Dodoma',
                                              'Morogoro', 'Pwani', 'Ruvuma', 
                                              'Tanga', 'Lindi', 'Mtwara'))
install.packages("RColorBrewer")
library(RColorBrewer)
set3_palette <- brewer.pal(n = 11, name = "Set3")


region_colors <- c(
  'Pwani' = set3_palette[1],
  'Morogoro' = set3_palette[2],
  'Kigoma' = set3_palette[3],
  'Tanga' = set3_palette[4],
  'Ruvuma' = set3_palette[5],
  'Kagera' = set3_palette[6],
  'Mtwara' = set3_palette[7],
  'Katavi' = set3_palette[8],
  'Dodoma' = set3_palette[9],
  'Lindi' = set3_palette[10],
  'Mwanza' = set3_palette[11]
)

#do this for froh and nroh
comb <- ggplot(summary_df_comb, aes(x = admin1_name, y = n_segs, fill = admin1_name)) +
  geom_boxplot() +
  xlab('') +
  ylab("Number of Segments") +  # Add y-axis label for clarity
  scale_fill_manual(values = region_colors) +  # Apply the specified colors
  theme_classic() +
  ggtitle("") +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 18), 
    axis.title.x = element_text(face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = '25'),  # Bold y-axis title
    axis.text.x = element_text(face = "bold", size = '22', angle = 90, hjust = 1),  # Bold x-axis text and angled
    legend.position = "none",  # Remove legend if not needed
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = 'white', color = NA),  # Set background to white for contrast
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  # Add black border around the plot
  )
# Display the plot
print(comb)

ggsave("ROH_boxplot_combined.tiff", plot = comb, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("ROH_boxplot_combined.png", plot = comb, width = 10, height = 8, dpi = 300, bg = "white")

#joint distribution

plot_roh_comb <- ggplot(summary_df_comb, aes(x = f_roh, y = n_segs, colour = admin1_name)) +
  geom_point(size = 6) +  # Increase point size
  ylab(bquote("Count ROH")) +
  xlab("Fraction of the genome in ROH") +
  ggtitle("") +
  scale_colour_manual(values = region_colors) +  # Apply the specified colors
  labs(colour = "Regions") +  # Add legend title
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 25),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 25),  # Bold y-axis title
    legend.title = element_text(face = "bold", size = 25),  # Bold legend title and increase font size
    legend.text = element_text(face = "bold", size = 25),  # Bold legend text and increase font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = 'white', color = NA),  # Set background to white for contrast
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  # Add black border around the plot
  )

# Display the plot
print(plot_roh_comb)


ggsave("ROH_combined.tiff", plot = plot_roh_comb, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("ROH_combined.png", plot = plot_roh_comb, width = 10, height = 8, dpi = 300, bg = "white")

#plot this as a 3 panel figure, order by geographic location and make colours consistent with the TZ funestus analysis notebook


