setwd("C:/Users/weixuan/Desktop/Pixy_LD/Heterzygoisity_N163")
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(reshape2)
library(stringr)
library(dplyr)
library(grid)   # for the textGrob() function
library(ggpubr)
library(ggplot2)
library(scales)
library(reshape2)
library(pheatmap)
library("cowplot")
######################################################################
safe_shape_palette <- c("K. drynarioides" = 19, 
                        "K. cookei" = 15, 
                        "K. kauaiensis" = 17)
species.color <- RColorBrewer::brewer.pal(name="Set1", 3) %>%
  setNames(., c("K. cookei", "K. drynarioides", "K. kauaiensis"))

dot.color <- RColorBrewer::brewer.pal(name="Set2", 4)
######################################################################

heteperce <- read.table(list.files(pattern = "\\.het$")[1], header = TRUE) 

heteperce2 <- heteperce %>%
  filter(!str_detect(INDV, "HV_1a|HV_1b|HV_1B")) %>%
  mutate(group = str_extract(INDV, "^[^_]+")) %>% 
  mutate(group = if_else(grepl("Kd_HV", INDV), "Kc", group)) %>%
  mutate(SampleGroup = if_else(grepl("Kd_HV", INDV), "Kc_HV", group)) %>%
  mutate(group = recode(group,
                        "Kc" = "K. cookei",
                        "Kd" = "K. drynarioides",
                        "Kk" = "K. kauaiensis"))  %>%
  mutate(percentage = (N_SITES-O.HOM.)/N_SITES) 


heteperce_mean <- heteperce2 %>%
  add_count(group, name = "freq") %>% 
  group_by(group,freq) %>%
  dplyr::summarize(percent_mean = mean(percentage, na.rm=TRUE)) 

group_ordered <- with(heteperce_mean,  reorder(group,  percent_mean,  mean))
heteperce2$group <- factor(heteperce2$group, levels = levels(group_ordered))

heteperce_Plot <- ggplot(heteperce2, aes(x=group, y=percentage, fill = group)) + 
  geom_boxplot(show.legend = F, width = 0.7, outlier.shape = NA) + 
  geom_text(data = heteperce_mean,
            aes(x = group, y = 0.01,
            label = paste0(formatC(percent_mean * 100, format = "f", digits = 2), "%")),
            size = 4 ) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue", fill="blue") +
  geom_jitter(aes(colour = SampleGroup), size=3, alpha=0.9) +
  scale_fill_manual(values=species.color)+
  scale_color_manual(values = dot.color) +
  ylab("Proportion of heterozygosity sites") +
  xlab(NULL) +
  theme_bw() +
  theme(#legend.position = "none",
        axis.text.x = element_text(size = 10, face = "italic" ),
        panel.border = element_rect(colour = "black", fill=NA))

heteperce_Plot <- heteperce_Plot + guides(fill="none")

############################################################################

FIS_mean <- heteperce2 %>%
  add_count(group, name = "freq") %>% 
  group_by(group,freq) %>%
  dplyr::summarize(F_mean = mean(F, na.rm=TRUE)) 

Fis_Plot <- ggplot(heteperce2, aes(x=group, y= F, fill = group)) + 
  geom_boxplot(show.legend = F, width = 0.7, outlier.shape = NA) + 
  geom_text(data = FIS_mean,
            aes(x = group, y = 0.4,
                label = formatC(F_mean, digits = 2)),
            size = 4 ) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue", fill="blue") +
  geom_jitter(aes(colour = SampleGroup), size=3, alpha=0.9) +
  scale_fill_manual(values=species.color)+
  scale_color_manual(values = dot.color) +
  ylab(expression(F[IS])) +
  xlab(NULL) +
  theme_bw() +
  theme(#legend.position = "none",
        
        axis.text.x = element_text(size = 10, face = "italic" ),
        panel.border = element_rect(colour = "black", fill=NA))

Fis_Plot2 <- Fis_Plot + guides(fill="none")

#########################################################################
#########################################################################

pairwise.wilcox.test(heteperce2$percentage, heteperce2$group, p.adjust.method = "fdr")
pairwise.wilcox.test(heteperce2$F, heteperce2$group, p.adjust.method = "fdr")



#######################################################################
#######################################################################

safe_colorblind_palette <- c("Kd" = "#E69F00", 
                             "Kc" = "#56B4E9", 
                             "Kk" = "#009E73")

safe_shape_palette <- c("K. drynarioides" = 19, 
                        "K. cookei" = 15, 
                        "K. kauaiensis" = 17)
species.color <- RColorBrewer::brewer.pal(name="Set1", 3) %>%
  setNames(., c("K. cookei", "K. drynarioides", "K. kauaiensis"))

######################################################################
### Pi ###############################################################
#######################################################################

inp.pi <-read.table("KcKdKk_Kkn160.combined.genic.pi.txt",sep="\t",header=T)

### count average per species 
inp.pi.avg <- inp.pi %>%
  mutate(pop = recode(pop,
                      "Kc" = "K. cookei",
                      "Kd" = "K. drynarioides",
                      "Kk" = "K. kauaiensis"))  %>%
  group_by(pop, chromosome) %>% 
  dplyr::summarize(Avg_pi = mean(avg_pi, na.rm=TRUE)) 

### change the species order in the plot by average 
group_ordered <- with(inp.pi.avg,  reorder(pop,  Avg_pi,  mean))
inp.pi.avg$pop <- factor(inp.pi.avg$pop,levels = levels(group_ordered))

### change chromosome order
chrm_ordered <- levels(as.factor(inp.pi.avg$chromosome))[c(1,12,2:11)]
inp.pi.avg$chromosome <- factor(inp.pi.avg$chromosome,levels = chrm_ordered)


inp.pi.avg.mean <- inp.pi.avg %>%
  group_by(pop) %>% 
  dplyr::summarize(Avg = mean(Avg_pi, na.rm=TRUE)) %>%
  mutate(label = gsub("e", " %*% 10^", sprintf("%.2e", Avg)))%>%
  as.data.frame()

plot.inp.pi.avg <- ggplot() +
  geom_point(data = inp.pi.avg, aes(x=chromosome, y=Avg_pi, color = pop, shape = pop), size = 5, show.legend = F) +
  #geom_pointrange(data = inp.pi, aes(x=chromosome, y=avg_pi, color = pop, shape = pop, ymin=avg_pi-pi_sd, ymax=avg_pi+pi_sd)) +
  scale_color_manual(values=species.color)+
  scale_shape_manual(values=safe_shape_palette)+
  facet_grid(. ~ pop) +
  geom_hline(data = inp.pi.avg.mean, aes(yintercept = Avg), color = "blue", show.legend = F) +
  geom_text(data = inp.pi.avg.mean, 
            aes(x = "Kk_08", y = 0.0002, label = label), parse = TRUE, 
            inherit.aes = FALSE, size = 5) +
  ylab("Avg Pi / Chromosome")+
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10,face = "italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1, size = 10))

plot.inp.pi.avg


#######################################################################
#######################################################################
### Dxy & Fst #########################################################
#https://stackoverflow.com/questions/41764818/ggplot2-boxplots-with-points-and-fill-separation
#######################################################################

inp.dxy <-read.table("KcKdKk_Kkn160.combined.genic.dxy.txt",sep="\t",header=T)
inp.dxy$Group <- str_c(inp.dxy$pop1, "_", inp.dxy$pop2)

### calcuate average by each group
inp.dxy.avg <- inp.dxy %>%
  group_by(pop1, pop2, chromosome, Group) %>%
  dplyr::summarize(Mean = mean(avg_dxy, na.rm=TRUE)) %>%
  as.data.frame()

### change the species order in the plot by average 
group_ordered <- with(inp.dxy.avg,  reorder(Group,  Mean,  mean))
inp.dxy.avg$Group <- factor(inp.dxy.avg$Group,levels = levels(group_ordered))

### change chromosome order
chrm_ordered <- levels(as.factor(inp.pi.avg$chromosome))[c(1,12,2:11)]
inp.pi.avg$chromosome <- factor(inp.pi.avg$chromosome,levels = chrm_ordered)

### plot the pairwise Dxy between three species chrm by chrm
plot.inp.dxy.avg <- ggplot(inp.dxy.avg, aes(x= Mean, y= Group)) + 
  geom_boxplot( fill = "#56B4E9", alpha = 0.5, outlier.shape=NA, width = 0.23) +
  geom_point( color = "#E69F00", alpha = 0.3, size = 2) +
  stat_summary(fun=mean, geom="point", shape=18, size=3, color="blue", fill="blue") +
  #  ylab("Group")+
  xlab("Avg Dxy / Chromosome")+
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size = 10))


### make it more informative by adding the actual numbers below each plot
summary_df <- inp.dxy.avg %>% group_by(Group) %>% summarize(m=mean(Mean))
plot.inp.dxy.avg2 <- plot.inp.dxy.avg + geom_text(data=summary_df,
                                                  aes(x=m, label=formatC(m, format = "f", digits  = 5)),
                                                  color='blue', nudge_y = -0.3)


########################################################################
########################################################################
finalplot <- ggdraw() +
  draw_plot(plot.inp.pi.avg, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(Fis_Plot2  , x = 0, y = 0, width = 0.5, height = 0.5) +
  
  draw_plot(heteperce_Plot, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(plot.inp.dxy.avg2, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,  fontface = "bold",
                  x = c(0, 0, 0.5, 0.5), y = c(1, 0.5, 1 , 0.5))



pdf("Fig3_PixyKKref_onehavo_n160_nosubset4.pdf", width = 12, height = 9)
finalplot
dev.off()
