library(LEA)
library(tidyverse)
library(hues)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(grid)
library("ggplotify")

setwd(getwd())

sample.names <- read.table("samplename.txt")
colnames(sample.names)[1] <- "ind"

sample.names$population <- stringr::str_extract(sample.names$ind, "[^_]*_[^_]*")
sample.names$seq <- rownames(sample.names)
sample.names2 <- sample.names[order(sample.names$ind),] 
sample.names2$order <- 1:nrow(sample.names2) 
samplename3 <- sample.names2[order(as.numeric(sample.names2$seq), decreasing = F),]
samplename3$population <- gsub("MK.*", "MK", samplename3$population)
 
file2 <- list.files(pattern = "\\.snmfProject$")
project = load.snmfProject(file = file2)


kvalue <- read.csv("LEA4_crvalues.txt")
plotkvalue <- Map(function(x, y) as.data.frame(Q(project, K=x, run=y)), kvalue$K, kvalue$run)
plotkvalue2 <- Map(function(x, y) paste("K",x,"run",y, sep=""), kvalue$K, kvalue$run)

for (i in 1:length(plotkvalue)){
names(plotkvalue)[[i]] <- plotkvalue2[[i]]
}


#K3run1 <- as.data.frame(Q(project, K=3, run=1))
#K4run1 <- as.data.frame(Q(project, K=4, run=1))
#K5run6 <- as.data.frame(Q(project, K=5, run=6))
#K8run9 <- as.data.frame(Q(project, K=8, run=9))
#K9run7 <- as.data.frame(Q(project, K=9, run=7))
#K10run9 <- as.data.frame(Q(project, K=10, run=9))
#krun <- list("K = 3"=K3run1, "K = 4"=K4run1, "K = 5"=K5run6)

cbbPalette <- colorRampPalette( brewer.pal(n = 12, name = "Paired") )(15)

#c("#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")

plot_list = list()
plot_list2 = list()

for (i in 1:(length(plotkvalue)-7)) {
q_mat = data.frame()
q_mat = as.data.frame(plotkvalue[i])
colnames(q_mat) <- paste0("P", 1:length(colnames(q_mat)))
head(q_mat)

q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = samplename3$ind,
         accession = samplename3$population,
         order = samplename3$order)
q_df

q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 
q_df_long


q_df_prates <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_df_prates

# Plot the barplot
#my.colors<-c("tomato","lightblue", "olivedrab","gold", "#CC79A7", "#009E73") 

q_df_ordered <- q_df_long %>% 
  # assign the population assignment according to the max q value (ancestry proportion) and include the assignment probability of the population assignment
  group_by(individual) %>%
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>%
  # arrange the data set by the ancestry coefficients
  arrange(likely_assignment, assignment_prob) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_df_ordered

plotkrun <- q_df_prates %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop), show.legend = FALSE) +
  scale_fill_manual(values = cbbPalette) +
  #scale_fill_manual(values = brewer.pal(5, name = "Set3")) +
  #scale_fill_viridis_d() +
  labs(fill = "Populations") +
  ylab(paste0(names(plotkvalue[i])))+
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
		axis.title.y = element_text(),
        axis.text.x = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

plot_list[[i]] = plotkrun
plot_list2[[i]] = ggplotGrob(plotkrun)

}

####################

for (i in 9) {
q_mat = data.frame()
q_mat = as.data.frame(plotkvalue[i])
colnames(q_mat) <- paste0("P", 1:length(colnames(q_mat)))
head(q_mat)

q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = samplename3$ind,
         accession = samplename3$population,
         order = samplename3$order)
q_df

q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 
q_df_long


q_df_prates <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_df_prates

# Plot the barplot
#my.colors<-c("tomato","lightblue", "olivedrab","gold", "#CC79A7", "#009E73") 

q_df_ordered <- q_df_prates %>% 
  # assign the population assignment according to the max q value (ancestry proportion) and include the assignment probability of the population assignment
  group_by(individual) %>%
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>%
  # arrange the data set by the ancestry coefficients
  arrange(likely_assignment, assignment_prob) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  ungroup() %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_df_ordered

plotkrun <- q_df_prates %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop), show.legend = FALSE) +
  scale_fill_manual(values = cbbPalette) +
  #scale_fill_manual(values = brewer.pal(5, name = "Set3")) +
  #scale_fill_viridis_d() +
  labs(fill = "Populations") +
  ylab(paste0(names(plotkvalue[i])))+
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
		axis.title.y = element_text(),
        axis.text.x = element_text(angle = -90, hjust= 0.001, vjust=0.5),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

plot_list[[i]] = plotkrun
plot_list2[[i]] = ggplotGrob(plotkrun)

}

plot_list3 <- as.ggplot(do.call("rbind", plot_list2), size = "last") +  
  scale_fill_manual(values=cbbPalette)


pdf("LEA6_Krun_structure.pdf", width = 25, height = 35)
#for (i in 1:length(plotkvalue)) {print(plot_list[[i]])}
#ggarrange(plotlist=plot_list, ncol = 1, nrow = 20)
plot_list3
dev.off()

pdf("LEA7_Best_Krun_structure.pdf", width = 25, height = 5)
plot_list[[7]] + 
  scale_fill_manual(values=cbbPalette) +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
		axis.title.y = element_text(),
        axis.text.x = element_text(angle = -90, hjust= 0.001, vjust=0.5),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
dev.off()



save.image(file = "LEAstruture_n65.RData")
saveRDS(plot_list, file="plotlist_n65.RData")

#plot_list <- readRDS("plotlist.RData")


