library(rtracklayer)
library(tidyverse)

dnasa=rtracklayer::import.bw("External_data/ENCFF165GHP_DNase.bigWig", as = "GRanges")
lib=rtracklayer::import.bed("Library_data/res/library.bed")
strand(lib) <- "*"
ol=plyranges::find_overlaps( dnasa, lib)
seqlevels(lib)->seqlevels(dnasa)

lib=rtracklayer::import.bed("Library_data/res/library.bed")
df=as_tibble(ol) %>% select(-strand) %>% left_join(as_tibble(lib) %>% select(name, strand), by="name") %>% filter(str_detect(name, "^FP")) %>% 
  group_by(name) %>% mutate(range_order=((ifelse(strand=="-", row_number(-start), row_number(start))-1)*25), rel_score=score.x/sum(score.x))
df %>% ggplot(aes(range_order %>% as_factor(), score.x))+geom_violin(draw_quantiles = 0.5)+geom_smooth(aes(x=range_order %>% as_factor() %>% as.numeric()))
df %>% ggplot(aes(range_order %>% as_factor(), fill=score.x>1))+geom_bar(position="fill")
df_matrix <- df %>% ungroup() %>%  select(name, score.x, range_order) %>% 
  pivot_wider(names_from = range_order, values_from = score.x) %>%
  select(-name) %>% # Remove the name column before converting to matrix
  as.matrix()
rownames(df_matrix) <- unique(df$name)
heatmap_order <- hclust(dist(df_matrix))$order

df_matrix_ordered <- df_matrix[heatmap_order, ]
pheatmap::pheatmap(df_matrix_ordered, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Color scale
         cluster_rows = FALSE,  # We've already clustered
         cluster_cols = FALSE,  # Usually you'd cluster columns too, but I'm leaving it out based on the previous code
         show_rownames = TRUE,  # Show row names
         show_colnames = TRUE,  # Show column names
         scale = "none", # If you want to scale the data before plotting.
         main = "Heatmap of Scores (Clustered)") # Title

####
df_matrix <- df %>% ungroup() %>%  select(name, rel_score, range_order) %>% 
  pivot_wider(names_from = range_order, values_from = rel_score) %>%
  select(-name) %>% # Remove the name column before converting to matrix
  as.matrix()
rownames(df_matrix) <- unique(df$name)
#heatmap_order <- hclust(dist(df_matrix))$order

df_matrix_ordered <- df_matrix[heatmap_order, ]
pheatmap::pheatmap(log10(df_matrix_ordered), 
                   color = colorRampPalette(c("blue", "white", "red"))(100), # Color scale
                   cluster_rows = FALSE,  # We've already clustered
                   cluster_cols = FALSE,  # Usually you'd cluster columns too, but I'm leaving it out based on the previous code
                   show_rownames = TRUE,  # Show row names
                   show_colnames = TRUE,  # Show column names
                   scale = "none", # If you want to scale the data before plotting.
                   main = "Heatmap of Scores (Clustered)") # Title


# df_plot <- df  %>%  select(name, score.x, range_order)%>% # Create a separate data frame for plotting
#   pivot_wider(names_from = range_order, values_from = score.x) %>%
#   mutate(name = factor(name, levels = unique(name)[heatmap_order])) %>%
#   # Add a 'y_level' column for easier plotting
#   mutate(y_level = factor(1)) # Initialize with 1
# 
# Create the heatmap (Corrected and Improved)
# p <- ggplot(df_plot, aes(x = name, y = y_level, fill = `1`)) + # Use df_plot
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.ticks = element_blank()) +
#   labs(title = "Heatmap of Scores (Clustered)", fill = "Score")
# 
# # Dynamically add geom_tile layers (Corrected)
# for (i in 2:12) {
#   p <- p + geom_tile(data = df_plot %>% mutate(y_level = factor(i)), # Data for this layer
#                      aes(x = name, y = y_level, fill = .data[[as.character(i)]]))
# }
# 
# print(p) # Print the plot

data=read_tsv("Analysis/Tables/full_data.tsv")

df %>% group_by(name) %>% summarise(mean_occupancy_dnasa=mean(score.x)) %>% left_join(data,., by=c("seq_id"="name")) %>% 
  group_by(mean_sw, rep) %>% summarise(mean_occupancy_dnasa=mean(mean_occupancy_dnasa)) %>% ggplot(aes(mean_sw, mean_occupancy_dnasa))+
  geom_point(col="#AD343E")+
  geom_smooth(col="#08415c")+facet_wrap(~rep)+labs(x="Mean activity (binned by rank)", y="Mean nucleosome occupancy (DNase-seq)")+theme(text=element_text(size=20))+scale_y_log10()

df %>% group_by(name) %>% summarise(mean_occupancy_dnasa=mean(score.x)) %>% left_join(data,., by=c("seq_id"="name")) %>% 
  ggplot(aes(mean, mean_occupancy_dnasa))+
  geom_hex(bins=200)+
  geom_smooth(col="#08415c")+facet_wrap(~rep)+labs(x="Mean activity (binned by rank)", y="Mean nucleosome occupancy (DNase-seq)")+theme(text=element_text(size=20))
df %>% group_by(name) %>% summarise(mean_occupancy_dnasa=mean(score.x)) %>% left_join(data,., by=c("seq_id"="name")) %>% 
  ggplot(aes(mean_rank_sw, mean_occupancy_dnasa))+
  geom_hex()+
  geom_smooth(col="#08415c")+facet_wrap(~rep)+labs(x="Mean activity (binned by rank)", y="Mean nucleosome occupancy (DNase-seq)")+theme(text=element_text(size=20))

df %>% select(name, score.x, range_order) %>% mutate(range_order=paste0("range_", range_order)) %>% pivot_wider(names_from = range_order, values_from = score.x)  %>% left_join(data,., by=c("seq_id"="name")) %>% 
  group_by(rep, mean_sw) %>% summarise(across(starts_with("range"), mean)) %>% pivot_longer(starts_with("range")) %>% 
  ggplot(aes(mean_sw, value, col=name))+geom_point()+geom_smooth()+facet_wrap(~rep)
df %>% select(name, score.x, range_order) %>% mutate(range_order=paste0("range_", range_order)) %>% pivot_wider(names_from = range_order, values_from = score.x)  %>% left_join(data,., by=c("seq_id"="name")) %>% 
  group_by(rep, mean_sw, CGI) %>% summarise(across(starts_with("range"), mean)) %>% pivot_longer(starts_with("range")) %>% 
  ggplot(aes(mean_sw, value, col=name))+geom_point()+geom_smooth()+facet_wrap(~rep+CGI)

a=df %>% select(name, score.x, range_order) %>% mutate(range_order=paste0("range_", range_order)) %>% pivot_wider(names_from = range_order, values_from = score.x)  %>% left_join(data,., by=c("seq_id"="name")) %>% 
  group_by(rep, mean_sw) %>% summarise(across(starts_with("range"), mean)) %>% pivot_longer(starts_with("range")) %>% filter(mean_sw<=10) %>% group_by(rep, name) %>% summarise(zero=mean(value))
df %>% select(name, score.x, range_order) %>% mutate(range_order=paste0("range_", range_order)) %>% pivot_wider(names_from = range_order, values_from = score.x)  %>% left_join(data,., by=c("seq_id"="name")) %>% 
  group_by(rep, mean_sw) %>% summarise(across(starts_with("range"), mean)) %>% pivot_longer(starts_with("range")) %>% left_join(a) %>% mutate(value=value/zero) %>% 
  ggplot(aes(mean_sw, value, col=name))+geom_point()+geom_smooth()+facet_wrap(~rep)


df %>% group_by(name) %>% summarise(mean_occupancy_dnasa=mean(score.x, na.rm=T)) %>% left_join(data,., by=c("seq_id"="name"))  %>% 
  mutate(high_occupancy=mean_occupancy_dnasa>median(mean_occupancy_dnasa)) %>% 
  ggplot(aes(as.numeric(mean_sw), var_rank_sw, col=high_occupancy))+geom_point(size=0.2, alpha=0.5)+facet_wrap(~rep)+
  xlab("Mean rank group")+ylab("Variance rank (SW)")+
  geom_smooth()+labs(col="high_occupancy")+scale_color_manual(values = c("#08415c", "#AD343E"))+ggpubr::theme_pubr(base_size = 20)

df %>% select(name, score.x, range_order) %>% mutate(range_order=paste0("range_", range_order)) %>% pivot_wider(names_from = range_order, values_from = score.x)  %>% 
  left_join(data,., by=c("seq_id"="name")) %>% 
  mutate(high_occupancy=range_225>median(range_225)) %>% 
  ggplot(aes(as.numeric(mean_sw), var_rank_sw, col=high_occupancy))+geom_point(size=0.2, alpha=0.5)+facet_wrap(~rep)+
  xlab("Mean rank group")+ylab("Variance rank (SW)")+
  geom_smooth()+labs(col="high_occupancy")+scale_color_manual(values = c("#08415c", "#AD343E"))+ggpubr::theme_pubr(base_size = 20)
df %>% select(name, score.x, range_order) %>% mutate(range_order=paste0("range_", range_order)) %>% 
  pivot_wider(names_from = range_order, values_from = score.x)  %>% 
  left_join(data,., by=c("seq_id"="name")) %>% 
  group_by(rep, var_rank_sw) %>% summarise(across(starts_with("range"), mean)) %>% 
  pivot_longer(starts_with("range")) %>% 
  ggplot(aes(var_rank_sw, value, col=name))+geom_point()+geom_smooth()+facet_wrap(~rep)

df %>% select(name, score.x, range_order) %>% mutate(range_order=paste0("range_", range_order)) %>% pivot_wider(names_from = range_order, values_from = score.x)  %>% left_join(data,., by=c("seq_id"="name"))  %>%
  ungroup() %>% group_split(rep) %>% expand.grid(unique(paste0("range_", df$range_order))) %>% 
  map2(.$Var1, .$Var2,~roc_curve(.x, .y) %>% mutate(rep=unique(.x$rep), AUC=auc(.), range_order=.y)) %>% list_rbind()%>% 
  ggplot(aes(FPR,TPR, col=range_order, linetype=rep))+geom_line()+geom_abline(linetype="dashed")
  #geom_text(data=.%>% select(AUC, rep) %>% mutate(AUC=paste("AUC:", round(AUC,3)))%>% unique() %>% bind_cols(tibble(TPR=c(0.9,0.8), FPR=0.02)), aes(label=AUC), hjust="left")+ggtitle("TATA-box")+scale_color_manual(values=c("#0D2C54", "#AD343E"))+theme_pubr(base_size = 20)+labs(linetype="", col="")
