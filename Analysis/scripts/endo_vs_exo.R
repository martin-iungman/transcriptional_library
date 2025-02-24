library(tidyverse)
# Max CAGE expression along samples vs. exogenous activity (to better consider tissue specific promoters)
tissue=read_tsv("Analysis/Tables/tissue_CAGE_activity.tsv") %>% mutate(group="tissue")
cell=read_tsv("Analysis/Tables/primary_cell_CAGE_activity.tsv")%>% mutate(group="primary_cell")
df=rbind(tissue ,cell) %>% mutate(tpm=1E6*counts/libsize) %>% arrange(desc(tpm)) %>% group_by(name) %>% mutate(order=row_number(desc(tpm))) %>% filter(order==1) 
df%>% write_tsv("Analysis/Tables/maxcount_all_samples.tsv")

data=read_tsv("Analysis/Tables/full_data.tsv") %>% left_join(df, by=c("seq_id"="name")) %>% mutate(tpm=replace_na(tpm, 0)) 

data %>% group_by(rep) %>% mutate(mean=dense_rank(mean)) %>% 
  ggplot(aes(mean, tpm+0.5))+geom_point()+facet_wrap(~rep)+scale_y_log10()+labs(x="Ranked mean exogenous activity",y="TPM (max value along CAGE FANTOM5)")

data %>% group_by(rep, mean_sw) %>% mutate( tpm=mean(tpm)) %>% 
  ggplot(aes(mean_sw, tpm+0.5))+geom_point()+facet_wrap(~rep)+scale_y_log10()+labs(x="Ranked groups of mean exogenous activity (bin size=100)",y="median TPM (max value along CAGE FANTOM5)")

data  %>% mutate(`High tissue specificity`=gini_tissue>median(gini_tissue, na.rm=T)) %>% group_by(rep) %>% mutate(mean=dense_rank(mean)) %>% 
  ggplot(aes(mean, tpm, col=`High tissue specificity`))+geom_point()+facet_wrap(~rep)+scale_y_log10()+labs(x="Ranked mean exogenous activity",y="TPM (max value along CAGE FANTOM5)")

data  %>% mutate(`High tissue specificity`=gini_tissue>median(gini_tissue, na.rm=T)) %>% 
  group_by(rep, `High tissue specificity`) %>% arrange(mean) %>% mutate(mean=dense_rank(mean) %>% cut_width(width=100) %>% as.numeric()) %>% group_by(rep, `High tissue specificity`, mean) %>% 
  mutate( tpm=mean(tpm)) %>% 
  ggplot(aes(mean, tpm, col=`High tissue specificity`))+geom_point()+facet_wrap(~rep)+scale_y_log10()+labs(x="Ranked groups of mean exogenous activity (bin size=100)",y="mean TPM (max value along CAGE FANTOM5)")

data  %>% mutate(`High tissue specificity`=gini_tissue>median(gini_tissue, na.rm=T)) %>% 
  group_by(rep, `High tissue specificity`) %>% arrange(mean) %>% mutate(mean=dense_rank(mean) %>% cut_width(width=100) %>% as.numeric()) %>% group_by(rep, `High tissue specificity`, mean) %>% 
  mutate( tpm=median(tpm)) %>% 
  ggplot(aes(mean, tpm, col=`High tissue specificity`))+geom_point()+facet_wrap(~rep)+scale_y_log10()+labs(x="Ranked groups of mean exogenous activity (bin size=100)",y="median TPM (max value along CAGE FANTOM5)")

data %>% mutate(gini_tissue=cut_number(gini_tissue, n=3) %>% as.numeric()) %>% 
  group_by(rep, gini_tissue) %>% arrange(mean) %>% mutate(mean=dense_rank(mean) %>% cut_width(width=100) %>% as.numeric()) %>% group_by(rep, gini_tissue, mean) %>% 
  mutate( tpm=median(tpm)) %>% filter(!is.na(gini_tissue), gini_tissue!=2) %>% mutate(`Expression pattern`=ifelse(gini_tissue==1, "Housekeeping", "Tissue-specific")) %>% 
  ggplot(aes(mean, tpm, col=`Expression pattern`))+geom_point()+geom_smooth(method="lm")+
  facet_wrap(~rep)+scale_y_log10()+labs(x="Ranked groups of mean exogenous activity (bin size=100)",y="median TPM (max value along CAGE FANTOM5)")+
  scale_colour_manual(values=c("#1B8C8E", "#0D2C54"))+ggpubr::theme_pubr()
ggsave("Plots/Fig3/tissue_specificty_scatter.pdf")




lm_l=data  %>% mutate(`High tissue specificity`=gini_tissue>median(gini_tissue, na.rm=T)) %>% filter(!is.na(`High tissue specificity`)) %>% 
  group_by(rep, `High tissue specificity`) %>% arrange(mean) %>% mutate(mean=dense_rank(mean) %>% cut_width(width=100) %>% as.numeric(), tpm=log10(tpm+0.1), group=paste(rep, `High tissue specificity`, sep="_")) %>%  split(.$group) %>%  
  map(~lm(tpm~mean, data=.x)) 
lm_l %>% imap(~tibble(residuals=residuals(.x), group=.y) %>% separate(group, into=c("rep", "High tissue specificity"), sep="_")) %>% list_rbind() %>% 
  ggplot(aes(residuals, fill=`High tissue specificity`))+geom_density(alpha=0.5)+facet_wrap(~rep)

lm_l=data  %>% mutate(gini_tissue=cut_number(gini_tissue, n=3) %>% as.numeric()) %>%  filter(!is.na(gini_tissue)) %>% 
  group_by(rep, gini_tissue) %>% arrange(mean) %>% mutate(mean=dense_rank(mean) , tpm=log10(tpm+0.1), group=paste(rep, gini_tissue, sep="_")) %>%  split(.$group) %>%  
  map(~lm(tpm~mean, data=.x)) 
lm_l %>% imap(~tibble(residuals=residuals(.x), group=.y) %>% separate(group, into=c("rep", "Expression pattern"), sep="_")) %>% list_rbind()%>% filter(!is.na(`Expression pattern`), `Expression pattern`!=2) %>% mutate(`Expression pattern`=ifelse(`Expression pattern`==1, "Housekeeping", "Tissue-specific")) %>% 
  ggplot(aes(residuals, fill=`Expression pattern`))+geom_density(alpha=0.5)+facet_wrap(~rep)+scale_fill_manual(values=c("#1B8C8E", "#0D2C54"))+ggpubr::theme_pubr()
ggsave("Plots/Fig3/tissue_specificity_residuals.pdf")

lm_l %>% imap(~tibble(residuals=residuals(.x), group=.y, mean_rank=.x$model$mean) %>% separate(group, into=c("rep", "High tissue specificity"), sep="_")) %>% list_rbind() %>% group_by(rep,`High tissue specificity`) %>% mutate(mean_rank=cut_width(mean_rank, width=100)) %>% group_by(rep, mean_rank, `High tissue specificity`) %>% summarise(residuals=median(abs(residuals))) %>% filter(`High tissue specificity`!=2) %>% 
  ggplot(aes(mean_rank %>% as.numeric(), residuals, col=`High tissue specificity`))+geom_point()+facet_wrap(~rep) +geom_smooth()+labs(x="Mean activity", y="Median absolute residual")