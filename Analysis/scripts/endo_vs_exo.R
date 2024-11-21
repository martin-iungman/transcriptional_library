library(tidyverse)
# Max CAGE expression along samples vs. exogenous activity (to better consider tissue specific promoters)
tissue=read_tsv("Tables/tissue_CAGE_activity.tsv") %>% mutate(group="tissue")
cell=read_tsv("Tables/primary_cell_CAGE_activity.tsv")%>% mutate(group="primary_cell")
df=rbind(tissue %>% rename(sample=tissue),cell) %>% mutate(tpm=1E6*counts/libsize) %>% arrange(desc(tpm)) %>% group_by(name) %>% mutate(order=row_number(desc(tpm))) %>% filter(order==1) 
df%>% write_tsv("Tables/maxcount_all_samples.tsv")

data=read_tsv("Tables/full_data.tsv") %>% left_join(df, by=c("seq_id"="name")) %>% mutate(tpm=replace_na(tpm, 0)) 

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


lm_l=data  %>% mutate(`High tissue specificity`=gini_tissue>median(gini_tissue, na.rm=T)) %>% filter(!is.na(`High tissue specificity`)) %>% 
  group_by(rep, `High tissue specificity`) %>% arrange(mean) %>% mutate(mean=dense_rank(mean) %>% cut_width(width=100) %>% as.numeric(), tpm=log10(tpm+0.1), group=paste(rep, `High tissue specificity`, sep="_")) %>%  split(.$group) %>%  
  map(~lm(tpm~mean, data=.x)) 
lm_l %>% imap(~tibble(residuals=residuals(.x), group=.y) %>% separate(group, into=c("rep", "High tissue specificity"), sep="_")) %>% list_rbind() %>% 
  ggplot(aes(residuals, fill=`High tissue specificity`))+geom_density(alpha=0.5)+facet_wrap(~rep)
