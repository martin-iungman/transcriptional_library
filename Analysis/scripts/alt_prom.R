library(rtracklayer)
library(tidyverse)
inputFiles = map(c("tissue", "primary_cell"), ~list.files(paste0("External_data/FANTOM5/hg38_",.x), full.names = T,pattern=".nobarcode.ctss.bed.gz$"))
lib_gr=import.bed("External_data/EPD/human38_epdnew.bed")
cell_ont=map(c("tissue", "primary_cell"), ~readxl::read_xls(paste0("External_data/FANTOM5/",.x,"_ontology_FANTOM5.xlsx")) %>% dplyr::rename("sample_id"=`Sample ID`))
libsize_tissues=map(c("tissue", "primary_cell"),~read_tsv(paste0("External_data/FANTOM5/library_size_", .x,"_by_facets.tsv")))


sample_id=map(inputFiles, ~.x%>%str_extract("CNhs.{5}"))
cell_ont=map2(inputFiles, sample_id, ~tibble(filename=.x, sample_id=.y)) %>% map2(cell_ont, ~.x%>%inner_join(.y, by=c("sample_id")))
cell_ont_list=cell_ont %>% list_rbind() %>%group_split(`Facet ontology term`)
names(cell_ont_list)<-group_keys(list_rbind(cell_ont)%>%group_by(`Facet ontology term`))%>%as_vector()

load_cage=function(df){
  ls=map(df$filename, ~import.bed(.x)%>% #load bed files
           plyranges::join_overlap_left_directed(lib_gr,.)%>% #overlap with the library
           values()%>%as_tibble()%>%group_by(name.x)%>%
           dplyr::summarise(counts=sum(score.y)%>%replace_na(0)))%>% #sum the counts overlapping each promoter
    list_rbind()%>%group_by(name.x)%>%dplyr::summarise(counts=sum(counts)) #join all the samples
  return(ls)
}

ls=read_rds("Analysis/Tables/ls.rds")
#ls=cell_ont_list%>%imap(~load_cage(.x)%>%mutate(sample=.y))%>%list_rbind()%>%dplyr::rename(name=name.x)
ls=ls%>%left_join(list_rbind(libsize_tissues), by=c("sample"="Facet ontology term"))
ls=mutate(ls, tpm=counts*1E6/libsize)
ls$gene_sym=str_remove(ls$name, "_.{1,2}$")

tot_counts=ls %>% group_by(name, gene_sym) %>% summarise(tpm=sum(tpm)) %>% 
  group_by(gene_sym) %>% mutate(max_tpm=tpm==max(tpm), tpm_gene=sum(tpm), perc_tpm_gene=tpm/tpm_gene)
unique_prom=tot_counts %>% filter(perc_tpm_gene>0.99, tpm_gene>5)
non_expressed=tot_counts %>% filter(tpm_gene<5) 

ls2=ls %>% group_by(gene_sym, sample) %>% mutate(tpm_gene_sample=sum(tpm), perc_tpm_gene_sample=tpm/tpm_gene_sample, is_highest_in_sample=tpm==max(tpm)) %>% ungroup()
ls2=ls2 %>% filter(tpm_gene_sample>0) %>% group_by(name) %>% summarise(highest_in_N_samples=sum(is_highest_in_sample)) %>% inner_join(ls2,.)
ls2=ls2 %>% ungroup() %>% group_by(gene_sym,sample) %>% mutate(main_N_samples=highest_in_N_samples/sum(highest_in_N_samples)) 
main_prom=ls2 %>% filter(main_N_samples>0.5) %>% left_join(tot_counts %>% select(-tpm), by=c("name", "gene_sym")) %>% filter(perc_tpm_gene>0.5) %>% filter(!gene_sym%in%c(unique_prom$gene_sym, non_expressed$gene_sym)) # %>% select(name, gene_sym) %>% unique()

wide_df=ls2 %>% select(name, gene_sym, tpm, sample) %>% inner_join(main_prom %>% select(name, gene_sym, tpm, sample) %>% unique() %>% rename(main_name=name, main_tpm=tpm), by=c("sample","gene_sym")) %>% filter(name!=main_name)
cor_data_pearson=wide_df %>% ungroup() %>% filter(tpm+main_tpm>0) %>% group_split(name) %>% map(~summarise(.x, cor=cor(tpm, main_tpm, method = "pearson"), n_samples=n(), name=unique(name), main_name=unique(main_name))) %>% list_rbind()
cor_data_spearman=wide_df %>% ungroup() %>% filter(tpm+main_tpm>0) %>% group_split(name) %>% map(~summarise(.x, cor=cor(tpm, main_tpm, method = "spearman"), n_samples=n(), name=unique(name), main_name=unique(main_name))) %>% list_rbind()
cor_data=left_join(cor_data_pearson, cor_data_spearman, by=c("n_samples", "name", "main_name"), suffix = c("_pearson", "_spearman")) 
cor_data %>% ggplot(aes(cor_pearson, cor_spearman, col=n_samples))+geom_point()
wide_df=wide_df %>% left_join(cor_data)
write_tsv(wide_df, "Analysis/Tables/prom_alt_cor.tsv")

#example of highly correlated promoters
wide_df %>% filter(name=="NABP1_2") %>% ggplot(aes(tpm, main_tpm))+geom_point()+scale_x_log10()+scale_y_log10()

#example of highly anti-correlated promoters
wide_df %>% filter(name=="ACY3_1") %>% ggplot(aes(tpm, main_tpm))+geom_point()+scale_x_log10()+scale_y_log10()


ls2  %>%  filter(gene_sym=="SH3YL1") %>% filter(tpm_gene_sample>0) %>% mutate(promid=paste0("prom", name %>% str_extract("_.{1,2}$"))) %>% select(gene_sym, promid, sample, tpm) %>% pivot_wider(names_from = "promid", values_from = "tpm") %>% 
  ggplot(aes(prom_1, prom_2)) +geom_point()+geom_abline(linetype="dashed")+scale_x_log10()+scale_y_log10()
  

data=read_tsv("Analysis/Tables/full_data.tsv")
data %>% mutate(main_prom=str_remove(seq_id, "^FP.{6}_")%in%unique(wide_df$main_name)) %>% summarise(sum(main_prom)) 
data %>% mutate(main_prom=str_remove(seq_id, "^FP.{6}_")%in%unique(wide_df$main_name)) %>% 
  group_by(mean_sw,rep) %>% summarise(main_prom=sum(main_prom)/n()) %>% 
  ggplot(aes(mean_sw,main_prom))+geom_point(col="#14AFB2")+geom_smooth(col="#216869")+theme(text=element_text(size=20))+
  facet_wrap(~rep)+ylim(0,1)+labs(x="Mean activity (binned by rank)",y="Proportion of promoters")+ggtitle("Main promoters")
data %>% mutate(main_prom=str_remove(seq_id, "^FP.{6}_")%in%unique(wide_df$main_name)) %>% 
  ggplot(aes(as.numeric(mean_sw), var_rank_sw, col=main_prom))+geom_point(size=0.2, alpha=0.5)+
  facet_wrap(~rep)+xlab("Mean rank group")+ylab("Variance rank (SW)")+
  geom_smooth()+labs(col="Unique promoters")+scale_color_manual(values=c("#0D2C54", "#AD343E"))+theme_pubr(base_size = 20)
data %>% mutate(main_prom=str_remove(seq_id, "^FP.{6}_")%in%unique(wide_df$main_name)) %>% ungroup() %>% group_split(rep) %>% 
  map(~roc_curve(.x, "main_prom") %>% mutate(rep=unique(.x$rep), AUC=auc(.))) %>% list_rbind()%>% 
  ggplot(aes(FPR,TPR, col=rep))+geom_line()+geom_abline(linetype="dashed")+
  geom_text(data=.%>% select(AUC, rep) %>% mutate(AUC=paste("AUC:", round(AUC,3)))%>% unique() %>% 
              bind_cols(tibble(TPR=c(0.9,0.8), FPR=0.02)), aes(label=AUC), hjust="left")+
  ggtitle("Main promoters ROC curve")+scale_color_manual(values=c("#0D2C54", "#AD343E"))+
  theme_pubr(base_size = 20)+labs(linetype="", col="")

data %>% mutate(main_prom=str_remove(seq_id, "^FP.{6}_")%in%unique(unique_prom$name)) %>% 
  group_by(mean_sw,rep) %>% summarise(main_prom=sum(main_prom)/n()) %>% 
  ggplot(aes(mean_sw,main_prom))+geom_point(col="#14AFB2")+geom_smooth(col="#216869")+theme(text=element_text(size=20))+
  facet_wrap(~rep)+ylim(0,1)+labs(x="Mean activity (binned by rank)",y="Proportion of promoters")+ggtitle("Unique promoters")
data %>% mutate(main_prom=str_remove(seq_id, "^FP.{6}_")%in%unique(unique_prom$name)) %>% 
  ggplot(aes(as.numeric(mean_sw), var_rank_sw, col=main_prom))+geom_point(size=0.2, alpha=0.5)+
  facet_wrap(~rep)+xlab("Mean rank group")+ylab("Variance rank (SW)")+
  geom_smooth()+labs(col="Unique promoters")+scale_color_manual(values=c("#0D2C54", "#AD343E"))+theme_pubr(base_size = 20)
data %>% mutate(main_prom=str_remove(seq_id, "^FP.{6}_")%in%unique(unique_prom$name)) %>% ungroup() %>% group_split(rep) %>% 
  map(~roc_curve(.x, "main_prom") %>% mutate(rep=unique(.x$rep), AUC=auc(.))) %>% list_rbind()%>% 
  ggplot(aes(FPR,TPR, col=rep))+geom_line()+geom_abline(linetype="dashed")+
  geom_text(data=.%>% select(AUC, rep) %>% mutate(AUC=paste("AUC:", round(AUC,3)))%>% unique() %>% 
              bind_cols(tibble(TPR=c(0.9,0.8), FPR=0.02)), aes(label=AUC), hjust="left")+
  ggtitle("Unique promoters ROC curve")+scale_color_manual(values=c("#0D2C54", "#AD343E"))+
  theme_pubr(base_size = 20)+labs(linetype="", col="")
