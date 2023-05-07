library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(phytools)
library(hdp)

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

mut_sigs_theme=theme(strip.text.x = element_text(size=6,margin = margin(0.6,0,0.6,0, "mm")),
                     strip.text.y=element_text(size=6,margin = margin(0,0.6,0,0.6, "mm")),
                     axis.text.x = element_text(size=3),
                     axis.text.y=element_text(size=5),
                     axis.title.x=element_text(size=7),
                     axis.title.y=element_text(size=7),
                     axis.ticks.x=element_blank(),
                     axis.ticks.y=element_line(size=0.25),
                     legend.text = element_text(size=5),
                     legend.title = element_text(size=7),
                     legend.key.size=unit(2.5,"mm"),
                     strip.background = element_rect(size=0.25),
                     panel.grid.major = element_line(size=0.25),
                     panel.border = element_rect(size=0.25))

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch119/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
source(paste0(root_dir,"/data/HSCT_functions.R"))
setwd(tree_mut_dir); source("treemut.R");setwd(root_dir)
vcf_header_path="~/R_work/Phylogeny_of_foetal_haematopoiesis/Data/vcfHeader.txt"
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired")
names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]
names(DorR_cols)<-c("D","R")
  
#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/sample_metadata_full.Rds"))
trees_list<-readRDS(paste0(root_dir,"/data/tree_lists.Rds"))
details_list<-readRDS(paste0(root_dir,"/data/details_lists.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_list)) {assign(x,details_list[[x]])}

#### MUTATIONAL SIGNATURE ANALYSIS USING HDP
mut_example_multi=readRDS(paste0(HDP_folder,"/HDP_multi_chain.Rdata"))
mutations=read.table(paste0(HDP_folder,"/trinuc_mut_mat.txt"))
key_table=read.table(paste0(HDP_folder,"/key_table.txt"))
sig_profiles=mut_mat_HDP_comp(mut_example_multi,plot=T)

#Altered plotting function to allow y axis to be free
plot_96_profile=function (mut_matrix, colors = NA, ymax = 0.2, condensed = FALSE) 
{
  freq <- full_context <- substitution <- context <- NULL
  if (MutationalPatterns:::.is_na(colors)) {
    colors <- MutationalPatterns:::COLORS6
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }
  norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
  tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>% 
    dplyr::mutate(substitution = stringr::str_replace(full_context, 
                                                      "\\w\\[(.*)\\]\\w", "\\1"), context = stringr::str_replace(full_context, 
                                                                                                                 "\\[.*\\]", "\\.")) %>% dplyr::select(-full_context) %>% 
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", 
                        values_to = "freq") %>% dplyr::mutate(sample = factor(sample, 
                                                                              levels = unique(sample)))
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  plot <- ggplot(data = tb, aes(x = context, y = freq, fill = substitution, width = width)) +
    geom_bar(stat = "identity", colour = "black", size = 0.2) +
    scale_fill_manual(values = colors) +
    facet_grid(sample ~ substitution,scales="free_y") +
    ylab("Relative contribution") +
    guides(fill = "none") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12, vjust = 1),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 9),
          panel.grid.major.x = element_blank(),
          panel.spacing.x = unit(spacing, "lines"))
  return(plot)
}

#Use this to plot the main four signatures
sigs_plot<-plot_96_profile(sig_profiles[,2:5],condensed=T)+mut_sigs_theme
ggsave(filename = paste0(plots_dir,"extracted_sigs_plot.pdf"),plot = sigs_plot,width=4,height=2.5)

#Compare signatures to known signatures
pheatmap::pheatmap(MutationalPatterns::cos_sim_matrix(sig_profiles,get_known_signatures()),
                   cluster_rows = F,
                   cluster_cols = F)

cosmic_sigs<-MutationalPatterns::get_known_signatures()
MutationalPatterns::cos_sim_matrix(sig_profiles,cosmic_sigs)%>%
  apply(1, function(x) {idxs<-order(x,decreasing = T)[1:2];return(data.frame(
                                                              COSMIC_sig1=colnames(cosmic_sigs)[idxs[1]],
                                                              cosine_sim1=x[idxs[1]],
                                                              COSMIC_sig2=colnames(cosmic_sigs)[idxs[2]],
                                                              cosine_sim2=x[idxs[2]]))})

#Plot signature contributions in each sample: relative & absolute
sig_cols=rev(RColorBrewer::brewer.pal(5,"Paired"))
plot.mutsigs.relative<-sample_metadata%>%
  filter(!is.na(N0))%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  dplyr::select(Sample,DorR,Pair_new,N0,N1,N2,N3,N4)%>%
  gather(-Sample,-DorR,-Pair_new,key = "Signature",value="Proportion")%>%
  mutate(Signature=factor(Signature,levels=rev(paste0("N",0:4))))%>%
  ggplot(aes(x=Sample,y=Proportion,fill=Signature))+
  geom_bar(stat="identity",col=NA,size=0.05)+
  geom_tile(aes(x=Sample,y=-0.05,col=factor(DorR),height=0.0125),size=2,fill=NA,inherit.aes = F)+
  scale_fill_manual(values = sig_cols)+
  scale_color_manual(values = DorR_cols)+
  facet_grid(cols=vars(Pair_new),scale="free",space = "free")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  my_theme+
  theme(legend.key.size = unit(3,"mm"),
        axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  labs(fill="Signature",col="Donor or\nRecipient")
ggsave(filename = paste0(plots_dir,"plot.mutsigs.relative.pdf"),plot = plot.mutsigs.relative,width=7,height=2.5)

plot.mutsigs.absolute<-sample_metadata%>%
  filter(!is.na(N0))%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  dplyr::select(Sample,DorR,Pair_new,N0_abs,N1_abs,N2_abs,N3_abs,N4_abs)%>%
  mutate(N4_abs=sapply(N4_abs,function(x) min(x,500)))%>% #To aid visualization limit total contribution of `N4 (APOBEC) to 500 mutations`- affects 4 samples
  mutate(N1_abs=sapply(N1_abs,function(x) min(x,2000)))%>% #To aid visualization limit total contribution of `N1 (HSPC sig) to 2000 mutations`- affects 4 samples
  gather(-Sample,-Pair_new,-DorR,key = "Signature",value="Number_of_mutations")%>%
  mutate(Signature=gsub("_abs","",Signature))%>%
  mutate(Signature=factor(Signature,levels=rev(paste0("N",0:4))))%>%
  ggplot(aes(x=Sample,y=Number_of_mutations,fill=factor(Signature)))+
  geom_bar(stat="identity",col=NA,size=0.05)+
  geom_tile(aes(x=Sample,y=-40,col=factor(DorR),height=15),fill=NA,size=2,inherit.aes = F)+
  scale_fill_manual(values = sig_cols)+
  scale_color_manual(values = DorR_cols)+
  facet_grid(cols=vars(Pair_new),scales="free",space = "free")+
  scale_y_continuous(breaks=seq(0,3000,400),limits=c(-50,3000))+
  my_theme+
  theme(legend.key.size = unit(3,"mm"),
        axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  labs(fill="Signature",col="Donor or\nRecipient",y="Absolute mutation burden")
ggsave(filename = paste0(plots_dir,"plot.mutsigs.absolute.pdf"),plot = plot.mutsigs.absolute,width=7,height=2.5)

plot.SBS31.muts<-sample_metadata%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  filter(!is.na(N3_abs))%>%
  ggplot(aes(x=DorR,y=N3_abs,col=DorR))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1,height=0,alpha=0.25,size=0.1)+
  facet_grid(cols=vars(Pair_new))+
  scale_color_manual(values=DorR_cols)+
  labs(y="Mutations assigned to\nplatinum-associated signature",x="Donor (D) or Recipient (R)")+
  theme_bw()+
  my_theme+
  theme(legend.position = "none")
ggsave(filename = paste0(plots_dir,"plot.SBS31.muts.pdf"),plot = plot.SBS31.muts,width=5,height=2)

#Plot the proportions of samples with APOBEC activation by DorR and Pair
n_samples_df<-sample_metadata%>%
  dplyr::filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  dplyr::count(Pair_new,DorR,name="Total")

sample_metadata%>%
  dplyr::filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  group_by(DorR)%>%
  summarise(APOBEC_pos=sum(N4_abs>15),Total=n(),APOBEC_pos_prop=sum(N4_abs>15)/n())

plot.APOBEC.samples<-sample_metadata%>%
  dplyr::filter(Sample%in%unlist(all.trees.cc.nodups) & Sample!="Ancestral")%>%
  filter(!is.na(N4_abs) & N4_abs>15)%>%
  dplyr::count(Pair_new,DorR,name="n_APOBEC")%>%
  tidyr::complete(Pair_new,DorR,fill=list(n_APOBEC=0))%>%
  left_join(n_samples_df)%>%
  mutate(prop_APOBEC=n_APOBEC/Total)%>%
  ggplot(aes(x=Pair_new,y=prop_APOBEC,fill=DorR))+
  geom_bar(stat="identity",position="dodge",col="black")+
  scale_y_continuous(labels=scales::label_percent(accuracy = 1))+
  scale_fill_manual(values=DorR_cols)+
  labs(y="Percentage of samples with\nAPOBEC-associated signature",x="",fill="Donor (D) or\nRecipient (R)")+
  theme_bw()+
  my_theme+
  theme(legend.key.size = unit(3,"mm"))
ggsave(filename = paste0(plots_dir,"plot.APOBEC.samples.pdf"),plot = plot.APOBEC.samples,width=3.5,height=2)


plot.APOBEC.mutationburdens<-sample_metadata%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  dplyr::filter(Sample%in%unlist(all.trees.cc.nodups) & Sample!="Ancestral")%>%
  filter(!is.na(N4_abs) & N4_abs>15)%>%
  ggplot(aes(x=Pair_new,y=N4_abs,col=Pair_new,shape=DorR))+
  geom_jitter(width=0.1,height=0,alpha=0.6)+
  scale_y_log10(breaks=c(25,50,100,250,500,1000,2500,5000,10000))+
  scale_color_manual(values=Pair_cols)+
  labs(x="ID",y="Number of APOBEC-associated\nmutations per sample",col="Pair",shape="Donor (D) or\nRecipient (R)")+
  theme_bw()+
  my_theme+
  theme(legend.key.size = unit(2.2,"mm"))
ggsave(filename = paste0(plots_dir,"plot.APOBEC.mutationburdens.pdf"),plot = plot.APOBEC.mutationburdens,width=3.5,height=2)
