#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest","pheatmap")
bioconductor_packages=c("clusterProfiler")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
}

#========================================#
# Set the ggplot2 themes for plotting ####
#========================================#

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 7),
                axis.title = element_text(size=8),
                axis.line = element_line(linewidth = 0.4),
                axis.ticks = element_line(linewidth = 0.3),
                legend.text = element_text(size=6),
                legend.title = element_text(size=8),
                strip.text = element_text(size=7),
                strip.background = element_rect(fill="lightgray",linewidth = 0.4),
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
                     strip.background = element_rect(linewidth =0.25),
                     panel.grid.major = element_line(size=0.25),
                     panel.border = element_rect(linewidth=0.25))

#========================================#
# Set the root directory and read in the necessary files ####
#========================================#

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
source(paste0(root_dir,"/data/HSCT_functions.R"))
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")
vcf_header_path=paste0(root_dir,"/data/reference_files/vcfHeader.txt") #A dummy vcf header that can be used

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=Pair_metadata$Pair_new;names(new_pair_names)<-Pair_metadata$Pair

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

#Generate information regarding loss-of-Y in male samples from X and Y coverage data
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

#Read in mutational signature extraction data
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)


#========================================#
# MUTATIONAL SIGNATURE ANALYSIS USING HDP ####
#========================================#

mut_example_multi=readRDS(paste0(HDP_folder,"/HDP_multi_chain.Rdata"))
mutations=read.table(paste0(HDP_folder,"/trinuc_mut_mat.txt"))
key_table=read.table(paste0(HDP_folder,"/key_table.txt"))
sig_profiles=mut_mat_HDP_comp(mut_example_multi,plot=T)

# Altered plotting function to allow y axis to be free
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

## Generate Extended Data Fig. 4a ----
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


## Generate Extended Data Fig. 4b ----
new_sig_names=c("SBS2 (APOBEC)","SBS31 (Platinum)","SBS1-like","HSPC Signature","Unattributed")
names(new_sig_names)=rev(paste0("N",0:4))

plot.mutsigs.absolute<-sample_metadata%>%
  filter(!is.na(N0))%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  dplyr::select(Sample,DorR,Pair_new,N0_abs,N1_abs,N2_abs,N3_abs,N4_abs)%>%
  mutate(N4_abs=sapply(N4_abs,function(x) min(x,500)))%>% #To aid visualization limit total contribution of `N4 (APOBEC) to 500 mutations`- affects 4 samples
  mutate(N1_abs=sapply(N1_abs,function(x) min(x,2000)))%>% #To aid visualization limit total contribution of `N1 (HSPC sig) to 2000 mutations`- affects 4 samples
  gather(-Sample,-Pair_new,-DorR,key = "Signature",value="Number_of_mutations")%>%
  mutate(Signature=gsub("_abs","",Signature))%>%
  mutate(Signature=factor(new_sig_names[Signature],levels=new_sig_names))%>%
  ggplot(aes(x=Sample,y=Number_of_mutations,fill=factor(Signature)))+
  geom_bar(stat="identity",col=NA,size=0.05)+
  geom_tile(aes(x=Sample,y=-40,col=factor(DorR),height=15),fill=NA,size=2,inherit.aes = F)+
  scale_fill_manual(values = sig_cols)+
  scale_color_manual(values = DorR_cols)+
  guides(col=guide_legend(order = 2),fill=guide_legend(order=1))+
  facet_grid(cols=vars(Pair_new),scales="free",space = "free")+
  scale_y_continuous(breaks=seq(0,3000,400),limits=c(-50,3000))+
  my_theme+
  theme(legend.key.size = unit(3,"mm"),
        axis.text.x = element_blank(),
        axis.line.x =element_blank(),
        axis.ticks.x = element_blank())+
  labs(fill="Signature",col="",y="Absolute mutation burden")
ggsave(filename = paste0(plots_dir,"ExtDatFig4b.pdf"),plot = plot.mutsigs.absolute,width=7,height=2.5)

#========================================#
# ANALYSIS OF PLATINUM SIGNATURE ####
#========================================#

## Generate Extended Data Fig. 4c ----
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
ggsave(filename = paste0(plots_dir,"ExtDatFig4c.pdf"),plot = plot.SBS31.muts,width=5,height=2)


#========================================#
# ANALYSIS OF APOBEC SIGNATURE ####
#========================================#

#Plot the proportions of samples with APOBEC activation by DorR and Pair
n_samples_df<-sample_metadata%>%
  dplyr::filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  dplyr::count(Pair_new,DorR,name="Total")

sample_metadata%>%
  dplyr::filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  group_by(DorR)%>%
  summarise(APOBEC_pos=sum(N4_abs>15),Total=n(),APOBEC_pos_prop=sum(N4_abs>15)/n())

## Generate Extended Data Fig. 4d ----
plot.APOBEC.samples<-sample_metadata%>%
  dplyr::filter(Sample%in%unlist(all.trees.cc.nodups) & Sample!="Ancestral")%>%
  filter(!is.na(N4_abs) & N4_abs>15)%>%
  dplyr::count(Pair_new,DorR,name="n_APOBEC")%>%
  tidyr::complete(Pair_new,DorR,fill=list(n_APOBEC=0))%>%
  left_join(n_samples_df)%>%
  mutate(prop_APOBEC=n_APOBEC/Total)%>%
  ggplot(aes(x=Pair_new,y=prop_APOBEC,fill=DorR))+
  geom_bar(stat="identity",position="dodge",col="black",linewidth=0.3)+
  scale_y_continuous(labels=scales::label_percent(accuracy = 1))+
  scale_fill_manual(values=DorR_cols,labels=c("Donor","Recipient"))+
  labs(y="Percentage of samples with\nAPOBEC-associated signature",x="",fill="")+
  theme_bw()+
  my_theme+
  theme(legend.key.size = unit(3,"mm"),axis.text.x = element_text(angle=90))
ggsave(filename = paste0(plots_dir,"ExtDatFig4d.pdf"),plot = plot.APOBEC.samples,width=3.5,height=2.5)

## Generate Extended Data Fig. 4e ----
plot.APOBEC.mutationburdens<-sample_metadata%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  dplyr::filter(Sample%in%unlist(all.trees.cc.nodups) & Sample!="Ancestral")%>%
  filter(!is.na(N4_abs) & N4_abs>15)%>%
  ggplot(aes(x=Pair_new,y=N4_abs,col=Pair_new,shape=DorR))+
  geom_jitter(width=0.1,height=0,alpha=0.6)+
  scale_y_log10(breaks=c(25,50,100,250,500,1000,2500,5000,10000))+
  scale_color_manual(values=Pair_cols,guide="none")+
  labs(x="ID",y="Number of APOBEC-associated\nmutations per sample",col="Pair",shape="Donor (D) or\nRecipient (R)")+
  theme_bw()+
  my_theme+
  theme(legend.key.size = unit(2.2,"mm"))
ggsave(filename = paste0(plots_dir,"plot.APOBEC.mutationburdens.pdf"),plot = plot.APOBEC.mutationburdens,width=3.5,height=2)
