#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest")
bioconductor_packages=c()

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
trees_list<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/details_lists.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

#Generate information regarding loss-of-Y in male samples from X and Y coverage data
LOY_files=list.files(path=paste0(root_dir,"/data/SV_and_CNA_data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv(paste0(root_dir,"/data/SV_and_CNA_data/Copy_number_changes.csv"))

#Read in mutational signature extraction data
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)

#========================================#
# ANALYSIS OF COPY NUMBER ALTERRATIONS (CNAs) ####
#========================================#
#Plot of copy number changes by individual

## Generate Extended Data Fig. 7a ----
CNA_autosomal_plot<-Y_loss_df%>%
  dplyr::rename("Sample"=id)%>%
  mutate(Copy_number_change=ifelse(loss_of_Y=="YES","ChrY_loss",NA))%>%
  left_join(sample_metadata%>%dplyr::select(Sample,DorR,ID))%>%
  dplyr::select(Sample,Copy_number_change,DorR,"Pair"=ID)%>%
  bind_rows(CN_change_df)%>%
  filter(!is.na(Copy_number_change))%>%
  left_join(Pair_metadata,by="Pair")%>%
  filter(!Copy_number_change%in%c("ChrY_loss" ,"ChrX_del"))%>% #Can remove the LOY to get better resolution on the other changes
  ggplot(aes(x=DorR,y=1,fill=stringr::str_wrap(factor(Copy_number_change),width=22)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.15)+
  facet_grid(cols=vars(factor(Pair_new,levels=paste0("Pair_",1:10))),drop=F)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  my_theme+
  theme(legend.key.size = unit(3,"mm"))+
  labs(x="Donor (D) or Recipient (R)",y="Number of HSPCs",fill="")
ggsave(filename = paste0(plots_dir,"ExtDatFig7a.pdf"),CNA_autosomal_plot,device = "pdf",width=7,height = 2)

#========================================#
# ANALYSIS OF STRUCTURAL VARIANTS (SVs) ####
#========================================#

SV_df<-read_csv(paste0(root_dir,"/data/SVs_combined.csv"))%>%dplyr::rename("Sample"=sample)
SV_df$length<-sapply(1:nrow(SV_df),function(i) ifelse(SV_df$`SV type`[i]%in%c("DEL","INV"),SV_df$pos2_start[i]-SV_df$pos1_start[i],NA))

#Number of samples with SVs
SV_df%>%dplyr::filter(Sample%in%final_sample_list)%>%group_by(Sample)%>%summarise(n=n())%>%nrow()

SV_summary_df<-SV_df%>%
  dplyr::filter(Sample%in%final_sample_list)%>%
  group_by(Sample)%>%
  summarise(n=n(),DorR=DorR[1],types=paste(unique(`SV type`),collapse=","))

SV_summary_df%>%
  group_by(types)%>%
  summarise(n=n())

SV_summary_df%>%
  group_by(DorR)%>%
  summarise(n=n())
prop.test(c(7,16),n=c(1262,1562))

SV_plot<-SV_summary_df%>%
  left_join(sample_metadata%>%dplyr::select(Sample,ID))%>%
  left_join(Pair_metadata,by=c("ID"="Pair"))%>%
  ggplot(aes(x=DorR,y=1,fill=stringr::str_wrap(factor(types),width=22)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.15)+
  facet_grid(cols=vars(factor(Pair_new,levels=paste0("Pair_",1:10))),drop=F)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  my_theme+
  labs(x="Donor (D) or Recipient (R)",y="Samples",fill="Structural\nvariant type")
ggsave(filename = paste0(plots_dir,"ExtDatFig7b.pdf"),SV_plot,device = "pdf",width=5.5,height = 2)


#========================================#
# CALCULATE COMBINED RISK OF SV OR CNA ####
#========================================#

SV_or_CNA_samples<-CN_change_df%>%filter(Copy_number_change!="ChrX_del")%>%
  dplyr::bind_rows(SV_summary_df)%>%pull(Sample)

## Generate Extended Data Fig. 7c ----
SV_or_CNA_prop_plot<-sample_metadata%>%
  filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  mutate(SV_or_CNA=ifelse(Sample%in%SV_or_CNA_samples,"Abnormal","Normal"))%>%
  group_by(SV_or_CNA,DorR)%>%
  summarise(n=n())%>%
  pivot_wider(names_from="SV_or_CNA",values_from="n")%>%
  mutate(abnormal_prop=Abnormal/(Abnormal+Normal),
         CI_low=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[1])}),
         CI_high=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[2])}))%>%
  ggplot(aes(x=DorR,col=DorR,y=100*abnormal_prop,ymin=100*CI_low,ymax=100*CI_high))+
  geom_point()+
  geom_errorbar(width=0.2)+
  scale_color_manual(values=DorR_cols)+
  theme_bw()+
  my_theme+
  scale_y_continuous(limits=c(0,3))+
  theme(legend.position = "none")+
  labs(x="Donor or Recipient",y="Proportion of HSPCs with\nSV or autosomal CNA (%)")
ggsave(filename = paste0(plots_dir,"ExtDatFig7c.pdf"),SV_or_CNA_prop_plot,device = "pdf",width=1.5,height = 2)

prop.test(c(9,26),c(1252,1536))

#Plot of copy number changes by individual
## Generate Extended Data Fig. 7d ----
CNA_sex_plot<-Y_loss_df%>%
  dplyr::rename("Sample"=id)%>%
  mutate(Copy_number_change=ifelse(loss_of_Y=="YES","ChrY_loss",NA))%>%
  left_join(sample_metadata%>%dplyr::select(Sample,DorR,ID))%>%
  dplyr::select(Sample,Copy_number_change,DorR,"Pair"=ID)%>%
  bind_rows(CN_change_df)%>%
  filter(!is.na(Copy_number_change))%>%
  left_join(Pair_metadata,by="Pair")%>%
  filter(Copy_number_change%in%c("ChrY_loss","ChrX_del"))%>% #Can remove the LOY to get better resolution on the other changes
  ggplot(aes(x=DorR,y=1,fill=stringr::str_wrap(factor(Copy_number_change),width=22)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.1)+
  facet_grid(cols=vars(factor(Pair_new,levels=paste0("Pair_",1:10))),drop=F)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  my_theme+
  labs(x="Donor or Recipient",y="Number of HSPCs",fill="Copy number\nabnormality")
ggsave(filename = paste0(plots_dir,"ExtDatFig7d.pdf"),CNA_sex_plot,device = "pdf",width=7,height = 2)