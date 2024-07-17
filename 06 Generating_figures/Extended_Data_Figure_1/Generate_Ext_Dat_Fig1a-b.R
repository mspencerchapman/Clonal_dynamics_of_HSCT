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
# Set the ggplot2 theme for plotting ####
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

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT" #Change this to the directory where you have cloned the github
source(paste0(root_dir,"/data/HSCT_functions.R"))
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")

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
details_list<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/details_lists.Rds"))

## Generate information regarding loss-of-Y in male samples from X and Y coverage data ----
LOY_files=list.files(path=paste0(root_dir,"/data/SV_and_CNA_data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

## Read in the spreadsheet listing other copy number changes ----
CN_change_df=read.csv(paste0(root_dir,"/data/SV_and_CNA_data/Copy_number_changes.csv"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_list)) {assign(x,details_list[[x]])}

#========================================#
# COVERAGE ANALYSIS ####
#========================================#

## Generate Extended Data Fig. 1a ----
Coverage_stats<-sample_metadata%>%
  filter(sample_status=="PASS")%>%
  group_by(Pair_new)%>%
  dplyr::summarise(mean_cov=mean(Coverage),median_cov=median(Coverage),sd_cov=sd(Coverage))

Coverage_histograms<-sample_metadata%>%
  filter(!is.na(Coverage))%>%
  filter(sample_status=="PASS")%>%
  ggplot(aes(x=Coverage))+
  geom_histogram(fill="lightblue",col="black",size=0.3)+
  geom_vline(data=Coverage_stats,aes(xintercept=mean_cov),col="red",linetype=2)+
  geom_text(data=Coverage_stats,aes(label=paste("µ =",round(mean_cov,digits=1),"x")),x=20,y=50,size=2,col="red")+
  #geom_vline(xintercept=4,col="black",linetype=1)+
  #geom_rect(xmin=0,xmax=4,ymin=0,ymax=80,fill="7a8baf#30")+
  theme_classic()+
  my_theme+
  facet_wrap(~Pair_new,nrow=2)+
  theme(strip.text.x=element_text(margin = unit(c(0.6,1,0.6,1),"mm")))+
  labs(x="Average genome-wide coverage",y="Count of colonies")

ggsave(filename=paste0(plots_dir,"ExtDatFig1a.pdf"),Coverage_histograms,width=6,height=2)

#========================================#
# COLONY OUTCOME SUMMARY ####
#========================================#
colony_outcomes<-c("PASS","Low coverage","Non-clonal","Different individual","Duplicate")
colony_outcome_cols=c("#33A02C","#A6CEE3","#1F78B4","#FB9A99","grey50")
names(colony_outcome_cols)<-colony_outcomes

## Generate Extended Data Fig. 1b ----
colony_outcome_summary_plot<-sample_metadata%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  mutate(sample_status=factor(sample_status,levels=rev(colony_outcomes)))%>%
  filter(sample_status!="Not Sequenced")%>%
  ggplot(aes(x=DorR,fill=sample_status))+
  geom_bar(col="black",size=0.1)+
  theme_bw()+
  my_theme+
  theme(strip.text.y=element_text(angle=0,size=6))+
  facet_grid(rows=vars(Pair_new))+
  scale_fill_manual(values=colony_outcome_cols)+
  coord_flip()+
  labs(x="Donor or Recipient",y="Number of colonies",fill="Sample outcome")

ggsave(filename = paste0(plots_dir,"ExtDatFig1b.pdf"),colony_outcome_summary_plot,width =4,height=3)