#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(phytools)
library(rsimpop)
library(abc)
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","tidyr","ape","dichromat","abc","stringr","readr","phytools","data.table","pheatmap")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
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

## Read in Pair metadata data frame----
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=Pair_metadata$Pair_new;names(new_pair_names)<-Pair_metadata$Pair

## Define colour themes for the Pairs & DorR----
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

## Read in other data objects ----
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))

#========================================#
# IMPORT PARAMETERS AND SUMMARY STATISTICS FOR THE SIMULATIONS & PERFORM ABC ####
#========================================#

## Summary statistics and parameters from simulations are saved in this folder
## However, to see the simulations used to generate these statistics, review the script in /root/simulation_scripts/HSCT_simulation_rsimpop_farm_Tcell_differences.R
setwd(paste0(root_dir,"/data/ABC_simulation_results/Tcell_fraction_simulations/"))

Tcell_frac_pairs=paste("Pair",3:10,sep="_")

#Now extract the parameters/ summary statistics
all_stats<-lapply(Pair_metadata%>%pull(Pair),function(pair){
  cat(pair,sep="\n")
  
  sim.files=list.files(path=pair,pattern=".Rds",full.names = T)
  params_file=paste0("all_params_",pair,".Rds")
  sumstats_file=paste0("all_sumstats_",pair,".Rds")
  
  if(file.exists(params_file) & file.exists(sumstats_file)){
    cat("Reading in saved files",sep="\n\n")
    params<-readRDS(params_file)
    sumstats<-readRDS(sumstats_file)
  } else {
    cat("Importing individual results files and concatenating",sep="\n\n")
    sim.output<-lapply(sim.files,function(file) {if(which(sim.files==file)%%100 == 0) {print(which(sim.files==file))};readRDS(file)})
    
    params<-lapply(sim.output,function(res) unlist(res$params))%>%dplyr::bind_rows()%>%mutate(idx=1:length(sim.output))
    sumstats<-lapply(1:length(sim.output),function(i) {
      if(i%%100 == 0) {print(i)}
      df<-sim.output[[i]]$ratios%>%mutate(idx=i)
      return(df)
      })%>%dplyr::bind_rows()
    
    saveRDS(params,file=params_file)
    saveRDS(sumstats,file=sumstats_file)
  }
  
  df<-left_join(params,sumstats,by="idx")%>%mutate(Pair=pair)
  return(df)
})%>%bind_rows()

## Generate Extended Data Fig. 13c ----
Tcell_fraction_simulation_plot<-all_stats%>%
  mutate(Pair=factor(new_pair_names[Pair],levels=paste("Pair",1:10,sep="_")))%>%
  filter(Pair%in%paste("Pair",3:10,sep="_") & Myeloid>0.05)%>%
  ggplot(aes(x=T_cell_clone_duration,y=`T/M_ratio`))+
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method="lm",col="black",fullrange=T)+
  theme_bw()+
  my_theme+
  facet_wrap(~Pair,nrow=2)+
  geom_hline(yintercept=1,linetype=2)+
  annotate(geom="rect", xmin = 8, xmax = 15, ymin = 0, ymax = 1.2,alpha = .2)+
  scale_x_continuous(limits=c(0,20))+
  scale_y_continuous(limits=c(0,1.2),breaks=seq(0,2,0.2))+
  labs(x="Simulated lifespan of T cell clones (years)",
       y="Expanded clone T-cell fraction/\nExpanded clone myeloid fraction",
       fill="Count of simulations")

ggsave(filename = paste0(plots_dir,"ExtDatFig13c.pdf"),Tcell_fraction_simulation_plot,device = "pdf",width=7,height = 3)

