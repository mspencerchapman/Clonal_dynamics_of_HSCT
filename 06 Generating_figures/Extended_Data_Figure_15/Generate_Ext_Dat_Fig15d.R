#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr")
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

if(!require("dndscv", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("im3sanger/dndscv")
  library("dndscv",character.only=T,quietly = T, warn.conflicts = F)
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

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
details_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_list)) {assign(x,details_list[[x]])}

#========================================#
# Define custom functions for this analysis ####
#========================================#


get_DR_ids=function(tree){
  PD_IDs=unique(substr(tree$tip.label[-which(tree$tip.label=="Ancestral")],1,8))
  #Get number elements only
  PD_numbers=readr::parse_number(PD_IDs)
  names(PD_IDs)<-sapply(PD_numbers,function(n) ifelse(n%%2==0,"donor_ID","recip_ID"))
  return(PD_IDs)
}

#Wrapper function to apply function directly on the details matrix
dndscv_on_details=function(details,id="this_sample",outp=1,max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf,...) {
  if(!"sampleID"%in%colnames(details)) {
    details$sampleID=id
  }
  muts=details[,c("sampleID","Chrom","Pos","Ref","Alt")]
  colnames(muts)<-c("sampleID","chr","pos","ref","alt")
  geneindels = NULL #Fixes a weird bug in dndscv
  dndscvout=dndscv(muts,outp=outp,max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,max_coding_muts_per_sample = max_coding_muts_per_sample,...)
  return(dndscvout)
}

#Calculate the number of mutations inferred to be under selection
calculate_number_of_muts_under_selection=function(dndscvout) {
  n_nonsynonymous<-dndscvout$annotmuts%>%filter(impact!="Synonymous")%>%nrow()
  dNdS_values_with_CI=as.numeric(dndscvout$globaldnds["wall",c("mle","cilow","cihigh")])
  n_drivers=round(sapply(dNdS_values_with_CI,function(x) {n_nonsynonymous*(x-1)/x}))
  names(n_drivers)=c("mle","CI_low","CI_high")
  return(n_drivers)
}

#========================================#
# RUN dNdScV ####
#========================================#
# Run across combined mutations sets across all donors & all recipients

#Apply over separated DONOR and RECIPIENT mutant sets

#1. Separate out a 'donor only' set of details matrices
details_D<-dplyr::bind_rows(Map(f=function(details,Pair_ID,tree) {
  details$sampleID=Pair_ID
  donor_ID<-get_DR_ids(tree)['donor_ID']
  donor_nodes=tree$edge[,2][sapply(tree$edge[,2],function(node) any(grepl(donor_ID,getTips(tree,node))))] #Define branches where at least one descendant tip is donor-derived
  return(details%>%dplyr::filter(node%in%donor_nodes))
},Pair_ID=names(all.muts),details=all.muts,tree=all.trees))

#2. Separate out a 'recipient only' set of details matrices
details_R<-dplyr::bind_rows(Map(f=function(details,Pair_ID,tree) {
  details$sampleID=Pair_ID
  recip_ID<-get_DR_ids(tree)['recip_ID']
  recip_nodes=tree$edge[,2][sapply(tree$edge[,2],function(node) any(grepl(recip_ID,getTips(tree,node))))] #Define branches where at least one descendant tip is recipient-derived
  return(details%>%dplyr::filter(node%in%recip_nodes))
},Pair_ID=names(all.muts),details=all.muts,tree=all.trees))

# Apply the dndscv wrapper function over each
dndscvout_D<-dndscv_on_details(details = details_D,outp=3)
dndscvout_R<-dndscv_on_details(details = details_R,outp=3)

#How many drivers are implied in the two datasets? Can compare this to the number of known drivers.
calculate_number_of_muts_under_selection(dndscvout_D)
calculate_number_of_muts_under_selection(dndscvout_R)

#Combine the results into a list, and then bind into a dataframe for easy plotting
dndscvout_list=list("Donor_only"=dndscvout_D$globaldnds,
                    "Recipient_only"=dndscvout_R$globaldnds)

dndscv_res_df<-Map(df=dndscvout_list,pair=names(dndscvout_list),function(df,pair) {
  as.data.frame(df)%>%
    mutate(pair=pair,.before=1)
})%>%dplyr::bind_rows()

#For visualizing, expand the names of the categories
types=c("all mutations","missense","nonsense","splice site","truncating")
names(types)=c("wall","wmis","wnon","wspl","wtru")

## Generate Fig. 15d ----
p.D_vs_R<-dndscv_res_df%>%
  dplyr::filter(grepl("_only",pair))%>%
  mutate(name=types[name])%>%
  mutate(pair=gsub("_"," ",gsub("_combined","",pair)))%>%
  ggplot(aes(x=pair,y=mle,ymin=cilow,ymax=cihigh))+
  geom_point(size=0.7)+
  geom_errorbar(width=0.2,linewidth=0.3)+
  scale_y_continuous(breaks=seq(0.5,2,0.1))+
  facet_grid(rows=vars(name))+
  geom_hline(yintercept=1,col="red")+
  labs(x="",y="dN/dS ratio")+
  coord_flip()+
  theme_bw()+
  my_theme+
  theme(strip.text.y = element_text(angle=0))

ggsave(filename = paste0(plots_dir,"ExtDatFig15d.pdf"),p.D_vs_R,width =3.5,height=2.5)
