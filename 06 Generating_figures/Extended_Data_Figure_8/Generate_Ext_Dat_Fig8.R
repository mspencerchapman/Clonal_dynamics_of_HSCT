###-----------------------TRANSPLANT BOTTLENECK SIMULATION TREES-----------------------
#This is a simulation of transplant bottlenecks of various sizes

#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ape","phangorn")

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

if(!require("rsimpop", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/rsimpop")
  library("rsimpop",character.only=T,quietly = T, warn.conflicts = F)
}

#========================================#
# Set the root directory and read in the necessary files ####
#========================================#

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
setwd(root_dir)
source(paste0(root_dir,"/data/HSCT_functions.R"))
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")
vcf_header_path=paste0(root_dir,"/data/reference_files/vcfHeader.txt") #A dummy vcf header that can be used

## Read in Pair metadata data frame----
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=Pair_metadata$Pair_new;names(new_pair_names)<-Pair_metadata$Pair

#========================================#
# Set parameters for the simulation ####
#========================================#

n_sampled=300 #The number of 'colonies' sequenced
bottleneck_sizes=c(1e2,1e3,1e4) #Size of the different bottlenecks
time_of_transplant=20
time_of_sampling=30

#========================================#
# RUN THE SIMULATIONS ####
#========================================#

#Early development - faster division rate for the first year
cfg=getDefaultConfig(target_pop_size  = 5e5,ndriver = 1,basefit = 0.2,rate = 0.1)
print(cfg)
sp=sim_pop(NULL,params=list(n_sim_days=365*1,b_stop_at_pop_size=0),cfg=cfg)
plot(sp)
fulltree1=get_tree_from_simpop(sp)

#Steady state - slower division rate from 1 year up time of transplant (decrease from 0.1 to 0.002/ day)
cfg$compartment$rate[2]<-0.002
sp2=sim_pop(fulltree1,params=list(n_sim_days=365*time_of_transplant,b_stop_at_pop_size=0),cfg=cfg)
sp2=combine_simpops(sp,sp2)
plot(sp2)
fulltree2=get_tree_from_simpop(sp2)

#Now simulate the bottleneck with 3 different bottleneck sizes
sim_trees<-lapply(bottleneck_sizes,function(n_engrafted){
  sampledtree1=get_subsampled_tree(fulltree2,n_engrafted)
  
  sp3=sim_pop(sampledtree1,params=list(n_sim_days=365*time_of_sampling,b_stop_at_pop_size=0),cfg=cfg)
  sp3=combine_simpops(sp2,sp3)
  plot(sp3)
  fulltree3=get_tree_from_simpop(sp3)
  sampledtree2=get_subsampled_tree(fulltree3,n_sampled)
  sampledtree2m=get_elapsed_time_tree(sampledtree2,mutrateperdivision=0.5,backgroundrate=16/365)
  tree=plot_tree(sampledtree2m,cex.label = 0)
  return(list(tree=sampledtree2,et_tree=sampledtree2m))
})

#========================================#
# VISUALIZE RESULTS ####
#========================================#

#Create ultrametric trees
sim_trees<-lapply(sim_trees,function(list){
  tree.ultra<-make.ultrametric.tree(list$et_tree)
  tree.ultra$edge.length[is.infinite(tree.ultra$edge.length)]<-1
  tree.ultra$edge.length<-tree.ultra$edge.length*mean(get_mut_burden(list$et_tree))
  list$tree.ultra<-tree.ultra
  return(list)
})

## Generate Extended Data Fig. 8 ----
pdf(file=paste0(plots_dir,"ExtDatFig8.pdf"),width=7,height=6)
par(mfrow=c(length(bottleneck_sizes),1))
temp=mapply(list=sim_trees,bottleneck_size=bottleneck_sizes,function(list,bottleneck_size) plot_tree(list$tree.ultra,cex.label=0,plot_axis=F,title = paste0("Engrafting HSCs = ",scales::comma(bottleneck_size))))
dev.off()

