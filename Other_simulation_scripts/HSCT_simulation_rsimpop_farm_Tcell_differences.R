#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("optparse", "ape","dplyr","tidyr","ggplot2","stringr","readr","phytools")

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

#----------------------------------
# Specify options to parse from the command line-----------
#----------------------------------

option_list = list(
  make_option(c("-n", "--n_of_pair"), action="store", default=1, type='numeric', help="index of pair forsimulation"),
  make_option(c("-i", "--increase_late_pop"), action="store_true", default=FALSE, type='logical', help="option to increase the population in last 10 years"),
  make_option(c("-p", "--pop_in_final_10years"), action="store", default=5e5, type='numeric', help="population size in last 10 years"),
  make_option(c("-s", "--sim"), action="store", default=100, type='numeric', help="number of simulation per run of the script"),
  make_option(c("-o", "--output_dir"), action="store", default='.', type='character', help="output directory")
)

#Parse the inputted options
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))
print(opt)
j=opt$n
nsim=opt$s
output_dir=opt$o
increase_late_pop=opt$i
HSC_pop_size_final_10yrs=opt$p

if(!exists("opt")){
  j<-1
  nsim<-1
  increase_late_pop=F
  HSC_pop_size_final_10yrs=5e5
  output_dir="."
}

#----------------------------------
# Set working directories-----------
#----------------------------------

root_dir<-"/lustre/scratch126/casm/team154pc/ms56/Clonal_dynamics_of_HSCT/" #Change this to wherever you have cloned the github directory

tree_folder=paste0(root_dir,"/data/trees_no_dups/")
posterior_parameters_file=paste0(root_dir,"/data/reference_files/driver_parameter_posterior_sample.txt")
HSC_pop_posteriors_file<-paste0(root_dir,"/data/reference_files/HSC_population_posterior_sample.txt")
R_function_files=paste0(root_dir,"/data/HSCT_functions.R")
source(R_function_files)
setwd(root_dir)
visualize=F

###-----------------SPECIFIC FUNCTIONS FOR THIS SCRIPT-----------------
get_ltt = function(tree,time_points) {
  nodeheights <- nodeHeights(tree)
  ltt_tree = sapply(time_points, function(x) {
    sum(nodeheights[,1] < x & !nodeheights[,2] < x)
  })
  return(ltt_tree)
}

get_cutoff_branches=function(tree,cut_off) {
  heights=nodeHeights(tree)
  cutoff_branches=tree$edge[,2][heights[,1]<cut_off & heights[,2]>=cut_off]
  return(cutoff_branches)
}

get_coalescences = function(ltt) {
  coals=sapply(2:length(ltt), function(i) {return(ltt[i]-ltt[i-1])})
  return(coals)
}

coals_within_time_window_per_clone=function(tree,define_clone_height=5,time_points) {
  clone_nodes=get_cutoff_branches(tree,cut_off = define_clone_height)
  nodes_within_time_points=tree$edge[,1][nodeHeights(tree)[,1]>time_points[1] & nodeHeights(tree)[,1]<time_points[2]]
  
  coals_within_time_window_by_clone=sapply(clone_nodes,function(node) {
    clone_daughters<-get_all_node_children(node,tree)
    coals_within_time_points=intersect(clone_daughters,nodes_within_time_points)
    return(length(coals_within_time_points))
  })
  return(coals_within_time_window_by_clone)
}

get_DR_ids=function(tree){
  tree=ape::drop.tip(tree,"Ancestral")
  PD_IDs=unique(substr(tree$tip.label,1,8))
  #Get number elements only
  PD_numbers=readr::parse_number(PD_IDs)
  names(PD_IDs)<-sapply(PD_numbers,function(n) ifelse(n%%2==0,"donor_ID","recip_ID"))
  return(PD_IDs)
}

get_expanded_clade_nodes=function(tree,height_cut_off=100,min_clonal_fraction=0.02,min_samples=1){
  nodeheights=nodeHeights(tree)
  
  #This pulls out nodes that fulfill on the criteria: branches cross the cut-off & contain the minimum proportion of samples
  nodes=tree$edge[,2][nodeheights[,1] < height_cut_off &
                        !nodeheights[,2] < height_cut_off &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))/length(tree$tip.label)})>min_clonal_fraction &
                        sapply(tree$edge[,2],function(node) {length(getTips(tree,node))})>=min_samples]
  df=data.frame(nodes=nodes,n_samples=sapply(nodes,function(node) {length(getTips(tree,node))}),MRCA_time=sapply(nodes,function(node) {nodeheight(tree,node)}),clonal_fraction=sapply(nodes,function(node) {length(getTips(tree,node))/length(tree$tip.label)}))
  return(df)
}

#Bulk totals for other expansions
find_latest_acquisition_node=function(tree,pos_samples){
  #Get list of ancestral nodes for all samples
  ancestral_nodes_list=lapply(pos_samples,function(Sample) {
    get_ancestral_nodes(node = which(tree$tip.label==Sample),edge=tree$edge)
    })
  #Find nodes that are ancestral to all the samples
  common_nodes=Reduce(intersect,ancestral_nodes_list)
  #Which of these is the most recent (i.e. has the maximum node height)
  nodeheights<-sapply(common_nodes,function(node) nodeheight(tree = tree,node = node))
  MRCA_node<-common_nodes[which.max(nodeheights)]
  return(MRCA_node)
}

#Function to extract clonal fractions from the total population based on a subsampled tree
#Need to make sure that the tip labels match
extract_bulk_cell_fractions=function(nodes,sub_pop,full_pop,states) {
  #Get a vector of the 'states' of the all tree tips (i.e. which compartment)
  full_pop_tip_states=full_pop$state[full_pop$edge[,2]<=length(full_pop$tip.label)]
  
  #Now go through each node from the sub_pop tree & find the full population fraction of the MRCA of that clade
  res<-lapply(nodes,function(node) {
    cat(node,sep = "\n")
    raw_tips<-stringr::str_split(getTips(sub_pop,node),pattern = "_",simplify = T)[,1] #Get the sample names from the clade
    bulk_node=find_latest_acquisition_node(tree=full_pop,pos_samples = raw_tips) #Get the MRCA of those samples from the full population
    clade_states=full_pop_tip_states[which(full_pop$tip.label%in%getTips(full_pop,node=bulk_node))] #Get the states of enclosed samples
    bulk_fracs=as.data.frame(table(clade_states)[as.character(states)]/table(full_pop_tip_states)[as.character(states)])%>%
      dplyr::mutate(clade_states=states)%>%
      tidyr::replace_na(replace = list(Freq=0))%>%
      pivot_wider(names_from = "clade_states",values_from="Freq",names_prefix="state_")
  })%>%dplyr::bind_rows()%>%
    mutate(node=nodes,.before=1)
  
  return(res)
}

#Do updated 'get_subsampled_tree' function that maintains original tip labels (needed to extract bulk clonal fractions)
get_subsampled_tree2=function (tree, N, tips = tree$edge[c(which(tree$state == 0 & 
                                                                   tree$edge[, 2] <= length(tree$tip.label)), sample(which(tree$state != 
                                                                                                                             0 & tree$edge[, 2] <= length(tree$tip.label)), N)), 2]) 
{
  N = length(tips)
  tip.labels<-tree$tip.label[sort(tips)]
  
  #Create the tree with the tips kept in the same order
  tree_same_order<-keep.tip(tree,tip=tip.labels)
  tmp = rsimpop:::C_subsample_pop(tree, sort(tips))
  tmp$tip.label = sprintf("s%d", 1:N)
  class(tmp) = c("simpop", "phylo")
  
  #Now reorder the tips
  node_trans<-all.equal(tree_same_order,tmp,use.tip.label=F,index.return = T)
  node_trans<-node_trans[node_trans[,2]<=length(tree_same_order$tip.label),]
  new_tips<-tip.labels[node_trans[,2]]
  tmp$tip.label<-new_tips
  
  #tmp$tip.label = tip.labels
  
  checkValidPhylo(tmp)
  tmp$is_combined = tree$is_combined
  tmp
}


###-----------------Pair metadata & trees-----------------

exp_nos<-c(11,13,21,24,25,28,31,38,40,41)
Pair_metadata=data.frame(Pair=paste0("Pair",exp_nos),
                         Pair_new=factor(c("Pair_9","Pair_7","Pair_5","Pair_1","Pair_4","Pair_10","Pair_6","Pair_8","Pair_3","Pair_2"),levels=paste0("Pair_",1:10)),
                         Age=c(74.8,65.5,64.5,34.2,58.4,79.9,65.2,65.8, 51.9,42.4),
                         Age_at_transplant=c(66,36,50,18,47,63,43,35,35,30),
                         MNC_dose=c(2.66,4.17,16.28,NA,10.94,13.9,NA,4.05,2.48,15.01),
                         CD34_dose=c(1.56,NA,7.1,7.9,8.97,2.4,NA,NA,NA,4.51),
                         stem_cell_source=c("BM","BM","PBSC","PBSC","PBSC","PBSC","BM","BM","BM","PBSC"),
                         conditioning=c("MAC","MAC","RIC","MAC","RIC","RIC","MAC","MAC","MAC","MAC"))

#Import the trees (the ones with duplicates removed)
tree_paths=list.files(tree_folder,pattern="vaf_post_mix_post_dup.tree",full.names = T)
all.trees<-lapply(exp_nos,function(exp_no){
  pair_tree_paths=grep(paste0("Pair",exp_no),tree_paths,value = T)
  if(length(pair_tree_paths)==1) {
    tree<-read.tree(pair_tree_paths)
  } else if(length(pair_tree_paths)>1){
    a_j_path<-grep("a_j",pair_tree_paths,value=T)
    tree<-read.tree(a_j_path)
  } else {
    stop(cat(paste("No tree for Pair",exp_no),sep="\n"))
  }
  return(tree)
  })
names(all.trees)<-paste0("Pair",exp_nos)

#Read in the posteriors from E. Mitchell et al, 2022 to use as the prior
param_posterior<-read.delim(posterior_parameters_file)
HSC_pop_posteriors<-read.delim(HSC_pop_posteriors_file,header = F)

###-----------------Set up parameters for this simulation-----------------

#Set the HSCT parameters
PairID=paste0("Pair",j)
age_of_donor_at_HSCT=Pair_metadata$Age_at_transplant[Pair_metadata$Pair==PairID]
age_of_donor_at_sampling=Pair_metadata$Age[Pair_metadata$Pair==PairID]
T_cell_time=age_of_donor_at_sampling-10

if(T_cell_time<age_of_donor_at_HSCT) {
  stop(return("T cells are older than time of HSCT"))
}

tree.data<-all.trees[[PairID]]
DR_ids<-get_DR_ids(tree.data)

no_of_WGS_recip=sum(grepl(DR_ids['recip_ID'],tree.data$tip.label))
no_of_WGS_donor=sum(grepl(DR_ids['donor_ID'],tree.data$tip.label))

###-----------------SIMULATION APPROACH-----------------
#Grown tree to age of transplant (driver process simulation)
#Create new cell compartment of transplanted (recipient) cells
#Grow each compartment independently to HSC population size
system(paste("mkdir -p",paste0(output_dir,"/tree_files")))
cat(paste("Running",nsim,"simulations for",PairID),sep="\n")

for(sim in 1:nsim) {
  cat(paste("Running simulation number",sim),sep="\n")
  
  uid=ids::random_id()
  
  #Set the fixed HSCT parameters
  age_10yrs_pre_sampling=max(age_of_donor_at_sampling-10,age_of_donor_at_HSCT)
  HSC_symmetric_division_rate=1/(2*365)
  
  #Set the variable parameters
  HSC_pop_size=round(sample(HSC_pop_posteriors$V1,1))
  log10_HSCT_bottleneck=runif(1,min=2.7,max=4.7)
  HSCT_bottleneck=round(10^log10_HSCT_bottleneck)
  
  #Driver mutation parameters distribution parameters (from posterior of E. Mitchell et al, 2022)
  param_idx=sample(1:nrow(param_posterior),1)
  number_drivers_per_year = param_posterior$number_drivers_per_year[param_idx]
  gamma_shape = param_posterior$gamma_shape[param_idx]
  gamma_rate = param_posterior$gamma_rate[param_idx]
  fitness_threshold = 0.05
  dpcpd=min(1e-5,number_drivers_per_year/(365*HSC_pop_size)) #Number of drivers per cell per year cannot be more than 1e-5 in RSimpop
  
  ##Function to generate gamma distribution based fitness
  genGammaFitness=function(fitness_threshold,shape,rate){
    function() rtrunc(n=1,a=fitness_threshold, b=Inf,"gamma",shape=shape,rate=rate)
  }
  fitnessGammaFn=genGammaFitness(fitness_threshold=fitness_threshold,shape = gamma_shape, rate=gamma_rate)
  
  ##Function to generate exp distribution based fitness - this isn't actually used currently in this simulation.
  genExpFitness=function(fitness_threshold,rate){
    function() rtrunc(n=1,a=fitness_threshold, b=Inf,"exp",rate=rate)
  }
  fitnessExpFn=genExpFitness(fitness_threshold=0.08,rate=40)
  
  ##Put parameters into a single list (to save later)
  params=list(age_of_donor_at_HSCT=age_of_donor_at_HSCT,
              age_of_donor_at_sampling=age_of_donor_at_sampling,
              no_of_WGS_recip=no_of_WGS_recip,
              no_of_WGS_donor=no_of_WGS_donor,
              HSC_symmetric_division_rate=HSC_symmetric_division_rate,
              HSC_pop_size=HSC_pop_size,
              HSCT_bottleneck=HSCT_bottleneck,
              number_drivers_per_year=number_drivers_per_year,
              gamma_shape=gamma_shape,
              gamma_rate=gamma_rate,
              fitness_threshold=fitness_threshold,
              dpcpd=dpcpd)
  
  print(params)
  
  ###START ACTUAL SIMULATION
  #Grow population to the age of HSCT
  cat("Starting simulation",sep="\n")
  cat("Growing population to time of HSCT",sep="\n")
  dps2=run_driver_process_sim(initial_division_rate=0.1,
                              final_division_rate = HSC_symmetric_division_rate,
                              target_pop_size = HSC_pop_size,
                              nyears = age_of_donor_at_HSCT,
                              fitnessGen=fitnessGammaFn,
                              drivers_per_cell_per_day = dpcpd)
  
  #Select number of cells to be transplanted into the recipient as a new compartment "Transplanted"
  #These can then grow separately to the HSC population size
  cfg=dps2$cfg
  cfg=addCellCompartment(cfg,population = HSC_pop_size,rate=cfg$compartment$rate[2],ndriver=1,descr="Transplanted",basefit = 0.1)
  
  cat("Specifying transplanted cells",sep="\n")
  post_HSCT<-addDifferentiationEvents(tree=dps2,cfg = cfg,2,nEvent = HSCT_bottleneck)
  post_HSCT$cfg$info=rbind(post_HSCT$cfg$info,post_HSCT$cfg$info %>% filter(val==1 & id>0) %>% mutate(val=2))
  
  #Grow & age donor and recipient compartments (driver process simulation)
  cat("Growing and ageing donor & recipient compartments until time of T-cell origin",sep="\n")
  pop_Tcell_origin<-run_driver_process_sim(simpop=post_HSCT,
                                    target_pop_size = HSC_pop_size,
                                    final_division_rate = HSC_symmetric_division_rate,
                                    nyears = T_cell_time,
                                    fitnessGen=fitnessGammaFn,
                                    #fitnessGen=fitnessExpFn,
                                    drivers_per_cell_per_day = dpcpd)
  
  #Select cells for WGS: need to have the right number of donor & recipient cells, as well as the outgroup
  donorT_WGS=sample(which(pop_Tcell_origin$state[which(pop_Tcell_origin$edge[,2]%in%1:length(pop_Tcell_origin$tip.label))]==1),size = no_of_WGS_donor,replace = F)
  recipT_WGS=sample(which(pop_Tcell_origin$state[which(pop_Tcell_origin$edge[,2]%in%1:length(pop_Tcell_origin$tip.label))]==2),size = no_of_WGS_recip,replace = F)
  outgroupT=which(pop_Tcell_origin$state[which(pop_Tcell_origin$edge[,2]%in%1:length(pop_Tcell_origin$tip.label))]==0)
  
  #EXTRACT THE SEQUENCED TREE
  DR_treeT<-get_subsampled_tree2(pop_Tcell_origin,tips = c(outgroupT,donorT_WGS,recipT_WGS))
  
  #Get a vector of the "states" of the tips (i.e. whether they are from the donor or recipient)
  tip_states=DR_treeT$state[DR_treeT$edge[,2]%in%1:length(DR_treeT$tip.label)]
  
  #Use this to relabel the tips to include either 'donor' or 'recipient' in their labels
  DR_treeT$tip.label[tip_states==1]<-paste(DR_treeT$tip.label[tip_states==1],"donor",sep="_")
  DR_treeT$tip.label[tip_states==2]<-paste(DR_treeT$tip.label[tip_states==2],"recip",sep="_")
  
  DR_treeT_m<-get_elapsed_time_tree(DR_treeT,mutrateperdivision = 1.8,backgroundrate = 15/365)
  
  #Convert this into an ultrametric 'time tree' for assessing coalescences
  DR_treeT_m_ultra<-make.ultrametric.tree(DR_treeT_m)
  DR_treeT_m_ultra$edge.length[1]<-0 #Set length of outgroup branch to 0
  DR_treeT_m_ultra$edge.length<-DR_treeT_m_ultra$edge.length*age_of_donor_at_sampling #This assumes linear accumulation through life
  DR_treeT_m_ultra$coords<-NULL
  
  #Visualize the combined tree
  if(visualize){
    DR_treeT_m_ultra=plot_tree(DR_treeT_m_ultra,cex.label=0)
    temp=add_annotation(DR_treeT_m_ultra,
                        annot_function=plot_sharing_info,
                        donor_ID="donor",
                        recip_ID="recip",
                        sharing_cols=c("black", "#11a0aa80", "#c8256580")
    )
    hm=matrix(c("white",c("#11a0aa80", "#c8256580")[tip_states]),nrow=1,ncol=length(DR_treeT_m_ultra$tip.label))
    colnames(hm)<-DR_treeT_m_ultra$tip.label;rownames(hm)<-"D or R"
    add_heatmap(DR_treeT_m_ultra,heatmap = hm,heatvals = c("#11a0aa80", "#c8256580"))
    rect(xleft = 0,xright = length(DR_treeT_m_ultra$tip.label),ybottom = (age_of_donor_at_sampling-age_of_donor_at_HSCT-5),ytop = (age_of_donor_at_sampling-age_of_donor_at_HSCT+5),col="#D3D3D380",border = NA)
  }
  
  cat("Growing and ageing donor & recipient compartments until time of sampling",sep="\n")
  pop_final<-run_driver_process_sim(simpop=pop_Tcell_origin,
                                    target_pop_size = HSC_pop_size,
                                    final_division_rate = HSC_symmetric_division_rate,
                                    nyears = age_of_donor_at_sampling,
                                    fitnessGen=fitnessGammaFn,
                                    #fitnessGen=fitnessExpFn,
                                    drivers_per_cell_per_day = dpcpd)
  
  #Select cells for WGS: need to have the right number of donor & recipient cells, as well as the outgroup
  donor_WGS=sample(which(pop_final$state[which(pop_final$edge[,2]%in%1:length(pop_final$tip.label))]==1),size = no_of_WGS_donor,replace = F)
  recip_WGS=sample(which(pop_final$state[which(pop_final$edge[,2]%in%1:length(pop_final$tip.label))]==2),size = no_of_WGS_recip,replace = F)
  outgroup=which(pop_final$state[which(pop_final$edge[,2]%in%1:length(pop_final$tip.label))]==0)
  
  #EXTRACT THE SEQUENCED TREE
  DR_tree<-get_subsampled_tree2(pop_final,tips = c(outgroup,donor_WGS,recip_WGS))
  
  #Get a vector of the "states" of the tips (i.e. whether they are from the donor or recipient)
  tip_states=DR_tree$state[DR_tree$edge[,2]%in%1:length(DR_tree$tip.label)]
  
  #Use this to relabel the tips to include either 'donor' or 'recipient' in their labels
  DR_tree$tip.label[tip_states==1]<-paste(DR_tree$tip.label[tip_states==1],"donor",sep="_")
  DR_tree$tip.label[tip_states==2]<-paste(DR_tree$tip.label[tip_states==2],"recip",sep="_")
  
  if(visualize) {
    plot_tree_events(get_elapsed_time_tree(DR_tree),cex.label=0)
  }
  
  #Convert this into a 'mutation tree' using an approximate division/ time associated mutation rate - the precise values
  #are not important as the tree is made ultrametric for assessing coalescences
  DR_tree_m<-get_elapsed_time_tree(DR_tree,mutrateperdivision = 1.8,backgroundrate = 15/365)
  
  #Convert this into an ultrametric 'time tree' for assessing coalescences
  DR_tree_m_ultra<-make.ultrametric.tree(DR_tree_m)
  DR_tree_m_ultra$edge.length[1]<-0 #Set length of outgroup branch to 0
  DR_tree_m_ultra$edge.length<-DR_tree_m_ultra$edge.length*age_of_donor_at_sampling #This assumes linear accumulation through life
  DR_tree_m_ultra$coords<-NULL
  
  #Visualize the combined tree
  if(visualize){
    DR_tree_m_ultra=plot_tree(DR_tree_m_ultra,cex.label=0)
    temp=add_annotation(DR_tree_m_ultra,
                        annot_function=plot_sharing_info,
                        donor_ID="donor",
                        recip_ID="recip",
                        sharing_cols=c("black", "#11a0aa80", "#c8256580")
    )
    hm=matrix(c("white",c("#11a0aa80", "#c8256580")[tip_states]),nrow=1,ncol=length(DR_tree_m$tip.label))
    colnames(hm)<-DR_tree_m$tip.label;rownames(hm)<-"D or R"
    add_heatmap(DR_tree_m_ultra,heatmap = hm,heatvals = c("#11a0aa80", "#c8256580"))
    rect(xleft = 0,xright = length(DR_tree_m_ultra$tip.label),ybottom = (age_of_donor_at_sampling-age_of_donor_at_HSCT-5),ytop = (age_of_donor_at_sampling-age_of_donor_at_HSCT+5),col="#D3D3D380",border = NA)
  }
  
  #Extract the separated donor & recipient trees
  D_tree_m<-drop.tip(DR_tree_m,tip = grep("recip",DR_tree_m$tip.label,value=T));D_tree_m$coords<-NULL
  R_tree_m<-drop.tip(DR_tree_m,tip = grep("donor",DR_tree_m$tip.label,value=T));R_tree_m$coords<-NULL
  D_tree_m_ultra<-drop.tip(DR_tree_m_ultra,tip = grep("recip",DR_tree_m_ultra$tip.label,value=T));D_tree_m_ultra$coords<-NULL
  R_tree_m_ultra<-drop.tip(DR_tree_m_ultra,tip = grep("donor",DR_tree_m_ultra$tip.label,value=T));R_tree_m_ultra$coords<-NULL
  
  #Visualize the separated donor/ recipient trees
  if(visualize){
    par(mfrow=c(2,1))
    zz=plot_tree(D_tree_m_ultra,cex.label=0,default_edge_color = "#11a0aa80")
    rect(xleft = 0,xright = length(D_tree_m_ultra$tip.label),ybottom = (age_of_donor_at_sampling-age_of_donor_at_HSCT-5),ytop = (age_of_donor_at_sampling-age_of_donor_at_HSCT+5),col="#D3D3D380",border = NA)
    yy=plot_tree(R_tree_m_ultra,cex.label=0,default_edge_color = "#c8256580")
    rect(xleft = 0,xright = length(R_tree_m_ultra$tip.label),ybottom = (age_of_donor_at_sampling-age_of_donor_at_HSCT-5),ytop = (age_of_donor_at_sampling-age_of_donor_at_HSCT+5),col="#D3D3D380",border = NA)
  }
  
  ##EXTRACT SUMMARY STATISTICS
  #(1) 3 LARGEST CLADES (Donor & Recipient)
  largest_clades_D<-get_expanded_clade_nodes(D_tree_m,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    pull(n_samples)%>%sort(decreasing = T)%>%.[1:3]
  largest_clades_R<-get_expanded_clade_nodes(R_tree_m,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    pull(n_samples)%>%sort(decreasing = T)%>%.[1:3]
  
  #(2) Number of singletons (Donor & Recipient)
  n_singletons_D<-get_expanded_clade_nodes(D_tree_m,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    dplyr::filter(n_samples==1)%>%nrow(.)
  n_singletons_R<-get_expanded_clade_nodes(R_tree_m,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
    dplyr::filter(n_samples==1)%>%nrow(.)
  
  #(3) Peri-HSCT LTT/ coalescences (Combined, Donor, Recipient)
  peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
  DR_ltt_peri<-get_ltt(tree=DR_tree_m_ultra,time_points=peri_HSCT_time_points)
  D_ltt_peri<-get_ltt(tree=D_tree_m_ultra,time_points=peri_HSCT_time_points)
  R_ltt_peri<-get_ltt(tree=R_tree_m_ultra,time_points=peri_HSCT_time_points)
  DR_coals_peri<-get_coalescences(DR_ltt_peri)
  D_coals_peri<-get_coalescences(D_ltt_peri)
  R_coals_peri<-get_coalescences(R_ltt_peri)
  
  #(4) Pre-HSCT LTT/ coalescences (Combined, Donor, Recipient)
  pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
  DR_ltt_pre<-get_ltt(tree=DR_tree_m_ultra,time_points=pre_HSCT_time_points)
  D_ltt_pre<-get_ltt(tree=D_tree_m_ultra,time_points=pre_HSCT_time_points)
  R_ltt_pre<-get_ltt(tree=R_tree_m_ultra,time_points=pre_HSCT_time_points)
  DR_coals_pre<-get_coalescences(DR_ltt_pre)
  D_coals_pre<-get_coalescences(D_ltt_pre)
  R_coals_pre<-get_coalescences(R_ltt_pre)
  
  #(5) Maximum peri-transplant coalescences within single recipient clade (this is a new stat)
  R_max_coals_in_single_clade=max(coals_within_time_window_per_clone(tree=R_tree_m_ultra,define_clone_height = 5,time_points = peri_HSCT_time_points))
  
  #(6) Bulk totals for expanded populations to look at clonal shifts
  # These are only a possible metric when there are expansions (defined in the same way as for the data)
  expanded_df=get_expanded_clade_nodes(DR_tree_m,min_clonal_fraction=0.001,height_cut_off=50,min_samples=3)
  
  if(nrow(expanded_df)>0){
    expanded_df<-expanded_df%>%
      mutate(D_count=sapply(nodes,function(node) {sum(grepl("donor",getTips(DR_tree_m,node)))}),
             R_count=sapply(nodes,function(node) {sum(grepl("recip",getTips(DR_tree_m,node)))}))%>%
      mutate(D_phylo=D_count/no_of_WGS_donor,
             R_phylo=R_count/no_of_WGS_recip)%>%
      left_join(DR_tree$events%>%filter(driverid>0)%>%dplyr::select(node,driverid),by=c("nodes"="node"))
    
    for(i in 1:nrow(expanded_df)) {
      if(is.na(expanded_df$driverid[i])) {
        ancestrals<-get_ancestral_nodes(expanded_df$nodes[i],DR_tree_m$edge)
        children<-get_all_node_children(expanded_df$nodes[i],DR_tree_m)
        clone_drivers<-DR_tree_m$events%>%filter(driverid>0 & node%in%c(ancestrals,children))%>%pull(driverid)
        expanded_df$driverid[i]<-clone_drivers[1]
      }
    }
    
    D_total=sum(pop_final$cfg$info$population[pop_final$cfg$info$val==1])
    R_total=sum(pop_final$cfg$info$population[pop_final$cfg$info$val==2])
    expanded_df<-expanded_df%>%
      left_join(pop_final$cfg$info%>%filter(val==1)%>%dplyr::select(id,population),by=c("driverid"="id"))%>%
      mutate(D_bulk_by_driver=population/D_total)%>%
      dplyr::select(-population)%>%
      left_join(pop_final$cfg$info%>%filter(val==2)%>%dplyr::select(id,population),by=c("driverid"="id"))%>%
      mutate(R_bulk_by_driver=population/R_total)%>%
      dplyr::select(-population)
    
    DTcell_total=sum(pop_Tcell_origin$cfg$info$population[pop_Tcell_origin$cfg$info$val==1])
    RTcell_total=sum(pop_Tcell_origin$cfg$info$population[pop_Tcell_origin$cfg$info$val==2])
    expanded_df<-expanded_df%>%
      left_join(pop_Tcell_origin$cfg$info%>%filter(val==1)%>%dplyr::select(id,population),by=c("driverid"="id"))%>%
      mutate(D_Tbulk_by_driver=population/DTcell_total)%>%
      dplyr::select(-population)%>%
      left_join(pop_Tcell_origin$cfg$info%>%filter(val==2)%>%dplyr::select(id,population),by=c("driverid"="id"))%>%
      mutate(R_Tbulk_by_driver=population/RTcell_total)%>%
      dplyr::select(-population)
    
    cell_fracs_plot<-expanded_df%>%dplyr::select(driverid,D_bulk_by_driver,R_bulk_by_driver,D_Tbulk_by_driver,R_Tbulk_by_driver)%>%
      filter(!is.na(driverid))%>%
      gather(-driverid,key="sample",value="clonal_frac")%>%
      mutate(DorR=substr(sample,1,1))%>%
      mutate(cell_type=ifelse(substr(sample,3,3)=="T","T-cells","Myeloid"))%>%
      ggplot(aes(x=cell_type,y=clonal_frac,fill=factor(driverid)))+
      geom_bar(stat = "identity",position="stack")+
      facet_grid(~DorR)+
      theme_bw()
    
    ratios<-expanded_df%>%dplyr::select(driverid,D_bulk_by_driver,R_bulk_by_driver,D_Tbulk_by_driver,R_Tbulk_by_driver)%>%
      filter(!is.na(driverid))%>%
      gather(-driverid,key="sample",value="clonal_frac")%>%
      mutate(DorR=substr(sample,1,1))%>%
      mutate(cell_type=ifelse(substr(sample,3,3)=="T","T-cells","Myeloid"))%>%
      group_by(DorR,cell_type)%>%
      summarise(expanded_frac=sum(clonal_frac,na.rm = T))%>%
      pivot_wider(id_cols="DorR",names_from="cell_type",values_from="expanded_frac")%>%
      mutate(`T/M_ratio`=`T-cells`/Myeloid)
  }
  
  ##--COMBINE SUMSTATS INTO SINGLE VECTOR--

  
  #Save the key model parameters/ outputs - including the tree & expanded clone df that could be used to get different summary stats
  res=list(PairID=PairID,
           uid=uid,
           params=params,
           plot=cell_fracs_plot,
           ratios=ratios
           )
  
  cat("Saving simulation output",sep="\n\n")
  saveRDS(res,file=paste0(output_dir,"/","res_",PairID,"_",uid,".Rds"))
}
