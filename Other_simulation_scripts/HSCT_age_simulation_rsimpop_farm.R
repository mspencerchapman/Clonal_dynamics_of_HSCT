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

#Can set options manually here if running manually as a test run
if(!exists("opt")){
  j<-1
  nsim<-1
  increase_late_pop=F
  HSC_pop_size_final_10yrs=5e5
  output_dir="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Age/output" #Change this to wherever you would like to put the model output
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

###-----------------SPECIFIC BESPOKE FUNCTIONS FOR THIS SCRIPT-----------------
get_ltt = function(tree,time_points) {
  nodeheights <- nodeHeights(tree)
  ltt_tree = sapply(time_points, function(x) {
    sum(nodeheights[,1] < x & !nodeheights[,2] < x)
  })
  return(ltt_tree)
}

get_coalescences = function(ltt) {
  coals=sapply(2:length(ltt), function(i) {return(ltt[i]-ltt[i-1])})
  return(coals)
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
tree_folder=paste0(root_dir,"data/trees_no_dups/")
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


colony_numbers_df<-Map(tree=all.trees,pair=names(all.trees),function(tree,pair) {
  DR_ids<-get_DR_ids(tree)
  
  no_of_WGS_recip=sum(grepl(DR_ids['recip_ID'],tree$tip.label))
  no_of_WGS_donor=sum(grepl(DR_ids['donor_ID'],tree$tip.label))
  
  return(data.frame(PairID=pair,Donor=no_of_WGS_donor,Recip=no_of_WGS_recip))
})%>%dplyr::bind_rows()%>%
  gather(-PairID,key="D_or_R",value="nWGS")

#Read in the posterior from E. Mitchell et al, 2022 to use as the prior
param_posterior<-read.delim(posterior_parameters_file)

###-----------------Set up parameters for this simulation-----------------

#Set the HSCT parameters
age=runif(1,min=20,max=100)

###-----------------SIMULATION APPROACH-----------------
#Grown tree to age of transplant (driver process simulation)
#Create new cell compartment of transplanted (recipient) cells
#Grow each compartment independently to HSC population size
system(paste("mkdir -p",output_dir))
cat(paste("Running",nsim," AGE simulations"),sep="\n")

for(sim in 1:nsim) {
  cat(paste("Running simulation number",sim),sep="\n")
  
  uid=ids::random_id()
  
  #Choose age to run to
  age=runif(1,min=20,max=100)

  #Set the variable parameters
  HSC_pop_size=round(runif(1,53500,250000))
  HSC_symmetric_division_rate=1/(2*365)
  
  #Driver mutation parameters distribution parameters (from posterior of E. Mitchell et al, 2022)
  param_idx=sample(1:nrow(param_posterior),1)
  number_drivers_per_year = param_posterior$number_drivers_per_year[param_idx]
  gamma_shape = param_posterior$gamma_shape[param_idx]
  gamma_rate = param_posterior$gamma_rate[param_idx]
  fitness_threshold = 0.05
  dpcpd=min(1e-5,number_drivers_per_year/(365*HSC_pop_size)) #Number of drivers per cell per year cannot be more than 1e-5
  
  ##Function to generate gamma distribution based fitness
  genGammaFitness=function(fitness_threshold,shape,rate){
    function() rtrunc(n=1,a=fitness_threshold, b=Inf,"gamma",shape=shape,rate=rate)
  }
  fitnessGammaFn=genGammaFitness(fitness_threshold=fitness_threshold,shape = gamma_shape, rate=gamma_rate)
  
  ##Function to generate exp distribution based fitness
  genExpFitness=function(fitness_threshold,rate){
    function() rtrunc(n=1,a=fitness_threshold, b=Inf,"exp",rate=rate)
  }
  fitnessExpFn=genExpFitness(fitness_threshold=0.08,rate=40)
  
  ##Put parameters into a single list (to save later)
  params=list(age=age,
              HSC_symmetric_division_rate=HSC_symmetric_division_rate,
              HSC_pop_size=HSC_pop_size,
              number_drivers_per_year=number_drivers_per_year,
              gamma_shape=gamma_shape,
              gamma_rate=gamma_rate,
              fitness_threshold=fitness_threshold,
              dpcpd=dpcpd)
  
  print(params)
  
  ###START ACTUAL SIMULATION
  #Grow population to the age of HSCT
  cat("Starting simulation",sep="\n")
  pop_final=run_driver_process_sim(initial_division_rate=0.1,
                              final_division_rate = HSC_symmetric_division_rate,
                              target_pop_size = HSC_pop_size,
                              nyears = age,
                              fitnessGen=fitnessGammaFn,
                              drivers_per_cell_per_day = dpcpd)
  

  #VISUALIZE RESULTS
  ss_by_tree<-lapply(colony_numbers_df$nWGS,function(nWGS) {
    
    #Subsample the tree with the right number of colonies
    tree<-get_subsampled_tree2(pop_final,N = nWGS)
    
    #Create the mutation tree/ ultrametric tree for extracting the summary statistics
    tree_m<-get_elapsed_time_tree(tree,mutrateperdivision = 1.8,backgroundrate = 15/365)
    tree_m_ultra<-make.ultrametric.tree(tree_m);tree_m_ultra$edge.length[1]<-0
    
    ##EXTRACT SUMMARY STATISTICS
    #(1) 3 LARGEST CLADES
    largest_clades<-get_expanded_clade_nodes(tree_m,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
      pull(n_samples)%>%sort(decreasing = T)%>%.[1:3]
    
    #(2) Number of singletons
    n_singletons<-get_expanded_clade_nodes(tree_m,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
      dplyr::filter(n_samples==1)%>%nrow(.)
    
    #(3) Peri-HSCT LTT/ coalescences
    ltt<-get_ltt(tree=tree_m_ultra,time_points=c(0.2,0.4,0.6))
    coals<-get_coalescences(ltt)
    
    #(4) Bulk totals for expanded populations to look at clonal shifts
    expanded_df=get_expanded_clade_nodes(tree_m,min_clonal_fraction=0.001,height_cut_off=50,min_samples=3)
    
    if(nrow(expanded_df)>0){
      expanded_df<-expanded_df%>%
        mutate(count=sapply(nodes,function(node) {length(getTips(tree_m,node))}))%>%
        mutate(prop=count/nWGS)
      
      prop_from_expansions=sum(expanded_df$prop)
      
    } else {
      prop_from_expansions=0
    }
    
    ##--COMBINE SUMSTATS INTO SINGLE VECTOR--
    ss_names=c(paste("largest_clades",1:3,sep="_"),
               "n_singletons",
               paste("coals",1:2,sep="_"),
               "prop_from_expansions")
    
    ss_comb=c(largest_clades,
              n_singletons,
              coals,
              prop_from_expansions)
    
    names(ss_comb)<-ss_names
    
    #Save the key model parameters/ outputs - including the tree & expanded clone df that could be used to get different summary stats
    res=list(uid=uid,
             params=params,
             sumstats=ss_comb,
             tree=tree_m,
             expanded_clone_df=expanded_df)
    return(res)
  })
  
  cat("Saving output",sep="\n\n")
  saveRDS(ss_by_tree,file=paste0(output_dir,"/","res_",uid,".Rds"))
}
