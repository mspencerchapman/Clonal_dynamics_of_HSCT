## Run on farm5
args = commandArgs(TRUE)
args

source_directory = toString(args[1])
output_directory = toString(args[2])
output_trees_directory = toString(args[3])
prior_sample_drift_file = toString(args[4])
prior_sample_drivers_file = toString(args[5])
pair_ID=as.integer(args[6])
model_name = toString(args[7])
n_sims_per_job=as.integer(args[8])
sim_job_index=as.integer(args[9])
n_jobs_accumulated=as.integer(args[10])

####################################################################################################
#!/software/R-4.1.0/bin/Rscript
library(optparse)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(phytools)
library(rsimpop)

####################################################################################################
#
# Set seed for RNG;
# This ensures that for this "job" (identified by the "sim_job_index")
# within the current "job array",
# the RNG is set to a unique seed.
#

n_jobs_accumulated = n_jobs_accumulated + sim_job_index
set.seed( n_jobs_accumulated )

#
# Generate a random integer: run > n_sims_per_job;
# This ensures that each simulation i (= 1, 2, ..., n_sims_per_job)
# will have a unique seed = run + i
#

# run <- round(runif( n=1, min= 1e6, max = 2e6), digits = 0)
run <- round(runif( n=1, min= 1e6, max = 2e6), digits = 0)
cat( "0.0. sim_prior_transplant.R main: run = ", run, "\n", sep = "\t" )
# readLines("stdin",n=1)

####################################################################################################
# This directory is for output saved in my output format!

# Create the output directory:
if( !file.exists( output_directory ) )
{
	dir.create( output_directory )
	
}

# Create the output directory:
if( !file.exists( output_trees_directory ) )
{
	dir.create( output_trees_directory )
	
}

####################################################################################################
#
# Create (donor-specific) directory paths;
#

# PairID=paste0("Pair",pair_index)
donor_ID = paste0("Pair",pair_ID) # kjd

cat( "pair_ID = ", pair_ID, "\n", sep = "\t" ) # kjd
cat( "donor_ID = ", donor_ID, "\n", sep = "\t" ) # kjd

####################################################################################################

nsim = n_sims_per_job
# increase_late_pop=F
# HSC_pop_size_final_10yrs=5e5

####################################################################################################
# 
# Create data frame "sim_table" columns:
#

if( model_name == "m1" )
{
  
  # param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "fitnessGammaFn", "dpcpd" )
  param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
  
}else if( model_name == "m2" )
{
  
  # param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
  param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd", "gamma_shape_engraftment", "gamma_rate_engraftment" )
  
}else if( model_name == "m3" )
{
  
  # param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
  param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd", "exaggeration_ratio", "exaggeration_prop" )
  
}

####################################################################################################
####################################################################################################
#
# Create list "stat_names_vect_list";
#
  
  stat_set_name_vect = c( "original", "pre_interval_divided", "peri_interval_divided", "peri_interval_narrower", "peri_interval_wider" )
  stat_set_name_list = as.list( stat_set_name_vect )
  
  stat_names_vect_list = list()
  
  original_stat_names_vect = c(paste("largest_clade_D",1:3,sep="_"),
                          paste("largest_clade_R",1:3,sep="_"),
                          "n_singletons_D",
                          "n_singletons_R",
                          "DR_coals_pre",
                          "D_coals_pre",
                          "R_coals_pre",
                          "DR_coals_peri",
                          "D_coals_peri",
                          "R_coals_peri",
                          "R_max_coals_in_single_clade",
                          "mean_abs_log2FC",
                          "median_abs_log2FC",
                          "max_abs_log2FC",
                          "mean_abs_change",
                          "median_abs_change",
                          "max_abs_change")
  
  new_list = list( original_stat_names_vect )
  names( new_list ) <- "original"
  
  stat_names_vect_list = append( stat_names_vect_list, new_list )
  
  
  pre_interval_divided_stat_names_vect = c(paste("largest_clade_D",1:3,sep="_"),
                          paste("largest_clade_R",1:3,sep="_"),
                          "n_singletons_D",
                          "n_singletons_R",
                          "DR_coals_pre_1", # new
                          "DR_coals_pre_2", # new
                          "D_coals_pre_1", # new
                          "D_coals_pre_2", # new
                          "R_coals_pre_1", # new
                          "R_coals_pre_2", # new
                          "DR_coals_peri",
                          "D_coals_peri",
                          "R_coals_peri",
                          "R_max_coals_in_single_clade",
                          "mean_abs_log2FC",
                          "median_abs_log2FC",
                          "max_abs_log2FC",
                          "mean_abs_change",
                          "median_abs_change",
                          "max_abs_change")
	
  new_list = list( pre_interval_divided_stat_names_vect )
  names( new_list ) <- "pre_interval_divided"
  
  stat_names_vect_list = append( stat_names_vect_list, new_list )
  
  
  peri_interval_divided_stat_names_vect = c(paste("largest_clade_D",1:3,sep="_"),
                          paste("largest_clade_R",1:3,sep="_"),
                          "n_singletons_D",
                          "n_singletons_R",
                          "DR_coals_pre",
                          "D_coals_pre",
                          "R_coals_pre",
                          "DR_coals_peri_1", # new
                          "DR_coals_peri_2", # new
                          "D_coals_peri_1", # new
                          "D_coals_peri_2", # new
                          "R_coals_peri_1", # new
                          "R_coals_peri_2", # new
                          "R_max_coals_in_single_clade",
                          "mean_abs_log2FC",
                          "median_abs_log2FC",
                          "max_abs_log2FC",
                          "mean_abs_change",
                          "median_abs_change",
                          "max_abs_change")
	
  
  new_list = list( peri_interval_divided_stat_names_vect )
  names( new_list ) <- "peri_interval_divided"
  
  stat_names_vect_list = append( stat_names_vect_list, new_list )
  
  
  new_list = list( original_stat_names_vect )
  names( new_list ) <- "peri_interval_narrower"
  
  stat_names_vect_list = append( stat_names_vect_list, new_list )
  
  
  new_list = list( original_stat_names_vect )
  names( new_list ) <- "peri_interval_wider"
  
  stat_names_vect_list = append( stat_names_vect_list, new_list )

cat( "1.0. sim_prior_transplant.R main: names( stat_names_vect_list ) = ", names( stat_names_vect_list ), "\n", sep = "\t" )
cat( "1.0. sim_prior_transplant.R main: length( stat_names_vect_list ) = ", length( stat_names_vect_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)

####################################################################################################
####################################################################################################
#
# More function defs:
#
####################################################################################################
####################################################################################################

make_sim_table <- function( n_sims, param_names_vect, stat_names_vect )
{
	
	sim_col_names_vect = c( "seed", param_names_vect, stat_names_vect )
	sim_table = as.data.frame( matrix( NA, nrow = n_sims, ncol = length( sim_col_names_vect ) ) )
	names( sim_table ) <- sim_col_names_vect
	
	return( sim_table )
	
}

####################################################################################################

update_sim_table <- function( sim_table, row_index, selected_col_names_vect, selected_info_vect )
{
	
	if( length( selected_info_vect ) == length( selected_col_names_vect ) )
	{
		sim_table[ row_index, selected_col_names_vect ] <- selected_info_vect
		
	}
	
	return( sim_table )
	
}

####################################################################################################

write_table_row_to_outfile <- function( sim, sim_table, sim_outfile )
{
	
	if( sim == 1 )
	{
		write.table( sim_table[ sim, ], file = sim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
		# write.table( sim_table[ sim, ], file = sim_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
		
	}else
	{
		# write.table( sim_table[ sim, ], file = sim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
		write.table( sim_table[ sim, ], file = sim_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
		
	}
	
}

####################################################################################################

####################################################################################################
# This function computes the vector of summary statistics from the input tree.

calc_sim_stats_vect <- function( stat_set_name, ss_names, DR_tree_m, D_tree_m, R_tree_m, DR_tree_m_ultra, D_tree_m_ultra, R_tree_m_ultra, largest_clades_D, largest_clades_R, n_singletons_D, n_singletons_R, age_of_donor_at_HSCT, pop_final, no_of_WGS_recip, no_of_WGS_donor )
{
	
  #(3) Peri-HSCT LTT/ coalescences (Combined, Donor, Recipient)
  
  # peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
  if( stat_set_name == "original" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
  
  }else if( stat_set_name == "peri_interval_narrower" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-2.5,x+2.5)))
  
  }else if( stat_set_name == "peri_interval_wider"  )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-7.5,x+7.5)))
  
  }else if( stat_set_name == "pre_interval_divided" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5))) # same as "original"!
  
  }else if( stat_set_name == "peri_interval_divided" )
  {
	peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x,x+5))) # notice extra point (at "x" = age_of_donor_at_HSCT) in vector c(x-5,x,x+5);
  
  }
  
  DR_ltt_peri<-get_ltt(tree=DR_tree_m_ultra,time_points=peri_HSCT_time_points)
  D_ltt_peri<-get_ltt(tree=D_tree_m_ultra,time_points=peri_HSCT_time_points)
  R_ltt_peri<-get_ltt(tree=R_tree_m_ultra,time_points=peri_HSCT_time_points)
  DR_coals_peri<-get_coalescences(DR_ltt_peri)
  D_coals_peri<-get_coalescences(D_ltt_peri)
  R_coals_peri<-get_coalescences(R_ltt_peri)
  
  #(4) Pre-HSCT LTT/ coalescences (Combined, Donor, Recipient)
  
  # pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
  if( stat_set_name == "original" )
  {
	pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
  
  }else if( stat_set_name == "peri_interval_narrower" )
  {
	pre_HSCT_time_points=c(5,peri_HSCT_time_points[1]) # same as "original"; except that peri_HSCT_time_points[1] has been moved up!
  
  }else if( stat_set_name == "peri_interval_wider"  )
  {
	pre_HSCT_time_points=c(5,peri_HSCT_time_points[1]) # same as "original"; except that peri_HSCT_time_points[1] has been moved down!
  
  }else if( stat_set_name == "pre_interval_divided" )
  {
	mid_point = 5 + ((peri_HSCT_time_points[1] - 5)/2)
	pre_HSCT_time_points=c(5,mid_point,peri_HSCT_time_points[1]) # notice extra point ("mid_point");
  
  }else if( stat_set_name == "peri_interval_divided" )
  {
	pre_HSCT_time_points=c(5,peri_HSCT_time_points[1]) # same as "original"!
  
  }
  
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
  
# expanded_df = matrix( NA, nrow = 0, ncol = 1 ) # kjd: The variable "pop_final" is a hidden variable;  "pop_final" is NOT a summary statistic of the tree;
  
  if(nrow(expanded_df)>0){
    expanded_df<-expanded_df%>%
      mutate(D_count=sapply(nodes,function(node) {sum(grepl("donor",getTips(DR_tree_m,node)))}),
             R_count=sapply(nodes,function(node) {sum(grepl("recip",getTips(DR_tree_m,node)))}))%>%
      mutate(D_phylo=D_count/no_of_WGS_donor,
             R_phylo=R_count/no_of_WGS_recip)%>%
      left_join(DR_tree$events%>%filter(driverid>0)%>%dplyr::select(node,driverid),by=c("nodes"="node"))
    
    D_total=sum(pop_final$cfg$info$population[pop_final$cfg$info$val==1])
    R_total=sum(pop_final$cfg$info$population[pop_final$cfg$info$val==2])
    expanded_df<-expanded_df%>%
      left_join(pop_final$cfg$info%>%filter(val==1)%>%dplyr::select(id,population),by=c("driverid"="id"))%>%
      mutate(D_bulk_by_driver=population/D_total)%>%
      dplyr::select(-population)%>%
      left_join(pop_final$cfg$info%>%filter(val==2)%>%dplyr::select(id,population),by=c("driverid"="id"))%>%
      mutate(R_bulk_by_driver=population/R_total)%>%
      dplyr::select(-population)
    
    #Make sure each node only has one entry in the table - summarise information by node
    all_driver_nodes<-expanded_df%>%pull(nodes)%>%unique()
    expanded_df<-lapply(all_driver_nodes,function(driver_node) {
      node_driver_df<-expanded_df%>%filter(nodes==driver_node)
      data.frame(nodes=driver_node,
                 n_samples=node_driver_df$n_samples[1],
                 MRCA_time=node_driver_df$MRCA_time[1],
                 clonal_fraction=node_driver_df$clonal_fraction[1],
                 D_count=node_driver_df$D_count[1],
                 R_count=node_driver_df$R_count[1],
                 D_phylo=node_driver_df$D_phylo[1],
                 R_phylo=node_driver_df$R_phylo[1],
                 driverid=paste0(node_driver_df$driverid,collapse=","),
                 D_bulk_by_driver=mean(node_driver_df$D_bulk_by_driver,na.rm=T),
                 R_bulk_by_driver=mean(node_driver_df$R_bulk_by_driver,na.rm=T))
    })%>%dplyr::bind_rows()
    
    #Get the bulk fractions from all nodes (the )
    bulk_fractions<-extract_bulk_cell_fractions(nodes=expanded_df$nodes,sub_pop = DR_tree,full_pop=pop_final,states=1:2)%>%
      dplyr::rename("D_bulk"=state_1,"R_bulk"=state_2)
    
    expanded_df<-left_join(expanded_df,bulk_fractions,by=c("nodes"="node"))%>%
      mutate(log2FC=log2(R_bulk/D_bulk))%>%
      mutate(abs_change=R_bulk-D_bulk)
    
    mean_abs_log2FC=mean(abs(expanded_df$log2FC))
    median_abs_log2FC=median(abs(expanded_df$log2FC))
    max_abs_log2FC=max(abs(expanded_df$log2FC))
    mean_abs_change=mean(abs(expanded_df$abs_change))
    median_abs_change=median(abs(expanded_df$abs_change))
    max_abs_change=max(abs(expanded_df$abs_change))
    
  } else {
    mean_abs_log2FC=NA
    median_abs_log2FC=NA
    max_abs_log2FC=NA
    mean_abs_change=NA
    median_abs_change=NA
    max_abs_change=NA
  }
  
  ##--COMBINE SUMSTATS INTO SINGLE VECTOR--
#   ss_names=c(paste("largest_clade_D",1:3,sep="_"),
#              paste("largest_clade_R",1:3,sep="_"),
#              "n_singletons_D",
#              "n_singletons_R",
#              paste(c("DR","D","R"),"coals_pre",sep="_"),
#              paste(c("DR","D","R"),"coals_peri",sep="_"),
#              "R_max_coals_in_single_clade",
#              "mean_abs_log2FC",
#              "median_abs_log2FC",
#              "max_abs_log2FC",
#              "mean_abs_change",
#              "median_abs_change",
#              "max_abs_change")
  
  if( stat_set_name %in% c( "original", "peri_interval_narrower", "peri_interval_wider" ) )
  {
	ss_comb=c(largest_clades_D,
            largest_clades_R,
            n_singletons_D,
            n_singletons_R,
            DR_coals_pre,
            D_coals_pre,
            R_coals_pre,
            DR_coals_peri,
            D_coals_peri,
            R_coals_peri,
            R_max_coals_in_single_clade,
            mean_abs_log2FC,
            median_abs_log2FC,
            max_abs_log2FC,
            mean_abs_change,
            median_abs_change,
            max_abs_change)
	  names(ss_comb)<-ss_names
	  
  }else if( stat_set_name == "pre_interval_divided" )
  {
	
	DR_coals_pre_1 = DR_coals_pre[ 1 ] # new
    DR_coals_pre_2 = DR_coals_pre[ 2 ] # new
    D_coals_pre_1 = D_coals_pre[ 1 ] # new
    D_coals_pre_2 = D_coals_pre[ 2 ] # new
    R_coals_pre_1 = R_coals_pre[ 1 ] # new
    R_coals_pre_2 = R_coals_pre[ 2 ] # new
    
	ss_comb=c(largest_clades_D,
            largest_clades_R,
            n_singletons_D,
            n_singletons_R,
            DR_coals_pre_1, # new
            DR_coals_pre_2, # new
            D_coals_pre_1, # new
            D_coals_pre_2, # new
            R_coals_pre_1, # new
            R_coals_pre_2, # new
            DR_coals_peri,
            D_coals_peri,
            R_coals_peri,
            R_max_coals_in_single_clade,
            mean_abs_log2FC,
            median_abs_log2FC,
            max_abs_log2FC,
            mean_abs_change,
            median_abs_change,
            max_abs_change)
	  names(ss_comb)<-ss_names
	  
  }else if( stat_set_name == "peri_interval_divided" )
  {
	
	DR_coals_peri_1 = DR_coals_peri[ 1 ] # new
    DR_coals_peri_2 = DR_coals_peri[ 2 ] # new
    D_coals_peri_1 = D_coals_peri[ 1 ] # new
    D_coals_peri_2 = D_coals_peri[ 2 ] # new
    R_coals_peri_1 = R_coals_peri[ 1 ] # new
    R_coals_peri_2 = R_coals_peri[ 2 ] # new
	
	ss_comb=c(largest_clades_D,
            largest_clades_R,
            n_singletons_D,
            n_singletons_R,
            DR_coals_pre,
            D_coals_pre,
            R_coals_pre,
            DR_coals_peri_1, # new
            DR_coals_peri_2, # new
            D_coals_peri_1, # new
            D_coals_peri_2, # new
            R_coals_peri_1, # new
            R_coals_peri_2, # new
            R_max_coals_in_single_clade,
            mean_abs_log2FC,
            median_abs_log2FC,
            max_abs_log2FC,
            mean_abs_change,
            median_abs_change,
            max_abs_change)
	  names(ss_comb)<-ss_names
	  
  }
  
	# return( stats_vect )
	return( ss_comb )
	
}

####################################################################################################

####################################################################################################
####################################################################################################
# 
# Create list (of data frames) "sim_table_list";
#

sim_table_list = Map( function(x){ make_sim_table( nsim, param_names, x ) }, stat_names_vect_list )
# names( sim_table_list ) <- stat_set_name_list
names( sim_table_list ) <- names( stat_names_vect_list )

cat( "2.0. sim_prior_transplant.R main: names( sim_table_list ) = ", names( sim_table_list ), "\n", sep = "\t" )
cat( "2.0. sim_prior_transplant.R main: length( sim_table_list ) = ", length( sim_table_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)

####################################################################################################
####################################################################################################
# 
# Create list "sim_tree_list";
#

sim_tree_list = list()

####################################################################################################

# my_working_directory<-"/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/"
R_functions_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/my_functions","/lustre/scratch126/casm/team154pc/ms56/my_functions")
tree_mut_dir=ifelse(Sys.info()["sysname"] == "Darwin","~/R_work/treemut","/lustre/scratch126/casm/team154pc/ms56/fetal_HSC/treemut")
R_function_files=list.files(R_functions_dir,pattern=".R",full.names = T)
sapply(R_function_files[-2],source)
setwd(tree_mut_dir); source("treemut.R") # ;setwd(my_working_directory)
visualize=F

####################################################################################################

my_working_directory = output_directory
setwd(my_working_directory)

####################################################################################################

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

####################################################################################################

###-----------------SPECIFIC FUNCTIONS FOR MODEl "m3"-----------------

## reset coefficients - for specified compartment
resetCoeff=function(tree,
                    driver_clones,
                    selection_change_ratio,
                    target_compartment=2,b.verbose=FALSE,b.check=TRUE){
  cat("Changing by ratio=",selection_change_ratio,"\n")
  #browser()
  # First take a copy of the drivers data.frame and then rescale the target drivers fitness
  drivers=tree$cfg$drivers %>% filter(fitness>0)
  drivers=drivers %>% mutate(fitness2=ifelse(driver %in% driver_clones,fitness*selection_change_ratio,fitness))
  
  # Get info augmented with the the NODE id of the most recent driver acquisition of each compartment
  # Note the info ID column corresponds to the driverid of the most recent driver and sub-compartment fitness is
  # the sum of the fitness of the ancestral drivers.
  inf=tree$cfg$info %>% filter(val==target_compartment) %>% # Filter to just the target compartment 
    left_join(tree$events %>% filter(driverid>0)  %>% select(-value),by=c("id"="driverid")) #Select just driver events (i.e. not differentiation events etc)
  
  # Identify a fitness and node with each driver event - note we can ignore the "val" column in cfg$drivers and for drivers we can also ignore the value column in events. 
  driver_events=tree$events %>% filter(driverid>0) %>% left_join(drivers,by=c("driverid"="driver"))
  
  # sub-compartment fitness is the sum of the fitness of the ancestral drivers.
  inf$chk_fitness=sapply(1:length(inf$node),function(i){
    node=inf$node[i]
    pnodes=rsimpop:::get_parents(node,tree$edge)
    alldrivers=driver_events %>% filter(node %in% pnodes)
    if(dim(alldrivers)[1]==0){
      0
    }else{
      sum(alldrivers$fitness2,na.rm = TRUE)
    }
  })
  changes=inf %>% filter(abs(chk_fitness-fitness)>1e-6) %>% head()
  if(abs(selection_change_ratio-1)<1e-6 && dim(changes)[1]>0 && b.check){
    stop("Inconsistent Selection Coefficient Reconstruction!")
  }
  if(b.verbose){
    cat("The compartments with altered fitness (head)\n")
    print(changes)
  }
  # browser()
  # Finally override the target compartment fitness in the original cfg$info dataframe.  
  # Note that the original cfg$drivers is left unchanged.
  tree$cfg$info=tree$cfg$info %>% left_join(inf %>% select(val,id,chk_fitness)) %>% 
    mutate(fitness=ifelse(is.na(chk_fitness),fitness,chk_fitness)) %>% select(-chk_fitness)
  tree
}

####################################################################################################

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
# tree_folder="/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/filtering_runs2/tree_files"
tree_folder="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/filtering_runs2/tree_files"
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

####################################################################################################
# 
# Read files:
# "prior_sample_drift_file";
# "prior_sample_drivers_file";
#

# prior_drift_table <- read.table( file = prior_sample_drift_file, sep = "\t", header=TRUE, stringsAsFactors = FALSE )
HSC_pop_posteriors <- read.table( file = prior_sample_drift_file, sep = "\t", header=TRUE, stringsAsFactors = FALSE )

# prior_drivers_table <- read.table( file = prior_sample_drivers_file, sep = "\t", header=TRUE, stringsAsFactors = FALSE )
param_posterior <- read.table( file = prior_sample_drivers_file, sep = "\t", header=TRUE, stringsAsFactors = FALSE )

####################################################################################################

###-----------------Set up parameters for this simulation-----------------

PairID=paste0("Pair",pair_ID)

cat( "pair_ID = ", pair_ID, "\n", sep = "\t" )
cat( "PairID = ", PairID, "\n", sep = "\t" )

age_of_donor_at_HSCT=Pair_metadata$Age_at_transplant[Pair_metadata$Pair==PairID]
age_of_donor_at_sampling=Pair_metadata$Age[Pair_metadata$Pair==PairID]

tree.data<-all.trees[[PairID]]
DR_ids<-get_DR_ids(tree.data)

donor_id<-DR_ids['donor_ID'] # kjd
recip_id<-DR_ids['recip_ID'] # kjd

cat( "donor_id = ", donor_id, "\n", sep = "\t" ) # kjd
cat( "recip_id = ", recip_id, "\n", sep = "\t" ) # kjd


no_of_WGS_recip=sum(grepl(DR_ids['recip_ID'],tree.data$tip.label))
no_of_WGS_donor=sum(grepl(DR_ids['donor_ID'],tree.data$tip.label))

####################################################################################################
# 
# Create outfiles;
#

sim_tree_Rds_outfile = paste( output_trees_directory, "/", "sim_trees_", "job_index_", sim_job_index, "_run_", run, ".Rds", sep = "" )

####################################################################################################

output_dir_list = lapply( stat_set_name_list, function(x){ paste( output_directory, x, "/", sep="" ) } )

outfile_list = lapply( output_dir_list, function(x){ paste( x, "sim_", "job_index_", sim_job_index, "_run_", run, ".txt", sep="" ) } )

####################################################################################################
# 
# Do sims;
#

# for( sim in 1:n_sims_per_job ) {
for(sim in 1:nsim) {
  cat(paste("Running simulation number",sim),sep="\n")
  
  # Re-set seed in RNG (in R and in C++);
  initSimPop(run+sim,bForce=TRUE)
  # set.seed(run+i)
  
  uid=ids::random_id()
  
  
  #Set the fixed HSCT parameters
  age_10yrs_pre_sampling=max(age_of_donor_at_sampling-10,age_of_donor_at_HSCT)
  HSC_symmetric_division_rate=1/(2*365)
  
  #Set the variable parameters
  # HSC_pop_size=round(sample(HSC_pop_posteriors$V1,1)) # kjd
  HSC_pop_size=round(sample(HSC_pop_posteriors$target_pop_size,1)) # kjd
  log10_HSCT_bottleneck=runif(1,min=2.7,max=4.7)
  HSCT_bottleneck=round(10^log10_HSCT_bottleneck)
  
  while(HSCT_bottleneck>HSC_pop_size) {
    cat("Value for bottleneck is greater than the total HSC population size... Reselecting parameters.",sep="\n\n")
    #Set the variable parameters
    # HSC_pop_size=round(sample(HSC_pop_posteriors$V1,1)) # kjd
    HSC_pop_size=round(sample(HSC_pop_posteriors$target_pop_size,1)) # kjd
    log10_HSCT_bottleneck=runif(1,min=2.7,max=4.7)
    HSCT_bottleneck=round(10^log10_HSCT_bottleneck)
  }
  
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
  
  
  ####################################################################################################
  # 
  # Create (named) vector "sim_params_vect" columns:
  #
  
  sim_params_vect = rep( NA, length( param_names ) )
  names( sim_params_vect ) <- param_names
  
  ####################################################################################################
  
  if( model_name == "m1" )
  {
	
	# param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
	
	# sim_table[ sim, param_names ] <- params[[ param_names ]] # kjd
	sim_params_vect[ "HSC_pop_size" ] = HSC_pop_size # kjd
	sim_params_vect[ "HSCT_bottleneck" ] = HSCT_bottleneck # kjd
	sim_params_vect[ "number_drivers_per_year" ] = number_drivers_per_year # kjd
	sim_params_vect[ "gamma_shape" ] = gamma_shape # kjd
	sim_params_vect[ "gamma_rate" ] = gamma_rate # kjd
	# sim_params_vect[ "fitnessGammaFn" ] = fitnessGammaFn # kjd
	sim_params_vect[ "dpcpd" ] = dpcpd # kjd
	
	
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
	
	
  }else if( model_name == "m2" )
  {
	
	# param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
	# param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd", "gamma_shape_engraftment", "gamma_rate_engraftment" )
	
	gamma_shape_engraftment = runif( 1, min=0.1, max=1.5 )
	gamma_rate_engraftment = runif( 1, min=0.1, max=1.5 )
	
	# sim_table[ sim, param_names ] <- params[[ param_names ]] # kjd
	sim_params_vect[ "HSC_pop_size" ] = HSC_pop_size # kjd
	sim_params_vect[ "HSCT_bottleneck" ] = HSCT_bottleneck # kjd
	sim_params_vect[ "number_drivers_per_year" ] = number_drivers_per_year # kjd
	sim_params_vect[ "gamma_shape" ] = gamma_shape # kjd
	sim_params_vect[ "gamma_rate" ] = gamma_rate # kjd
	# sim_params_vect[ "fitnessGammaFn" ] = fitnessGammaFn # kjd
	sim_params_vect[ "dpcpd" ] = dpcpd # kjd
	
	sim_params_vect[ "gamma_shape_engraftment" ] = gamma_shape_engraftment # kjd
	sim_params_vect[ "gamma_rate_engraftment" ] = gamma_rate_engraftment # kjd
	
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
	
  }else if( model_name == "m3" )
  {
	
	# param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd" )
	# param_names=c( "HSC_pop_size", "HSCT_bottleneck", "number_drivers_per_year", "gamma_shape", "gamma_rate", "dpcpd", "exaggeration_ratio", "exaggeration_prop" )
	
	#Set the selection parameters
	exaggeration_ratio=runif(1,min=1.5,max=6)
	exaggeration_prop=runif(1,min=0.1,max=0.3)
	
	
	# sim_table[ sim, param_names ] <- params[[ param_names ]] # kjd
	sim_params_vect[ "HSC_pop_size" ] = HSC_pop_size # kjd
	sim_params_vect[ "HSCT_bottleneck" ] = HSCT_bottleneck # kjd
	sim_params_vect[ "number_drivers_per_year" ] = number_drivers_per_year # kjd
	sim_params_vect[ "gamma_shape" ] = gamma_shape # kjd
	sim_params_vect[ "gamma_rate" ] = gamma_rate # kjd
	# sim_params_vect[ "fitnessGammaFn" ] = fitnessGammaFn # kjd
	sim_params_vect[ "dpcpd" ] = dpcpd # kjd
	
	sim_params_vect[ "exaggeration_ratio" ] = exaggeration_ratio # kjd
	sim_params_vect[ "exaggeration_prop" ] = exaggeration_prop # kjd
	
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
              dpcpd=dpcpd,
              exaggeration_ratio=exaggeration_ratio,
              exaggeration_prop=exaggeration_prop)
	
	
  }
  
  print(params)
  
  ####################################################################################################
  
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
  
  ####################################################################################################
  
  if( model_name == "m1" )
  {
  	
  #Select number of cells to be transplanted into the recipient as a new compartment "Transplanted"
  #These can then grow separately to the HSC population size
  cfg=dps2$cfg
  cfg=addCellCompartment(cfg,population = HSC_pop_size,rate=cfg$compartment$rate[2],ndriver=1,descr="Transplanted",basefit = 0.1)
  
  cat("Specifying transplanted cells",sep="\n")
  post_HSCT<-addDifferentiationEvents(tree=dps2,cfg = cfg,2,nEvent = HSCT_bottleneck)
  post_HSCT$cfg$info=rbind(post_HSCT$cfg$info,post_HSCT$cfg$info %>% filter(val==1 & id>0) %>% mutate(val=2))
  
  #Grow & age donor and recipient compartments (driver process simulation)
  cat("Growing and ageing donor & recipient compartments",sep="\n")
  pop_final<-run_driver_process_sim(simpop=post_HSCT,
                                    target_pop_size = HSC_pop_size,
                                    final_division_rate = HSC_symmetric_division_rate,
                                    nyears = age_of_donor_at_sampling,
                                    fitnessGen=fitnessGammaFn,
                                    #fitnessGen=fitnessExpFn,
                                    drivers_per_cell_per_day = dpcpd)
  
  ####################################################################################################
  
  }else if( model_name == "m2" )
  {
  
  #Select number of cells to be transplanted into the recipient as a new compartment "Transplanted"
  #These can then grow separately to the HSC population size
  cfg=dps2$cfg
  clones<-dps2$driverid[dps2$edge[,2]%in%1:dps2$ntips]
  
  #Assign a clone-specific "engraftment fitness" separate to the general fitness
  #The parameters set here are empirical, but designed to show that engraftment selection can mimic the features seen in the trees
  # engraftment_fitness=sapply(1:length(unique(clones)),function(i) genGammaFitness(fitness_threshold=0.05,shape = 0.5, rate=0.5)())
  # names(engraftment_fitness)<-sort(unique(clones))
  # engraftment_fitness['0']<-quantile(engraftment_fitness,0.33) #Make the germline clone equal to the 1/3 quantile - i.e. driver mutations are twice as likely to cause an engraftment advantage than a disadvantage relative to wild-type
  
  engraftment_fitness=sapply( 1:length( unique(clones) ), function(i) genGammaFitness(fitness_threshold=0.05,shape = gamma_shape_engraftment, rate = gamma_rate_engraftment )()) # Notice that there is now a prior on the parmaters (shape = gamma_shape_engraftment, rate = gamma_rate_engraftment) of this gamma distribution of fitness effects;
  names(engraftment_fitness)<-sort(unique(clones))
  engraftment_fitness['0']<-quantile(engraftment_fitness,0.33) #Make the germline clone equal to the 1/3 quantile - i.e. driver mutations are twice as likely to cause an engraftment advantage than a disadvantage relative to wild-type
  
  
  #Choose the tips based on these 'engraftment fitness' values - i.e. sample from the tips vector based on
  dps2$engraftmentFitness<-engraftment_fitness[as.character(clones)]
  idx<-sample(dps2$edge[,2][dps2$edge[,2]%in%1:dps2$ntips],size=HSCT_bottleneck,replace = F,prob = dps2$engraftmentFitness)
  idx<-idx[idx!=1] #Remove outgroup if this has been selected

  #Now specifiy the new compartment
  cfg=addCellCompartment(cfg,population = HSC_pop_size,rate=cfg$compartment$rate[2],ndriver=1,descr="Transplanted",basefit = 0.1)
  
  ######UPDATE THE 'addDifferentiationEvents' function to allow engraftment-related selection
  cat("Specifying transplanted cells",sep="\n")
  post_HSCT<-addDifferentiationEvents(tree=dps2,cfg = cfg,2,idx=idx,nEvent = HSCT_bottleneck)
  post_HSCT$cfg$info=rbind(post_HSCT$cfg$info,post_HSCT$cfg$info %>% filter(val==1 & id>0) %>% mutate(val=2))
  
  #Grow & age donor and recipient compartments (driver process simulation)
  cat("Growing and ageing donor & recipient compartments",sep="\n")
  pop_final<-run_driver_process_sim(simpop=post_HSCT,
                                    target_pop_size = HSC_pop_size,
                                    final_division_rate = HSC_symmetric_division_rate,
                                    nyears = age_of_donor_at_sampling,
                                    fitnessGen=fitnessGammaFn,
                                    #fitnessGen=fitnessExpFn,
                                    drivers_per_cell_per_day = dpcpd)
  
  ####################################################################################################
  
  }else if( model_name == "m3" )
  {
  
    if(packageVersion("rsimpop")<"2.2.7"){
      stop("Please use at least rsimpop 2.2.7 for this analysis!")
    }
    #Select number of cells to be transplanted into the recipient as a new compartment "Transplanted"
    #These can then grow separately to the HSC population size
    cfg=dps2$cfg
    cfg=addCellCompartment(cfg,population = HSC_pop_size,rate=cfg$compartment$rate[2],ndriver=1,descr="Transplanted",basefit = 0.1)
    
    cat("Specifying transplanted cells",sep="\n")
    post_HSCT<-addDifferentiationEvents(tree=dps2,cfg = cfg,2,nEvent = HSCT_bottleneck)
    post_HSCT$cfg$info=rbind(post_HSCT$cfg$info,post_HSCT$cfg$info %>% filter(val==1 & id>0) %>% mutate(val=2))
    
    ##ALTER/ EXAGGERATE SELECTION CO-EFFICIENTS FOR THE FOLLOWING 5 YEARS (This is the Type 2 selection bit)
    recipient_driver_ids<-post_HSCT$cfg$drivers%>%filter(val==1 & driver>0)%>%pull(driver)
    exaggerated_postengraftment_selection_clones<-sample(recipient_driver_ids,
                                                         size = round(exaggeration_prop*length(recipient_driver_ids)),
                                                         replace=F)
    
    #Cycle through the drivers that will have enhanced postengraftment selection
    # Reset the coefficients just in the recipient val=2 compartment
    change_compartment=2 ## The recipient compartment
    post_HSCT=resetCoeff(post_HSCT,
                         exaggerated_postengraftment_selection_clones,
                         exaggeration_ratio,
                         target_compartment = change_compartment)
    #Grow & age donor and recipient compartments for 5 years with exaggerate selection
    cat("Growing and ageing donor & recipient compartments for 5 years with exaggerated selection",sep="\n")
    pop_post5years<-run_driver_process_sim(simpop=post_HSCT,
                                           target_pop_size = HSC_pop_size,
                                           final_division_rate = HSC_symmetric_division_rate,
                                           nyears = age_of_donor_at_HSCT+5,
                                           fitnessGen=fitnessGammaFn,
                                           #fitnessGen=fitnessExpFn,
                                           drivers_per_cell_per_day = dpcpd,
                                           b_verbose=0)
    
    #Return selective coefficients to what they were (recognize these clones either by having very high fitness, or having a selected id)
    # Actually this just involves recalculating using the stored driver fitness as these haven't been changed.
    # So we're *not* going to divide through by exaggeration_ratio
    pop_post5years=resetCoeff(pop_post5years,exaggerated_postengraftment_selection_clones,1,target_compartment = change_compartment,b.check=FALSE)
    
    
    #post_HSCT$cfg$info%>%filter(val>0 & id>0)%>%ggplot(aes(x=fitness,y=population))+geom_point()+facet_grid(~val)
    #Grow & age donor and recipient compartments (driver process simulation)
    cat("Final growing and ageing donor & recipient compartments",sep="\n")
    pop_final<-run_driver_process_sim(simpop=pop_post5years,
                                      target_pop_size = HSC_pop_size,
                                      final_division_rate = HSC_symmetric_division_rate,
                                      nyears = age_of_donor_at_sampling,
                                      fitnessGen=fitnessGammaFn,
                                      #fitnessGen=fitnessExpFn,
                                      drivers_per_cell_per_day = dpcpd,b_verbose=0)
    
  
  } # if( model_name == "m3" )
  
  ####################################################################################################
  
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
  
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  
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
  
  
  # sim_stats_vect_list = Map( function(x,y){ calc_sim_stats_vect( x, y, DR_tree_m, D_tree_m, R_tree_m, DR_tree_m_ultra, D_tree_m_ultra, R_tree_m_ultra, largest_clades_D, largest_clades_R, n_singletons_D, n_singletons_R, age_of_donor_at_HSCT, pop_final, no_of_WGS_recip, no_of_WGS_donor ) }, stat_set_name_list, stat_names_vect_list )
  stats_result = try( sim_stats_vect_list <- Map( function(x,y){ calc_sim_stats_vect( x, y, DR_tree_m, D_tree_m, R_tree_m, DR_tree_m_ultra, D_tree_m_ultra, R_tree_m_ultra, largest_clades_D, largest_clades_R, n_singletons_D, n_singletons_R, age_of_donor_at_HSCT, pop_final, no_of_WGS_recip, no_of_WGS_donor ) }, stat_set_name_list, stat_names_vect_list ) )
  if( class( stats_result ) != "try-error" )
  {
	
	#
	# update sim RNG seed info;
	#
	
	sim_table_list = Map( function(x){ update_sim_table( x, sim, c( "seed" ), c( run+sim ) ) }, sim_table_list )
	
cat( "3.0. sim_prior_transplant.R main: names( sim_table_list ) = ", names( sim_table_list ), "\n", sep = "\t" )
cat( "3.0. sim_prior_transplant.R main: length( sim_table_list ) = ", length( sim_table_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	
	#
	# update sim parameter info;
	#
	
	sim_table_list = Map( function(x){ update_sim_table( x, sim, param_names, sim_params_vect ) }, sim_table_list )
	
cat( "3.1. sim_prior_transplant.R main: names( sim_table_list ) = ", names( sim_table_list ), "\n", sep = "\t" )
cat( "3.1. sim_prior_transplant.R main: length( sim_table_list ) = ", length( sim_table_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	
	#
	# update sim stats info;
	#
	
	sim_table_list = Map( function(x,y,z){ update_sim_table( x, sim, y, z ) }, sim_table_list, stat_names_vect_list, sim_stats_vect_list )
	
cat( "3.2. sim_prior_transplant.R main: names( sim_table_list ) = ", names( sim_table_list ), "\n", sep = "\t" )
cat( "3.2. sim_prior_transplant.R main: length( sim_table_list ) = ", length( sim_table_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	
  }
  
  ####################################################################################################
  
  #
  # write sim stats to outfile;
  #
  
  Map( function(x,y){ write_table_row_to_outfile( sim, x, y ) }, sim_table_list, outfile_list )
  
cat( "3.3. sim_prior_transplant.R main: names( sim_table_list ) = ", names( sim_table_list ), "\n", sep = "\t" )
cat( "3.3. sim_prior_transplant.R main: length( sim_table_list ) = ", length( sim_table_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)
  
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  
  #
  # update list "sim_tree_list":
  # insert new element "tree";
  # insert new element name "sim_seed_string";
  #
  
  uid_string <- paste( "uid_", uid, sep = "" )
 
  sim_tree_info = list( seed = run+sim, params_vect = sim_params_vect, tree = DR_tree )
  
  new_tree_info_list = list( sim_tree_info )
  names( new_tree_info_list ) <- uid_string
  # names( new_tree_info_list ) <- sim_seed_string
  
  sim_tree_list = append( sim_tree_list, new_tree_info_list ) # kjd 15/09/2023
  
  #
  # write sim trees to (.Rds) outfile;
  # (over-write previous version of .Rds outfile)
  #
  
  # saveRDS( sim_tree_list, file = sim_tree_Rds_outfile )
  saveRDS( sim_tree_list, file = sim_tree_Rds_outfile ) # kjd 15/09/2023
  
  ####################################################################################################
  
} # for(sim in 1:nsim)
# } # for( sim in 1:n_sims_per_job )

###############################################################################################

cat( "done!", "\n" )














