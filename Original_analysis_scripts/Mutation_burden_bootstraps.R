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




##BOOTSTRAP MUTATION BURDENS
pair="Pair40"
all_pairs=paste0("Pair",c(11,13,21,24,25,28,31,38,40,41))

bootstraps_res=list()
for(pair in all_pairs){
  cat(pair,sep="\n")
  tree=all.trees.adj2[[pair]]
  plot.phylo(tree,show.tip.label=F,direction="downwards")
  
  donor_id<-get_DR_ids(tree)['donor_ID']
  recip_id<-get_DR_ids(tree)['recip_ID']
  
  age_of_donor_at_HSCT<-Pair_metadata$Age_at_transplant[Pair_metadata$Pair==pair]
  age_of_donor_at_sampling<-Pair_metadata$Age[Pair_metadata$Pair==pair]
  
  
  #Model the mutation burden by a negative binomial to estimate the overdispersion parameter, theta
  tree_no_ancestral=drop.tip(tree,"Ancestral")
  mutation_burdens=nodeHeights(tree_no_ancestral)[tree_no_ancestral$edge[,2]%in%1:length(tree_no_ancestral$tip.label),2]
  mut_burden_df=data.frame(mutation_burden=round(mutation_burdens))
  nb.mod=glm.nb(mutation_burden~1,data=mut_burden_df)
  pois.mod=glm(mutation_burden~1,data=mut_burden_df,family = "poisson")
  theta_est=summary(nb.mod)$theta
  
  #Confirm the timing stats from the original tree
  #Make this tree ultrametric in the same way
  mean_mutation_burden=mean(get_mut_burden(tree))
  tree.ultra<-make.ultrametric.tree(tree)
  tree.ultra$edge.length=tree.ultra$edge.length*mean_mutation_burden
  tree.ultra$edge.length[which(tree$edge[,2]==which(tree$tip.label=="Ancestral"))]<-0 #Set ancestral branch length back to 0
  
  #(3) Peri-HSCT LTT/ coalescences
  tree.time<-tree.ultra
  tree.time$edge.length<-tree.time$edge.length*(age_of_donor_at_sampling/median(get_mut_burden(tree)))
  D_tree.time<-keep.tip(tree.time,tip=grep(donor_id,tree.time$tip.label))
  R_tree.time<-keep.tip(tree.time,tip=grep(recip_id,tree.time$tip.label))
  
  peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
  DR_ltt_peri<-get_ltt(tree=tree.time,time_points=peri_HSCT_time_points)
  D_ltt_peri<-get_ltt(tree=D_tree.time,time_points=peri_HSCT_time_points)
  R_ltt_peri<-get_ltt(tree=R_tree.time,time_points=peri_HSCT_time_points)
  DR_coals_peri<-get_coalescences(DR_ltt_peri)
  D_coals_peri<-get_coalescences(D_ltt_peri)
  R_coals_peri<-get_coalescences(R_ltt_peri)
  
  #(4) Pre-HSCT LTT/ coalescences
  pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
  DR_ltt_pre<-get_ltt(tree=tree.time,time_points=pre_HSCT_time_points)
  D_ltt_pre<-get_ltt(tree=D_tree.time,time_points=pre_HSCT_time_points)
  R_ltt_pre<-get_ltt(tree=R_tree.time,time_points=pre_HSCT_time_points)
  DR_coals_pre<-get_coalescences(DR_ltt_pre)
  D_coals_pre<-get_coalescences(D_ltt_pre)
  R_coals_pre<-get_coalescences(R_ltt_pre)
  
  data_stats=data.frame(Pair=pair,
                        DR_coals_peri=DR_coals_peri,
                        D_coals_peri=D_coals_peri,
                        R_coals_peri=R_coals_peri,
                        DR_coals_pre=DR_coals_pre,
                        D_coals_pre=D_coals_pre,
                        R_coals_pre=R_coals_pre)%>%
    gather(-Pair,key="Stat",value = "N_coals")
  
  ##Now do the bootstraps
  out_list=lapply(1:100,function(k) {
    cat(k,sep="\n")
    tree.boot=tree
    model="nb"
    if(model=="poisson"){
      tree.boot$edge.length=sapply(tree$edge.length,function(n) {rpois(1,lambda=n)})
    } else if(model=="nb") {
      tree.boot$edge.length=sapply(tree$edge.length,function(n) {rnegbin(1,mu=n,theta = theta_est)})
    }
    
    #Make this tree ultrametric in the same way
    mean_mutation_burden=mean(get_mut_burden(tree.boot))
    tree.boot.ultra<-make.ultrametric.tree(tree.boot)
    tree.boot.ultra$edge.length=tree.boot.ultra$edge.length*mean_mutation_burden
    tree.boot.ultra$edge.length[which(tree$edge[,2]==which(tree$tip.label=="Ancestral"))]<-0 #Set ancestral branch length back to 0
    
    #(3) Peri-HSCT LTT/ coalescences
    tree.time<-tree.boot.ultra
    tree.time$edge.length<-tree.time$edge.length*(age_of_donor_at_sampling/median(get_mut_burden(tree.boot)))
    D_tree.time<-keep.tip(tree.time,tip=grep(donor_id,tree.time$tip.label))
    R_tree.time<-keep.tip(tree.time,tip=grep(recip_id,tree.time$tip.label))
    
    peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))
    DR_ltt_peri<-get_ltt(tree=tree.time,time_points=peri_HSCT_time_points)
    D_ltt_peri<-get_ltt(tree=D_tree.time,time_points=peri_HSCT_time_points)
    R_ltt_peri<-get_ltt(tree=R_tree.time,time_points=peri_HSCT_time_points)
    DR_coals_peri<-get_coalescences(DR_ltt_peri)
    D_coals_peri<-get_coalescences(D_ltt_peri)
    R_coals_peri<-get_coalescences(R_ltt_peri)
    
    #(4) Pre-HSCT LTT/ coalescences
    pre_HSCT_time_points=c(5,peri_HSCT_time_points[1])
    DR_ltt_pre<-get_ltt(tree=tree.time,time_points=pre_HSCT_time_points)
    D_ltt_pre<-get_ltt(tree=D_tree.time,time_points=pre_HSCT_time_points)
    R_ltt_pre<-get_ltt(tree=R_tree.time,time_points=pre_HSCT_time_points)
    DR_coals_pre<-get_coalescences(DR_ltt_pre)
    D_coals_pre<-get_coalescences(D_ltt_pre)
    R_coals_pre<-get_coalescences(R_ltt_pre)
    
    sumstats_df<-data.frame(Pair=pair,
                            DR_coals_peri=DR_coals_peri,
                            D_coals_peri=D_coals_peri,
                            R_coals_peri=R_coals_peri,
                            DR_coals_pre=DR_coals_pre,
                            D_coals_pre=D_coals_pre,
                            R_coals_pre=R_coals_pre)
    return(list(sumstats_df=sumstats_df,
                tree=tree.time))
  })
  out_df=dplyr::bind_rows(lapply(out_list,function(list) list$sumstats_df))
  bootstraps_res[[pair]]<-list(data_stats=data_stats,out_df=out_df)
}

# #Plot the data summary stats and the bootstraps
# gather(out_df,-Pair,key="Stat",value = "N_coals")%>%
#   ggplot(aes(x=N_coals,fill=Stat))+
#   geom_histogram(binwidth = 1,col="gray",size=0.2)+
#   geom_vline(aes(xintercept = N_coals,col=Stat),linetype=2,data=data_stats)+
#   theme_classic()+
#   my_theme+
#   labs(x="Number of coalescences")

all_simstats_df<-lapply(bootstraps_res,function(list) {
  sim_stats<-gather(list$out_df,-Pair,key="Stat",value = "N_coals")%>%
    filter(!grepl("DR",Stat))%>%
    mutate(type="simulation")%>%
    mutate(Stat=gsub("_","\n",gsub("coals","coalescences",gsub("D","Donor",gsub("R","Recipient",gsub("peri","Peri-HCT",gsub("pre","Pre-HCT",Stat)))))))
  return(sim_stats)
})%>%bind_rows()%>%
  left_join(Pair_metadata)

all_datastats_df<-lapply(bootstraps_res,function(list) {
  data_stats<-list$data_stats%>%filter(!grepl("DR",Stat))%>%
    mutate(type="data")%>%
    mutate(Stat=gsub("_","\n",gsub("coals","coalescences",gsub("D","Donor",gsub("R","Recipient",gsub("peri","Peri-HCT",gsub("pre","Pre-HCT",Stat)))))))
  return(data_stats)
})%>%bind_rows()%>%
  left_join(Pair_metadata)

sumstat_confidence_intervals<-all_simstats_df%>%
  filter(!Pair%in%c("Pair41","Pair24"))%>%
  ggplot(aes(x=Stat,y=N_coals))+
  geom_boxplot(width=0.6,linewidth=0.4,outlier.shape = NA)+
  geom_jitter(alpha=0.1,size=0.5,width=0.1,height=0.5)+
  geom_point(aes(x= Stat,y=N_coals),col="red",shape=4,data=all_datastats_df%>%filter(!Pair%in%c("Pair41","Pair24")))+
  theme_classic()+
  facet_wrap(~Pair_new,nrow=2)+
  my_theme+
  labs(x="ABC summary statistic",y="Number of coalescences in bin")

ggsave(plot=sumstat_confidence_intervals,filename = paste0(plots_dir,"sumstat_confidence_intervals.pdf"),width=9,height=5)

get_internal_node_heights=function(tree,exclude_before=NULL) {
  nh<-nodeHeights(tree)
  h<-nh[!tree$edge[,2]%in%1:length(tree$tip.label),2]
  if(!is.null(exclude_before)) {
    h<-h[!h<exclude_before]
  }
  return(h)
}
bind_rows(data.frame(type="Donor",
                     heights=get_internal_node_heights(D_tree.time,exclude_before=2)),
          data.frame(type="Recipient",
                     heights=get_internal_node_heights(R_tree.time,exclude_before=2))
)%>%ggplot(aes(x=type,y=heights,col=type))+
  geom_jitter(height=0,width=0.1,alpha=0.6)+
  scale_y_continuous(limits = c(0,age_of_donor_at_sampling))+
  geom_hline(yintercept=c(age_of_donor_at_HSCT-5,age_of_donor_at_HSCT+5),linetype=2)+
  theme_classic()+
  my_theme+
  theme(legend.position = "none")

