#Functions required for script
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

get_DR_ids=function(tree){
  tree=ape::drop.tip(tree,"Ancestral")
  PD_IDs=unique(substr(tree$tip.label,1,8))
  #Get number elements only
  PD_numbers=readr::parse_number(PD_IDs)
  names(PD_IDs)<-sapply(PD_numbers,function(n) ifelse(n%%2==0,"donor_ID","recip_ID"))
  return(PD_IDs)
}

dndscv_on_details=function(details,id="this_sample",outp=1,max_muts_per_gene_per_sample = Inf,max_coding_muts_per_sample = Inf,...) {
  require(dndscv)
  if(!"sampleID"%in%colnames(details)) {
    details$sampleID=id
  }
  muts=details[,c("sampleID","Chrom","Pos","Ref","Alt")]
  colnames(muts)<-c("sampleID","chr","pos","ref","alt")
  
  dndscvout=dndscv(muts,outp=outp,max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,max_coding_muts_per_sample = max_coding_muts_per_sample,...)
  return(dndscvout)
}

remove_samples_from_tree_and_update_details=function(remove_samples,tree,details){
  if(!any(remove_samples%in%tree$tip.label)){
    cat(paste("None of the specified samples (",paste0(remove_samples,collapse=","),") are in the tree. Returning original tree and details objects."),sep = "\n")
    stop(return(list(details=details,tree=tree)))
  } else {
    remove_samples<-remove_samples[which(remove_samples%in%tree$tip.label)]
    cat(paste("Sample",remove_samples,"found in tree and will be removed"),sep="\n")
  }
  
  sample_nodes=which(tree$tip.label%in%remove_samples)
  
  #Drop samples from the tree to create the new tree
  new_tree<-drop.tip(tree,tip=remove_samples)
  
  #Match up new nodes to the old nodes - create reference list, where each new node is named with the corresponding old node
  new_tree_clades=lapply(new_tree$edge[,2],function(node) getTips(new_tree,node))
  names(new_tree_clades)<-new_tree$edge[,2]
  
  old_nodes<-tree$edge[,2]
  new_nodes<-sapply(old_nodes,function(node) {
    old_clade_samples<-getTips(tree,node)
    if(all(old_clade_samples%in%remove_samples)) {
      return(NA)
    } else {
      new_clade_samples<-old_clade_samples[!old_clade_samples%in%remove_samples]
      new_idx<-which(sapply(new_tree_clades,function(clade) {setequal(clade,new_clade_samples)}))
      new_node<-new_tree$edge[,2][new_idx]
      return(new_node)
    }
  })
  names(new_nodes)<-old_nodes
  
  details$node<-new_nodes[as.character(details$node)]
  
  #Remove muts from branches that now have all tips in remove samples i.e. (removed branches)
  cat(paste(sum(is.na(details$node)),"mutations being removed from details matrix."),sep = "\n")
  details_new<-details[!is.na(details$node),]
  
  return(list(details=details_new,tree=new_tree))
}

##Now convert mutations to age
age_from_muts=function(n_muts,tree,sampling_age,birth_muts=50) {
  sampling_age_muts=mean(get_mut_burden(tree))
  if(n_muts<=birth_muts) {
    return(0)
  } else {
    age=sampling_age * (n_muts-birth_muts)/(sampling_age_muts-birth_muts)
    return(age)
  }
}
muts_from_age=function(age,tree,sampling_age,birth_muts=50) {
  sampling_age_muts=mean(get_mut_burden(drop.tip(tree,"Ancestral")))
  n_muts=birth_muts + (age*(sampling_age_muts-birth_muts)/sampling_age)
  return(n_muts)
}

get_mut_vafs=function(SampleID,COMB_mats,tree) {
  sample_nodes=get_ancestral_nodes(which(tree$tip.label==SampleID),edge = tree$edge)
  mutations=COMB_mats$mat$mut_ref[COMB_mats$mat$node%in%sample_nodes]
  NR=COMB_mats$NR
  NV=COMB_mats$NV
  rownames(NV)=rownames(NR)=COMB_mats$mat$mut_ref
  mutation_NR<-NR[mutations,SampleID]
  mutation_NV<-NV[mutations,SampleID]
  mutation_vaf<-calculate_vaf(mutation_NV,mutation_NR)
  return(data.frame(SampleID=SampleID,mut_ref=mutations,NV=mutation_NV,NR=mutation_NR,vaf=mutation_vaf))
}

get_duplicate_sets=function(tree,mut_threshold){
  pseudo_terminal_nodes=sapply(tree$edge[,2][!tree$edge[,2]%in%1:length(tree$tip.label)],function(node) {
    node_height=nodeheight(tree = tree,node=node)
    samples=getTips(tree=tree,node=node)
    sample_heights=nodeHeights(tree)[tree$edge[,2]%in%which(tree$tip.label%in%samples),2]
    
    if(all((sample_heights-node_height)<mut_threshold)){ #This is the method of determining the duplicates - may need to alter the "30" for low mutation burden samples
      return(node)
    }else{
      return(NA)
    }
  })
  pseudo_terminal_nodes=pseudo_terminal_nodes[!is.na(pseudo_terminal_nodes)]
  duplicate_samples=lapply(pseudo_terminal_nodes,function(node) getTips(tree,node=node))
  return(duplicate_samples)
}


get_DR_expanded_clades=function(tree,DorR,height_cut_off=50,min_samples,min_clonal_fraction=0.01){
  expanded_nodes=get_expanded_clade_nodes(tree=tree,height_cut_off = 50,min_samples = 2,min_clonal_fraction = 0)
  if(DorR=="D"){
    total_DR_samples=sum(grepl(get_DR_ids(tree)['donor_ID'],tree$tip.label))
  } else if(DorR=="R"){
    total_DR_samples=sum(grepl(get_DR_ids(tree)['recip_ID'],tree$tip.label))
  }
  DR_expanded_nodes=unlist(lapply(expanded_nodes$nodes,function(node) {
    node_samples=getTips(tree,node)
    if(DorR=="D"){
      DR_node_samples=sum(grepl(get_DR_ids(tree)['donor_ID'],node_samples))
    } else if(DorR=="R"){
      DR_node_samples=sum(grepl(get_DR_ids(tree)['recip_ID'],node_samples))
    }
    DR_clonal_frac=DR_node_samples/total_DR_samples
    if(DR_clonal_frac>=min_clonal_fraction & DR_node_samples>=min_samples){
      return(node)
    } else {
      return(NULL)
    }
  }))
  return(DR_expanded_nodes)
}

#Function to import the HDP data and convert into an exposures data frame
generate_exposures_df=function(HDP_multi_chain_RDS_path,trinuc_mut_mat_path,key_table_path){
  mut_example_multi=readRDS(HDP_multi_chain_RDS_path)
  mutations=read.table(trinuc_mut_mat_path)
  key_table=read.table(key_table_path)
  
  sample_remove=rownames(mutations)[rowSums(mutations)<50]
  mutations=mutations[!rownames(mutations)%in%sample_remove,]
  key_table=key_table[!key_table$Sample%in%sample_remove,]
  freq=nrow(mutations)
  
  dp_distn <- comp_dp_distn(mut_example_multi)
  ndp <- nrow(dp_distn$mean)
  ncomp <- ncol(dp_distn$mean)
  exposures <- t(dp_distn$mean[length(freq)+1+1:nrow(mutations),,drop=FALSE])
  colnames(exposures)=rownames(mutations)
  rownames(exposures)<-paste0("N",rownames(exposures))
  sigs=rownames(exposures)
  sig_profiles=mut_example_multi@comp_categ_distn$mean
  
  exposures_df<-as.data.frame(t(exposures),stringsAsFactors=F)%>%
    tibble::rownames_to_column("branch")%>%
    tidyr::separate(col="branch",into=c("node","exp_ID"),sep="_")%>%
    mutate(node=as.numeric(node))
  
  return(exposures_df)
}

mut_mat_HDP_comp=function(HDP_multi,ymax=0.2,plot=T){
  require(MutationalPatterns)
  sig_profiles=t(mut_example_multi@comp_categ_distn$mean)
  colnames(sig_profiles)<-paste0("N",0:(ncol(sig_profiles)-1))
  bases=c("A","C","G","T")
  subs=c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]")
  rownames(sig_profiles)<-paste0(rep(rep(bases,each=4),times=6),rep(subs,each=16),rep(bases,times=24))
  if(plot){
    plot_96_profile(sig_profiles,ymax=ymax,condensed = T) 
  }
  return(sig_profiles)
}


#Function to extract average signature contributions in each sample
#The tree should be the corrected, non-ultrametric tree, and the exposures df should have unique node numbers 
get_signatures_in_samples=function(tree,signature_names,exposures_df) {
  df<-dplyr::bind_cols(lapply(signature_names,function(signature) {
    sigs_in_samples=sapply(tree$tip.label,function(sample) {
      sample_node=which(tree$tip.label==sample)
      nodes_included=get_ancestral_nodes(node = sample_node,edge = tree$edge,exclude_root = T)
      branch_lengths=sapply(nodes_included,function(node) tree$edge.length[tree$edge[,2]==node])
      branch_sig_prop=sapply(nodes_included,function(node) {ifelse(node%in%exposures_df$node,sum(exposures_df[exposures_df$node==node,signature]),NA)})
      overall_contribution=weighted.mean(x=branch_sig_prop[!is.na(branch_sig_prop)],w = branch_lengths[!is.na(branch_sig_prop)])
      if(is.nan(overall_contribution)){overall_contribution<-0}
      return(overall_contribution)
    })
    return(sigs_in_samples)
  }))
  colnames(df)<-signature_names
  return(cbind(data.frame("Sample"=tree$tip.label),df))
}

plot_duplicate_plate_map=function(tree,sample_metadata,private_mut_threshold=NA,height_cut_off=NA,labels="set"){
  #Default mode is the 'look back from tips' threshold - if 'private_mut_threshold' is set, this is used
  if(!is.na(private_mut_threshold)){
    duplicate_clade_samples=get_duplicate_sets(tree,mut_threshold = private_mut_threshold)
  } else if(is.na(height_cut_off)){ #else looks for the 'height_cut_off' - and assumes
    duplicate_clades<-get_expanded_clade_nodes(tree,min_samples = 2,height_cut_off=height_cut_off,min_clonal_fraction = 0)
    duplicate_clade_samples<-lapply(duplicate_clades$nodes,function(node) getTips(tree,node=node))
  }
  if(length(duplicate_clade_samples)==0) {stop(print("No duplicates detected"))}
  dups_df<-Map(duplicate_clades=duplicate_clade_samples,dup_set=1:length(duplicate_clade_samples),f=function(duplicate_clades,dup_set){
    dup_set_metadata<-sample_metadata%>%
      dplyr::filter(Sample%in%duplicate_clades)%>%
      mutate(dup_set=dup_set)
    return(dup_set_metadata)
  })%>%
    dplyr::bind_rows()
  
  plot<-dups_df%>%
    mutate(column=factor(substr(position,1,1),levels=LETTERS[1:8]),
           row=factor(substr(position,2,3),levels=c("01","02","03","04","05","06","07","08","09","10","11","12")))%>%
    mutate(Sample=stringr::str_sub(Sample,start=10))%>%
    ggplot(aes(x=column,y=row))+
    geom_point(aes(col=factor(dup_set)),size=10,alpha=0.8)+
    facet_wrap(~plate_ID,drop = F)+
    scale_color_discrete(guide="none")+
    scale_x_discrete(drop=F)+
    scale_y_discrete(drop=F)+
    theme_bw()+
    labs(x="Plate column",y="Plate row",col="Duplicate set")
  if(labels=="set"){
    print(plot+geom_text(aes(label=dup_set),size=3))
  } else if(labels=="sample"){
    print(plot+geom_text(aes(label=Sample),size=3))
  }
}
