#Functions required for script
nodeHeights = function (tree, ...) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  nHeight <- function(tree) {
    tree <- reorder(tree)
    edge <- tree$edge
    el <- tree$edge.length
    res <- numeric(max(tree$edge))
    for (i in seq_len(nrow(edge))) res[edge[i, 2]] <- res[edge[i, 
                                                               1]] + el[i]
    res
  }
  nh <- nHeight(tree)
  return(matrix(nh[tree$edge], ncol = 2L) + ROOT)
}

nodeheight=function (tree, node, ...) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (node == (Ntip(tree) + 1)) 
    h <- 0
  else {
    a <- setdiff(c(getAncestors(tree, node), node), Ntip(tree) + 
                   1)
    h <- sum(tree$edge.length[sapply(a, function(x, e) which(e == 
                                                               x), e = tree$edge[, 2])])
  }
  h + ROOT
}

getAncestors=function (tree, node, type = c("all", "parent")) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  type <- type[1]
  if (type == "all") {
    aa <- vector()
    rt <- Ntip(tree) + 1
    currnode <- node
    while (currnode != rt) {
      currnode <- getAncestors(tree, currnode, "parent")
      aa <- c(aa, currnode)
    }
    return(aa)
  }
  else if (type == "parent") {
    aa <- tree$edge[which(tree$edge[, 2] == node), 1]
    return(aa)
  }
  else stop("do not recognize type")
}

getTips = function(tree,node) {
  require(ape)
  if(node <= length(tree$tip.label)) {
    daughters <- tree$tip.label[node]
  } else {
    daughters <- extract.clade(tree, node = node)$tip.label
  }
  return(daughters)
}

#Note - only apply the "get_edge_from_tree" option if starting from an SNV only tree. Function will assume that all the existing edge length is SNVs.
correct_edge_length = function(node, tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE, get_edge_from_tree=FALSE) {
  daughters <- getTips(tree = tree, node = node)
  #correct SNVs on edge, or set to 0 if want an indel only tree
  if(include_SNVs == TRUE) {
    if(get_edge_from_tree) {
      nSNV=tree$edge.length[tree$edge[,2]==node]
    } else {
      nSNV = sum(details$node == node & details$Mut_type == "SNV")
    }   		
    all_sens_SNVs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"SNV_sensitivity"]
    branch_SNV_sens = 1 - prod(1-all_sens_SNVs)  
    new_nSNV = nSNV/branch_SNV_sens
  } else {
    new_nSNV <- 0
  }
  #correct INDELs on edge, or set to 0 if want an SNV only tree
  if(include_indels == TRUE) {
    nINDEL = sum(details$node == node & details$Mut_type == "INDEL")
    all_sens_INDELs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"INDEL_sensitivity"]
    branch_INDEL_sens = 1 - prod(1-all_sens_INDELs)
    new_nINDEL = nINDEL/branch_INDEL_sens
  } else {
    new_nINDEL <- 0
  }
  new_edge_length = new_nSNV + new_nINDEL
  return(new_edge_length)
}

get_subset_tree = function(tree, details, v.field = "Mut_type", value = "SNV") {
  get_new_edge_length = function(node, tree, details,v.field,value) {
    sum(details$node == node & details[v.field] == value)
  }
  tree_subset = tree
  tree_subset$edge.length = sapply(tree$edge[,2], get_new_edge_length, tree = tree, details = details,v.field = v.field,value=value)
  return(tree_subset)
}

get_corrected_tree = function(tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE,get_edge_from_tree=FALSE) {
  tree_c = tree
  tree_c$edge.length = sapply(tree$edge[,2], correct_edge_length, tree = tree, details = details, sensitivity_df = sensitivity_df, include_indels = include_indels, include_SNVs=include_SNVs,get_edge_from_tree=get_edge_from_tree)
  return(tree_c)
}

get_mut_burden = function(tree) {
  mut_burden = nodeHeights(tree)[tree$edge[,2] %in% 1:length(tree$tip.label),2]
  return(mut_burden)
}

get_mut_burden_stats = function(tree) {
  mut_burden = get_mut_burden(tree)
  cat(paste("Mean mutation burden is", round(mean(mut_burden),digits = 1),"\n"))
  cat(paste("Range of mutation burden is", round(range(mut_burden)[1],digits = 1),"to",round(range(mut_burden)[2],digits = 1),"\n"))
  cat(paste("Standard deviation of mutation burden is", round(sd(mut_burden),digits = 1),"\n"))
}

#Function to calculate the absolute minimum number of clones by counting the number of times a parent node is shared, but
#a daughter node is recipient only (this would give you the number of extant transplanted clones if had full phylogeny)
get_minimum_clones=function(tree,donor_ID,recip_ID){
  shared_node_test=function(tree,node,donor_ID,recip_ID) {
    node_samples=getTips(tree,node)
    n_donor=sum(grepl(donor_ID,node_samples))
    n_recip=sum(grepl(recip_ID,node_samples))
    sharing_info=ifelse(n_donor>0&n_recip>0,"shared",ifelse(n_donor>0,"donor","recipient")) 
  }
  N=dim(tree$edge)[1]
  by_node=sapply(1:N,function(i) {
    node=tree$edge[i,2]
    sharing_info=shared_node_test(tree,node,donor_ID,recip_ID)
    if(sharing_info=="shared") {
      daughter_nodes=get_node_children(node,tree = tree)
      evidence_of_clone=0
      for(i in daughter_nodes) {
        daughter_sharing=shared_node_test(tree,node=i,donor_ID,recip_ID)
        if(daughter_sharing=="recipient") {
          evidence_of_clone=sum(1,evidence_of_clone)
        }
      }
    } else {
      evidence_of_clone=0
    }
    return(evidence_of_clone)
  })
  total_clones=sum(unlist(by_node))
  return(total_clones)
}

#Setup functions for AMOVA
# function to perform amova.
amova.fn <- function(distmat, groupnames, cell_key) {
  groupnums <- length(groupnames)
  dw <- c()
  cellnums <- c()
  for (i in 1:groupnums) {
    tgroup <- groupnames[i]
    tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type==tgroup])
    cellnums <- c(cellnums, length(tcells))
    dw <- c(dw, sum(distmat[tcells, tcells]))
  }
  tdist <- distmat[colnames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames], rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames]]
  dap <- (sum(tdist) - sum(dw))/2
  
  dfAP <- length(groupnames) - 1 
  dfWP <- sum(cellnums-1)
  N <- sum(cellnums)
  SSwp <- sum(dw/(cellnums*2))
  SSap <- sum(((dw + dap)/(2*N)) - (dw/(2*cellnums)))
  
  MSwp <- SSwp/dfWP
  MSap <- SSap/dfAP
  nc <- (N - (sum(cellnums^2)/N))/dfAP
  varwp <- MSwp
  varap <- (MSap - MSwp)/nc
  obsphi <- varap/(varwp + varap)
  return(obsphi)
}

# function to randomise sample labels and repeat
randamova.fn <- function(distmat, groupnames, cell_key) {
  
  # change added 2018.01.22: when randomizing, only include the part of the distance matrix that involves the cell types being considered.
  tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames])
  #
  
  randmat <- distmat[tcells, tcells]
  colnames(randmat) <- sample(colnames(randmat))
  rownames(randmat) <- colnames(randmat)
  randphi <- amova.fn(distmat=randmat, groupnames=groupnames, cell_key=cell_key)
  return(randphi)
}

# function tying it all in together
amovapval.fn <- function(distmat, groupnames, cell_key, iterations, plottitle) {
  # calculate observed
  obsphi <- amova.fn(distmat=distmat, groupnames=groupnames, cell_key=cell_key)
  # calculate null
  randphis <- sapply(1:iterations, function(cell) randamova.fn(distmat = distmat, groupnames=groupnames, cell_key=cell_key))
  # calculate pval
  pval <- length(which(randphis>obsphi))/length(randphis) 
  
  hist(randphis, col="grey", 100, main=plottitle, xlab="Phi statistic")
  abline(v=obsphi, col="red", lwd=2)
  legend("topright", legend=paste0("Observed\n p = ", signif(pval,digits=2)), lwd=2, col="red", bty="n")
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

set_cedge=function(parent,tree){
  for(child in get_node_children(parent,tree)){
    #cat("\nsetting parent - child ",parent,child,"\n")
    child.idx=which(tree$edge[,2]==child)
    parent.idx=which(tree$edge[,2]==parent)
    if(length(parent.idx)==0){
      pedge=0
    }else{
      pedge=tree$cedge[parent.idx]
    }
    #cat("before:",length(which(tree$cedge>0)),"\n")
    tree$cedge[child.idx]=pedge+tree$edge.length[child.idx];
    #cat("after:",length(which(tree$cedge>0)),"\n")
    #cat(parent,child,":",tree$cedge)
    tree=set_cedge(child,tree)
  }
  tree
}

##Not very efficient -- recursively calculates height as average of children's height.
get_height=function(tree,node){
  #if(is.null(tree$height)){
  #  tree$height=rep(NA,1:length(tree$edge.length))
  #}
  N=length(tree$tip.label)
  if(node<=N){
    return(node)#tree$height[which(tree$edge[,2]==node)]=node
  }else{
    children=get_node_children(node,tree)
    return(mean(sapply(children,function(child) get_height(tree,child))))
  }
}

set_height=function(tree){
  tree$height_end=sapply(tree$edge[,2],function(i) get_height(tree,i))
  tree$height_start=tree$height[match(tree$edge[,1],tree$edge[,2])]
  N=length(tree$tip.label)
  root=N+1
  idx=which(tree$edge[,1]==root)
  tree$height_start[idx]=mean(tree$height_end[idx])
  tree
}

##horizontal elbow
elbow=function(x0,x1,y0,y1,...){ 
  arrows(x0=x0,y0=y0,y1=y1,length=0,...)
  arrows(x0=x0,x1=x1,y0=y1,length=0,...)
}
##vertical elbow
elbowv=function(x0,x1,y0,y1,...){ 
  #browser()
  arrows(x0=x0,x1=x1,y0=y0,length=0,...)
  arrows(x0=x1,y0=y0,y1=y1,length=0,...)
  
}

set_tree_coords=function(atree){
  ##get the cumulative distance from the root.
  tree=atree
  tree$cedge=rep(0,length(tree$edge.length))
  N=length(tree$tip.label)
  root=N+1
  tree=set_cedge(root,tree)
  tt=set_height(tree)
  atree$coords=data.frame(a0=tt$cedge-tt$edge.length,a1=tt$cedge,
                          b0=tt$height_start,b1=tt$height_end,stringsAsFactors = FALSE)
  if(!is.null(atree$color)){
    atree$coords$color=atree$color
  }
  atree
}

plot_tree=function(tree,direction="down",cex.label=5,offset=0,plot_axis=T,title=NULL,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="darkgrey",ymax=NULL,cex.terminal.dots=0,vspace.reserve=0,cex.axis=1,cex.title=1){
  par(mar=c(1, 1, 1, 3) + 0.1)
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ymax=max(coords$a1)*1.05
    }
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  if(is.null(bars)){
    ymin=0-ymax*0.05-vspace.reserve*ymax
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(ymin,ymax),xlab="",ylab="")
  }else{
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=max(coords$a1)
    ##browser()
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
      text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,srt=270,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=c("darkgrey", "blueviolet","deeppink")[Y_loss], cex=cex.terminal.dots,pch=15)
    }
  }
  tree$direction=direction
  tree$top=top
  #scale =10
  scales=c(0,0.5,2,5,10,100,200,500,1000,2000,5000)
  scale=scales[max(which(ymax/4>scales))]
  #scale=scales[max(which(ymax/2>=scales))]
  #browser()
  cat("scale=",scale,"\n")
  if(plot_axis){
    axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2, cex.axis = cex.axis, lwd = 1, lwd.ticks = 1, col = "black") 
  }
  if(!is.null(title)) {
    text(x = 0,y=ymax,pos = 4,labels = title,cex=cex.title)
  }
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("Mutation number",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,border = "black",lwd=0.25,col="darkred")
    
  }
  tree$ymax=ymax
  tree$vspace.reserve=vspace.reserve
  tree
}

## Add heatmap with additional information under tree
add_heatmap=function(tree,heatmap,heatvals=NULL,border="white",cex.label=2){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.01*ymax
  gap=tree$vspace.reserve/dim(heatmap)[1]
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-0.05*ymax
    #bot=top-(0.05/dim(heatmap)[1])*ymax
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd = 0.25)
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      text(labels[i],x=-0.5,y=0.5*(top+bot),pos = 2,cex = cex.label)
    }
    top=bot
  }
  tree
}


plot_tree_old=function(tree,direction="down",cex.label=1,offset=0,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="darkgrey",ymax=NULL,cex.terminal.dots=0,plot_axis=TRUE){
  par(mar=c(1, 1, 1, 3) + 0.1)
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ymax=max(coords$a1)*1.05
    }
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  if(is.null(bars)){
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.05),ymax),xlab="",ylab="")
  }else{
    
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=max(coords$a1)
    ##browser()
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
      text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,srt=90,pos = 1,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=tipcol,cex=cex.terminal.dots,pch=19)
    }
  }
  tree$direction=direction
  tree$top=top
  scales=c(0,1,10,50,100,200,500,1000,2000,5000)
  #scale=scales[max(which(ymax/2>scales))]
  scale=scales[max(which(ymax/5>=scales))]
  #browser()
  cat("scale=",scale,"\n")
  if(plot_axis) {
    axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2)
  }
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("%s Muts",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,col = "grey")
    
  }
  tree
}


get_all_node_children=function(node,tree){
  children=tree$edge[which(tree$edge[,1]==node),2]
  offspring=children
  for(child in children){
    offspring=c(offspring,get_all_node_children(child,tree))
  }
  offspring
}
get_node_children=function(node,tree){
  tree$edge[which(tree$edge[,1]==node),2]
}

get_samples_in_clade=function(node,tree){
  if(node<=length(tree$tip.label)){
    return(tree$tip.label[node])
  }
  tree$tip.label[intersect(get_all_node_children(node,tree),1:length(tree$tip.label))]
}

get_y_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$top-tree$coords[idx,c("a1","a0")])
}

get_x_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$coords[idx,c("b1","b0")])
}

require("RColorBrewer")
get_idx_for_node=function(details,node){
  which(details$node==node)
}

get_edge_info=function(tree,details,node){
  y=get_y_range(tree,node)
  x=get_x_range(tree,node)
  idx=get_idx_for_node(details,node)
  samples=get_samples_in_clade(node,tree)
  list(yb=y[1],yt=y[2],x=x[1],xm=x[2],idx.in.details=idx,samples=samples)
}

add_annotation=function(tree,details=NULL,matrices=NULL,annot_function,...){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(tree,details,matrices,tree$edge[i,2],...))
}

add_binary_proportion=function(tree,##<< enhanced phylo returned from plot_tree
                               details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                               matrices,##<< a list of matrices parallel to details with columns named by tip labels
                               node,
                               bfield,
                               b.add.line=TRUE,
                               b.add.text=FALSE,
                               ...
){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  
  bdat=details[[bfield]][info$idx]
  if(is.null(bdat) || class(bdat)!="logical"){
    stop("Error in provided bfield (does it exist and is it boolean?)")
  }
  pass=sum(bdat,na.rm=TRUE)
  fail=sum(!bdat,na.rm=TRUE)
  tot=pass+fail
  ycross=info$yb+(fail/tot)*(info$yt-info$yb)
  ##Could add in a third category NA
  #missing=sum(is.na(bdat))
  if(b.add.line){
    arrows(y0=info$yb,y1=ycross,x0=info$x,length = 0,col="black",lend=1,...)
    arrows(y0=ycross,y1=info$yt,x0=info$x,length = 0,col="red",lend=1,...)
  }
  if(b.add.text){
    text(y=ycross,x=info$x,label=pass,pos = 4,offset = 0,...)
  }
}

add_simple_labels=function(tree,##<< enhanced phylo returned from plot_tree
                           details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                           matrices,##<< a list of matrices parallel to details with columns named by tip labels
                           node,##<< Node (see details)
                           query.field,##<< Name of column in details to query against
                           query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                           label.field,##<< Name of column in details specifying the label text.
                           cex.label=1,
                           b.add.label=TRUE,
                           b.add.marker=TRUE,
                           ... ##<< paremeters for points (not color)
){
  info=get_edge_info(tree,details,node)
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>50){
    stop("too many variants to annotate")
  }
  query.value=details[[query.field]][idx]
  idx.match=match(query.value,query.allowed.df$value)
  cols=query.allowed.df$col[idx.match]
  pch=query.allowed.df$pch[idx.match]
  
  vlabels=details[[label.field]][idx]
  ## spread out
  N=length(idx)
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  if(N>0){
    yd=info$yt-info$yb
    if(N==1){
      y=0.5*(info$yb+info$yt)
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    if(b.add.marker){
      points(rep(info$x,N),y,col=cols,pch=pch,...)
    }
    if(b.add.label){
      text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}

add_vaf=function(tree,##<< enhanced phylo returned from plot_tree
                 details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                 matrices,##<< a list of matrices parallel to details with columns named by tip labels
                 node,
                 samples=NULL,
                 b.plot.bars=TRUE,
                 lwd.rect=1,
                 min.depth=1,
                 vc.field,
                 vc.df,
                 filter.on=NULL,
                 ...
){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  #browser()
  ##Can do additional filter e.g. missense 
  #if(!is.null(filter.on)){
  #info$idx=info$idx,3)#info$idx.in.details[which(details$VC=="missense")]
  # }
  
  if(length(info$idx)==0){
    return(NULL)
  }
  if(b.plot.bars){
    plotF=plotBars
  }else{
    plotF=plotPie
  }
  if(is.null(samples)){
    samples=info$samples
  }
  cat(length(info$idx),"\n")
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
  }else{
    df=data.frame(mtr=matrices$mtr[info$idx,samples],
                  dep=matrices$dep[info$idx,samples],stringsAsFactors = FALSE)
  }
  df=cbind(df,details[info$idx,])
  df=df[which(df$dep>=min.depth),]
  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  N=dim(df)[1]
  if(N==0){
    return(df)
  }
  
  
  df=df[order(df$vaf),]
  yd=info$yt-info$yb
  
  
  if(N==1){
    y=0.5*(info$yb+info$yt)
    width=yd
  }else{
    y1=seq(info$yb,info$yt,length.out = N+2)
    #Don't use the ends..
    y=y1[2:(N+1)]
    width=y[2]-y[1]
  }
  
  if(!b.plot.bars){
    r=0.8  ##r>0.5 will cause likely overlap problems
  }else{
    r=0.4
  }
  #arrows(x0=info$x-w,x1=info$x-w+df$vaf*2*w,y0=y,lend=1,length=0,col="black")
  #arrows(x0=info$x-w+df$vaf*2*w,x1=info$x+w,y0=y,lend=1,length=0,col="grey")
  for(i in 1:N){
    vaf=min(df$vaf[i],0.999)
    if(is.na(vaf)){
      plotF(x=info$x,y = y[i],radius=r,col=c("lightgray","lightgray"),prop=c(0,1),border="lightgray",width=width)
    }else{
      plotF(x=info$x,y = y[i],radius = r,col=c("black","white"),prop = c(vaf,1-vaf),width = width)
    }
  }
  if( !b.plot.bars){
    
    return(df)
  }
  ##Now check to see if we need to highlight
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  min.mean.vaf=0.45
  z1=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
  z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
  max.p.value=max(z1$p.value,z2$p.value)
  txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
  if(z1$p.value>0.05 & z2$p.value<0.05) {
    border.color="green"
  } else if(max.p.value<0.05){
    if(max.p.value<0.05/dim(tree$edge)[1]){
      border.color="red"
    }else{
      border.color="blue"
    }
  }else{
    border.color="darkgrey"
  }
  
  
  rect(xleft=info$x-r,xright=info$x+r,ybottom=y[1]-width/2,ytop=y[N]+width/2,border=border.color,lwd=lwd.rect)
  if(border.color!="darkgrey"){
    text(txt,x=info$x,y=y[1]+0.3*(y[N]-y[1]),col="black",cex=0.6)
  }
  arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)
  
  #df
  
}


plotBars=function(x,y,radius,col,prop,border="black",width=1){
  #cat(prop,"\n")
  if(width<2){
    arrows(x0 = x-radius,y0=y,x1=x-radius+2*radius*prop[1],col="darkgrey",lend=2,length=0)
    arrows(x0 = x-radius+2*radius*prop[1],y0=y,x1=x-radius+2*radius,col=rgb(0.98,0.98,0.98),lend=2,length=0)
  }else{
    rect(xleft = x-radius,xright =x-radius+2*radius*prop[1],ybottom = y-width/2,ytop=y+width/2,border = NA,col="darkgrey")
    rect(xleft =  x-radius+2*radius*prop[1],xright =x+radius,ybottom = y-width/2,ytop=y+width/2,border = NA,col=rgb(0.98,0.98,0.98))
  }
  1
}

plotPie=function(x,y,radius,col,prop,llwd=0.5,border="black",width=NA){
  lims=par("usr")
  as=dev.size()
  asr=as[1]/as[2]
  yscale=asr*(lims[4]-lims[3])/(lims[2]-lims[1])
  prop=prop/sum(prop)
  cutpoint=c(0,cumsum(prop)*2*pi)
  
  N=2*pi/0.05
  n=ceiling(N*diff(cutpoint)/(2*pi))
  d=diff(cutpoint)/n
  if(length(prop)>1){
    for(i in 2:length(cutpoint)){
      polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              border=border,col=col[i-1],lwd=llwd)
    }
  }else{
    i=2
    polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            border=border,col=col[1],lwd=llwd)
  }
  yscale
}

plot_tree_vaf=function(tree,details,matrices,samples=NULL){
  tree=plot_tree(tree)
  res=add_annotation(tree,details,matrices,
                     function(tree,details,matrices,node){
                       add_vaf(tree,details,matrices,node,samples=samples,b.plot.bars = FALSE)
                     }
  )
  ##Post process res to add legend
  
}


plot_tree_labels_genes=function(tree,details,query.field="GENE",label.field="GENE",genes=c("JAK2","CBL","TET2","DNMT3A"),cex.label=1){
  qdf=data.frame(value=genes,col=rainbow(length(genes)),pch=19)
  plot_tree_labels(tree,details,
                   query.allowed.df = qdf,
                   query.field=query.field,
                   label.field=label.field,
                   cex.label=cex.label)
}



##Gets unique colour pch combos and returns in dataframe with columns "col" and "pch"
get_color_pch_df=function(n){
  pch.list=c(18,17,16,15,0:6)
  if(n>length(pch.list)*8){
    stop("Too many colours requested")
  }
  cols=rep(RColorBrewer::brewer.pal(8,"Set1"),times=length(pch.list))
  pch=rep(pch.list,each=8)
  data.frame(col=cols,pch=pch,stringsAsFactors = FALSE)[1:n,]
  
}

get_qdf=function(values){
  if(length(values)>length(unique(values))){
    stop("get_qdf: please provide values without duplication")
  }
  cbind(data.frame(value=values,stringsAsFactors = FALSE),
        get_color_pch_df(length(values)))
}

plot_tree_labels_consequence=function(tree,details,consequences,
                                      query.allowed.df=get_qdf(consequences),
                                      query.field="VC",
                                      label.field="GENE",
                                      cex.label=1){
  ##qdf=get_qdf(consequences)
  plot_tree_labels(tree,details,
                   query.allowed.df = query.allowed.df,
                   query.field=query.field,
                   label.field=label.field,
                   cex.label=cex.label)
}

#Useful function from online for drawing background boxes to the labels (then used in the "add_simple_labels_line" function)
boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA, 
                    border.bg = NA, adj = NULL, pos = NULL, offset = 0.5, 
                    padding = c(0.5, 0.5), cex = 1, font = graphics::par('font')){
  
  ## The Character expansion factro to be used:
  theCex <- graphics::par('cex')*cex
  
  ## Is y provided:
  if (missing(y)) y <- x
  
  ## Recycle coords if necessary:    
  if (length(x) != length(y)){
    lx <- length(x)
    ly <- length(y)
    if (lx > ly){
      y <- rep(y, ceiling(lx/ly))[1:lx]           
    } else {
      x <- rep(x, ceiling(ly/lx))[1:ly]
    }       
  }
  
  ## Width and height of text
  textHeight <- graphics::strheight(labels, cex = theCex, font = font)
  textWidth <- graphics::strwidth(labels, cex = theCex, font = font)
  
  ## Width of one character:
  charWidth <- graphics::strwidth("e", cex = theCex, font = font)
  
  ## Is 'adj' of length 1 or 2?
  if (!is.null(adj)){
    if (length(adj == 1)){
      adj <- c(adj[1], 0.5)            
    }        
  } else {
    adj <- c(0.5, 0.5)
  }
  
  ## Is 'pos' specified?
  if (!is.null(pos)){
    if (pos == 1){
      adj <- c(0.5, 1)
      offsetVec <- c(0, -offset*charWidth)
    } else if (pos == 2){
      adj <- c(1, 0.5)
      offsetVec <- c(-offset*charWidth, 0)
    } else if (pos == 3){
      adj <- c(0.5, 0)
      offsetVec <- c(0, offset*charWidth)
    } else if (pos == 4){
      adj <- c(0, 0.5)
      offsetVec <- c(offset*charWidth, 0)
    } else {
      stop('Invalid argument pos')
    }       
  } else {
    offsetVec <- c(0, 0)
  }
  
  ## Padding for boxes:
  if (length(padding) == 1){
    padding <- c(padding[1], padding[1])
  }
  
  ## Midpoints for text:
  xMid <- x + (-adj[1] + 1/2)*textWidth + offsetVec[1]
  yMid <- y + (-adj[2] + 1/2)*textHeight + offsetVec[2]
  
  ## Draw rectangles:
  rectWidth <- textWidth + 2*padding[1]*charWidth
  rectHeight <- textHeight + 2*padding[2]*charWidth    
  graphics::rect(xleft = xMid - rectWidth/2, 
                 ybottom = yMid - rectHeight/2, 
                 xright = xMid + rectWidth/2, 
                 ytop = yMid + rectHeight/2,
                 col = col.bg, border = border.bg)
  
  ## Place the text:
  graphics::text(xMid, yMid, labels, col = col.text, cex = theCex, font = font, 
                 adj = c(0.5, 0.5))    
  
  ## Return value:
  if (length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                yMid + rectHeight/2))
  } else {
    invisible(cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                    yMid + rectHeight/2))
  }    
}

add_simple_labels_line=function(tree,##<< enhanced phylo returned from plot_tree
                                details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                                matrices,##<< a list of matrices parallel to details with columns named by tip labels
                                node,##<< Node (see details)
                                query.field,##<< Name of column in details to query against
                                query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                                label.field,##<< Name of column in details specifying the label text.
                                cex.label=2,
                                b.add.label=TRUE,
                                b.add.marker=TRUE,
                                lty=1,
                                lwd=3,
                                ... ##<< paremeters for points (not color)
){
  info=get_edge_info(tree,details,node)
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>1){
    print("Some branches have multiple variants, plotting only the first variant of each branch")
    query.value=details[[query.field]][idx]
    vlabels=paste0(details[[label.field]][idx],collapse="\n")
  } else {
    query.value=details[[query.field]][idx]
    idx.match=match(query.value,query.allowed.df$value)
    cols=query.allowed.df$col[idx.match]
    vlabels=details[[label.field]][idx]
  }
  
  ## spread out
  N=ifelse(length(idx)>0,1,0)
  ##Vertical offset so that labels sit slightly above the markers.
  if(N>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col="red",lend=1,lwd=lwd,lty=lty,...)
    if(b.add.label){
      #text(rep(info$x,N),y=info$yb+0.5*(info$yt-info$yb),labels = vlabels,pos = 2,offset = 0.25,cex=cex.label)
      boxtext(info$x-1,info$yb+runif(1,min=0.35*(info$yt-info$yb),max=0.65*(info$yt-info$yb)),col.bg="white",border.bg="black",padding = c(0.5, 5),labels = vlabels,pos=2,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}


plot_tree_labels=function(tree,details,
                          query.field="VC",
                          type="label",
                          query.allowed.df=data.frame(value=c("nonsense","frameshift"),
                                                      col=c("red","black"),pch=c(17,18)
                          ),
                          label.field="GENE",
                          cex.label=1,
                          lty=1,
                          lwd=1){
  if(type=="label") {res=add_annotation(tree,
                                        details,list(),
                                        function(tree,details,matrices,node){
                                          add_simple_labels(tree,details,matrices,node,
                                                            query.field =query.field,
                                                            query.allowed.df = query.allowed.df,
                                                            label.field = label.field,
                                                            cex.label =cex.label)})
  with(query.allowed.df,legend("topleft",legend=value,col=col,pch=pch))
  }
  if(type=="line") {res=add_annotation(tree,
                                       details,list(),
                                       function(tree,details,matrices,node){
                                         add_simple_labels_line(tree,details,matrices,node,
                                                                query.field =query.field,
                                                                query.allowed.df = query.allowed.df,
                                                                label.field = label.field,
                                                                cex.label=cex.label,
                                                                lty=lty,
                                                                lwd=lwd)})}
}


plot_tree_vaf=function(tree,details,matrices,samples=NULL,b.plot.bars =TRUE,filter.on=NULL){
  res=add_annotation(tree,details,matrices,
                     function(tree,details,matrices,node){
                       add_vaf(tree,details,matrices,node,samples=samples,b.plot.bars = b.plot.bars,filter.on=filter.on)})
}

#add_vaf_bar
#add_vof_pie
#add_label

#add_var_col function which colours each mutation line according to a numeric variable (scaled between 0 and 1)
library(dichromat)
add_var_col=function(tree, ##<< enhanced phylo returned from plot_tree
                     details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                     matrices,##<< a list of matrices parallel to details with columns named by tip labels
                     node,
                     var_field,
                     pval_based=FALSE,
                     b.add.line=TRUE,
                     colours = c("black","green","red"),
                     scale_muts_to_branch=TRUE,
                     ...){
  
  #Define the col.scale from the colours vector
  require(dichromat)
  colfunc = colorRampPalette(colours)
  col.scale = colfunc(101)
  
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  muts_on_edge=length(info$idx.in.details)
  edge_length=tree$edge.length[tree$edge[,2]==node]
  
  if(muts_on_edge > 0 & edge_length>0) {
    if(var_field == "vaf") {
      NV_vec = rowSums(matrices$mtr[info$idx.in.details,info$samples,drop =FALSE])
      NR_vec = rowSums(matrices$dep[info$idx.in.details,info$samples,drop = FALSE])
      chroms = unlist(lapply(strsplit(names(NV_vec),split = "-"),function(x) return(x[1])))
      if(pval_based) {
        bdat=sapply(1:length(info$idx.in.details), function(i) binom.test(NV_vec[i],NR_vec[i],p = ifelse(chroms[i] %in% c("X","Y"),0.95,0.5), alternative = "two.sided")$p.value)
      } else {
        bdat = NV_vec/NR_vec
        bdat[chroms %in% c("X","Y")] <- bdat[chroms %in% c("X","Y")]/2
      }
    } else {
      bdat=details[[var_field]][info$idx]
      if(is.null(bdat) || class(bdat)!="numeric"){
        stop("Error in provided bfield (does it exist and is it numeric?)")
      }
    }
    bdat = sort(bdat, decreasing = TRUE)
    if(scale_muts_to_branch) {
      mut_unit_of_edge=edge_length/muts_on_edge
    } else {
      mut_unit_of_edge=1
    }
    ##Could add in a third category NA
    #missing=sum(is.na(bdat))
    if(b.add.line){
      y0_next = info$yt
      for(i in 1:muts_on_edge) {
        arrows(y0=y0_next,y1=(y0_next - mut_unit_of_edge),x0=info$x,length = 0,col=col.scale[ceiling(100*bdat[i])],lend=1,...)
        y0_next = y0_next - mut_unit_of_edge
      }
    }
  }
}



highlight_nodes=function(tree,details,matrices,node,nodes,...) {
  info=get_edge_info(tree,details,node=node)
  if(node %in% nodes){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col="red",lend=1,...)
  }
}

#The "plot_node_number" function
plot_node_number = function(tree,details,matrices,node,cex=0.4) {
  info=get_edge_info(tree,details,node)
  text(info$x,info$yb,node,cex = cex,col="black",font=2)
}

#Plot the tip point colour according to whether is donor or recipient
plot_d_or_r_tip_point = function(sample,tree,details,donor_ID,recip_ID,cols=c("dark green","red")) {
  node=which(tree$tip.label==sample)
  info=get_edge_info(tree,details,node)
  tip_col=ifelse(grepl(donor_ID,sample),cols[1],cols[2])
  points(x=info$x,y=info$yb,type="p",pch=20,bg=tip_col,col=tip_col)
}

plot_category_tip_point = function(sample_ID,tree,details=NULL,cat_df,cat_name="cat",cols=RColorBrewer::brewer.pal(8,"Set1"),col="black",...) {
  all_categories<-cat_df%>%pull(cat_name)%>%unique()
  
  if(!all(all_categories%in%names(cols))) {
    cols=cols[1:length(cat_df%>%pull(cat_name)%>%unique())]
    names(cols)<-sort(all_categories)
  }
  
  node=which(tree$tip.label==sample_ID)
  info=get_edge_info(tree,details,node)
  tip_col=cols[cat_df%>%filter(sample==sample_ID)%>%pull(cat_name)]
  points(x=info$x,y=info$yb,type="p",pch=21,bg=tip_col,col=col,...)
}

plot_postGT_tree=function(tree,details,matrices,node,highlight="post",sharing_cols=c("black","gray92"),cat_df){  #sharing_cols is a vector of colours for "shared", "donor only" and "recipient only" branches.
  info=get_edge_info(tree,details,node=node)
  if(highlight=="pre"){
    n=sum(cat_df%>%dplyr::filter(Sample%in%info$samples)%>%pull(Time_point)==0)
  } else if(highlight=="post"){
    n=sum(cat_df%>%dplyr::filter(Sample%in%info$samples)%>%pull(Time_point)>0)
  }
  sharing_info=ifelse(n>0,"highlight","lowlight")
  names(sharing_cols)=c("highlight","lowlight")
  #if(length(tree$edge.length[tree$edge[,2]==node])>0){
  arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=sharing_cols[sharing_info],lend=1)
  #}
}

plot_sharing_info=function(tree,details,matrices,node,donor_ID,recip_ID,sharing_cols=c("black","dark green","red"),...){  #sharing_cols is a vector of colours for "shared", "donor only" and "recipient only" branches.
  info=get_edge_info(tree,details,node=node)
  n_donor=sum(grepl(donor_ID,info$samples))
  n_recip=sum(grepl(recip_ID,info$samples))
  sharing_info=ifelse(n_donor>0&n_recip>0,"shared",ifelse(n_donor>0,"donor","recipient"))
  names(sharing_cols)=c("shared","donor","recipient")
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=sharing_cols[sharing_info],lend=1,lwd=ifelse(sharing_info=="shared",0.8,0.4),...)
  }
}

plot_sharing_multiple=function(tree,details,matrices,node,sharing_cols=c("black","dark green","red","orange","brown","green"),...){  #sharing_cols is a vector of colours for "shared", "donor only" and "recipient only" branches.
  categories=unique(tree$tip.label)
  if(length(categories)>5) {stop("Too many tip label categories")}
  info=get_edge_info(tree,details,node=node)
  if(length(unique(info$samples))>1) {
    sharing_info <- "shared"
  } else {
    sharing_info <- unique(info$samples)
  }
  sharing_cols=sharing_cols[1:(1+length(categories))]  	
  names(sharing_cols)=c("shared",categories)
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=sharing_cols[sharing_info],lend=1,...)
  }
}

plot_mut_vaf_by_branch=function(tree,
                                details,
                                matrices,
                                node,
                                mut,
                                colours=c("black","green","red"),
                                cex=0.4,
                                #show_pval=FALSE,
                                ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  colfunc = colorRampPalette(colours)
  col.scale = colfunc(101)
  
  #Get the allocated node number
  allocated_node=details$node[details$mut_ref==mut]
  expected_samples=get_edge_info(tree,details,allocated_node)$samples
  
  #Get the vaf
  info=get_edge_info(tree,details,node=node)
  samples=info$samples
  variant_reads=sum(matrices$NV[mut,samples])
  total_reads=sum(matrices$NR[mut,samples])
  vaf=variant_reads/total_reads
  
  #Plot the vaf on the tree using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=col.scale[ceiling(100*vaf)],lend=1,...)
  }
  
  #Print the total depth for that branch at the node
  if(variant_reads>0|node %in% which(tree$tip.label %in% expected_samples))
    text(info$x,info$yb,paste0(variant_reads,"/",total_reads),srt=90,cex = cex,col="black",font=2)
  
  #Print the mutation name and log10(pvalue)
  text(1,1,pos=4,mut)
  #if(show_pval&"pval"%in%colnames(details)){text(1,50,pos=4,paste0("Log10 p-value for allocated node:",round(log10(details$pval[details$mut_ref==mut]))))} 
}

plot_MAV_mut=function(tree,
                      details,
                      matrices,
                      node,
                      lesion_node=NA,
                      mut1,
                      mut2=NULL,
                      colours=c("dark gray","red","blue"),
                      cex=0.4,
                      #show_pval=FALSE,
                      ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  mut_colfunc = colorRampPalette(colours[-1])
  mut_colscale = mut_colfunc(11)
  
  #Get the tips within the lesion node
  if(!is.na(lesion_node)) {
    expected_samples=getTips(tree,lesion_node)
  }
  
  #Get the vaf
  info=get_edge_info(tree,details,node=node)
  samples=info$samples
  mut1_reads=sum(matrices$NV[mut1,samples])
  if(is.null(mut2)){mut2_reads=0}else{mut2_reads=sum(matrices$NV[mut2,samples])}
  variant_reads=mut1_reads+mut2_reads
  total_reads=variant_reads+(sum(matrices$NR[mut1,samples]) - mut1_reads)
  wt_reads=total_reads-variant_reads
  
  if((mut1_reads+mut2_reads)!=0) {
    if(is.null(mut2)) {
      base_col=colours[2]
    } else {
      base_col=mut_colscale[1+round(10*mut1_reads/(mut1_reads+mut2_reads),digits=0)] #How red or blue should the "mut" element of the colour be
    }
    final_colfunc=colorRampPalette(c("light gray",base_col))
    final_colscale=final_colfunc(101)
    branch_col=final_colscale[1+round(100*(mut1_reads+mut2_reads)/total_reads)] #Now how "concentrated" should the mut colour be
  } else {
    branch_col="light gray"
  }  
  
  #Plot the branches using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=branch_col,lend=1,...)
  }
  
  #Print the total depth for that branch at the node
  if(!is.na(lesion_node)){
    if(node %in% which(tree$tip.label %in% expected_samples)) {
      text(info$x,info$yb,paste0(variant_reads,"/",total_reads),srt=90,cex = cex,col="black",font=2)
    }
  }
  
  #Print the mutation name
  if(node==1) { #Do for node==1 so that only prints the name once
    if(is.null(mut2)) {
      text(1,1,pos=4,mut1)
    } else {
      text(1,1,pos=4,paste0(mut1,"/",strsplit(x=mut2,split="-")[[1]][4]))
    }
  }
}

confirm_PVV_phylogeny=function(tree,
                               details,
                               matrices,
                               node,
                               mut,
                               PVV_mut=NULL,
                               lesion_node,
                               colours=c("black","green","red"),
                               cex=0.4,
                               #show_pval=FALSE,
                               ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  colfunc = colorRampPalette(colours)
  col.scale = colfunc(101)
  
  #Get the samples within the lesion node
  lesion_node_samples=getTips(tree,lesion_node)
  
  #Get the vaf
  info=get_edge_info(tree,details,node=node)
  samples=info$samples
  variant_reads=sum(matrices$NV[mut,samples])
  total_reads=sum(matrices$NR[mut,samples])
  vaf=variant_reads/total_reads
  
  #Plot the vaf on the tree using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=col.scale[ceiling(100*vaf)],lend=1,...)
  }
  
  #Print the total depth for that branch at the node
  if(variant_reads>0|node %in% which(tree$tip.label %in% lesion_node_samples))
    text(info$x,info$yb,paste0(variant_reads,"/",total_reads),srt=90,cex = cex,col="black",font=2)
  
  #Print the mutation name and log10(pvalue)
  if(node==1) {text(1,1,pos=4,paste("Confirmation of phylogeny for",PVV_mut,":", mut))}
  #if(show_pval&"pval"%in%colnames(details)){text(1,50,pos=4,paste0("Log10 p-value for allocated node:",round(log10(details$pval[details$mut_ref==mut]))))} 
}


###FUNCTIONS TO GO WITH THE APE "PLOT.PHYLO" FUNCTION

highlight_groups=function(tree,group1,group2,cols=c("#17698E","#17A258")) {
  any_group=c(group1,group2) #make combined list of samples assigned to at least one group
  
  #Iterate through all the nodes & work if (1) all the tips are in group1 (2) all are in group2 (3) mixture of both
  col_vec=sapply(tree$edge[,2],function(node) {
    node_daughters=getTips(tree,node)
    if(any(node_daughters%in%any_group)) {
      node_daughters<-node_daughters[node_daughters%in%any_group] #only include those tips that are assigned to a specific group
      if(all(node_daughters%in%group1)) {
        col<-cols[1]
      } else if(all(node_daughters%in%group2)) {
        col<-cols[2]
      } else {
        col<-"black" #If has tips in both groups, colour black
      }
    } else {
      col<-"lightgrey"
    }
  })
  return(col_vec)
}

highlight_samples=function(tree,samples) {
  tips=which(tree$tip.label%in%samples)
  edge_cols=sapply(tree$edge[,2],function(node) ifelse(node%in%tips,"red","black"))
  return(edge_cols)
}


##FUNCTION FOR THE SIMULATED TREES
drivers_per_sample=function(tree){
  if(is.null(tree$events)){stop(print("Need tree with events matrix recording driver information."))}
  driver_nodes=tree$events$node[tree$events$value==1 & tree$events$driverid>0]
  n_drivers<-sapply(tree$tip.label[-1],function(tip) { #Don't include the 'ancestral tip'
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    return(sum(driver_nodes%in%tip_nodes))
  })
  driver_ids<-sapply(tree$tip.label[-1],function(tip) { #Don't include the 'ancestral tip'
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    driver_ids<-tree$events%>%dplyr::filter(driverid!=0 & node %in% tip_nodes)%>%pull(driverid)
    return(paste0(driver_ids,collapse=","))
  })
  return(data.frame(Sample=tree$tip.label[-1],n_drivers=n_drivers,driver_ids=ifelse(nchar(driver_ids)==0,NA,driver_ids)))
}

##FUNCTION FOR THE DATA TREES
drivers_per_sample_data=function(tree,details){
  require(dplyr)
  if(is.null(details$is.driver)){stop(print("Need variable 'is.driver' in the details matrix"))}
  driver_nodes=details%>%dplyr::filter(is.driver==1)%>%pull(node)
  n_drivers<-sapply(tree$tip.label,function(tip) {
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    return(sum(driver_nodes%in%tip_nodes))
  })
  driver_ids<-sapply(tree$tip.label,function(tip) {
    tip_nodes<-get_ancestral_nodes(which(tree$tip.label==tip),edge = tree$edge)
    driver_ids<-details%>%dplyr::filter(node %in% tip_nodes & is.driver==1)%>%pull(variant_ID)
    return(paste0(driver_ids,collapse=","))
  })
  out_df<-data.frame(Sample=tree$tip.label,n_drivers=n_drivers,driver_ids=ifelse(nchar(driver_ids)==0,NA,driver_ids))
  out_df<-out_df%>%dplyr::filter(Sample!="Ancestral")
  return(out_df)
}

find.distance <- function(tree, from, to) {
  path <- nodepath(tree, from, to)
  res <- 0
  for (j in 2:length(path)) {
    index <- which(tree$edge[,2] == path[j] & tree$edge[,1] == path[j-1], arr.ind = TRUE)
    res <- res + tree$edge.length[index]
  }
  return(res)
}

length.normalise <- function(orig.tree, new.tree, curr.node, remaining.stick) {
  curr.node.children <- unlist(Descendants(orig.tree, curr.node, "children"))
  
  for (j in curr.node.children) {
    index <- which(orig.tree$edge[,2] == j & orig.tree$edge[,1] == curr.node, arr.ind = TRUE)
    
    if (j %in% orig.tree$tip.label) {
      new.tree$edge.length[index] <- remaining.stick
    } else {
      curr.node.tips <- unlist(Descendants(orig.tree, j, "tips"))
      curr.dist <- find.distance(orig.tree, curr.node, j)
      if (curr.dist == 0) {curr.dist <- 0.01} # So that no edge lengths are zero
      desc.lengths <- sapply(curr.node.tips, FUN = find.distance, tree = orig.tree, from = curr.node)
      new.tree$edge.length[index] <- remaining.stick * curr.dist / mean(desc.lengths)
      shorter.stick <- remaining.stick - new.tree$edge.length[index]
      
      # Call as recursive function
      new.tree <- length.normalise(orig.tree, new.tree, j, shorter.stick)
    }
  }
  return(new.tree)
} 


make.ultrametric.tree <- function(tree) {
  root.number <- length(tree$tip.label) + 1
  ultra.tree <- length.normalise(tree, tree, root.number, 1)
  return(ultra.tree)
}

generate.bespoke.plots <- function(tree) {
  # Generate different versions of ultrametric tree
  tree.bespoke <- make.ultrametric.tree(tree)
  
  plot(tree.bespoke, show.tip.label = FALSE)
  
  # Generate population size trajectories
  tree.BNPR.bespoke <- BNPR(tree.bespoke)
  plot_BNPR(tree.BNPR.bespoke)
  
}

# Throw away function to compare methods
generate.diagnostic.plots <- function(tree) {
  #Make the tree dichotomous - required for phylodyn
  tree<- multi2di(tree)
  
  # Generate different versions of ultrametric tree
  tree.nnls <- force.ultrametric(tree, method="nnls")
  tree.extend <- force.ultrametric(tree, method="extend")
  tree.bespoke <- make.ultrametric.tree(tree)
  
  par(mfrow=c(2,2))
  plot(tree, show.tip.label = FALSE, main="Unadjusted")
  plot(tree.nnls, show.tip.label = FALSE, main="NNLS")
  plot(tree.extend, show.tip.label = FALSE, main="Extend")
  plot(tree.bespoke, show.tip.label = FALSE, main="Bespoke")
  
  # Generate population size trajectories
  tree.BNPR <- BNPR(tree)
  tree.BNPR.nnls <- BNPR(tree.nnls)
  tree.BNPR.extend <- BNPR(tree.extend)
  tree.BNPR.bespoke <- BNPR(tree.bespoke)
  
  plot_BNPR(tree.BNPR, col = "red", main="Population trajectory: unadjusted", xlab = "Time", ylab = "Relative population size")
  plot_BNPR(tree.BNPR.nnls, col = "red", main="Population trajectory: NNLS method", xlab = "Time", ylab = "Relative population size")
  plot_BNPR(tree.BNPR.extend, col = "red", main="Population trajectory: extend method", xlab = "Time", ylab = "Relative population size")
  plot_BNPR(tree.BNPR.bespoke, col = "red", main="Population trajectory: custom method", xlab = "Time", ylab = "Relative population size")
  
}

designUltra <- function(tree, sparse = TRUE) {
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorder(tree, "postorder")
  leri <- allChildren(tree)
  bp <- bip(tree)
  n <- length(tree$tip.label)
  l <- tree$Nnode
  nodes <- integer(l)
  k <- 1L
  u <- numeric(n * (n - 1) / 2)
  v <- numeric(n * (n - 1) / 2)
  m <- 1L
  for (i in seq_along(leri)) {
    if (length(leri[[i]]) > 1) {
      if (length(leri[[i]]) == 2) ind <- getIndex(bp[[leri[[i]][1] ]],
                                                  bp[[leri[[i]][2] ]], n)
      else {
        ind <- NULL
        le <- leri[[i]]
        nl <- length(le)
        for (j in 1:(nl - 1)) ind <- c(ind, getIndex(bp[[le[j] ]],
                                                     unlist(bp[ le[(j + 1):nl] ]), n))
      }
      li <- length(ind)
      v[m:(m + li - 1)] <- k
      u[m:(m + li - 1)] <- ind
      nodes[k] <- i
      m <- m + li
      k <- k + 1L
    }
  }
  if (sparse) X <- sparseMatrix(i = u, j = v, x = 2L)
  else {
    X <- matrix(0L, n * (n - 1) / 2, l)
    X[cbind(u, v)] <- 2L
  }
  colnames(X) <- nodes
  attr(X, "nodes") <- nodes
  X
}

find_latest_acquisition_node=function(tree,pos_samples){
  #Get list of ancestral nodes for all samples
  ancestral_nodes_list=lapply(pos_samples,function(Sample) {
    get_ancestral_nodes(node = which(tree$tip.label==Sample),edge=tree$edge,exclude_root = F)
  })
  #Find nodes that are ancestral to all the samples
  common_nodes=Reduce(intersect,ancestral_nodes_list)
  #Which of these is the most recent (i.e. has the maximum node height)
  nodeheights<-sapply(common_nodes,function(node) nodeheight(tree = tree,node = node))
  MRCA_node<-common_nodes[which.max(nodeheights)]
  return(MRCA_node)
}

#Nick's "get_ancestral_nodes" function. Needed for the validation functions above
get_ancestral_nodes= function(node,edge,exclude_root=TRUE){
  idx=which(edge[,2]==node)
  parents=node ##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=edge[idx,1]
    parents=c(parents,parent)
    #This finds the parent of the current parent - thus navigating up to the root.
    idx=which(edge[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)] ##The last node is the root.
  }else{
    parents
  }
}


#The "squash tree" function to cut the tree at any given node height.  Tree structure is maintained, but edge lengths are shortened.
squash_tree=function(tree,cut_off=50) {
  tree$edge.length[nodeHeights(tree)[,1]>=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
  idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
  new_edge_lengths=cut_off - nodeHeights(tree)[idxs_to_squash,1] #work-out the edge lengths that these should be such that they finish at the cut-off
  tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
  return(tree)
}
