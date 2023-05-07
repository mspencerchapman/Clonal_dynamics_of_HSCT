###TARGETED_SEQ_PLOTTING
##Contains my own functions for plotting the targeted seqeuncing data on the tree
##Mainly contains a custom function for plotting how far down the tree mutations shared by both donor & recipient go

add_var_col=function(tree, ##<< enhanced phylo returned from plot_tree
                     details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                     node,
                     var_field,
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
    bdat=details[[var_field]][info$idx]
    if(is.null(bdat) || class(bdat)!="numeric"){
      stop("Error in provided bfield (does it exist and is it numeric?)")
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

####
Gibbs_targ_seq_plots=function(SampleID,
                              tree,
                              details_targ,
                              pair_cell_fracs,
                              scale_muts_to_branch=TRUE,
                              colour.scale=c("#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#BD0026"),
                              log_min=-5,
                              vaf_lwd=5,
                              title=NULL) {
  
  posterior_cell_fracs.SampleID <-pair_cell_fracs[[SampleID]]
  
  get_median_cellfracs=function(post.df) {
    require(dplyr)
    median_cell_fracs=post.df%>%
      dplyr::select(-Node_assignment,-mutation_ID)%>%
      as.matrix()%>%
      apply(1,median)
    return(data.frame(mutation_ID=post.df$mutation_ID,median_cell_frac=median_cell_fracs))
  }
  
  details_targ_full<-left_join(details_targ,
                               get_median_cellfracs(post.df=posterior_cell_fracs.SampleID),
                               by=c("mut_ref"="mutation_ID"))
  
  #Generate the rescaled log cell fraction for plotting with contrast
  details_targ_full$log_median_cell_frac<-log(details_targ_full$median_cell_frac)
  details_targ_full$log_median_cell_frac[details_targ_full$log_median_cell_frac<log_min]<-log_min
  details_targ_full$log_median_cell_frac=plotrix::rescale(details_targ_full$log_median_cell_frac,newrange = c(0,1))
  
  ##Generate the plot
  tree=plot_tree(tree,cex.label=F)
  add_annotation(tree=tree,
                 details=details_targ_full,
                 matrices=NULL,
                 annot_function=function(tree,details,matrices,node) {
                   add_var_col(tree,
                               details,
                               node,
                               var_field = "log_median_cell_frac",
                               lwd = 3,
                               colours=colour.scale,
                               scale_muts_to_branch=scale_muts_to_branch)
                 }
  )
}

Gibbs_targ_seq_plots(SampleID="PD45792e",
                     tree=all.trees.cc.nodups[["Pair11"]],
                     details_targ=all_targeted_res[["Pair11"]]$details_targ,
                     pair_cell_fracs=posterior_cell_fracs[["Pair11"]],
                     scale_muts_to_branch=TRUE,
                     colour.scale=c("#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#BD0026"),
                     log_min=-5,
                     vaf_lwd=5,
                     title=NULL)

#Plot the VAFS for donor myeloid cells (monocytes, which are available for all patients)
donor_monocyte_ids<-bulk_smry_all%>%filter(individual_type=="Donor" & cell_type=="Monocytes" & time_point==0)%>%dplyr::select(Pair,tissueID)%>%filter(!duplicated(.))
pdf(paste0(plots_dir,"donor_targseq_trees.pdf"),width=15,height=10)
sapply(1:nrow(donor_monocyte_ids),function(i) {
  PairID=donor_monocyte_ids$Pair[i]
  Gibbs_targ_seq_plots(SampleID=donor_monocyte_ids$tissueID[i],
                       tree=all.trees.cc.nodups[[PairID]],
                       details_targ=all_targeted_res[[PairID]]$details_targ,
                       pair_cell_fracs=posterior_cell_fracs[[PairID]],
                       scale_muts_to_branch=TRUE,
                       colour.scale=c("#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#BD0026"),
                       log_min=-6,
                       vaf_lwd=5,
                       title=paste(PairID,"donor monocytes"))
})
dev.off()

plot_min_of_recip_or_donor=function(pair,
                                    tree,
                                    details_targ,
                                    Cell_type,
                                    log_min) {
  require(dplyr)
  donor_id=bulk_smry_all%>%filter(Pair==pair & time_point==0 & individual_type=="Donor" & cell_type==Cell_type)%>%pull(tissueID)%>%unique()
  recip_id=bulk_smry_all%>%filter(Pair==pair & time_point==0 & individual_type=="Recipient" & cell_type==Cell_type)%>%pull(tissueID)%>%unique()
  
  get_median_cellfracs=function(post.df) {
    require(dplyr)
    median_cell_fracs=post.df%>%
      dplyr::select(-Node_assignment,-mutation_ID)%>%
      as.matrix()%>%
      apply(1,median)
    return(median_cell_fracs)
  }
  
  pmin=pmin(get_median_cellfracs(posterior_cell_fracs[[pair]][[donor_id]]),
            get_median_cellfracs(posterior_cell_fracs[[pair]][[recip_id]]))
  
  min_median_cellfrac=data.frame(mutation_ID=posterior_cell_fracs[[pair]][[donor_id]]$mutation_ID,
                                 median_cell_frac=pmin)
  
  details_targ_full<-left_join(details_targ,
                               min_median_cellfrac,
                               by=c("mut_ref"="mutation_ID"))
  
  #Generate the rescaled log cell fraction for plotting with contrast
  details_targ_full$log_median_cell_frac<-log(details_targ_full$median_cell_frac)
  details_targ_full$log_median_cell_frac[details_targ_full$log_median_cell_frac<log_min]<-log_min
  details_targ_full$log_median_cell_frac=plotrix::rescale(details_targ_full$log_median_cell_frac,newrange = c(0,1))
  
  ##Generate the plot
  tree=plot_tree(tree,cex.label=F,plot_axis = F)
  add_annotation(tree=tree,
                 details=details_targ_full,
                 matrices=NULL,
                 annot_function=function(tree,details,matrices,node) {
                   add_var_col(tree,
                               details,
                               node,
                               var_field = "log_median_cell_frac",
                               lwd = 3,
                               colours=c("#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#BD0026"),
                               scale_muts_to_branch=T)
                 }
  )
}

ABC.trees<-readRDS(paste0(root_dir,"/data/trees_for_ABC.Rds"))

###Now plot the minimum values across pairs - consider anything of 10^-6 or less as absent (minimum sensitivity)
pdf(paste0(root_dir,"/plots/Minimum_targeted_fraction.pdf"),width=15,height=10)
temp=lapply(names(all_targeted_res),function(pair) {
  #pair="Pair28"
  pair_age_at_transplant<-Pair_metadata%>%dplyr::filter(Pair==pair)%>%pull(Age_at_transplant)
  pair_age<-Pair_metadata%>%dplyr::filter(Pair==pair)%>%pull(Age)
  
  #Set up the axis for plotting
  binwidth=ifelse(pair_age>50,10,5)
  max_age_to_plot=binwidth*ceiling(pair_age/binwidth)
  axis_at=c(0,sapply(seq(0,max_age_to_plot,by=binwidth),muts_from_age,ABC.trees[[pair]],pair_age))
  labels_at=c("Zyg.","Birth",seq(binwidth,max_age_to_plot,by=binwidth))
  samples_to_drop<-setdiff(all.trees.cc.nodups[[pair]]$tip.label,ABC.trees[[pair]]$tip.label)
  
  if(length(samples_to_drop)>0) {
    output=remove_samples_from_tree_and_update_details(remove_samples=samples_to_drop,
                                                       tree = all.trees.cc.nodups[[pair]],
                                                       details=all_targeted_res[[pair]]$details_targ)
    details_targ=output$details
  } else {
    #tree = all.trees.cc.nodups[[pair]]
    details_targ=all_targeted_res[[pair]]$details_targ
  }
  
  temp=plot_min_of_recip_or_donor(pair=pair,
                                  tree=ABC.trees[[pair]],
                                  details_targ=details_targ,
                                  Cell_type="Monocytes",
                                  log_min=-6)
  
  axis(side=4,at=mean(get_mut_burden(drop.tip(ABC.trees[[pair]],"Ancestral")))-axis_at,labels = labels_at,las=2,cex.axis=0.7)
  transplant_time_median=muts_from_age(pair_age_at_transplant,ABC.trees[[pair]],sampling_age=pair_age)
  CI_lower=transplant_time_median-2*sqrt(transplant_time_median)
  CI_upper=transplant_time_median+2*sqrt(transplant_time_median)
  rect(xleft = -1,
       xright=1+length(ABC.trees[[pair]]$tip.label),
       ybottom=mean(get_mut_burden(ABC.trees[[pair]]))-CI_upper,
       ytop=mean(get_mut_burden(ABC.trees[[pair]]))-CI_lower,
       col=rgb(0.1,0.1,0.1,alpha=0.3),border = NA)
  
  
})
dev.off()
