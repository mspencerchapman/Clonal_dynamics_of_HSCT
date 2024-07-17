#Read in the posteriors from E. Mitchell et al, 2022 to use as the parameter distribution
param_posterior<-read.delim(ifelse(Sys.info()["sysname"] == "Darwin",paste0(root_dir,"/data/reference_files/driver_parameter_posterior_sample.txt"),"/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/posterior_sample.txt"))
these_params=colnames(param_posterior)
HSC_pop_posteriors<-read.delim(ifelse(Sys.info()["sysname"] == "Darwin",paste0(root_dir,"/data/reference_files/HSC_population_posterior_sample.txt"),"/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/ABC_models/ABC_Apr2022/KX001_KX002_combined_Nt_posterior.txt"),header = T)%>%
  dplyr::mutate("Parameter"="HSC_population_size",.before=1)%>%
  dplyr::rename("Value"=2)

old_params=c("HSC_population_size",these_params)
labels=c("HSC population size","Number of drivers\nper year","Gamma distribution\nrate parameter","Gamma distribution\nshape parameter")
names(labels)<-old_params

#Plot of input parameter distributions into the model
input_params_plot<-param_posterior%>%
  tidyr::pivot_longer(cols=all_of(these_params),names_to="Parameter",values_to = "Value")%>%
  dplyr::bind_rows(HSC_pop_posteriors)%>%
  mutate(Parameter=labels[Parameter])%>%
  ggplot(aes(x=Value))+
  geom_density(fill="lightblue")+
  facet_wrap(~factor(Parameter,levels=labels),scales="free",nrow=1)+
  theme_bw()+
  my_theme

#Plot of the prior distribution of the size of the HSCT bottleneck (i.e. the number of engrafting cells)
log10_HSCT_bottleneck=runif(2e4,min=2.7,max=4.7)
Engrafting_HSCs_prior_plot<-data.frame(Parameter="Number of engrafting HSCs \n(Prior distribution)",Value=10^log10_HSCT_bottleneck)%>%
  ggplot(aes(x=Value))+
  geom_density(fill="plum1")+
  facet_wrap(~Parameter,scales="free",nrow=1)+
  scale_x_log10(limits=c(200,1e5),breaks=c(500,5000,50000),labels=scales::label_comma())+
  scale_y_continuous(limits=c(0,1))+
  theme_bw()+
  my_theme

library(gridExtra)
comb_plot<-gridExtra::arrangeGrob(grobs=list(Engrafting_HSCs_prior_plot,input_params_plot),widths=c(1,3),nrow=1)
plot(comb_plot)

ggsave(filename=paste0(plots_dir,"ABC_parameters.pdf"),comb_plot,width=7,height=1.4)

##ILLUSTRATE THE SUMMARY STATISTICS USED IN THE ABC
#Bunch of functions to plot trees in a way to illustrate the summary statistics
highlight_nodes=function(tree,nodes,col="black",cex=0.5) {
  temp=sapply(nodes,function(node) {
    y=get_y_range(tree,node)
    x=get_x_range(tree,node)
    points(x=x[1],y=y[1],pch=19,col=col,cex=cex)
  })
  return(NULL)
}

highlight_nodes_in_window=function(tree,time_points,col="black",cex=0.5) {
  nodeheights <- nodeHeights(tree)
  nodes_in_window=tree$edge[,2][nodeheights[,2]>time_points[1] & nodeheights[,2]<time_points[2]]
  temp=highlight_nodes(tree = tree,nodes=nodes_in_window,col=col,cex=cex)
  cat(paste(length(nodes_in_window),"nodes within this window"),sep="\n")
  return(length(nodes_in_window))
}

annotate_coalescence_sumstats=function(tree,age_of_donor_at_HSCT,age_of_donor_at_sampling,title=NULL,cols=c("black","blue")) {
  #Define the time point sets
  peri_HSCT_time_points=c(age_of_donor_at_HSCT-5,age_of_donor_at_HSCT+5)
  pre_HSCT_time_points=c(5,age_of_donor_at_HSCT-5)
  
  #Plot these
  tree=plot_tree(tree,cex.label=0,title=title)
  temp=Map(time_points=list(pre_HSCT_time_points,peri_HSCT_time_points),
      col=cols,
      label=c("Pre-transplant\nwindow","Peri-transplant\nwindow"),function(time_points,col,label) {
        #Peri-HSCT window
        rect(xleft = -1,
             xright=1+length(tree$tip.label),
             ybottom=age_of_donor_at_sampling-time_points[2],
             ytop=age_of_donor_at_sampling-time_points[1]-0.1,
             col=rgb(0.1,0.1,0.1,alpha=0.1),border = col)
        temp=highlight_nodes_in_window(tree=tree,time_points = time_points,col=col)
        text(x=-15,
             y=(age_of_donor_at_sampling-time_points[2])+0.5*((age_of_donor_at_sampling-time_points[1])-(age_of_donor_at_sampling-time_points[2])),
             labels = label,srt=90)
  })
  return(NULL)
}

highlight_max_coals_in_single_clade_within_time_window=function(tree,title=title,define_clone_height=5,time_points) {
  #Plot the tree
  tree=plot_tree(tree,cex.label=0,title=title)
  expanded_clones=get_expanded_clade_nodes(tree,height_cut_off = define_clone_height,min_samples = 2)
  nodes_within_time_points=tree$edge[,1][nodeHeights(tree)[,1]>time_points[1] & nodeHeights(tree)[,1]<time_points[2]]
  
  #Work out the relevant clade
  expanded_clones$coals_within_time_window_by_clone=sapply(expanded_clones$nodes,function(node) {
    clone_daughters<-get_all_node_children(node,tree)
    coals_within_time_points=intersect(clone_daughters,nodes_within_time_points)
    return(length(coals_within_time_points))
  })
  clone_with_max_window_coals<-expanded_clones$nodes[which.max(expanded_clones$coals_within_time_window_by_clone)]
  clone_daughters<-get_all_node_children(clone_with_max_window_coals,tree)
  coals_within_time_points=intersect(clone_daughters,nodes_within_time_points)
  cat(paste(length(coals_within_time_points)," is the maximum number of nodes within the time window in a single clade"),sep="\n")
  
  #Highlight the relevant clade
  highlight_clade(tree = tree,
                  node = clone_with_max_window_coals,
                  col="darkblue")
  
  #Highlight the time window
  rect(xleft = -1,
       xright=1+length(tree$tip.label),
       ybottom=age_of_donor_at_sampling-time_points[2],
       ytop=age_of_donor_at_sampling-time_points[1]-0.1,
       col=rgb(0.1,0.1,0.1,alpha=0.1),border = col)
  
  #Highlight the relevant coalescences
  highlight_nodes(tree = tree,nodes = coals_within_time_points,col = "darkorange")
  
}

highlight_clade=function(tree,node,col="blue") {
  all_children=get_all_node_children(node,tree)
  temp=sapply(all_children,function(node) {
    y=get_y_range(tree,node)
    x=get_x_range(tree,node)
    arrows(x0=x[1],y0=y[1],x1=x[1],y1=y[2],col=col,length=0)
    arrows(x0=x[1],y0=y[2],x1=x[2],y1=y[2],col=col,length=0)
  })
  return(NULL)
}

highlight_largest_clades=function(tree,n_clades,cols=c("blue","darkorange","magenta"),height_cut_off=50,title=NULL) {
  require(dplyr)
  tree=plot_tree(tree,cex.label=0,title=title)
  clades_df<-get_expanded_clade_nodes(tree,height_cut_off = height_cut_off,min_clonal_fraction = 0.001,min_samples=1)%>%
    arrange(desc(n_samples))%>%slice_head(n=n_clades)
  biggest_clade_sizes<-clades_df%>%pull(n_samples)
  cat(paste("Largest clade sizes:",biggest_clade_sizes,collapse=","))
  biggest_clade_nodes<-clades_df%>%pull(nodes)
  temp=Map(node=biggest_clade_nodes,col=cols[1:n_clades],function(node,col){
    temp=highlight_clade(tree=tree,node=node,col=col)
  })
}

highlight_singletons=function(tree,height_cut_off=50,col="blue",title=NULL) {
  tree=plot_tree(tree,cex.label=0,title=title)
  singleton_nodes<-get_expanded_clade_nodes(tree,height_cut_off = height_cut_off,min_clonal_fraction = 0.001,min_samples=1)%>%
    dplyr::filter(n_samples==1)%>%pull(nodes)
  cat(paste(length(singleton_nodes),"singletons"),sep="\n")
  temp=sapply(singleton_nodes,function(node) {
    y=get_y_range(tree,node)
    x=get_x_range(tree,node)
    arrows(x0=x[1],y0=y[1],x1=x[1],y1=y[2],col=col,length=0)
  })
  return(NULL)
}

### Plot the tree coalescences
pdf(file=paste0(plots_dir,"Sumstat_plots.pdf"),width=10,height=7)
par(mfrow=c(3,3))
Map(tree=list(D_tree.time,R_tree.time),title=c("Donor only","Recipient only"),function(tree,title) {
  annotate_coalescence_sumstats(tree=tree,
                                age_of_donor_at_HSCT = age_of_donor_at_HSCT,
                                age_of_donor_at_sampling = age_of_donor_at_sampling,
                                title=title)
})

### Now plot the largest clades of each
#par(mfrow=c(1,3))
Map(tree=list(D_tree,R_tree),title=c("Donor only","Recipient only"),function(tree,title) {
  highlight_largest_clades(tree=tree,n_clades=3,title=title)
})

### Plot the tree singletons
#par(mfrow=c(1,3))
Map(tree=list(D_tree,R_tree),title=c("Donor only","Recipient only"),function(tree,title) {
  highlight_singletons(tree=tree,height_cut_off = 50,title=title,col="black")
})

### Maximum coalescences within single recipient clade
highlight_max_coals_in_single_clade_within_time_window(tree=R_tree.time,time_points=peri_HSCT_time_points,title="Recipient only")
dev.off()



