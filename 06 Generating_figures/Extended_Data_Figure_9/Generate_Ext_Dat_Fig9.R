#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest","phangorn")
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
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
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
vcf_header_path=paste0(root_dir,"/data/reference_files/vcfHeader.txt") #A dummy vcf header that can be used

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
new_pair_names=Pair_metadata$Pair_new;names(new_pair_names)<-Pair_metadata$Pair

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
ABC.trees<-readRDS(paste0(root_dir,"/data/trees_for_ABC.Rds"))

#========================================#
# Generate 'Priors' plots for the figure ####
#========================================#

#Read in the posteriors from E. Mitchell et al, 2022 to use as the parameter distribution
param_posterior<-read.delim(paste0(root_dir,"/data/reference_files/posterior_sample.txt"))
these_params=colnames(param_posterior)
HSC_pop_posteriors<-read.delim(paste0(root_dir,"/data/reference_files/driver_parameter_posterior_sample.txt"),header = F)%>%
  dplyr::mutate("Parameter"="HSC_population_size",.before=1)%>%
  dplyr::rename("Value"=2)

#Rename the parameters
old_params=c("HSC_population_size",these_params)
labels=c("HSC population size","Number of drivers\nper year","Gamma distribution\nrate parameter","Gamma distribution\nshape parameter")
names(labels)<-old_params

#Plot of input parameter distributions into the model
input_params_plot<-param_posterior%>%
  pivot_longer(cols=all_of(these_params),names_to="Parameter",values_to = "Value")%>%
  bind_rows(HSC_pop_posteriors)%>%
  mutate(Parameter=labels[Parameter])%>%
  ggplot(aes(x=Value))+
  geom_density(fill="lightblue")+
  facet_wrap(~factor(Parameter,levels=labels),scales="free",nrow=1)+
  theme_bw()+
  my_theme

#Plot of the prior distribution of the size of the HSCT bottleneck (i.e. the number of engrafting cells)
log10_HSCT_bottleneck=runif(1e5,min=2.7,max=4.7)
Engrafting_HSCs_prior_plot<-data.frame(Parameter="Number of engrafting HSCs \n(Prior distribution)",Value=10^log10_HSCT_bottleneck)%>%
  ggplot(aes(x=Value))+
  geom_density(fill="plum1")+
  facet_wrap(~Parameter,scales="free",nrow=1)+
  scale_x_log10(limits=c(200,1e5),breaks=c(500,5000,50000),labels=scales::label_comma())+
  scale_y_continuous(limits=c(0,1))+
  theme_bw()+
  my_theme

comb_plot<-gridExtra::arrangeGrob(grobs=list(Engrafting_HSCs_prior_plot,input_params_plot),widths=c(1,3),nrow=1)
plot(comb_plot)

ggsave(filename=paste0(plots_dir,"ExtDatFig9_priors.pdf"),comb_plot,width=7,height=1.4)


#========================================#
# Generate 'Summary statistics' plots for the figure ####
#========================================#

pair="Pair21" #Pair21 (Pair_5) is used as the example in the plot
tree=ABC.trees[[pair]]
donor_id<-get_DR_ids(tree)['donor_ID']
recip_id<-get_DR_ids(tree)['recip_ID']

age_of_donor_at_HSCT<-Pair_metadata$Age_at_transplant[Pair_metadata$Pair==pair]
age_of_donor_at_sampling<-Pair_metadata$Age[Pair_metadata$Pair==pair]

D_tree<-keep.tip(tree,tip=grep(donor_id,tree$tip.label))
R_tree<-keep.tip(tree,tip=grep(recip_id,tree$tip.label))

largest_clades_D<-get_expanded_clade_nodes(D_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
  pull(n_samples)%>%sort(decreasing = T)%>%.[1:3]
largest_clades_R<-get_expanded_clade_nodes(R_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
  pull(n_samples)%>%sort(decreasing = T)%>%.[1:3]

#(2) Number of singletons
n_singletons_D<-get_expanded_clade_nodes(D_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
  dplyr::filter(n_samples==1)%>%nrow(.)
n_singletons_R<-get_expanded_clade_nodes(R_tree,height_cut_off = 50,min_clonal_fraction = 0.001,min_samples=1)%>%
  dplyr::filter(n_samples==1)%>%nrow(.)

#(3) Peri-HSCT LTT/ coalescences
tree.time<-tree
tree.time$edge.length<-tree$edge.length*(age_of_donor_at_sampling/median(get_mut_burden(tree)))
D_tree.time<-keep.tip(tree.time,tip=grep(donor_id,tree$tip.label))
R_tree.time<-keep.tip(tree.time,tip=grep(recip_id,tree$tip.label))
peri_HSCT_time_points=unlist(lapply(age_of_donor_at_HSCT,function(x) c(x-5,x+5)))

## These are to ILLUSTRATE THE SUMMARY STATISTICS USED IN THE ABC
#Bunch of functions to plot trees in a way to illustrate the summary statistics
highlight_nodes=function(tree,nodes,col="black",cex=0.5) {
  temp=sapply(nodes,function(node) {
    y=get_y_range(tree,node)
    x=get_x_range(tree,node)
    points(x=x[1],y=y[1],pch=19,col=col,cex=cex)
  })
  return(NULL)
}

highlight_nodes_in_window=function(tree,time_points,col="black",cex=0.3) {
  nodeheights <- nodeHeights(tree)
  nodes_in_window=tree$edge[,2][nodeheights[,2]>time_points[1] & nodeheights[,2]<time_points[2]]
  temp=highlight_nodes(tree = tree,nodes=nodes_in_window,col=col,cex=cex)
  cat(paste(length(nodes_in_window),"nodes within this window"),sep="\n")
  return(length(nodes_in_window))
}

annotate_coalescence_sumstats=function(tree,age_of_donor_at_HSCT,age_of_donor_at_sampling,title=NULL,labels=T,cols=c("black","blue")) {
  #Define the time point sets
  peri_HSCT_time_points=c(age_of_donor_at_HSCT-5,age_of_donor_at_HSCT+5)
  pre_HSCT_time_points=c(5,age_of_donor_at_HSCT-5)
  
  #Plot these
  tree=plot_tree(tree,cex.label=0,title=title,lwd=0.2)
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
        if(labels) {
          text(x=-15,
               y=(age_of_donor_at_sampling-time_points[2])+0.5*((age_of_donor_at_sampling-time_points[1])-(age_of_donor_at_sampling-time_points[2])),
               labels = label,srt=90)
        }
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
       col=rgb(0.1,0.1,0.1,alpha=0.1),border = "gray90")
  
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
  tree=plot_tree(tree,cex.label=0,title=title,lwd=0.2)
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
  tree=plot_tree(tree,cex.label=0,title=title,lwd = 0.2)
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

### Make the plot ----
pdf(file=paste0(plots_dir,"ExtDatFig9_sumstats.pdf"),width=10,height=2)
par(mfrow=c(1,4))
annotate_coalescence_sumstats(tree=R_tree.time,
                                age_of_donor_at_HSCT = age_of_donor_at_HSCT,
                                age_of_donor_at_sampling = age_of_donor_at_sampling,
                                title="Recipient only",
                              labels=F)

### Now plot the largest clades of each
highlight_largest_clades(tree=R_tree,n_clades=3,title="Recipient only")


### Plot the tree singletons
highlight_singletons(tree=R_tree,height_cut_off = 50,title="Recipient only",col="black")

### Maximum coalescences within single recipient clade
highlight_max_coals_in_single_clade_within_time_window(tree=R_tree.time,time_points=peri_HSCT_time_points,title="Recipient only")
dev.off()



