#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest")
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
if(!require("hdp", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/hdp", build_vignettes = F)
  library("hdp",character.only=T,quietly = T, warn.conflicts = F)
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
trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

#Generate information regarding loss-of-Y in male samples from X and Y coverage data
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

#Read in mutational signature extraction data
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)


#========================================#
# MUTATION BURDEN ANALYSIS ####
#========================================#
#Generate list of all samples in the final trees to include in this analysis
final_sample_list=unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label))

#Filter the metadata to include only these samples, and add in the Donor age information
mut_burden_df<-sample_metadata%>%
  dplyr::filter(Sample%in%final_sample_list)%>%
  mutate(Donor_age=sapply(ID,function(pair) as.numeric(Pair_metadata$Age[Pair_metadata$Pair==pair])[1]))

#LME regression to infer the age relationship & intercept for DONOR samples only
lme.D<-lmerTest::lmer(SNV_burden~Donor_age+(1|ID),data = mut_burden_df%>%dplyr::filter(DorR=="D"))
summary(lme.D)
confint(lme.D)

## Generate Extended Data Fig. 6a ----
donor_mut_burden_plot<-mut_burden_df%>%
  filter(DorR=="D")%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=Donor_age,y=SNV_burden,col=Pair_new))+
  geom_point(alpha=0.2,size=0.5)+
  scale_y_continuous(limits=c(0,1750),breaks=seq(0,1750,200))+
  scale_x_continuous(limits=c(0,81),breaks=seq(0,80,10))+
  geom_smooth(formula='y~x',col="black",method="lm",lwd=0.5,fullrange=T,se = F)+
  theme_classic()+
  my_theme+
  theme(legend.key.size=unit(3,"mm"))+
  scale_color_manual(values=Pair_cols)+
  guides(col=guide_legend(override.aes = list(alpha=1,size=1)))+
  annotate("text",size=2,x=20,y=1500,label=paste("y = ",round(lme.D@beta[1],0),"+ (",round(lme.D@beta[2],1),"x Donor age )"))+
  #annotate("text",size=2,x=20,y=1250,label=paste("R^2 = ",round(summary(lm.D)$r.squared,2)))+
  labs(x="Donor age (years)",y="SNV burden",col="")

ggsave(filename = paste0(plots_dir,"ExtDatFig6a.pdf"),donor_mut_burden_plot,width = 3,height=2)

## Generate Extended Data Fig. 6b ----
#Two features can alter estimates of differences between donors/ recipient in ways that we are not interested in
#1. OUTLIER SAMPLES
# Single high mutation burden samples may result from APOBEC activation or unusually numbers of in vitro mutations but
#2. CLONAL EXPANSIONS
# Clonal expansions have many samples for which the mutation burdens are non independent (e.g. a cell with a high mutation burden by chance may go on to produce a large high mutation burden clone)
# This is not accounted for unless the phylogenetic relationships are taken into account.
# A simple way to adjust for this is to exclude samples that are part of a clonal expansion, keeping only single samples from each clone
exclude_outliers=function(x,n_sd=2,return_idx=F,reverse=F) {
  mu<-mean(x)
  sd<-sd(x)
  idxs<-if(reverse) {which(x<(mu-n_sd*sd)|x>(mu+n_sd*sd))} else {which(x>(mu-n_sd*sd)&x<(mu+n_sd*sd))}
  if(return_idx) {
    return(idxs)
  } else {
    return(x[idxs])
  }
}

#This function returns any nodes later than the cutoff time that encompass ≥2 samples
get_expanded_clades=function(tree,cut_off=100) {
  nodeHeights<-nodeHeights(tree)
  res<-nodeHeights[,2]>cut_off & !tree$edge[,2]%in%1:length(tree$tip.label)
  return(tree$edge[,2][res])  
}

#Produce list of expanded clades. These can then be filtered for the analysis to see if this is distorting things i.e. test only singletons
expanded_clade_samples<-unlist(lapply(Pair_metadata$Pair,function(test_pair) {
  tree<-all.trees[[test_pair]]
  unique(unlist(lapply(get_expanded_clades(tree),function(node) getTips(tree,node)[-1]))) #Remove all but one sample from each clone
}))

#Produce list of mut burden outliers. These can then be filtered as may represent residual mixed colonies
outlier_samples<-unlist(lapply(Pair_metadata$Pair,function(test_pair) {
  Pair_df<-mut_burden_df%>%filter(ID==test_pair)
  Pair_outliers<-Pair_df$Sample[exclude_outliers(Pair_df$SNV_burden_adj1,n_sd=2.5,return_idx = T,reverse=T)]
}))

#LME regression to infer the age relationship & intercept
lme.combined.adj2<-lmerTest::lmer(SNV_burden_adj1~Donor_age+DorR+(1|ID),data = mut_burden_df%>%dplyr::filter(Sample%in%unlist(all.trees.cc.nodups) & !Sample%in%expanded_clade_samples))
summary(lme.combined.adj2)
confint(lme.combined.adj2)

#Function to compare donor vs recipient mutation burdens of all the pairs
#This applies a simple t-test
pval_comparisons=function(mut_burden_df,test_field="SNV_burden_adj",test_type="two.sided") {
  test_samples=unique(mut_burden_df$ID)
  results<-lapply(test_samples,function(test_pair) {
    donor_burdens<-mut_burden_df%>%
      filter(ID==test_pair&DorR=="D")%>%
      pull(!!as.symbol(test_field))
    recip_burdens<-mut_burden_df%>%
      filter(ID==test_pair&DorR=="R")%>%
      pull(!!as.symbol(test_field))
    res<-t.test(x=recip_burdens,y=donor_burdens,alternative=test_type)
    return(data.frame(pval=res$p.value,donor_mean=res$estimate[2],recip_mean=res$estimate[1],diff_mean=res$estimate[1]-res$estimate[2],diff_lower_CI=res$conf.int[1],diff_upper_CI=res$conf.int[2]))
  })
  return(cbind(data.frame(Pair=test_samples),dplyr::bind_rows(results))%>%tibble::remove_rownames())
}

#Apply the function - removing (1) all but one sample from each clonal expansion and (2) outlier samples
difference_df<-pval_comparisons(mut_burden_df%>%filter(!Sample%in%expanded_clade_samples & !Sample %in%outlier_samples),test_field="SNV_burden_adj1",test_type = "two.sided")

#Plot the results
D_vs_R_difference_CI_plot<-difference_df%>%
  left_join(Pair_metadata)%>%
  ggplot(aes(x=Pair_new,
             y=diff_mean,
             ymin=diff_lower_CI,
             ymax=diff_upper_CI))+
  geom_point()+
  geom_errorbar(width=0.25)+
  theme_classic()+
  my_theme+
  geom_hline(yintercept = 0)+
  theme(axis.line.x = element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_text(angle=90))+
  labs(x="",y="Increase in recipient mean\nmutation burden from HSCT\n(95% confidence interval)")

ggsave(filename = paste0(plots_dir,"ExtDatFig6b.pdf"),D_vs_R_difference_CI_plot,width =3.5,height=2)

#========================================#
# DRIVER ANALYSIS ####
#========================================#
# Create a driver dataframe with key driver info
driver_df<-Map(details=all.muts,pair=names(all.muts),function(details,pair) return(details%>%filter(Decision%in%c("Possible","Oncogenic"))%>%mutate(Pair=pair)))%>%bind_rows()
gene_order<-driver_df%>%group_by(Gene)%>%summarise(n=n())%>%arrange(desc(n))%>%pull(Gene)

## Generate Extended Data Fig. 6c ----
driver_gene_numbers<-driver_df%>%
  mutate(Gene=factor(Gene,levels=gene_order))%>%
  mutate(mut_type=ifelse(grepl("splice",Type),"Splice variant",ifelse(grepl("stop_gained|frameshift",Type),"Truncating variant",ifelse(grepl("inframe",Type),"Inframe deletion","Missense variant"))))%>%
  dplyr::count(Gene,mut_type)%>%
  ggplot(aes(y=n,x=Gene,fill=mut_type))+
  geom_bar(stat="identity",col="black",linewidth=0.2)+
  guides(fill=guide_legend(position="inside"))+
  theme_classic()+
  my_theme+
  theme(axis.text.x=element_text(angle=90,face = "italic"),legend.key.size = unit(3,"mm"),legend.position.inside = c(0.7,0.7))+
  labs(x="",y="# of unique mutations",fill="")
ggsave(filename=paste0(plots_dir,"ExtDatFig6c.pdf"),driver_gene_numbers,width=3,height=2)

## Generate Extended Data Fig. 6d ----
driver_gene_heatmap<-driver_df%>%
  dplyr::count(Pair,Gene)%>%
  tidyr::complete(Pair,Gene,fill=list(n=NA))%>%
  mutate(label=ifelse(!is.na(n),n,""))%>%
  left_join(Pair_metadata,by="Pair")%>%
  bind_rows(driver_df%>%
              mutate(Gene=factor(Gene,levels=gene_order))%>%
              dplyr::count(Pair,Gene)%>%
              tidyr::complete(Pair,Gene,fill=list(n=0))%>%
              mutate(label=ifelse(n>0,n,""))%>%
              left_join(Pair_metadata,by="Pair")%>%group_by(Pair_new)%>%summarise(Gene="Total",n=sum(n),label=as.character(sum(n))))%>%
  mutate(Pair_new=factor(Pair_new,levels=Pair_metadata%>%arrange(Age)%>%pull(Pair_new)))%>%
  mutate(Gene=factor(Gene,levels=c(gene_order,"Total")))%>%
  ggplot(aes(x=Gene,y=Pair_new,fill=n,label=label))+
  geom_tile()+
  scale_fill_gradientn(colours=c("white",brewer.pal(8,name="YlOrRd")),na.value = "gray90")+
  geom_text(size=1.5)+
  scale_x_discrete(position="top")+
  my_theme+
  theme(axis.text.x=element_text(angle=90,face="italic",hjust = 0),legend.position = "none",axis.line = element_blank())+
  labs(x="",y="")
ggsave(filename=paste0(plots_dir,"ExtDatFig6d.pdf"),driver_gene_heatmap,width=3,height=2)


#========================================#
# DRIVER TIMING ANALYSIS ####
#========================================#
# Create a driver dataframe with additional info about molecular time of acquisition
# Iterate through the Pairs, find the driver branches and find the molecular time of the top/ bottom of these branches
driver_timing_df<-Map(tree=all.trees.ultra,details=all.muts.nodups,pair=names(all.trees.ultra),function(tree,details,pair) {
  drivers_df<-details%>%filter(Decision%in%c("Oncogenic","Possible"))
  tree=plot_tree(tree,cex.label = 0)
  temp=add_annotation(tree=tree,details=details,annot_function=highlight_nodes,nodes=drivers_df$node)
  
  timing_df=lapply(1:nrow(drivers_df),function(i) {
    node=drivers_df$node[i]
    heights=nodeHeights(tree = tree)[which(tree$edge[,2]==node),]
    return(data.frame(mut_ref=drivers_df$mut_ref[i],min_molecular_time=heights[1],max_molecular_time=heights[2]))
  })%>%bind_rows()%>%
    right_join(drivers_df)%>%
    mutate(PairID=pair,.before=1)
  return(timing_df)
})%>%bind_rows()

## Generate Extended Data Fig. 6e ----
gene_cols=c("#aee39a", "#197959", "#4ad9e1", "#2f5672", "#24a5f7", "#6457d9", "#f3c5fa", "#8d2973", "#fd92fa", "#652cf6", "#e0079b", "#34f199", "#19a71f", "#bce333", "#6a9012", "#2cf52b", "#683d0d", "#fab899", "#ae301f", "#f47d0d", "#a27c59", "#f9bd3a", "#ff1c5d", "#f87574", "#ac82b4")
driver_mut_timing_plot<-driver_timing_df%>%
  mutate(Type=ifelse(Gene%in%gene_order[1:5],Gene,"Other"))%>%
  mutate(Type=factor(Type,levels=c(gene_order[1:5],"Other")))%>%
  mutate(Gene=factor(Gene,levels=gene_order))%>%
  arrange(Gene,max_molecular_time)%>%
  mutate(pos=1:nrow(.))%>%
  ggplot(aes(ymin=min_molecular_time,ymax=max_molecular_time,xmax=pos+0.4,xmin=pos-0.4,fill=factor(Gene)))+
  geom_rect()+
  theme_classic()+
  scale_fill_manual(values=gene_cols)+
  scale_y_reverse()+
  facet_grid(cols=vars(Type),scales="free",space="free")+
  my_theme+
  theme(strip.text.x = element_text(size=6,face="italic"),legend.text=element_text(face="italic"),axis.line.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.key.size = unit(2.5,"mm"))+
  labs(y="Molecular time of acquisition\n(# of mutations)",fill="")
ggsave(filename = paste0(plots_dir,"ExtDatFig6e.pdf"),driver_mut_timing_plot,width=7,height = 2.5)










