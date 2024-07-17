#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest","pheatmap","MASS")
bioconductor_packages=c("clusterProfiler","GenomicRanges","IRanges","Rsamtools","BSgenome","BSgenome.Hsapiens.UCSC.hg19","MutationalPatterns")

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
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome,character.only=TRUE)
options(stringsAsFactors = FALSE)

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

mut_sigs_theme=theme(strip.text.x = element_text(size=6,margin = margin(0.6,0,0.6,0, "mm")),
                     strip.text.y=element_text(size=6,margin = margin(0,0.6,0,0.6, "mm")),
                     axis.text.x = element_text(size=3),
                     axis.text.y=element_text(size=5),
                     axis.title.x=element_text(size=7),
                     axis.title.y=element_text(size=7),
                     axis.ticks.x=element_blank(),
                     axis.ticks.y=element_line(size=0.25),
                     legend.text = element_text(size=5),
                     legend.title = element_text(size=7),
                     legend.key.size=unit(2.5,"mm"),
                     strip.background = element_rect(linewidth =0.25),
                     panel.grid.major = element_line(linewidth=0.25),
                     panel.border = element_rect(linewidth=0.25))

#========================================#
# Set the root directory and read in the necessary files ####
#========================================#

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
genome_file=ifelse(Sys.info()['sysname']=="Darwin","~/R_work/reference_files/genome.fa","/nfs/cancer_ref02/human/GRCh37d5/genome.fa")
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
trees_list<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/details_lists.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

#Generate information regarding loss-of-Y in male samples from X and Y coverage data
LOY_files=list.files(path=paste0(root_dir,"/data/SV_and_CNA_data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv(paste0(root_dir,"/data/SV_and_CNA_data/Copy_number_changes.csv"))

#Read in mutational signature extraction data
exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)


#========================================#
# WRITE VCFS OF APOBEC MUTATIONS####
#========================================#
#VCFs are used for testing for Extended Data Fig. 5c

## Find branches with >10% of mutations assigned to APOBEC (N4) ----
APOBEC_branches<-exposures_df%>%
  filter(N4>0.1)

## For each branch, get all branch mutations, annotate the trinucleotide context, and select only those from the specific APOBEC peaks ----
temp=lapply(1:nrow(APOBEC_branches),function(i) {
  pair<-APOBEC_branches$Pair[i]
  APOBEC_node<-APOBEC_branches$node[i]
  mutations<-all.muts[[pair]]%>%filter(node==APOBEC_node)%>%dplyr::select(Chrom,Pos,Ref,Alt)
  
  colnames(mutations) = c("chr","pos","ref","mut")
  mutations$pos=as.numeric(mutations$pos)
  mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
  mutations$trinuc_ref = as.vector(scanFa(genome_file, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                      as.numeric(mutations$pos)+1))))
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  mutations_APOBEC_only<-mutations%>%filter(sub=="C>T"&(trinuc_ref_py=="TCA"|trinuc_ref_py=="TCT"))
  mutations_APOBEC_only%>%dplyr::select(chr,pos,ref,mut)%>%dplyr::rename(Chrom=chr,Pos=pos,Ref=ref,Alt=mut)
  
  write.vcf(mutations_APOBEC_only%>%dplyr::select(chr,pos,ref,mut)%>%dplyr::rename(Chrom=chr,Pos=pos,Ref=ref,Alt=mut),vcf_path = paste0(root_dir,"/data/APOBEC_VCFs/APOBEC_",pair,"_",APOBEC_node,".vcf"),vcf_header_path = vcf_header_path)
})

#========================================#
# LOOK AT TIMING OF BRANCHES WITH APOBEC ####
#========================================#
# Update the exposure df with absolute mutation numbers of APOBEC &
# the minimum & maximum molecular time of the branch

exposures_df_abs<-Map(tree=all.trees,pair=names(all.trees),function(tree,pair) {
  
  #Make ultrametric versions of these trees
  mean_mutation_burden=mean(get_mut_burden(tree))
  tree.ultra<-make.ultrametric.tree(tree)
  tree.ultra$edge.length=tree.ultra$edge.length*mean_mutation_burden
  tree.ultra$edge.length[which(tree$edge[,2]==which(tree$tip.label=="Ancestral"))]<-0
  
  #Update the exposures df to reflect absolute mutation numbers
  exposures_df_ind<-exposures_df%>%filter(Pair==pair)
  exposures_df_ind_abs=dplyr::bind_rows(lapply(1:nrow(exposures_df_ind),function(i) {
    branch_length=tree$edge.length[tree$edge[,2]==exposures_df_ind$node[i]]
    x=exposures_df_ind[i,]
    x[,3:7]<-x[,3:7]*branch_length
    return(x)
  }))
  
  #
  tree.ultra.heights<-nodeHeights(tree = tree.ultra)
  exposures_df_ind_abs<-bind_cols(exposures_df_ind_abs,dplyr::bind_rows(lapply(exposures_df_ind_abs$node,function(node) {
    as.data.frame(t(tree.ultra.heights[tree.ultra$edge[,2]==node,]))
  })))
  
  exposures_df_ind_abs<-exposures_df_ind_abs%>%
    dplyr::rename("Min_height"=V1,"Max_height"=V2)
  
  #
  return(exposures_df_ind_abs)
})%>%dplyr::bind_rows()

#Estimate the molecular times of the HCT
get_molecular_time_of_HCT=function(mol_time_at_sampling,
                                   age_of_donor_at_sampling,
                                   age_of_donor_at_HCT,
                                   mutations_at_birth=60) {
  mutations_per_year=(mol_time_at_sampling-mutations_at_birth)/age_of_donor_at_sampling
  return(mutations_at_birth+(mutations_per_year*age_of_donor_at_HCT))
}

#Add this info to the Pair metadata dataframe
Pair_metadata$mol_time_at_HCT=sapply(Pair_metadata$Pair,function(pair) {
  get_molecular_time_of_HCT(mol_time_at_sampling = exposures_df_abs%>%filter(Pair==pair)%>%pull(Max_height)%>%max(),
                            age_of_donor_at_sampling=Pair_metadata%>%filter(Pair==pair)%>%pull(Age),
                            age_of_donor_at_HCT=Pair_metadata%>%filter(Pair==pair)%>%pull(Age_at_transplant),
  )
})

## Generate Extended Data Fig. 5a ---- 
APOBEC_timing_plot<-exposures_df_abs%>%
  filter(N4>10)%>%
  mutate(idx=row_number())%>%
  left_join(Pair_metadata)%>%
  ggplot(aes(xmin=idx-0.4,xmax=idx+0.4,ymin=Min_height,ymax=Max_height))+
  geom_rect(aes(fill=Pair_new))+
  scale_fill_manual(values=Pair_cols)+
  facet_grid(~Pair_new,scales = "free",space="free")+
  geom_hline(aes(yintercept=mol_time_at_HCT),data=Pair_metadata%>%filter(!Pair_new%in%c("Pair_4","Pair_8")),linetype=2)+
  theme_classic()+
  my_theme+
  labs(x="Branches with >10 APOBEC mutations",y="Molecular time")+
  theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x=element_blank(),strip.text.x = element_text(angle=90))

ggsave(filename = paste0(plots_dir,"ExtDatFig5a.pdf"),APOBEC_timing_plot,width=6,height=3)

#========================================#
# ADD SOME ADDITIONAL SAMPLE-LEVEL METADATA ####
#========================================#

# Enhance the metadata table with (1) if the sample contains a driver, (2) the sample is part of a clonal expansion, (3) whether the recipient received TBI
temp=Map(pair=names(all.trees.ultra),details=all.muts.nodups,tree=all.trees.ultra,function(pair,details,tree) {
  individual_metadata=sample_metadata%>%filter(ID==pair & sample_status=="PASS")
  
  #Add the "expansion" or "no expansion" info
  individual_metadata$clonal_expansion=sapply(individual_metadata$Sample,function(sample) {
    sample_anc=getAncestors(tree,which(tree$tip.label==sample),type="parent")
    anc_height=nodeheight(tree,sample_anc)
    ifelse(anc_height>100,"YES","NO")
  })
  
  #Add the "driver" or "no driver" info
  driver_nodes<-details%>%filter(Decision=="Oncogenic"|Decision=="Possible")%>%pull(node)
  driver_samples<-unique(unlist(lapply(driver_nodes,function(node) getTips(tree,node))))
  individual_metadata$driver=sapply(individual_metadata$Sample,function(sample) {
    ifelse(sample%in%driver_samples,"YES","NO")
  })
  
  return(individual_metadata)
})%>%dplyr::bind_rows()

Pair_metadata$TBI=c("NO","YES","NO","NO","NO","NO","NO","YES","YES","YES")
dat_full<-left_join(temp%>%filter(sample_status=="PASS"),Pair_metadata%>%dplyr::select(-Pair),by="Pair_new")


#========================================#
# RUN LINEAR MIXED EFFECTS MODELS FOR FACTORS THAT MAY INFLUENCE APOBEC ACTIVATION ####
#========================================#

#Test using a linear mixed effects model - define APOBEC+ as having more than 15 mutations assigned to APOBEC signature
lme1=lmerTest::lmer(N4_abs>15~DorR+(1|ID),data=dat_full)
summary(lme1)
confint(lme1)

#Now add in whether the sample is part of a clonal expansion (with DoR as a random effect)
lme2.1=lmerTest::lmer(N4_abs>15~clonal_expansion+(1|ID)+(1|DorR),data=dat_full)
summary(lme2.1)
confint(lme2.1)

#Now add in whether the sample has a driver mutation (with DoR as a random effect)
lme2.2=lmerTest::lmer(N4_abs>15~driver+(1|ID)+(1|DorR),data=dat_full)
summary(lme2.2)
confint(lme2.2)

#Now assess whether there is a link with conditioning type (with DoR as a random effect)
lme3.1=lmerTest::lmer(N4_abs>15~conditioning+(1|ID)+(1|DorR),data=dat_full)
summary(lme3.1)
confint(lme3.1)

#Now assess whether there is a link with stem cell source (with DoR as a random effect)
lme3.2=lmerTest::lmer(N4_abs>15~stem_cell_source+(1|ID)+(1|DorR),data=dat_full)
summary(lme3.2)
confint(lme3.2)

#Now assess whether there is a link with TBI exposure (with DoR as a random effect)
lme3.3=lmerTest::lmer(N4_abs>15~TBI+(1|ID)+(1|DorR),data=dat_full)
summary(lme3.3)
confint(lme3.3)

#Now assess whether there is a link with donor sex (with DoR as a random effect)
lme3.4=lmerTest::lmer(N4_abs>15~donor_sex+(1|ID)+(1|DorR),data=dat_full)
summary(lme3.4)
confint(lme3.4)

# Now combine lme models into a list & name them, for plotting purposes --------
all_lmes=list(lme1,lme2.1,lme2.2,lme3.1,lme3.2,lme3.3,lme3.4)
lme_test=c("Donor vs Recipient",
           "Clonal expansion vs Singleton",
           "Driver vs No driver",
           "MAC vs RIC conditioning",
           "BM vs PBSC stem cell source",
           "TBI vs No TBI",
           "Male vs female donor")

# Summarise the individual models with confidence intervals and p-values using the lmerTest package ----
summary=bind_rows(lapply(all_lmes,function(lme) {
  CI_table=confint(lme)
  idx=which(!grepl(".sig",rownames(CI_table)) & !grepl("(Intercept)",rownames(CI_table)))
  return(CI_table[idx,])
  }))
summary$test=lme_test

summary2=bind_rows(lapply(all_lmes,function(lme) {
  lme_coefs=summary(lme)$coefficients
  idx=which(!grepl(".sig",rownames(lme_coefs)) & !grepl("(Intercept)",rownames(lme_coefs)))
  return(lme_coefs[idx,])
}))
summary2$test=lme_test

## Generate Extended Data Fig. 5b ---- 
pvalue_plot<-summary2%>%
  ggplot(aes(y=`Pr(>|t|)`,x=factor(test,levels = rev(lme_test))))+
  geom_point(size=0.75)+
  scale_y_log10()+
  geom_hline(yintercept=0.05,linetype=2)+
  theme_classic()+
  my_theme+
  labs(x="Comparison",y="P value")+
  coord_flip()

ggsave(filename = paste0(plots_dir,"ExtDatFig5b.pdf"),pvalue_plot,width=3,height=1.5)


