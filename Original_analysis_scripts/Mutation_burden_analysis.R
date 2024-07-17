#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
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

#----------------------------------
# Set the ggplot2 theme for plotting
#----------------------------------

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

#----------------------------------
# Set the root directory and read in the necessary files
#----------------------------------

root_dir<-"~/R_work/Clonal_dynamics_of_HSCT"
source(paste0(root_dir,"/data/HSCT_functions.R"))
plots_dir=paste0(root_dir,"/plots/")
HDP_folder=paste0(root_dir,"/data/HDP")

#Read in Pair metadata data frame
Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))

#Define colour themes for the Pairs & DorR
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired"); names(Pair_cols)<-levels(Pair_metadata$Pair_new)
DorR_cols<-RColorBrewer::brewer.pal(8,"Dark2")[1:2]; names(DorR_cols)<-c("D","R")

#Read in other data objects
sample_metadata<-readRDS(paste0(root_dir,"/data/metadata_files/sample_metadata_full.Rds"))
trees_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/tree_lists.Rds"))
details_list<-readRDS(paste0(root_dir,"/data/trees_and_muts_files/details_lists.Rds"))
LOY_files=list.files(path=paste0(root_dir,"/data/LOY_files"),pattern="meanCoverage",full.names = T)
male_PDIDs<-c("PD45792","PD45793","PD45794","PD45795")
Y_loss_df=dplyr::bind_rows(lapply(LOY_files,read.delim))%>%
  mutate(donor=substr(id,1,7))%>%
  mutate(loss_of_Y=ifelse(!donor%in%male_PDIDs,NA,ifelse(y/x<0.15,"YES","NO")))

#Read in the spreadsheet listing other copy number changes
CN_change_df=read.csv(paste0(root_dir,"/data/Copy_number_changes.csv"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_list)) {assign(x,details_list[[x]])}

#----------------------------------
# COVERAGE ANALYSIS
#----------------------------------

Coverage_stats<-sample_metadata%>%
  filter(sample_status=="PASS")%>%
  group_by(Pair_new)%>%
  dplyr::summarise(mean_cov=mean(Coverage),median_cov=median(Coverage),sd_cov=sd(Coverage))

Coverage_histograms<-sample_metadata%>%
  filter(!is.na(Coverage))%>%
  filter(sample_status=="PASS")%>%
  ggplot(aes(x=Coverage))+
  geom_histogram(fill="lightblue",col="black",size=0.3)+
  geom_vline(data=Coverage_stats,aes(xintercept=mean_cov),col="red",linetype=2)+
  geom_text(data=Coverage_stats,aes(label=paste("µ =",round(mean_cov,digits=1),"x")),x=20,y=50,size=2,col="red")+
  #geom_vline(xintercept=4,col="black",linetype=1)+
  #geom_rect(xmin=0,xmax=4,ymin=0,ymax=80,fill="7a8baf#30")+
  theme_classic()+
  my_theme+
  facet_wrap(~Pair_new,nrow=2)+
  theme(strip.text.x=element_text(margin = unit(c(0.6,1,0.6,1),"mm")))

ggsave(filename=paste0(plots_dir,"Coverage_histograms.pdf"),Coverage_histograms,width=6,height=2)

#----------------------------------
# MUTATION BURDEN ANALYSIS
#----------------------------------

final_sample_list=unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label))
mut_burden_df<-sample_metadata%>%
  dplyr::filter(Sample%in%final_sample_list)%>%
  mutate(Donor_age=sapply(ID,function(pair) as.numeric(Pair_metadata$Age[Pair_metadata$Pair==pair])[1]))

#LME regression to infer the age relationship & intercept for DONOR samples only
lme.D<-lmerTest::lmer(SNV_burden~Donor_age+(1|ID),data = mut_burden_df%>%dplyr::filter(DorR=="D"))
summary(lme.D)
confint(lme.D)


#LME regression to infer the age relationship & intercept
lme.combined<-lmerTest::lmer(SNV_burden~Donor_age+DorR+(1|ID),data = mut_burden_df)
summary(lme.combined)
confint(lme.combined)

#LME regression to infer the age relationship & intercept
lme.combined.adj<-lmerTest::lmer(SNV_burden_adj1~Donor_age+DorR+(1|ID),data = mut_burden_df)
summary(lme.combined.adj)
confint(lme.combined.adj)

#Now two further adjustments
#1. exclude outlier samples
#2. exclude samples that are part of a clonal expansion (non-independent measures)
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
  unique(unlist(lapply(get_expanded_clades(tree),function(node) getTips(tree,node)[-1])))
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

#Plot + linear regression of donor mutation burdens
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
  labs(x="Donor age",y="SNV burden",col="Pair")

ggsave(filename = paste0(plots_dir,"donor_mut_burden_plot.pdf"),donor_mut_burden_plot,width = 3,height=2)

#Violin plot of corrected SNV burdens: donor vs recipient
D_vs_R_mut_burden_comparison_plot<-mut_burden_df%>%
  mutate(DorR=factor(DorR,levels=c("D","R")))%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=DorR,y=SNV_burden,col=DorR))+
  geom_violin(fill="white")+
  geom_jitter(alpha=0.25,size=0.25)+
  facet_grid(~Pair_new)+
  scale_y_continuous(limits=c(0,2000))+
  scale_color_manual(values=DorR_cols)+
  theme_bw()+
  my_theme+
  theme(legend.position = "none",strip.text.x = element_text(size=6))+
  labs(x="Donor (D) or Recipient (R)",y="Somatic SNV burden")

ggsave(filename = paste0(plots_dir,"DvsR_mutburden_comparison_plot.pdf"),D_vs_R_mut_burden_comparison_plot,width = 4,height=2)

#Violin plot of adjusted SNV burdens (APOBEC/ platinum removed): donor vs recipient
mut_burden_df%>%
  mutate(DorR=factor(DorR,levels=c("D","R")))%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  ggplot(aes(x=DorR,y=SNV_burden_adj1,col=DorR))+
  geom_violin(fill="white")+
  geom_jitter(alpha=0.25)+
  facet_grid(~Pair_new)+
  scale_y_continuous(limits=c(0,2000))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="Donor (D) or Recipient (R)",y="Somatic SNV burden")

#Violin plot comparison of donors vs recipients
D_vs_R_mut_burden_comparison_adjusted_plot<-mut_burden_df%>%
  mutate(DorR=factor(DorR,levels=c("D","R")))%>%
  mutate(Pair_new=factor(new_pair_names[ID],levels=paste("Pair",1:nrow(Pair_metadata),sep = "_")))%>%
  filter(!Sample%in%expanded_clade_samples)%>%
  filter(!Sample %in%outlier_samples)%>%
  ggplot(aes(x=DorR,y=SNV_burden_adj1,fill=DorR,col=DorR))+
  geom_violin(fill="white")+
  geom_jitter(alpha=0.25,size=0.25)+
  facet_grid(~Pair_new)+
  scale_y_continuous(limits=c(0,2000))+
  scale_color_manual(values=DorR_cols)+
  theme_bw()+
  my_theme+
  theme(legend.position = "none",strip.text.x = element_text(size=6))+
  labs(x="Donor (D) or Recipient (R)",y="Somatic SNV burden (adjusted)")

ggsave(filename = paste0(plots_dir,"DvsR_mutburden_comparison_adjusted_plot.pdf"),D_vs_R_mut_burden_comparison_adjusted_plot,width = 4,height=2)

#Function to compare donor vs recipient mutation burdens of all the pairs
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

#Apply the function
pval_comparisons(mut_burden_df,test_field="SNV_burden_adj1",test_type = "two.sided")
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
  theme(axis.line.x = element_blank(),axis.ticks.x=element_blank())+
  labs(x="Pair",y="Increase in recipient mean\nmutation burden from HSCT\n(95% confidence interval)")

ggsave(filename = paste0(plots_dir,"DvsR_difference_CI_plot.pdf"),D_vs_R_difference_CI_plot,width =3.5,height=2)

#Compare signatures for Donor & Recipient branches
exp_no=25
tree<-all.trees[[paste0("Pair",exp_no)]]
exposures_df%>%
  filter(Pair==paste0("Pair",exp_no))%>%
  mutate(branch_type=sapply(node,function(node) {
    tips<-getTips(tree,node)
    tip_cats<-ifelse(as.numeric(substr(tips,3,7))%%2==0,"D","R")
    if(all(tip_cats=="D")) {return("D")} else if(all(tip_cats=="R")) {return("R")} else {return("M")}}))%>%
  mutate(branch_type=factor(branch_type,levels = c("M","D","R")))%>%
  gather(key="Mutation_signature",value="Proportion",-node,-Pair,-branch_type)%>%
  ggplot(aes(x=branch_type,y=Proportion,col=Mutation_signature))+
  #geom_point(alpha=0.5)+
  geom_boxplot()+
  facet_grid(~Mutation_signature)+
  theme_classic()


#----------------------------------
# ANALYSIS OF COPY NUMBER ALTERRATIONS (CNAs)
#----------------------------------

#Plot of copy number changes by individual
CNA_autosomal_plot<-Y_loss_df%>%
  dplyr::rename("Sample"=id)%>%
  mutate(Copy_number_change=ifelse(loss_of_Y=="YES","ChrY_loss",NA))%>%
  left_join(sample_metadata%>%select(Sample,DorR,ID))%>%
  dplyr::select(Sample,Copy_number_change,DorR,"Pair"=ID)%>%
  bind_rows(CN_change_df)%>%
  filter(!is.na(Copy_number_change))%>%
  left_join(Pair_metadata,by="Pair")%>%
  filter(!Copy_number_change%in%c("ChrY_loss" ,"ChrX_del"))%>% #Can remove the LOY to get better resolution on the other changes
  ggplot(aes(x=DorR,y=1,fill=stringr::str_wrap(factor(Copy_number_change),width=22)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.15)+
  facet_grid(cols=vars(factor(Pair_new,levels=paste0("Pair_",1:10))),drop=F)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  my_theme+
  theme(legend.key.size = unit(3,"mm"))+
  labs(x="Donor or Recipient",y="Number of HSPCs",fill="Copy number\nabnormality")
ggsave(filename = paste0(plots_dir,"CNA_autosomal_plot.pdf"),CNA_autosomal_plot,device = "pdf",width=7,height = 2)

#Plot of copy number changes by individual
CNA_sex_plot<-Y_loss_df%>%
  dplyr::rename("Sample"=id)%>%
  mutate(Copy_number_change=ifelse(loss_of_Y=="YES","ChrY_loss",NA))%>%
  left_join(sample_metadata%>%select(Sample,DorR,ID))%>%
  dplyr::select(Sample,Copy_number_change,DorR,"Pair"=ID)%>%
  bind_rows(CN_change_df)%>%
  filter(!is.na(Copy_number_change))%>%
  left_join(Pair_metadata,by="Pair")%>%
  filter(Copy_number_change%in%c("ChrY_loss","ChrX_del"))%>% #Can remove the LOY to get better resolution on the other changes
  ggplot(aes(x=DorR,y=1,fill=stringr::str_wrap(factor(Copy_number_change),width=22)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.1)+
  facet_grid(cols=vars(factor(Pair_new,levels=paste0("Pair_",1:10))),drop=F)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  my_theme+
  labs(x="Donor or Recipient",y="Number of HSPCs",fill="Copy number\nabnormality")
ggsave(filename = paste0(plots_dir,"CNA_sex_plot.pdf"),CNA_sex_plot,device = "pdf",width=7,height = 2)

#Proportion of colonies with an autosomal CNA in donors vs recipients (excluding LOY)
autosomal_CNA_prop_plot<-right_join(CN_change_df%>%filter(Copy_number_change!="ChrX_del"),sample_metadata%>%mutate(Pair=paste0("Pair",Pair)))%>%
  filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  replace_na(replace=list(Copy_number_change="Normal"))%>%
  mutate(CNA=ifelse(Copy_number_change=="Normal","Normal","Abnormal"))%>%
  group_by(CNA,DorR)%>%
  summarise(n=n())%>%
  pivot_wider(id_cols=c("CNA","DorR"),names_from="CNA",values_from="n")%>%
  mutate(abnormal_prop=Abnormal/(Abnormal+Normal),
         CI_low=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[1])}),
         CI_high=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[2])}))%>%
  ggplot(aes(x=DorR,col=DorR,y=100*abnormal_prop,ymin=100*CI_low,ymax=100*CI_high))+
  geom_point()+
  geom_errorbar(width=0.2)+
  scale_color_manual(values=DorR_cols)+
  theme_bw()+
  my_theme+
  scale_y_continuous(limits=c(0,2))+
  theme(legend.position = "none")+
  labs(x="Donor or Recipient",y="Proportion of HSPCs\nwith CNA (%)")
ggsave(filename = paste0(plots_dir,"autosomal_CNA_prop_plot.pdf"),autosomal_CNA_prop_plot,device = "pdf",width=2,height = 2.5)

#----------------------------------
# ANALYSIS OF STRUCTURAL VARIANTS (SVs)
#----------------------------------

SV_df<-read_csv(paste0(root_dir,"/data/SVs_combined.csv"))%>%dplyr::rename("Sample"=sample)
SV_df$length<-sapply(1:nrow(SV_df),function(i) ifelse(SV_df$`SV type`[i]%in%c("DEL","INV"),SV_df$pos2_start[i]-SV_df$pos1_start[i],NA))

#Number of samples with SVs
SV_df%>%dplyr::filter(Sample%in%final_sample_list)%>%group_by(Sample)%>%summarise(n=n())%>%nrow()

SV_summary_df<-SV_df%>%
  dplyr::filter(Sample%in%final_sample_list)%>%
  group_by(Sample)%>%
  summarise(n=n(),DorR=DorR[1],types=paste(unique(`SV type`),collapse=","))

SV_summary_df%>%
  group_by(types)%>%
  summarise(n=n())

SV_summary_df%>%
  group_by(DorR)%>%
  summarise(n=n())
prop.test(c(7,16),n=c(1262,1562))

SV_plot<-SV_summary_df%>%
  left_join(Pair_metadata,by="Pair")%>%
  ggplot(aes(x=DorR,y=1,fill=stringr::str_wrap(factor(types),width=22)))+
  geom_bar(stat="identity",position="stack",col="black",size=0.15)+
  facet_grid(cols=vars(factor(Pair_new,levels=paste0("Pair_",1:10))),drop=F)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  my_theme+
  labs(x="Donor or Recipient",y="Samples",fill="Structural\nvariant type")
ggsave(filename = paste0(plots_dir,"SV_plot.pdf"),SV_plot,device = "pdf",width=5.5,height = 2)

#Proportion of colonies with a SV in donors vs recipients
SV_prop_plot<-right_join(SV_summary_df%>%dplyr::select(-DorR,"SV"=n),sample_metadata%>%mutate(Pair=paste0("Pair",Pair)))%>%
  filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  replace_na(replace=list(SV="Normal"))%>%
  mutate(SV=ifelse(SV=="Normal","Normal","Abnormal"))%>%
  group_by(SV,DorR)%>%
  summarise(n=n())%>%
  pivot_wider(id_cols=c("SV","DorR"),names_from="SV",values_from="n")%>%
  mutate(abnormal_prop=Abnormal/(Abnormal+Normal),
         CI_low=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[1])}),
         CI_high=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[2])}))%>%
  ggplot(aes(x=DorR,col=DorR,y=100*abnormal_prop,ymin=100*CI_low,ymax=100*CI_high))+
  geom_point()+
  geom_errorbar(width=0.2)+
  scale_color_manual(values=DorR_cols)+
  theme_bw()+
  my_theme+
  scale_y_continuous(limits=c(0,2))+
  theme(legend.position = "none")+
  labs(x="Donor or Recipient",y="Proportion of HSPCs\nwith SV (%)")
ggsave(filename = paste0(plots_dir,"SV_prop_plot.pdf"),SV_prop_plot,device = "pdf",width=2,height = 2.5)

#----------------------------------
# CALCULATE COMBINED RISK OF SV OR CNA
#----------------------------------

SV_or_CNA_samples<-CN_change_df%>%filter(Copy_number_change!="ChrX_del")%>%
  dplyr::bind_rows(SV_summary_df)%>%pull(Sample)

SV_or_CNA_prop_plot<-sample_metadata%>%
  filter(Sample%in%unlist(lapply(all.trees.cc.nodups,function(tree) tree$tip.label)))%>%
  mutate(SV_or_CNA=ifelse(Sample%in%SV_or_CNA_samples,"Abnormal","Normal"))%>%
  group_by(SV_or_CNA,DorR)%>%
  summarise(n=n())%>%
  pivot_wider(id_cols=c("SV_or_CNA","DorR"),names_from="SV_or_CNA",values_from="n")%>%
  mutate(abnormal_prop=Abnormal/(Abnormal+Normal),
         CI_low=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[1])}),
         CI_high=mapply(x=Abnormal,n=(Abnormal+Normal),function(x,n){return(binom.test(x=x,n=n)$conf.int[2])}))%>%
  ggplot(aes(x=DorR,col=DorR,y=100*abnormal_prop,ymin=100*CI_low,ymax=100*CI_high))+
  geom_point()+
  geom_errorbar(width=0.2)+
  scale_color_manual(values=DorR_cols)+
  theme_bw()+
  my_theme+
  scale_y_continuous(limits=c(0,3))+
  theme(legend.position = "none")+
  labs(x="Donor or Recipient",y="Proportion of HSPCs with\nSV or autosomal CNA (%)")
ggsave(filename = paste0(plots_dir,"SV_or_CNA_prop_plot.pdf"),SV_or_CNA_prop_plot,device = "pdf",width=1.5,height = 2)

prop.test(c(9,26),c(1252,1536))
