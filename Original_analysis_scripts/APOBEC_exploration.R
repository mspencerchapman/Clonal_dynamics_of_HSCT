#ADD IN INFO ABOUT WHETHER
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


#Test using a linear mixed effects model - define APOBEC+ as having more than 15 mutations assigned to APOBEC signature
lme1=lmerTest::lmer(N4_abs>15~DorR+(1|ID),data=dat_full)
summary(lme1)
confint(lme1)

lme2.1=lmerTest::lmer(N4_abs>15~clonal_expansion+(1|ID)+(1|DorR),data=dat_full)
summary(lme2.1)
confint(lme2.1)

#Test using a linear mixed effects model
lme2.2=lmerTest::lmer(N4_abs>15~driver+(1|ID)+(1|DorR),data=dat_full)
summary(lme2.2)
confint(lme2.2)

####NOW LOOK TO LINK WITH CONDITIONING TYPE/ PATIENT TYPE/ AGE
lme3.1=lmerTest::lmer(N4_abs>15~conditioning+(1|ID)+(1|DorR),data=dat_full)
summary(lme3.1)
confint(lme3.1)

lme3.2=lmerTest::lmer(N4_abs>15~stem_cell_source+(1|ID)+(1|DorR),data=dat_full)
summary(lme3.2)
confint(lme3.2)

lme3.3=lmerTest::lmer(N4_abs>15~TBI+(1|ID)+(1|DorR),data=dat_full)
summary(lme3.3)
confint(lme3.3)

lme3.4=lmerTest::lmer(N4_abs>15~donor_sex+(1|ID)+(1|DorR),data=dat_full)
summary(lme3.4)
confint(lme3.4)

all_lmes=list(lme1,lme2.1,lme2.2,lme3.1,lme3.2,lme3.3,lme3.4)
lme_test=c("Donor vs Recipient",
           "Clonal expansion vs Singleton",
           "Driver vs No driver",
           "MAC vs RIC conditioning",
           "BM vs PBSC stem cell source",
           "TBI vs No TBI",
           "Male vs female donor")

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

CI_plot<-summary%>%
  ggplot(aes(xmin=`2.5 %`,xmax=`97.5 %`,y=test))+
  geom_errorbar(width=0.1)+
  geom_vline(xintercept=0,linetype=2)+
  theme_classic()+
  my_theme

pvalue_plot<-summary2%>%
  ggplot(aes(y=`Pr(>|t|)`,x=factor(test,levels = rev(lme_test))))+
  geom_point(size=0.75)+
  scale_y_log10()+
  geom_hline(yintercept=0.05,linetype=2)+
  theme_classic()+
  my_theme+
  labs(x="Comparison",y="P value")+
  coord_flip()

ggsave(filename = paste0(plots_dir,"APOBEC_associations_pvalue_plot.pdf"),pvalue_plot,width=3,height=1.5)


##Can do individual level comparisons using Chi-squared testing
#
temp2=dat_full%>%
  filter(DorR=="R")%>%
  group_by(Pair_new)%>%
  summarise(n_clonal_expansions=sum(clonal_expansion=="YES"),
            n_singleton=sum(clonal_expansion=="NO"),
            n_clonal_with_APOBEC=sum(clonal_expansion=="YES" & N4_abs>15),
            n_singleton_with_APOBEC=sum(clonal_expansion=="NO" & N4_abs>15))


for(i in 1:nrow(temp2)) {
  prop.res=prop.test(x=as.numeric(temp2[i,4:5]),n=as.numeric(temp2[i,2:3]))
  print(prop.res)
}


temp4=temp3%>%
  filter(DorR=="R")%>%
  group_by(Pair_new)%>%
  summarise(n_driver=sum(driver=="YES"),
            n_no_driver=sum(driver=="NO"),
            n_driver_with_APOBEC=sum(driver=="YES" & N4_abs>15),
            n_no_driver_with_APOBEC=sum(driver=="NO" & N4_abs>15))

for(i in 1:nrow(temp4)) {
  prop.res=prop.test(x=as.numeric(temp4[i,4:5]),n=as.numeric(temp4[i,2:3]))
  print(prop.res)
}

#Now do aggregating over all individuals.
#This is heavily skewed by Pair 10 who has the most drivers and the most APOBEC
#Therefore more appropriate to do as a mixed effects model where ID is a random effect

temp5=temp3%>%
  filter(DorR=="R")%>%
  summarise(n_driver=sum(driver=="YES"),
            n_no_driver=sum(driver=="NO"),
            n_driver_with_APOBEC=sum(driver=="YES" & N4_abs>15),
            n_no_driver_with_APOBEC=sum(driver=="NO" & N4_abs>15))

prop.res=prop.test(x=as.numeric(temp5[1,3:4]),n=as.numeric(temp5[1,1:2]))
print(prop.res)


### LOOK AT TIMING OF BRANCHES WITH APOBEC
#Update the exposure df with absolute mutation numbers of APOBEC &
#the minimum & maximum molecular time of the branch
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

Pair_metadata$mol_time_at_HCT=sapply(Pair_metadata$Pair,function(pair) {
  get_molecular_time_of_HCT(mol_time_at_sampling = exposures_df_abs%>%filter(Pair==pair)%>%pull(Max_height)%>%max(),
                            age_of_donor_at_sampling=Pair_metadata%>%filter(Pair==pair)%>%pull(Age),
                            age_of_donor_at_HCT=Pair_metadata%>%filter(Pair==pair)%>%pull(Age_at_transplant),
                            )
})

#Plot the possible timing of APOBEC exposure relative to HCT
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

ggsave(filename = paste0(plots_dir,"APOBEC_timing_plot.pdf"),APOBEC_timing_plot,width=6,height=3)

