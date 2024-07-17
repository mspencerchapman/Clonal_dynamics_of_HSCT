#========================================#
# Load packages (and install if they are not installed yet) ####
#========================================#
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","ape","dichromat","seqinr","stringr","readr","phytools","devtools","lmerTest","pheatmap","tidyr")
bioconductor_packages=c("clusterProfiler","MutationalPatterns","GenomicFeatures","BSgenome.Hsapiens.UCSC.hg19")

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

if(!require("nrmisc", character.only=T,quietly = T, warn.conflicts = F)){
  devtools::install_github("nicolaroberts/nrmisc", build_vignettes = F)
  library("nrmisc",character.only=T,quietly = T, warn.conflicts = F)
}

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
trees_list<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/tree_lists.Rds"))
details_lists<-readRDS(paste0(root_dir,"/data/tree_and_mutation_files/details_lists.Rds"))

#Extract objects from these lists in a 'for' loop
for(x in names(trees_list)) {assign(x,trees_list[[x]])}
for(x in names(details_lists)) {assign(x,details_lists[[x]])}

#========================================#
# Define custom functions for this script ####
#========================================#

#Now define the function for getting the properties for all these sites
# This is adapted from Nicola Roberts code
#If no 'control positions' are used, then the whole callable genome is the default control
get_hg19_properties <- function(gpos, location='farm',control_gpos=NULL,root='/lustre/scratch126/casm/team154pc/ms56/reference_files/nr3/'){
  require(GenomicFeatures)
  require(stringr)
  require("BSgenome.Hsapiens.UCSC.hg19")
  
  # all gpos elements must have width 1.
  if(any(width(gpos)!=1)){stop('All gpos elements must have width 1.')}
  # gpos must have mcol called "hist"
  if(!"hist" %in% colnames(mcols(gpos))) {stop("gpos must have column 'hist'")}
  
  gpos <- sort(gpos, ignore.strand=TRUE)
  
  # genome property files - list by property type
  in_dir <- file.path(root, 'results/GRanges/')
  f_in <- list.files(in_dir)
  f_in <- split(f_in, sapply(strsplit(f_in, "_"), `[`, 1))
  
  # don't include plain g_quadruplex_GR.RData, just use g4 subset with loops <=4
  f_in <- f_in[-which(names(f_in)=='g')]
  
  # load callable genome for ECDF calc
  load(file.path(root, 'results/callable_genome.RData'))
  gpos<-subsetByOverlaps(gpos,callable_genome)
  
  for (i in 1:length(f_in)){
    
    if (length(f_in[[i]])==1) {
      load(file.path(in_dir, f_in[[i]]))
      pname <- gsub(".RData", "", f_in[[i]])
      cat(pname,sep="\n")
      prop <- get(pname)
      rm(list=pname)
      
      suppressWarnings(seqinfo(prop) <- seqinfo(gpos))
      prop <- trim(prop)
      
      ovl <- findOverlaps(gpos, prop)
      
      # only use one metric
      if (ncol(mcols(prop))==1) j <- 1 else {
        if (pname %in% c("LAD_GR", "gene_GR")) {
          j <- which(names(mcols(prop))=="dens_1e6")
        } else if (pname %in% c("cruciform_inverted_rep_GR", "short_tandem_rep_GR")){
          j <- which(names(mcols(prop))=="dens_3e3")
        } else  {
          j <- which(names(mcols(prop))=="dist_log10")
        }
      }
      
      cname <- gsub(" ", "_", paste(gsub("_GR", "", pname), names(mcols(prop))[j]))
      mcols(gpos)[,cname] <- mcols(prop)[subjectHits(ovl),j]
      
      # get quantile value - callable overlap (force 1kb tiles)
      kbs <- unlist(tile(prop, width=1e3))
      ovl <- findOverlaps(kbs, prop)
      kbs$value <- mcols(prop)[subjectHits(ovl),j]
      
      if(is.null(control_gpos)) {
        prop <- subsetByOverlaps(kbs, callable_genome)
        p_ecdf <- ecdf(jitter(prop$value))
        mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))
        rm(prop, kbs, j, p_ecdf, ovl, cname, pname)
      } else {
        #New code to use the full mutation set as 'control' to calculate background quantiles
        overlaps <- findOverlaps(kbs, control_gpos)
        temp<-control_gpos
        temp$value[overlaps@to]<-kbs$value[overlaps@from]
        p_ecdf <- ecdf(jitter(temp$value))
        mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))
        rm(kbs, j, p_ecdf, ovl, overlaps, cname, pname,temp)
      }
      
    } else {
      for(k in seq_along(f_in[[i]])){
        fname <- f_in[[i]][k]
        pname <- gsub(".RData", "", fname)
        tname <- strsplit(pname, "_")[[1]][2]
        
        if (!tname %in% gpos$hist) next
        
        cat(pname,sep="\n")
        load(file.path(in_dir, fname))
        prop <- get(pname)
        suppressWarnings(seqinfo(prop) <- seqinfo(gpos))
        prop <- trim(prop)
        rm(list=pname)
        
        ovl <- findOverlaps(gpos[gpos$hist==tname], prop)
        cname <- names(f_in)[i]
        
        mcols(gpos)[gpos$hist==tname,cname] <- prop$value[subjectHits(ovl)]
        
        # get quantile value - callable overlap (force 1kb tiles)
        kbs <- unlist(tile(prop, width=1e3))
        ovl <- findOverlaps(kbs, prop)
        kbs$value <- mcols(prop)[subjectHits(ovl),1] #edited this to 1 as there is no 'j' in this section
        
        if(is.null(control_gpos)) {
          prop <- subsetByOverlaps(kbs, callable_genome)
          p_ecdf <- ecdf(jitter(prop$value))
          mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))
          rm(prop, kbs, p_ecdf, ovl, cname, pname)
        } else {
          #New code to use the full mutation set as 'control' to calculate background quantiles
          overlaps <- findOverlaps(kbs, control_gpos)
          temp<-control_gpos
          temp$value[overlaps@to]<-kbs$value[overlaps@from]
          p_ecdf <- ecdf(jitter(temp$value))
          mcols(gpos)[,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[,cname]))
          rm(kbs, p_ecdf, ovl, overlaps, cname, pname,temp)
        }
      }
    }
  }
  return(gpos)
}


plot_quantile_metrics=function(df,
                               metrics_to_plot="signif",
                               signif_levels=c(1e-10,1e-5,1e-2)) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  qs <- df[,union(1:6, grep('^q_', colnames(df)))]
  
  ##Calculate significance of each metric
  bin_data=function(num_vec,bins) {
    x=vector(mode="numeric",length=(length(bins)-1))
    for(i in 1:(length(bins)-1)) {
      x[i]<-sum(num_vec>bins[i] & num_vec<=bins[i+1])
    }
    return(x)
  }
  all_metrics=grep("^q_",colnames(qs),value = T)
  metric_res=data.frame(metric=all_metrics)
  metric_res$p.value=sapply(all_metrics,function(metric) {
    #cat(metric,sep="\n")
    nbin=10; bins=seq(0,1,1/nbin)
    dat_binned=bin_data(num_vec = qs[[metric]][!is.na(qs[[metric]])],bins=bins)
    chisq_res=chisq.test(x = dat_binned)
    return(chisq_res$p.value)
  })
  metric_res$q.value=p.adjust(p=metric_res$p.value,method="BH")
  metric_res$signif=ifelse(metric_res$q.value<signif_levels[1],"***",ifelse(metric_res$q.value<signif_levels[2],"**",ifelse(metric_res$q.value<signif_levels[3],"*","n.s.")))
  
  if(metrics_to_plot[1]=="signif") {
    significant_metrics=metric_res%>%filter(signif!="n.s.")%>%pull(metric)
  } else if(metrics_to_plot[1]=="all"){
    significant_metrics<-all_metrics
  } else {
    significant_metrics<-metrics_to_plot
  }
  p1<-qs[,-c(1:6)]%>%
    gather(key="Metric",value="Quantile")%>%
    filter(Metric%in%significant_metrics)%>%
    mutate(Metric=gsub("^q_","",Metric))%>%
    ggplot(aes(x=Quantile))+
    geom_density(fill="lightblue")+
    facet_wrap(~Metric,ncol=5)+
    theme_classic()
  
  return(p1)
}

combined_plot=function(df_list,metrics_to_plot=NULL) {
  library(dplyr)
  library(ggplot2)
  
  datasets=names(df_list)
  
  if(is.null(metrics_to_plot)) {
    all_quantile_metrics<-Map(df=df_list,dataset=datasets,function(df,dataset) {
      df<-df[,grep('^q_', colnames(df))]%>%
        gather(key="Metric",value="Quantile")%>%
        mutate(Metric=gsub("^q_","",Metric))%>%
        mutate(dataset=dataset,.before=1)
      return(df)
    })%>%dplyr::bind_rows()
  } else {
    all_quantile_metrics<-Map(df=df_list,dataset=datasets,function(df,dataset) {
      df<-df[,grep('^q_', colnames(df))]%>%
        gather(key="Metric",value="Quantile")%>%
        filter(Metric%in%metrics_to_plot)%>%
        mutate(Metric=gsub("^q_","",Metric))%>%
        mutate(dataset=dataset,.before=1)
      return(df)
    })%>%dplyr::bind_rows()
  }
  
  plot<-ggplot(all_quantile_metrics,aes(x=Quantile))+
    geom_density(fill="lightblue")+
    facet_grid(Metric~dataset)+
    scale_x_continuous(breaks=c(0,1))+
    theme_classic()+
    theme(strip.text.y = element_text(angle=0),
          axis.text.x=element_blank())
  
  plot
}

#========================================#
# IMPORT MUTATION LISTS AS GRANGES ####
#========================================#

##1. Create the GRanges object of all HCT mutations to use as a 'control'
all_HCT_mutations_GR<-GenomicRanges::makeGRangesFromDataFrame(df=dplyr::bind_rows(all.muts.nodups)%>%mutate(Chrom=paste0("chr",Chrom)),
                                        seqnames.field="Chrom",
                                        start.field="Pos",
                                        end.field="Pos")
all_HCT_mutations_GR$hist<-"Myeloid"

#2. Import all the APOBEC vcfs as GRangesList
vcf_dir=paste0(root_dir,"/data/APOBEC_VCFs/")
vcf_files=list.files(vcf_dir,pattern=".vcf$",full.names = T)
sample_names=lapply(stringr::str_split(list.files(vcf_dir,pattern=".vcf$"),pattern = "_"),`[`,2:3)
sample_names=sapply(sample_names,function(x) {
  paste0(x[1],"_node",gsub(".vcf","", x[2]))
  })
APOBEC_GR<-read_vcfs_as_granges(vcf_files = vcf_files,
                     sample_names = sample_names,
                     genome = "BSgenome.Hsapiens.UCSC.hg19",
                     change_seqnames = T,
                     type="snv",
                     predefined_dbs_mbs =T)

#Combine APOBEC mutations into a single range
APOBEC_GR_all=unlist(as(APOBEC_GR, "GRangesList"))
APOBEC_GR_all$hist<-"Myeloid"

#========================================#
# RUN THE FUNCTION ####
#========================================#

#Run function for APOBEC mutations, using all HCT control----
APOBEC_with_properties=get_hg19_properties(gpos=APOBEC_GR_all,control_gpos = all_HCT_mutations_GR,root=paste0(root_dir,"/data/genomic_loci_reference_files"))
APOBEC_with_properties_df=as.data.frame(APOBEC_with_properties)

#Run function for APOBEC mutations, using no control----
APOBEC_no_control_with_properties=get_hg19_properties(gpos=APOBEC_GR_all,root=paste0(root_dir,"/data/genomic_loci_reference_files"))
APOBEC_no_control_with_properties_df=as.data.frame(APOBEC_no_control_with_properties)

#Run function over all HCT mutations, using no control----
all_HCT_mutations_with_properties=get_hg19_properties(gpos=all_HCT_mutations_GR,root=paste0(root_dir,"/data/genomic_loci_reference_files"))
all_HCT_mutations_with_properties_df=as.data.frame(all_HCT_mutations_with_properties)

#========================================#
# VISUALIZE THE RESULTS ####
#========================================#

APOBEC_with_ctrl_plot<-plot_quantile_metrics(df=APOBEC_with_properties_df,
                                             metrics_to_plot = significant_metrics)
APOBEC_no_ctrl_plot<-plot_quantile_metrics(df=APOBEC_no_control_with_properties_df,
                                           metrics_to_plot = significant_metrics)
all_HCT_mutations_plot<-plot_quantile_metrics(df=all_HCT_mutations_with_properties_df,
                                              metrics_to_plot = significant_metrics)

# This is a list of the metrics that come up as significant
significant_metrics=c("q_gc_content_value",
                      "q_rep_timing_value_value",
                      "q_cruciform_inverted_rep_dens_3e3",
                      "q_gene_dens_1e6",
                      "q_telomere_dist_log10",
                      "q_seq_complexity_value",
                      "q_ALU_rep_dist_log10",
                      "q_recomb_rate_nearest_value")

## Generate Extended Data Fig.5c ----  
##Plot all together
all_dat=list("APOBEC_only\n(with_background_model)"=APOBEC_with_properties_df,
     "APOBEC_only\n(no_background_model)"=APOBEC_no_control_with_properties_df,
     "all_HCT_mutations"=all_HCT_mutations_with_properties_df)

p_comb<-combined_plot(all_dat,metrics_to_plot = significant_metrics)
ggsave(paste0(plots_dir,"ExtDatFig5c"),plot = p_comb,width=7,height=7)
