##Using Nicola Roberts' code on assessing genomic positions

#Get some data for testing -----
#
# load('/Volumes/nfs_home/data/PanCan_callable_genome_lustre/results/callable_genome.RData')
# library(nrmisc)
# gpos <- sample_GRanges(callable_genome, 100)
# remove(callable_genome)
# gpos$hist <- sample(c("Breast", "Lung", "Thyroid"), 100, replace=TRUE)

library(MutationalPatterns)
library(GenomicFeatures)
library(stringr)
library(dplyr)
library(tidyr)
library("BSgenome.Hsapiens.UCSC.hg19")
library(nrmisc)

##Create the GRanges object of all HCT mutations to use as a 'control'
all_HCT_mutations_path<-"/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/HCT_sites_per_branch.txt"
all_HCT_mutations<-read.delim(all_HCT_mutations_path,stringsAsFactors = F)%>%
  filter(Mut_type=="SNV")%>%
  mutate(Chrom=paste0("chr",Chrom))
all_HCT_mutations_GR<-GenomicRanges::makeGRangesFromDataFrame(df=all_HCT_mutations,
                                        seqnames.field="Chrom",
                                        start.field="Pos",
                                        end.field="Pos")
all_HCT_mutations_GR$hist<-"Myeloid"

#Import all the APOBEC vcfs as GRangesList
vcf_dir="/lustre/scratch126/casm/team154pc/ms56/Zur_HSCT/APOBEC_VCFs/"
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

#Combine into single range
APOBEC_GR_all=unlist(as(APOBEC_GR, "GRangesList"))
APOBEC_GR_all$hist<-"Myeloid"

##Now do the same for the blood PVVs
PVVs<-read.delim("/lustre/scratch126/casm/team154pc/ms56/lesion_segregation/mutations_filtered.tsv",stringsAsFactors = F)
PVVs<-PVVs%>%filter(Type=="PVV" & Class=="PASS" & cat=="Adult_HSPC")%>%
  tidyr::separate(col="Chrom_pos",into = c("Chrom","Pos"),sep="-")%>%
  mutate(Chrom=paste0("chr",Chrom))
PVVs_GR<-GenomicRanges::makeGRangesFromDataFrame(df=PVVs,
                                                              seqnames.field="Chrom",
                                                              start.field="Pos",
                                                              end.field="Pos")
PVVs_GR$hist<-"Myeloid"

#Run function over the gpos object
PVVs_with_properties=get_hg19_properties(gpos=PVVs_GR,control_gpos = all_HCT_mutations_GR)
#PVVs_with_properties=get_hg19_properties(gpos=PVVs_GR)libr
PVVs_with_properties_df=as.data.frame(PVVs_with_properties)
PVVs_plot<-plot_quantile_metrics(df=PVVs_with_properties_df,
                      metrics_to_plot = "signif",
                      signif_levels = c(1e-3,1e-2,0.1))

#Now define the function for getting the properties for all these sites
#If no 'control positions' are used, then the whole callable genome is the default control
get_hg19_properties <- function(gpos, location='farm',control_gpos=NULL){
  require(GenomicFeatures)
  require(stringr)
  require("BSgenome.Hsapiens.UCSC.hg19")
  
  # all gpos elements must have width 1.
  if(any(width(gpos)!=1)){stop('All gpos elements must have width 1.')}
  # gpos must have mcol called "hist"
  if(!"hist" %in% colnames(mcols(gpos))) {stop("gpos must have column 'hist'")}
  
  gpos <- sort(gpos, ignore.strand=TRUE)
  
  # check location value, and set root accordingly
  location <- tolower(location)
  if (!location %in% c('farm', 'local')){stop("location must be one of 'farm' or 'local'")}
  if (location=='farm'){
    root <- '/lustre/scratch126/casm/team154pc/ms56/reference_files/nr3/'
  } else {
    root <- '/Volumes/nfs_home/data/'
  }
  
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
        # kbs <- unlist(tile(prop, width=1e3))
        # ovl <- findOverlaps(kbs, prop)
        # kbs$value <- prop$value[subjectHits(ovl)]
        # prop <- subsetByOverlaps(kbs, callable_genome)
        # 
        # p_ecdf <- ecdf(jitter(prop$value))
        # mcols(gpos)[gpos$hist==tname,paste0('q_', cname)] <- p_ecdf(jitter(mcols(gpos)[gpos$hist==tname,cname]))
        # 
        # 
        # rm(prop, kbs, p_ecdf, ovl, cname, pname, tname, fname)
      }
    }
    
  }
  
  
  # # nucleosome occupancy - use RAW data (at 1e5 random callable pos)
  # rpos <- nrmisc::sample_GRanges(callable_genome, 1e5)
  # 
  # mnase <- import.bw(file.path(root, "PanCan_genome_properties_lustre/raw/ENCFF000VNN.bigWig"),
  #                    which=c(rpos, granges(gpos)))
  # mnase <- keepSeqlevels(mnase, seqlevels(gpos))
  # seqinfo(mnase) <- seqinfo(gpos)
  # fo <- findOverlaps(gpos, mnase)
  # gpos$nucleosome_occupancy_raw <- as.numeric(NA)
  # gpos$nucleosome_occupancy_raw[queryHits(fo)] <- mnase$score[subjectHits(fo)]
  # 
  # # get quantile value - the rpos from callable
  # p_ecdf <- ecdf(jitter(subsetByOverlaps(mnase, rpos)$score))
  # gpos$q_nucleosome_occupancy_raw <- p_ecdf(jitter(gpos$nucleosome_occupancy_raw))
  # 
  
  return(gpos)
  
}

#Run function over the gpos object
APOBEC_with_properties=get_hg19_properties(gpos=APOBEC_GR_all,control_gpos = all_HCT_mutations_GR)
APOBEC_with_properties_df=as.data.frame(APOBEC_with_properties)

#Run function over the gpos object
APOBEC_no_control_with_properties=get_hg19_properties(gpos=APOBEC_GR_all)
APOBEC_no_control_with_properties_df=as.data.frame(APOBEC_no_control_with_properties)

#Run function over the gpos object
all_HCT_mutations_with_properties=get_hg19_properties(gpos=all_HCT_mutations_GR)
all_HCT_mutations_with_properties_df=as.data.frame(all_HCT_mutations_with_properties)

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




APOBEC_with_ctrl_plot<-plot_quantile_metrics(df=APOBEC_with_properties_df,
                                             metrics_to_plot = significant_metrics)
APOBEC_no_ctrl_plot<-plot_quantile_metrics(df=APOBEC_no_control_with_properties_df,
                                           metrics_to_plot = significant_metrics)
all_HCT_mutations_plot<-plot_quantile_metrics(df=all_HCT_mutations_with_properties_df,
                                              metrics_to_plot = significant_metrics)



significant_metrics=c("q_gc_content_value",
                      "q_rep_timing_value_value",
                      "q_cruciform_inverted_rep_dens_3e3",
                      "q_gene_dens_1e6",
                      "q_telomere_dist_log10",
                      "q_seq_complexity_value",
                      "q_ALU_rep_dist_log10",
                      "q_recomb_rate_nearest_value")
  

##Plot all together
all_dat=list("APOBEC_only\n(with_background_mod.)"=APOBEC_with_properties_df,
     "APOBEC_only\n(no_background_mod.)"=APOBEC_no_control_with_properties_df,
     "all_HCT_mutations"=all_HCT_mutations_with_properties_df)

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

p_comb<-combined_plot(all_dat,metrics_to_plot = significant_metrics)
ggsave("significant_metrics_APOBEC_muts.pdf",plot = p_comb,width=7,height=7)

p_comb_all_metrics<-combined_plot(all_dat)
ggsave("all_metrics_APOBEC_muts.pdf",plot = p_comb_all_metrics,width=5,height=7)
