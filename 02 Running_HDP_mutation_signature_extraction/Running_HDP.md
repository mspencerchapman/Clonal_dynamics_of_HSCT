## Running HDP

HDP standards for 'Hierarchical Dirichlet Process' and is a method that can be applied to deconvolute multiple mutational profiles into component 'mutational signatures'. It is an alternative to non-negative matrix factorization (NMF) which is also commonly applied for this purpose. HDP has been implemented for mutational signature extraction by Nicola Roberts ( https://github.com/nicolaroberts/hdp ). However, there are some compatibility issues with linux Ubuntu v22, and therefore it has recently been forked and updated by Nick Williams ( https://github.com/NickWilliamsSanger/hdp ).

HDP has options to run multiple single MCMC chains in parallel, and to then combine the results.

Running HDP has 4 major steps: \
1. Extracting the 96-mutation profiles from each sample in the dataset. In our case, a 'sample' is a branch from a phylogenetic tree. Only branches with ≥50 mutations are included, as below that the profiles are too imprecise for inclusion. \
2. Running multiple HDP single chains on the dataset. \
3. Combine these results to create the 'HDP_multichain' \
4. Extracting the final signature (component) and exposure results from this.

### 1. Extracting the 96-mutation profiles from each sample

This is done using the script '01 HDP_input.R'
You will need to update the file paths to reflect your local genome file (this should be GRCh37, ensemble version)

### 2. Running multiple HDP single chains on the dataset

This is the code we use for submitting the hdp_single_chain.R script in parallel
You need to have run the 'HDP_input.R' script prior to this so that the trinucleotide contexts and key_table.txt files are in the directory
This is on an LSF farm; may need to alter depnding on your local compute farm setup
This typically takes ~1-2 hours per chain

```bash

ROOT_DIR="/PATH/TO/REPO/"#Set this to the location of the cloned github directory
HDP_SINGLE_CHAIN_SCRIPT="${ROOT_DIR}02 Running_HDP_mutation_signature_extraction/hdp_single_chain.R"
MEM=10000 #May need more memory depending on the number of mutations

#We typically run 20 chains
for n in $(seq 1 20); do
bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R"select[mem>${MEM}] rusage[mem=${MEM}]" -M${MEM} -J HDP_${n} Rscript $HDP_SINGLE_CHAIN_SCRIPT $n
done
```

### 3. Combining single chains into a single 'multi chain'

Run simple script in the HDP directory to combine these and create the 'HDP_multi_chain.Rds' output

```bash
ROOT_DIR="/PATH/TO/REPO/"#Set this to the location of the cloned github directory
cd ${ROOT_DIR}data/HDP

HDP_COMBINE_SCRIPT="${ROOT_DIR}02 Running_HDP_mutation_signature_extraction/hdp_combine_results.R"
Rscript $HDP_COMBINE_SCRIPT

```

### 4. Extracting the exposures and signature information

This is done using the 'generate_exposures_df' function (defined in the HSCT_functions.R script) using the paths to the HDP_multi_chain.R, trinucleotide_context.txt and key_table.txt objects outputted from the previous steps as arguments. This is in fact all done within the 'Compile_data_objects.R' and other scripts, but example code is shown below.

```r
root_dir="path/to/repo"
HDP_folder=paste0(root_dir,"/data/HDP")

generate_exposures_df=function(HDP_multi_chain_RDS_path,trinuc_mut_mat_path,key_table_path){
  require(hdp)
  require(dplyr)
  require(tibble)
  require(tidyr)
  
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


exposures_df=generate_exposures_df(HDP_multi_chain_RDS_path=paste0(HDP_folder,"/HDP_multi_chain.Rdata"),
                                   trinuc_mut_mat_path=paste0(HDP_folder,"/trinuc_mut_mat.txt"),
                                   key_table_path = paste0(HDP_folder,"/key_table.txt"))%>%dplyr::rename("Pair"=exp_ID)

#This function can generate the signature profiles
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

                                   
```


