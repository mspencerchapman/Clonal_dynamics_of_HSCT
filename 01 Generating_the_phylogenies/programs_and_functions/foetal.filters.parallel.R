#Function to subset the list of matrices
list_subset = function(list, select_vector) {
  for (i in 1:length(list)) {
    if(!is.null(dim(list[[i]]))) {
      list[[i]] = list[[i]][select_vector,]
    }
  }
  return(list)
}

#Function to create pval matrix based on the NV and NR matrix
pval_matrix = function(COMB_mats) {
  cat("Starting pval matrix generation\n")
  COMB_mats$NR[COMB_mats$NR == 0] <- 1 
  pval_mat <- matrix(0, nrow = nrow(COMB_mats$NV), ncol = ncol(COMB_mats$NV))
  if(COMB_mats$gender == "male") {
    for(i in 1:nrow(COMB_mats$NV)) {
    for (j in 1:ncol(COMB_mats$NR)) {
      if (!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {pval_mat[i,j] <- binom.test(COMB_mats$NV[i,j], COMB_mats$NR[i,j], p = 0.5, alternative = "less")$p.value}
      else {pval_mat[i,j] <- binom.test(COMB_mats$NV[i,j], COMB_mats$NR[i,j], p = 0.95, alternative = "less")$p.value}
    }
    if (i %% 1000 == 0) {print(i)}
    }
  } else if(COMB_mats$gender == "female") {
    for(i in 1:nrow(COMB_mats$NV)) {
      for (j in 1:ncol(COMB_mats$NR)) {
        pval_mat[i,j] <- binom.test(COMB_mats$NV[i,j], COMB_mats$NR[i,j], p = 0.5, alternative = "less")$p.value
      }
      if (i %% 1000 == 0) {print(i)}
    }
  }
  cat("Completed pval matrix generation\n")
  return(pval_mat)
}

#mat object needs the Chrom column with chromosome.  Returns the pval vector for each mutation.
germline.binomial.filter = function(COMB_mats){
  cat("Starting the germline binomial filter\n")
  XY_chromosomal = COMB_mats$mat$Chrom %in% c("X","Y")
  autosomal = !XY_chromosomal
  
  if(COMB_mats$gender=="female"){
    NV_vec = rowSums(COMB_mats$NV) #create vector of COMBINED variant reads across all samples
    NR_vec = rowSums(COMB_mats$NR) #create vector of COMBINED depth across all samples
    pval = rep(1,length(NV_vec))
    #For loop to calculate whether each one is likely to have come from a binomial distribution where true probability is 0.5 (i.e. would be expected for autosomal heterozygous germline variant)
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      } 
      if (n%%1000==0){
        print(n)
      }
    }
  }
  # IF MALE - need to separate off the XY chromosomes, as expected probability of germline variant is now close to 1.
  if(COMB_mats$gender=="male"){
    pval=rep(1,nrow(COMB_mats$NV))
    NV_vec = rowSums(COMB_mats$NV)[autosomal]
    NR_vec = rowSums(COMB_mats$NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    NV_vec = rowSums(COMB_mats$NV)[XY_chromosomal]
    NR_vec = rowSums(COMB_mats$NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  cat("Completed the germline binomial filter\n")
  return(pval)
}

#THE BETA-BINOMIAL FUNCTIONS
require(VGAM)
estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(COMB_mats){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  cat("Starting the beta-binomial filter\n")
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  rho_est = rep(NA,nrow(COMB_mats$NR))
  for (k in 1:nrow(COMB_mats$NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(COMB_mats$NV[k,]),
                                  NR_vec=as.numeric(COMB_mats$NR[k,]))
    if (k%%1000==0){
      print(k)
    }
  }
  cat("Completed the beta-binomial filter\n")
  return(rho_est)
}

#FILTER VARIANTS WITH LOW VAF AMONGST POSITIVE SAMPLES
#Filters mutations present in multiple samples that are (on aggregate) unlikely to have come from true somatic mutation distribution (i.e. 0.5 in auto, 1 in XY)
low_vaf_in_pos_samples_dp2 = function(COMB_mats, define_pos = 2) {
  cat("Starting the low_vaf_in_pos_samples_dp2 filter\n")
  if(all(c(nrow(COMB_mats$mat) == nrow(COMB_mats$NV),dim(COMB_mats$NV) == dim(COMB_mats$NR)))) {print("Input matrices are of correct dimensions")}
  pval=rep(0,nrow(COMB_mats$mat))
  if(COMB_mats$gender == "male") {
    for(k in 1:nrow(COMB_mats$mat)) {
      NV_vec <- COMB_mats$NV[k,]
      NR_vec <- COMB_mats$NR[k,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
        if (COMB_mats$mat$Chrom[k] %in% c("X","Y")) {
          pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.95, alt = "less")$p.value
        }
        else {
          pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        }
        if(k %% 1000 == 0) {print(k)}
      }
    }
  } else if(COMB_mats$gender == "female") {
    for(k in 1:nrow(COMB_mats$mat)) {
      NV_vec <- COMB_mats$NV[k,]
      NR_vec <- COMB_mats$NR[k,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
        pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        if(k %% 1000 == 0) {print(k)}
      }
    }
  }
  cat("Completed the low_vaf_in_pos_samples_dp2 filter\n")
  return(pval)
}


low_vaf_in_pos_samples_dp3 = function(COMB_mats, define_pos = 3) {
  cat("Starting the low_vaf_in_pos_samples_dp3 filter\n")
  if(all(c(nrow(COMB_mats$mat) == nrow(COMB_mats$NV),dim(COMB_mats$NV) == dim(COMB_mats$NR)))) {print("Input matrices are of correct dimensions")}
  pval=rep(0,nrow(COMB_mats$mat))
  if(COMB_mats$gender == "male") {
    for(k in 1:nrow(COMB_mats$mat)) {
      NV_vec <- COMB_mats$NV[k,]
      NR_vec <- COMB_mats$NR[k,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
        if (COMB_mats$mat$Chrom[k] %in% c("X","Y")) {
          pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.95, alt = "less")$p.value
        }
        else {
          pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        }
        if(k %% 1000 == 0) {print(k)}
      }
    }
  } else if(COMB_mats$gender == "female") {
    for(k in 1:nrow(COMB_mats$mat)) {
      NV_vec <- COMB_mats$NV[k,]
      NR_vec <- COMB_mats$NR[k,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
        pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        if(k %% 1000 == 0) {print(k)}
      }
    }
  }
  cat("Completed the low_vaf_in_pos_samples_dp3 filter\n")
  return(pval)
}

get_mean_depth = function(COMB_mats) {
  mean_depth = rowMeans(COMB_mats$NR)
  cat("Completed the mean depth filter\n")
  return(mean_depth)
}

get_max_depth_in_pos = function(COMB_mats) {
  cat("Starting the max_depth_in_pos filter\n")
  apply_max_depth_in_pos = function(i, mat_list) {
    if(!mat_list$mat$Chrom[i] %in% c("X","Y") | COMB_mats$gender == "female") {
      max_depth_in_pos_samples <- max(mat_list$NR[i,which(mat_list$NV[i,]>=min_variant_reads_auto)])
    } else {
      max_depth_in_pos_samples <- max(mat_list$NR[i,which(mat_list$NV[i,]>=min_variant_reads_xy)])
    }
    return(max_depth_in_pos_samples)
  }
  max_depth_in_pos_vec = sapply(1:nrow(COMB_mats$NV), apply_max_depth_in_pos, mat_list = COMB_mats)
  cat("Completed the max_depth_in_pos filter\n")
  return(max_depth_in_pos_vec)
}

get_max_pval_in_pos = function(COMB_mats) {
  apply_max_pval_in_pos = function(i, COMB_mats) {
    if(COMB_mats$gender == "male") {
      if(!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {
        max_pval_in_pos_samples <- max(COMB_mats$PVal[i,which(COMB_mats$NV[i,]>=min_variant_reads_auto)])
      } else {
        max_pval_in_pos_samples <- max(COMB_mats$PVal[i,which(COMB_mats$NV[i,]>=min_variant_reads_xy)])
      }
    } else if(COMB_mats$gender == "female") {
      max_pval_in_pos_samples <- max(COMB_mats$PVal[i,which(COMB_mats$NV[i,]>=min_variant_reads_auto)])
    }
    return(max_pval_in_pos_samples)  
  }
  max_pval_in_pos_vec = sapply(1:nrow(COMB_mats$NV), apply_max_pval_in_pos, COMB_mats = COMB_mats)
  return(max_pval_in_pos_vec)
}

get_max_vaf = function(COMB_mats) {
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  vaf_mat = COMB_mats$NV/COMB_mats$NR
  max_vaf = apply(vaf_mat,1,max)
  return(max_vaf)
}

#Removing columns with low coverage, or that are otherwise unwanted
remove_low_coverage_samples = function(COMB_mats,
                                       filter_params=NULL,
                                       min_sample_mean_cov,
                                       other_samples_to_remove = NULL,
                                       min_variant_reads_auto = 3,
                                       min_variant_reads_xy = 2) {
  mean_sample_cov = colMeans(COMB_mats$NR)
  low_cov_names = gsub(pattern = "_DEP", replacement = "", x = names(mean_sample_cov)[mean_sample_cov < min_sample_mean_cov])
  remove_cols <- which(gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats$NV)) %in% c(low_cov_names, other_samples_to_remove))
  if(length(remove_cols) > 0) {
    cat("Removing",length(remove_cols),"samples\n")
    COMB_mats$NV <- COMB_mats$NV[,-remove_cols]
    COMB_mats$NR <- COMB_mats$NR[,-remove_cols]
    COMB_mats$PVal <- COMB_mats$PVal[,-remove_cols]
    
    null_remove = rowSums(COMB_mats$NV >= min_variant_reads_auto|(COMB_mats$NV >= min_variant_reads_xy & COMB_mats$mat$Chrom %in% c("X","Y") & COMB_mats$gender == "male")) == 0
    cat(sum(null_remove),"mutations removed as no positives in any remaining samples.\n")
    COMB_mats = list_subset(COMB_mats, select_vector = !null_remove)
    if(!is.null(filter_params)) {
      filter_params = filter_params[!null_remove,]
      output = list(COMB_mats=COMB_mats,filter_params=filter_params)
    } else {
      output=COMB_mats
    }
  } else {
    cat("No samples removed\n")
    if(!is.null(filter_params)) {
      output = list(COMB_mats=COMB_mats,filter_params=filter_params)
    } else {
      output=COMB_mats
    }
  }
  return(output)
}

#Functions for filtering from the filter_params and COMB_mats object, setting the desired cut-offs
assess_mean_depth = function(i, COMB_mats, AUTO_low_depth_cutoff, AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff) {
  if(COMB_mats$gender == "male") {
    if(!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {
      ifelse((filter_params$mean_depth[i] >= AUTO_low_depth_cutoff & filter_params$mean_depth[i] <= AUTO_high_depth_cutoff), 1,0)
    } else {
      ifelse(filter_params$mean_depth[i] >= XY_low_depth_cutoff & filter_params$mean_depth[i] <= XY_high_depth_cutoff, 1,0)
    }
  } else if(COMB_mats$gender == "female") {
    ifelse((filter_params$mean_depth[i] >= AUTO_low_depth_cutoff & filter_params$mean_depth[i] <= AUTO_high_depth_cutoff), 1,0)
  }
}

assess_max_depth_in_pos = function(i, COMB_mats, min_depth_auto, min_depth_xy) {
  if(COMB_mats$gender == "male") {
    if(!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {
      ifelse(filter_params$max_depth_in_pos_samples[i] >= min_depth_auto, 1,0)
    } else {
      ifelse(filter_params$max_depth_in_pos_samples[i] >= min_depth_xy, 1,0)
    }
  } else if(COMB_mats$gender == "female") {
    ifelse(filter_params$max_depth_in_pos_samples[i] >= min_depth_auto, 1,0)
  }
}

assess_max_vaf = function(i, COMB_mats, min_vaf_auto, min_vaf_xy) {
  if(COMB_mats$gender == "male") {
    if(!COMB_mats$mat$Chrom[i] %in% c("X","Y")) {
      ifelse(filter_params$max_mut_vaf[i] >= min_vaf_auto, 1,0)
    } else {
      ifelse(filter_params$max_mut_vaf[i] >= min_vaf_xy, 1,0)
    }
  } else if(COMB_mats$gender == "female") {
    ifelse(filter_params$max_mut_vaf[i] >= min_vaf_auto, 1,0)
  }
}

get_filtered_mut_set = function(input_set_ID,
                                COMB_mats,
                                filter_params,
                                gender,
                                
                                ##PARAMETERS FOR FILTERING THE MUTATION SET. These arguments set the cut-offs for selecting the set of "true somatic mutations"
                                #NB. If parameter is set to NA, filter will not be applied.
                                retain_muts = NA, #vector of "mut_refs" (in format Chrom-Pos-Ref-Alt) to retain even if fail the filtering
                                exclude_muts = NA, #vector of "mut_refs" (in format Chrom-Pos-Ref-Alt) to exclude even if pass the filtering
                                germline_pval = -10,  
                                rho = 0.1, #Beta-binomial filter, rho value cut-off (i.e. rho must be > cut-off)
                                mean_depth=NA, #Numeric vector of length 4 in the order (1) AUTO min depth cutoff, (2) AUTO max depth cutoff, (3) XY min depth cutoff, (4) XY max depth cutoff
                                pval_dp2=NA, #Numeric vector of length 1. 
                                pval_dp3=0.01, #Numeric vector of length 1. 
                                min_depth = c(6,4), #Numeric vector of length 2 for (1) AUTO and (2) XY muts. Minimum depth that at least one positive sample must have.
                                min_pval_for_true_somatic = 0.1, #Numeric vector of length 1. At least one sample must have a p-value of >this cutoff for coming from a true somatic distribution.
                                min_vaf = c(0.2,0.8), #Numeric vector of length 2 for (1) AUTO and (2) XY muts. At least one sample must meet the minimum vaf for the mutation.
                                
                                #These arguments decide the thresholds for genotyping each sample. This may be less stringent.
                                min_variant_reads_SHARED = 2, #PARAMETERS FOR DECIDING GENOTYPE FOR EACH SAMPLE FOR EACH RETAINED MUTATION - should be equal to, or more relaxed than the above. At least one parameter must be used.
                                min_pval_for_true_somatic_SHARED = 0.05,
                                min_vaf_SHARED = c(0.2,0.7) #Numeric vector of length 2 for AUTO and XY muts
) {
  #List (1) the names of each of the initial filters that can be applied,
  #(2) the names of their input parameters,
  #(3) the variable name that must be present in the filter_params dataframe if filter is being used.
  #The elements of each must match
  filter_name=c("Germline filter","Beta-binomial filter","Mean depth filter","P-value within positive samples [defining positive as ≥2 variant reads] filter","P-value within positive samples [defining positive as ≥3 variant reads] filter","Minimum depth filter","Minimum p-value for true somatic mutation filter","Minimum VAF filter")
  set_parameter=list(germline_pval,rho,mean_depth,pval_dp2,pval_dp3,min_depth,min_pval_for_true_somatic,min_vaf)
  var_name = c("germline_pval","bb_rhoval","mean_depth","pval_within_pos_dp2","pval_within_pos_dp3","max_depth_in_pos_samples","max_pval_in_pos_sample","max_mut_vaf")
  
  #List the functions that operate on the input parameters to decide if the mutation is a pass (1) or fail (0) for each filter
  test_function=list(function(x) {ifelse(log10(filter_params$germline_pval)<x,1,0)},
                     function(x) {ifelse(filter_params$bb_rhoval > x, 1, 0)},
                     function(x) {sapply(1:nrow(COMB_mats$NV), assess_mean_depth, COMB_mats = COMB_mats, AUTO_low_depth_cutoff = x[1], AUTO_high_depth_cutoff = x[2], XY_low_depth_cutoff = x[3], XY_high_depth_cutoff=x[4])},
                     function(x) {ifelse(filter_params$pval_within_pos_dp2 > x, 1, 0)},
                     function(x) {ifelse(filter_params$pval_within_pos_dp3 > x, 1, 0)},
                     function(x) {sapply(1:nrow(COMB_mats$NV), assess_max_depth_in_pos, COMB_mats = COMB_mats, min_depth_auto = x[1], min_depth_xy = x[2])},
                     function(x) {ifelse(filter_params$max_pval_in_pos_sample > x, 1, 0)},
                     function(x) {sapply(1:nrow(COMB_mats$NV), assess_max_vaf, COMB_mats = COMB_mats, min_vaf_auto = x[1], min_vaf_xy = x[2])}
  )
  
  #CHECK THE INPUT DATA is all present for the set filters
  if(any(!c("NR","NV")%in%names(COMB_mats))|!all.equal(dim(COMB_mats$NR),dim(COMB_mats$NV))) {
    stop(return("COMB_mats object must contain matched NR and NV matrices of equal dimensions"))
  }
  
  if(dim(filter_params)[1]!=dim(COMB_mats$NV)[1]) {
    stop(return("The filter_params data frame and the COMB_mats$NV and NR matrices must contain the same number of mutations"))
  }
  #For any filter that has a set parameter, check that the corresponding variable is included in the filter_params data frame. If not, stop & return an error.
  for(i in 1:length(filter_name)) {
    if(!is.na(set_parameter[[i]][1]) & !any(names(filter_params)==var_name[i])) {
      stop(print(paste(filter_name[i],"needs a variable named",var_name[i],"in the filter_params data frame. Update the filter_params data frame or set the appropriate argument to NA")))
    }
  }
  #If using either 'pval for true somatic' parameter, need to include the PVal matrix in COMB_mats
  if((!is.na(min_pval_for_true_somatic_SHARED[1])|!is.na(min_pval_for_true_somatic[1]))&!any(names(COMB_mats)=="PVal")) {
    stop(return("If using the 'min_pval_for_true_somatic' or 'min_pval_for_true_somatic_SHARED' parameters, the COMB_mats list must contain a matrix called PVal"))
  }
  
  #Apply each of the filter functions, for those with NULL parameters, a vector of 1's is returned (i.e. all mutations "pass" the filter)
  out=mapply(FUN=function(param,test_function,var_name,filter_name) {
    if(is.na(param[1])) {
      return(rep(1,nrow(filter_params)))
    } else if(!any(names(filter_params)==var_name)){
      stop(return(paste(filter_name,"needs a variable named:",var_name,"in the filter_params data frame")))
    } else {
      return(test_function(param))
    }
  },
  param=set_parameter,
  test_function=test_function,
  var_name=var_name,
  filter_name=filter_name,
  SIMPLIFY = F)
  
  filter_pass=Reduce(cbind,out);rownames(filter_pass)<-rownames(filter_params);colnames(filter_pass)<-var_name #Combine the output & name the rows
  select_muts = (apply(filter_pass,1, function(x) all(x == 1))&!rownames(filter_pass)%in%exclude_muts)|rownames(filter_pass) %in% retain_muts #Select the mutations for output. These must pass ALL the applied filters.
  filter_code = apply(filter_pass, 1, paste, collapse = "-") #Save a "filter_code" vector. This can be used as a quick test for which filter is removing most mutations.
  COMB_mats.tree.build = list_subset(COMB_mats, select_vector = select_muts) #Subset matrices to include only the PASS mutations
  
  #Set the rownames to the mut_ref, and colnames to the sample names (without MTR or DEP)
  rownames(COMB_mats.tree.build$NV) = rownames(COMB_mats.tree.build$NR) = rownames(COMB_mats.tree.build$PVal) <- COMB_mats.tree.build$mat$mut_ref
  colnames = gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats.tree.build$NV))
  colnames(COMB_mats.tree.build$NV) = colnames(COMB_mats.tree.build$NR) = colnames(COMB_mats.tree.build$PVal) <- colnames
  
  #BUILD THE GENOTYPE MATRIX
  #Select the "definite positive" samples by setting genotype to 1
  #First build individual matrices, the same dimensions as the NV matrix, where each cell is 1 if it passes that criteria, or 0 if not. If criteria is NULL, set to 1.
  if(!is.na(min_variant_reads_SHARED[1])) {min_variant_reads_mat <- COMB_mats.tree.build$NV >= min_variant_reads_SHARED} else {min_variant_reads_mat <- 1}
  if(!is.na(min_pval_for_true_somatic_SHARED[1])) {min_pval_for_true_somatic_mat <- COMB_mats.tree.build$PVal > min_pval_for_true_somatic_SHARED} else {min_pval_for_true_somatic_mat <- 1}
  if(!is.na(min_vaf_SHARED[1]) & gender == "female") {
    depth_no_zero = COMB_mats.tree.build$NR
    depth_no_zero[depth_no_zero == 0] <- 1
    min_vaf_mat <- COMB_mats.tree.build$NV/depth_no_zero > min_vaf_SHARED[1]
  } else if(!is.na(min_vaf_SHARED[1]) & gender == "male") {
    min_vaf_mat = matrix(0, ncol = ncol(COMB_mats.tree.build$NV), nrow = nrow(COMB_mats.tree.build$NV))
    xy_muts = COMB_mats.tree.build$mat$Chrom %in% c("X","Y")
    depth_no_zero = COMB_mats.tree.build$NR
    depth_no_zero[depth_no_zero == 0] <- 1
    min_vaf_mat[xy_muts,] <- COMB_mats.tree.build$NV[xy_muts,]/depth_no_zero[xy_muts,] > min_vaf_SHARED[2]
    min_vaf_mat[!xy_muts,] <- COMB_mats.tree.build$NV[!xy_muts,]/depth_no_zero[!xy_muts,] > min_vaf_SHARED[1]
  } else {min_vaf_mat <- 1}
  
  COMB_mats.tree.build$Genotype_bin = min_variant_reads_mat * min_pval_for_true_somatic_mat * min_vaf_mat
  
  #Select the "not sure" samples by setting genotype to 0.5.  THIS IS THE ONLY SLIGHTLY OPAQUE BIT OF THIS FUNCTION - SET EMPIRICALLY FROM EXPERIMENTATION.
  COMB_mats.tree.build$Genotype_bin[COMB_mats.tree.build$NV > 0 & COMB_mats.tree.build$PVal > 0.01 & COMB_mats.tree.build$Genotype_bin != 1] <- 0.5 #If have any mutant reads, set as "?" as long as p-value > 0.01
  COMB_mats.tree.build$Genotype_bin[COMB_mats.tree.build$NV >= 3 & COMB_mats.tree.build$PVal > 0.001 & COMB_mats.tree.build$Genotype_bin != 1] <- 0.5 #If have high numbers of mutant reads, should set as "?" even if incompatible p-value (may be biased sequencing)
  COMB_mats.tree.build$Genotype_bin[(COMB_mats.tree.build$NV == 0) & (COMB_mats.tree.build$PVal > 0.05)] <- 0.5 #Essentially if inadequate depth to exclude mutation, even if no variant reads
  
  Genotype_shared_bin = COMB_mats.tree.build$Genotype_bin[rowSums(COMB_mats.tree.build$Genotype_bin == 1) > 1,]
  
  #Save the input parameters to a list
  params = list(input_set_ID = input_set_ID,
                gender = gender,
                retain_muts = retain_muts,
                exclude_muts = exclude_muts,
		germline_pval = germline_pval,
                rho = rho,
                mean_depth = mean_depth,
                pval_dp2= pval_dp2,
                pval_dp3= pval_dp3,
                min_depth = min_depth,
                min_pval_for_true_somatic = min_pval_for_true_somatic,
                min_vaf = min_vaf,
                min_variant_reads_SHARED = min_variant_reads_SHARED,
                min_pval_for_true_somatic_SHARED = min_pval_for_true_somatic_SHARED,
                min_vaf_SHARED=min_vaf_SHARED)
  
  #Save the summary stats of the run
  summary = data.frame(total_mutations = sum(select_muts),
                       total_SNVs = sum(COMB_mats.tree.build$mat$Mut_type == "SNV"),
                       total_INDELs = sum(COMB_mats.tree.build$mat$Mut_type == "INDEL"),
                       shared_mutations = nrow(Genotype_shared_bin),
                       shared_SNVs = sum(COMB_mats.tree.build$mat$mut_ref %in% rownames(Genotype_shared_bin) & COMB_mats.tree.build$mat$Mut_type == "SNV"),
                       shared_INDELs = sum(COMB_mats.tree.build$mat$mut_ref %in% rownames(Genotype_shared_bin) & COMB_mats.tree.build$mat$Mut_type == "INDEL"))
  #Extract the dummy dna_strings for tree building with MPBoot
  dna_strings = dna_strings_from_genotype(Genotype_shared_bin)
  
  #Print the run stats to the screen
  cat(summary$total_mutations,"total mutations\n")
  cat(summary$total_SNVs,"total SNVs\n")
  cat(summary$total_INDELs,"total INDELs\n")
  cat(summary$shared_mutations,"shared mutations\n")
  cat(summary$shared_SNVs,"shared SNVs\n")
  cat(summary$shared_INDELs,"shared INDELs\n")
  
  output = list(COMB_mats.tree.build = COMB_mats.tree.build, Genotype_shared_bin= Genotype_shared_bin, filter_code=filter_code, params = params, summary = summary, dna_strings = dna_strings)
  return(output)
}


#Function to create the dummy dna strings from the shared genotype matrix
dna_strings_from_genotype = function(genotype_mat) {
  Muts = rownames(genotype_mat)
  Ref = rep("W", length(Muts)) #W = Wild-type
  Alt = rep("V", length(Muts)) #V = Variant
  
  dna_strings = list()
  dna_strings[1] = paste(Ref, sep = "", collapse = "")
  
  for (k in 1:ncol(genotype_mat)){
    Mutations = Ref
    Mutations[genotype_mat[,k]==0.5] = '?'
    Mutations[genotype_mat[,k]==1] = Alt[genotype_mat[,k]==1]
    dna_string = paste(Mutations,sep="",collapse="")
    dna_strings[k+1]=dna_string
  }
  names(dna_strings)=c("Ancestral",colnames(genotype_mat))
  return(dna_strings)
}

#Function to create vcf files from the "mat" object
create_vcf_files = function(mat, select_vector = NULL) {
  if(is.null(select_vector)) {vcf_file = mat[,c("Chrom","Pos","Ref","Alt")]} else {vcf_file = mat[select_vector,c("Chrom","Pos","Ref","Alt")]}
  names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
  vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
  vcf_file = vcf_file[,c(1,2,8,3,4,7,6,5)]
  return(vcf_file)
}

#Function to extract the nodes from the first few divisions of the tree (set divisions - default = 2)
get_early_nodes = function(tree,divisions=2) {
  root_no = length(tree$tip.label)+1
  nodes=root_no
  early_nodes = NULL
  for(i in 1:divisions) {
    daughter_nodes = tree$edge[tree$edge[,1] %in% nodes,2]
    early_nodes=c(early_nodes,daughter_nodes)
    nodes<-daughter_nodes
  }
  return(early_nodes)
}

#Function to calculate peak VAFs from sample mutations (after applying germline and beta-binomial filters) to screen for mixed colonies
check_peak_vaf = function(sample, COMB_mats, filter_params,rho_cutoff=0.3) {
  colnames(COMB_mats$NV) <- gsub(pattern = "_MTR", replacement = "",x = colnames(COMB_mats$NV))
  dens <- density((COMB_mats$NV/COMB_mats$NR)[COMB_mats$NV[,sample] >=2 & !COMB_mats$mat$Chrom %in% c("X","Y") & log10(filter_params$germline_pval) <(-10) & filter_params$bb_rhoval > rho_cutoff ,sample])
  return(dens$x[which.max(dens$y)])
}

#Function to visualize VAF plots in similar way to the above
vaf_density_plot = function(sample, COMB_mats, filter_params,rho_cutoff=0.3) {
  colnames(COMB_mats$NV)=colnames(COMB_mats$NR)=gsub(pattern = "_MTR", replacement = "",x = colnames(COMB_mats$NV))
  sample_mean_depth = mean(COMB_mats$NR[,sample])
  dens <- density((COMB_mats$NV/COMB_mats$NR)[COMB_mats$NV[,sample] >=2 & !COMB_mats$mat$Chrom %in% c("X","Y") & log10(filter_params$germline_pval) <(-10) & filter_params$bb_rhoval > rho_cutoff ,sample])
  plot(dens, main = sample,xlim=c(-0.05,1.05))
  abline(v = dens$x[which.max(dens$y)])
  text(0.7, max(dens$y) - 0.2, paste("Peak VAF dens=",round(dens$x[which.max(dens$y)], digits = 2),"\nMean coverage=",round(sample_mean_depth,digits =2)), col = "red", cex = 0.7)
}

vaf_density_plot_final=function(sample,tree,COMB_mats){
  node <- which(tree$tip.label==sample)
  sample_muts <- COMB_mats$mat$mut_ref[COMB_mats$mat$node %in% get_ancestral_nodes(node,tree$edge)]
  COMB_mats$NR[COMB_mats$NR == 0] <- 1
  dens <- density((COMB_mats$NV/COMB_mats$NR)[sample_muts,sample])
  plot(dens,xlim = c(0,1),main=sample)
  abline(v = dens$x[which.max(dens$y)])
  text(0.7, max(dens$y) - 0.2, paste("Peak VAF dens=",round(dens$x[which.max(dens$y)], digits = 2)),col="red",cex = 0.7)
}



#Functions to import the cgpVAF output matrices neatly into R in the format that is used in my scripts
import_cgpvaf_output=function(cgpvaf_output_file,ref_ID="PDv37is") {
  mat<-read.delim(cgpvaf_output_file,stringsAsFactors = FALSE)
  mat<- mat[,!grepl(ref_ID,colnames(mat))] #Remove the reference sample columns
  NV<- mat[,grepl("_MTR", colnames(mat))]; NR <- mat[,grepl("_DEP",colnames(mat))]
  colnames(NR)=colnames(NV)=gsub(pattern = "_MTR",replacement = "", colnames(NV))
  mat$mut_ref=paste(mat$Chrom,mat$Pos,mat$Ref,mat$Alt,sep = "-")
  mat<-mat[,c("Chrom","Pos","Ref","Alt","mut_ref")]
  rownames(NV)=rownames(NR)=mat$mut_ref
  combined_mats=list(mat,NV,NR); names(combined_mats) <-c("mat","NV","NR")
  return(combined_mats)
}

import_cgpvaf_SNV_and_INDEL = function(SNV_output_file,INDEL_output_file=NULL,ref_ID="PDv37is") {
  #Import cgpVAF snp output file for the single-cell colonies, create the mut_ref column, and extract the mut and dep cols
  SNV_mats=import_cgpvaf_output(SNV_output_file)
  SNV_mats$mat$Mut_type="SNV"
  if(!is.null(INDEL_output_file)) {
    INDEL_mats = import_cgpvaf_output(INDEL_output_file,ref_ID = ref_ID)
    INDEL_mats$mat$Mut_type="INDEL"
    
    #Only include samples that are included in both SNV and INDEL cgpVAF output
    samples_in_both=intersect(colnames(SNV_mats$NV),colnames(INDEL_mats$NV))
    print(paste(length(samples_in_both),"samples in both SNV and INDEL cgpVAF output matrices, and will be combined"))
    SNV_mats$NV<-SNV_mats$NV[,samples_in_both]
    SNV_mats$NR<-SNV_mats$NR[,samples_in_both]
    INDEL_mats$NV<-INDEL_mats$NV[,samples_in_both]
    INDEL_mats$NR<-INDEL_mats$NR[,samples_in_both]

    combined_mats=mapply(SNV_mats,INDEL_mats,FUN = rbind) #Bind the indel and snp matrices together
    return(combined_mats)
  } else {
    return(SNV_mats)
  }
}

#Function to split up the imported output from VAGRENT in a neat way
split_vagrent_output = function(df,split_col,col_IDs = c("Gene","Transcript","RNA","CDS","Protein","Type","SO_codes")) {
  col = df[[split_col]]
  output = matrix(nrow = nrow(df), ncol = length(col_IDs))
  for(i in 1:length(col_IDs)) {
    output[,i] = str_split(col, pattern = "\\|", simplify = TRUE)[,i]
  }
  colnames(output) = col_IDs
  output<-as.data.frame(output,stringsAsFactors=F)
  return(output)
}

#Function to check for mutations that have been called as germline that are in fact absent in a clade
#Run using: res=check_for_false_germline_calls(tree,COMB_mats = COMB_mats, filter_params=filter_params)
is.snv=function(mut_ref) {
  sub=stringr::str_split(mut_ref,pattern = "-",simplify=T)[3:4]
  res=ifelse(nchar(sub[1])==1 & nchar(sub[2])==1,T,F)
  return(res)
}

check_for_false_germline_calls = function(tree,
                                          COMB_mats,
                                          filter_params,
                                          max_clade_prop=0.1, #the cutoff size (proportion of total samples included in clade) to test the clade for absent germline mutations. At >10% the germline filter is unlikely to wrongly remove mutations.
                                          SNVs_only=T, #Only re-add SNVs (indels are likely to be high frequency artefacts)
                                          CN_table=NULL
) {
  #Pull out the root clades
  get_root_clades=function(tree) {
    tree=di2multi(tree)
    ROOT=tree$edge[1,1]
    clades=tree$edge[tree$edge[,1]==ROOT,2]
    root_clade_samples=lapply(clades,function(node) return(getTips(tree,node)))
    return(root_clade_samples)
  }
  
  #Drop the ancestral tip if present, as this messes up the "get_root_clades" function
  tree.noancestral<-drop.tip(tree,"Ancestral")
  
  root_clades=get_root_clades(tree.noancestral)
  nsamp=length(tree.noancestral$tip.label)
  which_small=sapply(root_clades,length)<max_clade_prop*nsamp
  root_clades[!which_small]<-NULL
  
  if(length(root_clades)>0) {
    res<-lapply(root_clades,function(outlier_sample_group) {
      print(paste("Testing outlier group:",paste(outlier_sample_group,collapse=" ")))
      #Ensure all names are consistent
      rownames(COMB_mats$NR)=rownames(COMB_mats$NV)=rownames(COMB_mats$PVal)<-COMB_mats$mat$mut_ref
      colnames(COMB_mats$NR)=colnames(COMB_mats$PVal)=colnames(COMB_mats$NV)<-gsub("_MTR","",colnames(COMB_mats$NV))
      
      #Select mutations that were filtered by the germline filter
      germline_filtered=rownames(filter_params)[log10(filter_params$germline_pval)>(-10)]
      
      cat(paste("initially there are",length(germline_filtered),"germline filtered mutations"),sep="\n")
      #Exclude those that are subject to a deletion/ LOH in any of the outlier sample group
      if(!is.null(CN_table)){
        CN_table_restricted<-CN_table[CN_table$Sample%in%outlier_sample_group & CN_table$Type%in%c("LOH","DEL"),]
        if(nrow(CN_table_restricted)>0){
          cat(paste("Excluding mutations at copy number sites in group",paste(outlier_sample_group,collapse=" ")),sep="\n")
          for(k in 1:nrow(CN_table_restricted)){
            muts_df=data.frame(mut_ref=germline_filtered)
            muts_df$Chrom=str_split(germline_filtered,pattern = "-",simplify = T)[,1]
            muts_df$Pos=as.integer(str_split(germline_filtered,pattern = "-",simplify = T)[,2])
            exclude<-(muts_df$Chrom==CN_table_restricted$Chrom[k] & muts_df$Pos>CN_table_restricted$Pos_min[k] & muts_df$Pos<CN_table_restricted$Pos_max[k])
            print(sum(exclude))
            cat(germline_filtered[exclude],sep = "\n")
            germline_filtered<-germline_filtered[!exclude]
          }
        }
      }
      cat(paste("After filtering there are",length(germline_filtered),"germline filtered mutations"),sep="\n")
 
      #Aggregate counts across an individual outlier sample/ sample group
      NR_outlier=apply(COMB_mats$NR[germline_filtered,outlier_sample_group,drop=F],1,sum)
      NV_outlier=apply(COMB_mats$NV[germline_filtered,outlier_sample_group,drop=F],1,sum)
      
      outlier_pvals=mapply(FUN=function(NV,NR) {if(NR==0){return(1)}else{binom.test(NV,NR,alternative="less")$p.value}},NV=NV_outlier,NR=NR_outlier)
      outlier_pval.adj=p.adjust(outlier_pvals,method = "BH")
      #hist(log10(outlier_pvals),breaks=50,main="Unadjusted p-values for mutations being present in outlier group") #Review the p-value histogram - any clear low outliers?
      
      #Test for germline filtered mutations that are likely to be absent (with Bon-Ferroni correction for multiple testing)
      mut_refs<-germline_filtered[outlier_pval.adj<0.05 & NV_outlier==0]
      if(SNVs_only & length(mut_refs)>0) {
        mut_refs<-mut_refs[sapply(mut_refs,is.snv)]
      }
      if(length(mut_refs)>0) {
        print(paste(mut_refs,"is convincingly absent in this group"))
        return(mut_refs)
      } else {
        print("There are no mutations called as germline that are robustly absent in this outlier group, though this would relies on adequate coverage")
        return(NULL)
      }
    })
    return(res)
  } else {
    print(paste0("There are no clades from the root that include <",max_clade_prop*100,"% of samples"))
    return(NULL)
  }
}

# #This is Nick's version of the function: adds the ancestral tip at the beginning (i.e. as position 1)
# add_ancestral_outgroup=function(tree,outgroup_name="Ancestral"){
#   tmp=tree$edge
#   N=length(tree$tip.label)
#   ##Renumber what was root->max+1
#   ##Renumber the node with the sameid as new root as max+2
#   renamedroot=max(tmp+1)
#   tmp=ifelse(tmp==N+1,renamedroot,tmp)
#   ##tmp[which(tmp[,1]==(N+1)),1]=renamedroot
#   tmp=ifelse(tmp==N+2,renamedroot+1,tmp)
#   ##Increment tips by 1
#   tmp[,2]=ifelse(tmp[,2]<=N,tmp[,2]+1,tmp[,2])
#   
#   tree$edge=rbind(matrix(c(N+2,N+2,renamedroot,1),ncol=2,byrow  = FALSE),tmp)
#   tree$edge.length=c(0,0,tree$edge.length)
#   
#   tree$tip.label=c(outgroup_name,tree$tip.label)
#   tree$Nnode=tree$Nnode+1
#   mode(tree$Nnode)="integer"
#   mode(tree$edge)="integer"
#   tree
# }

#This version of the function adds the ancestral tip at the end
add_ancestral_outgroup=function(tree,outgroup_name="Ancestral"){
  tmp=tree$edge
  N=length(tree$tip.label)
  newroot=N+2
  renamedroot=N+3
  ancestral_tip=N+1
  tmp=ifelse(tmp>N,tmp+2,tmp)
  
  tree$edge=rbind(c(newroot,renamedroot),tmp,c(newroot,ancestral_tip))
  tree$edge.length=c(0,tree$edge.length,0)
  
  tree$tip.label=c(tree$tip.label,outgroup_name)
  tree$Nnode=tree$Nnode+1
  mode(tree$Nnode)="integer"
  mode(tree$edge)="integer"
  tree
}

assign_mutations_to_branches=function(tree,filtered_muts,keep_ancestral=T,create_multi_tree=T,p.error.value=0.01,treefit_pval_cutoff=1e-3) {
  tree=drop.tip(tree,"Ancestral")
  if(!keep_ancestral) {
    print("Assigning mutation without an ancestral branch")
    tree <- multi2di(tree)
    tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
    
    #ASSIGN MUTATIONS TO THE TREE USING THE TREE_MUT PACKAGE
    df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape
    
    #Get matrices in order, and run the main assignment functions
    mtr = filtered_muts$COMB_mats.tree.build$NV; mtr = as.matrix(mtr)
    depth = filtered_muts$COMB_mats.tree.build$NR; depth = as.matrix(depth)
    p.error = sapply(df$samples,function(x) ifelse(x=="Ancestral",1e-6,p.error.value))
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
    
  } else {
    print("Assigning mutation with an ancestral branch")
    tree <- multi2di(tree)
    tree <- add_ancestral_outgroup(tree) #Re add the ancestral outgroup after making tree dichotomous - avoids the random way that baseline polytomy is resolved
    tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
    
    #ASSIGN MUTATIONS TO THE TREE USING THE TREE_MUT PACKAGE
    df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape
    
    #Get matrices in order, and run the main assignment functions
    mtr = filtered_muts$COMB_mats.tree.build$NV; mtr$Ancestral=0; mtr = as.matrix(mtr)
    depth = filtered_muts$COMB_mats.tree.build$NR; depth$Ancestral=10; depth = as.matrix(depth)
    p.error = sapply(df$samples,function(x) ifelse(x=="Ancestral",1e-6,p.error.value))
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  }
  
  if(create_multi_tree){
    print("Converting to a multi-furcating tree structure")
    tree$edge.length <- res$df$df$edge_length #Assign edge lengths from the initial res object
    #Maintain the dichotomy with the ancestral branch
    if(keep_ancestral) {
      ROOT=tree$edge[1,1]
      current_length<-tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]
      new_length<-ifelse(current_length==0,1,current_length)
      tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]<-new_length
    }
    tree<-di2multi(tree) #Now make tree multifurcating
    df = reconstruct_genotype_summary(tree) #Define df (data frame) for new treeshape
    
    #Re-run the mutation assignment algorithm from the new tree
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  }
  
  tree$edge.length <- res$df$df$edge_length #Assign edge lengths from the most recent res object
  res$tree<-tree #Add the tree to the res object
  
  #See how many mutations are "poor fit"
  poor_fit = res$summary$pval < treefit_pval_cutoff  #See how many mutations don't have read counts that fit the tree very well
  print(paste(sum(poor_fit),"mutations do not have read counts that fit any tree branch well"))
  return(res)
}


