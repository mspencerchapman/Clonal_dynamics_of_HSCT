## THIS SCRIPT IS FOR IMPORTING THE CAVEMAN OUTPUT FOR FILTERING & PHYLOGENY BUILDING
# It is quite specific to the Wellcome Sanger Institute file structure
# Unlikely to be useful elsewhere

library(optparse)

option_list = list(
  make_option(c("-p", "--projects"), action="store", type='character', help="CanPipe project ID to import from"),
  make_option(c("-s", "--samples"), action="store", type='character', help="PD numbers of samples to import"),
  make_option(c("-o", "--output_dir"), action="store", default=".", type='character', help="Output directory for the filtered files"),
  make_option(c("-t", "--type"), action="store", default="SNVs", type='character', help="Whether to import caveman or pindel output files")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

## FIXED PATHS
filter_indels_path="/lustre/scratch126/casm/team154pc/ms56/my_programs/filter_indels_new.pl"
filter_pass_subs_path="/lustre/scratch126/casm/team154pc/ms56/my_programs/filter_pass_subs_new.pl"

##SET RUN_ID AND FILEPATH
project_IDs=unlist(strsplit(opt$p,split=","))
PD_IDs=unlist(strsplit(opt$s,split=","))
output_dir=opt$o
file_type=opt$t

print(paste0("Importing samples with PD numbers ",paste0(PD_IDs,collapse=" & "), " from project numbers ", paste0(project_IDs,collapse=" & ")))

#Define how many matching samples are in the projects
samples_by_project_list=lapply(project_IDs, function(project) {
  project_samples_list=lapply(PD_IDs,function(PD_ID) {
    all_project_samples=list.dirs(path=paste0("/nfs/cancer_ref01/nst_links/live/",project),full.names = F,recursive = F)
    PD_project_samples=grep(PD_ID,all_project_samples,value = T)
    return(PD_project_samples)
  })
  project_samples_vec=Reduce(c,project_samples_list)
  return(project_samples_vec)
})
all_matching_samples_in_projects=unlist(samples_by_project_list)
print(all_matching_samples_in_projects)
print(paste(length(all_matching_samples_in_projects),"matching samples found in the listed projects"))

#Define which files are ready to be imported by looking through each project for the specified PD numbers and file type
paths_by_project_list=lapply(project_IDs, function(project) {
  project_paths_list=lapply(PD_IDs,function(PD_ID) {
    all_project_samples=list.dirs(path=paste0("/nfs/cancer_ref01/nst_links/live/",project),full.names = F,recursive = F)
    PD_project_samples=grep(PD_ID,all_project_samples,value = T)
    PD_paths_vec=sapply(PD_project_samples,function(sample) {
      return(list.files(path=paste0("/nfs/cancer_ref01/nst_links/live/",project,"/",sample),pattern=ifelse(file_type=="SNVs","caveman_c.annot.vcf.gz$","pindel.annot.vcf.gz$"),full.names = T))
    })
    return(PD_paths_vec)
  })
  project_paths_vec=Reduce(c,project_paths_list)
  return(project_paths_vec)
})
ready=unlist(paths_by_project_list)
#print(ready)
print(paste(length(ready),"of these samples have annotated",ifelse(file_type=="SNVs","CAVEMAN","PINDEL"),"files ready for import"))

ready_split = strsplit(ready,split = "/",fixed=TRUE)
ready_samples=sapply(ready_split,function(x) x[7])
not_yet_annotated_samples<-all_matching_samples_in_projects[!all_matching_samples_in_projects%in%ready_samples]
print("The following samples are found but do not yet have annotated variant files")
print(not_yet_annotated_samples)
ready_projects=sapply(ready_split,function(x) x[6])
ready_df=data.frame(samples=as.character(ready_samples),projects=as.character(ready_projects))

#Now create a list of the files that have already been processed and have "pass_flags" output
already=list.files(output_dir,pattern = "*pass_flags")
already_split = strsplit(already,split = ".",fixed=TRUE)
already_samples=unlist(lapply(already_split,function(x) x[1]))
print(paste(length(already_samples),"of these samples have already been imported and filtered"))

#Use these two to create a list of only those samples that are ready, but haven't already been processed
new_df = ready_df[!ready_df$samples %in% already_samples,]

new_df$samples = as.character(new_df$samples)
new_df$projects = as.character(new_df$projects)

print(paste0("Found ",length(new_df$samples)," new annotated ",ifelse(file_type=="SNVs","CAVEMAN","PINDEL")," files that have not yet been imported and filtered"))

if(length(new_df$samples)>0) {
	#Use this to copy over the new_sample .gz files
	print("Copying files from the nst_links directory")
	if(file_type=="SNVs") {
		for(i in 1:nrow(new_df)) {
			file_path = paste0("/nfs/cancer_ref01/nst_links/live/",new_df$projects[i],"/",new_df$samples[i],"/",new_df$samples[i],".caveman_c.annot.vcf.gz")
			system(paste0("cp -u ",file_path," ",output_dir))
		}
	} else if (file_type=="indels") {
		for(i in 1:nrow(new_df)) {
			file_path = paste0("/nfs/cancer_ref01/nst_links/live/",new_df$projects[i],"/",new_df$samples[i],"/",new_df$samples[i],".pindel.annot.vcf.gz")
			system(paste0("cp ",file_path," ",output_dir))
		}
	}

#Write a list of the new_samples that have been imported
#writeLines(new_df$samples,"new_samples")

#Write a list of the vcf file names to use for the filter_subs.pl script
if(file_type=="SNVs") {
  new_samples_vcfs = paste0(new_df$samples,".caveman_c.annot.vcf")
} else if (file_type=="indels") {
  new_samples_vcfs = paste0(new_df$samples,".pindel.annot.vcf")
}

writeLines(new_samples_vcfs,paste0(output_dir,"/new_samples_vcfs")) #This is needed for the perl script

#Unzip the newly imported files
print("Starting to unzip newly imported files")
system("gunzip *.gz")

#Run the filter_pass_subs.pl script
print("Running script to filter out only the PASS variants")
current_wd<-getwd()
if(file_type=="SNVs") {
  setwd(output_dir)
  system(paste0("perl ",filter_pass_subs_path))
  setwd(current_wd)
} else if (file_type=="indels") {
  setwd(output_dir)
  system(paste0("perl ",filter_indels_path))
  setwd(current_wd)
}
}

if(file.exists(paste0(output_dir,"/new_samples_vcfs"))) {system(paste0("rm ",output_dir,"/new_samples_vcfs"))}
