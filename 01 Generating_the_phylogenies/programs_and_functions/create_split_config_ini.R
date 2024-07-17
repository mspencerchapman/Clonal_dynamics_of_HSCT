library(optparse)

option_list = list(
  make_option(c("-p", "--project"), action="store", type='character', help="Project ID of config file to use as template"),
  make_option(c("-v", "--vcf_folder"), action="store", type='character', help="Path of folder containing all vcf files to include"),
  make_option(c("-b", "--batch_size"), action="store", default=10,type='numeric', help="Number of samples to include per batch"),
  make_option(c("-s", "--samples_file"), action="store", default=NULL,type='character', help="Number of samples to include per batch")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

##SET RUN_ID AND FILEPATH
project_IDs=unlist(strsplit(opt$p,split=","))
vcf_folder=ifelse(is.null(opt$v),"../../",opt$v)
batch_size=opt$b
samples_path=opt$s

config_file = readLines(paste0(project_IDs[1],"_cgpVafConfig.ini"))
header=config_file[1]
footer=tail(config_file,7)
Ref_insilico_id=strsplit(config_file[2],split="_")[[1]][1]


if(!is.null(samples_path)){
	all_IDs=readLines(samples_path)
} else {
	pass_flag_files=list.files(path=vcf_folder,pattern = "*pass_flags")
	all_IDs=sapply(strsplit(pass_flag_files,split="\\."),function(x) x[1])
}


n_split=ceiling(length(all_IDs)/batch_size)
remainder=length(all_IDs)%%batch_size
cat_split=NULL
for(j in 1:n_split){
  if(j==n_split & remainder>0){
    split=NULL
    split[1]=paste0(Ref_insilico_id,"_UNM",j,"= <<EOT")
    start_index=1+((j-1)*batch_size)
    split[2:(1+remainder)]=all_IDs[start_index:length(all_IDs)]
    split[(2+remainder)]="EOT"
    cat_split=c(cat_split,split)
  } else {
    split=NULL
    split[1]=paste0(Ref_insilico_id,"_UNM",j,"= <<EOT")
    start_index<-(1+((j-1)*batch_size))
    split[2:11]<-all_IDs[start_index:(start_index+9)]
    split[12]="EOT"
    cat_split=c(cat_split,split)
  }
}

final_file=c(header,cat_split,footer)
writeLines(final_file,paste0(project_IDs[1],"_cgpVafConfig_split.ini"))
  
