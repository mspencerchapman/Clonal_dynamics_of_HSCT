using Phylo
using RCall
using Plots
using DataFrames
using Distributions
using Random

Random.seed!(28)
base_dir = "/Users/pc8/Documents/Parallel_sequencing/Normal_tissue/Blood/Stem_cell_transplantation/Scripts/" #UPDATE TO LOCAL DIRECTORY
include("$base_dir/src/Deep_seq_tree_GS.jl")

####################################################
# Get R data and trees into Julia
# First, call R and read RDS files in; break into components
R"""
    library(ape)
    baitset=2
    mut_matrices <- readRDS(paste0("data/targeted_sequencing_data/data_by_baitset/baitset",baitset,"/baitset",baitset,"_cgpvaf_matrices.RDS"))
    tree_info <- readRDS(paste0("data/targeted_sequencing_data/data_by_baitset/baitset",baitset,"/baitset",baitset,"_trees_and_muts.RDS"))
    tree38 <- tree_info[["Pair38"]][["tree"]]
    tree40 <- tree_info[["Pair40"]][["tree"]]
    tree41 <- tree_info[["Pair41"]][["tree"]]
    muts38 <- tree_info[["Pair38"]][["details"]]
    muts40 <- tree_info[["Pair40"]][["details"]]
    muts41 <- tree_info[["Pair41"]][["details"]]
    nv <- mut_matrices[["NV"]]
    nr <- mut_matrices[["NR"]]
    if(all(rownames(mut_matrices[["NV"]]) == rownames(mut_matrices[["NR"]])) &&
        all(colnames(mut_matrices[["NV"]]) == colnames(mut_matrices[["NR"]]))) {
        nv_rownames <- rownames(mut_matrices[["NV"]])
        nv_colnames <- colnames(mut_matrices[["NV"]])
    } else {
        stop("Names of matrices do not match")
    }
    bulk_smry=data.frame(Pair=rep(c("Pair38","Pair40","Pair41"),each=24),
                       PD_number=rep(c("PD45808","PD45809","PD45810","PD45811","PD45812","PD45813"),each=12),
                       sampleID=paste(paste0(rep(c("PD45808","PD45809","PD45810","PD45811","PD45812","PD45813"),each=12),rep(c("c","d","e"),each=4,times=6)),rep(c("lo0001","lo0002","lo0003","lo0004"),times=18),sep="_"),
                       tissueID=paste0(rep(c("PD45808","PD45809","PD45810","PD45811","PD45812","PD45813"),each=12),rep(c("c","d","e"),each=4,times=6)),
                       individual_type=rep(c("Donor","Recipient"),each=12,times=3),
                       cell_type=rep(c("B_cells","T_cells","Monocytes"),each=4,times=6))
"""

# Now transfer objects to Julia
@rget tree38 tree40 tree41 muts38 muts40 muts41 nv nr nv_rownames nv_colnames bulk_smry


####################################################
# Define the samples to work on and basic variables for the Gibbs sampler

muts_file_dict = Dict("Pair38" => muts38, "Pair40" => muts40, "Pair41" => muts41)
tree_file_dict = Dict("Pair38" => tree38, "Pair40" => tree40, "Pair41" => tree41)
pt_pairs = Dict("PD45808" => "PD45809", "PD45810" => "PD45811", "PD45812" => "PD45813",
                "PD45809" => "PD45808", "PD45811" => "PD45810", "PD45813" => "PD45812")
for curr_samp in unique(bulk_smry.tissueID)
    mkdir("$base_dir/output/$curr_samp")
    master_row = first(bulk_smry[bulk_smry.tissueID .== curr_samp,:])
    muts_file = muts_file_dict[master_row.Pair]
    patient_ID = String(master_row.PD_number)
    partner_ID = pt_pairs[patient_ID]
    cell_type = curr_samp
    tree = tree_file_dict[master_row.Pair]
    cell_name = master_row.cell_type
    indvdl_type = master_row.individual_type
    iter = 1000
    burn_in = 100
    thin = 10

    # Define a mapping for the nodes assigned in the mutation info DataFrame to the branch names of tree
    node_to_branch = Dict{Int64, LinkBranch}()
    for i in 1:nleaves(tree)
        node_to_branch[i] = getinbound(tree, getleaves(tree)[i])
    end
    for i in (nleaves(tree)+2):nnodes(tree) # Skipping root, which in R is always nleaves+1
        node_to_branch[i] = getinbound(tree, "Node $i")
    end

    ####################################################
    # Initialise basics of trees and samples
    # Which columns to use for the sample of interest and the seq error estimation
    test_cols = findall(occursin.(cell_type, nv_colnames))
    ctrl_cols = findall(.!(occursin.(patient_ID, nv_colnames) .|
                                occursin.(partner_ID, nv_colnames)))

    # Initialise the mutation fields and branch start and end VAFs for Gibbs sampler
    start_VAF = Dict{LinkBranch, Array{Float64, 1}}()
    end_VAF = Dict{LinkBranch, Array{Float64, 1}}()
    muts = Dict{LinkBranch, Array{Mutation, 1}}()

    for i in branchiter(tree)
        start_VAF[i] = isroot(tree, src(tree, i)) ? [1.0;] : [0.0;]
        end_VAF[i] = [0.0;]
        muts[i] = Array{Mutation, 1}()
        setbranchdata!(tree, i, "start_VAF", start_VAF[i])
        setbranchdata!(tree, i, "end_VAF", end_VAF[i])
        setbranchdata!(tree, i, "Muts", muts[i])
    end

    ####################################################
    # Map deep-sequenced mutations to branches
    for i in eachrow(muts_file)
        if i.mut_ref ∈ nv_rownames
            curr_row = findall(nv_rownames .== i.mut_ref)[1]
            obsV = sum(nv[curr_row, test_cols])
            obsR = sum(nr[curr_row, test_cols])
            ctrlV = sum(nv[curr_row, ctrl_cols])
            ctrlR = sum(nr[curr_row, ctrl_cols])
            seqerr = ctrlV + ctrlR < 20 ? 0.01 : (ctrlV + 0.5) / (ctrlV + ctrlR + 1)

            push!(muts[node_to_branch[i.node]],
                Mutation(i.mut_ref, obsV, obsR, obsV+obsR, seqerr, Array{Float64, 1}()))
        end
    end

    (split_tree, VAFs) = split_tree_by_mut(tree)
    temp = plot(split_tree, line_z=log10.(VAFs.+0.0001), linewidth=2, showtips=false,
        linecolor=:RdPu, title = "$patient_ID ($indvdl_type): $cell_name")
    savefig(temp, "$base_dir/output/$curr_samp/Raw_$(patient_ID)_$(indvdl_type)_$(cell_name).pdf")

    GS_out = deep_seq_GS(tree, start_VAF, end_VAF, muts, 20000, 10000, 100; scale_pm = 50)
    (split_tree, VAFs) = split_tree_by_mut(tree,
            FUN = x->median(x.posterior_VAF) < 1.5 / x.obs_depth ? 0 : median(x.posterior_VAF))
    temp = plot(split_tree, line_z=log10.(VAFs.+0.0001), linewidth=2, showtips=false,
               linecolor=:RdPu, title = "Estimated $patient_ID ($indvdl_type): $cell_name")
    savefig(temp, "$base_dir/output/$curr_samp/GS_$(patient_ID)_$(indvdl_type)_$(cell_name).pdf")

    trunc_tree = truncate_tree(tree, 50)
    (split_tree, VAFs) = split_tree_by_mut(trunc_tree,
               FUN = x->median(x.posterior_VAF) < 1.5 / x.obs_depth ? 0 : median(x.posterior_VAF))
    temp = plot(split_tree, line_z=log10.(VAFs.+0.0001), linewidth=2, showtips=false,
                  linecolor=:RdPu, title = "Estimated $patient_ID ($indvdl_type): $cell_name")
    savefig(temp, "$base_dir/output/$curr_samp/GS_$(patient_ID)_$(indvdl_type)_$(cell_name)_truncated.pdf")

    write_GS_output(tree, GS_out, "$base_dir/output/$curr_samp/", curr_samp, node_to_branch)

end
