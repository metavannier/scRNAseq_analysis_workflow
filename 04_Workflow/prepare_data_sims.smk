### a voir pour le shell

#----------------------------------------
# Data preparation for SIMS
#----------------------------------------

rule data_for_sims:
    input:
        seurat_output = expand(OUTPUTDIR + "02_seurat/seurat_output.txt"),
    
    output: 
        data_for_sims_output = expand(OUTPUTDIR + "03_sims/data_for_sims_output.txt"),

    params:
        # General
        r_script = config["reference_sims"]["r_script"],
        sample_id = config["reference_sims"]["sample_id"],
        # Arlotta
        arlotta_metadata = config["reference_sims"]["arlotta_metadata"],
        arlotta_matrix = config["reference_sims"]["arlotta_matrix"],
        arlotta_cells = config["reference_sims"]["arlotta_cells"],
        arlotta_features = config["reference_sims"]["arlotta_features"],
        pattern_to_keep = config["reference_sims"]["pattern_to_keep"],
    
    conda:
        CONTAINER + "sims.yaml"

    message:
        "Create matrix for sims (CSV format)"

    script:
        "/media/lea/data/inmed_dechevigny_scrnaseq_cortex/03_Script/arlotta.R"
        



# #----------------------------------------
# # Use the csv create in the precedent
# # rule to create an anndata file 
# # for SIMS
# #----------------------------------------
# rule anndata_for_sims:
#     input:
#         data_for_sims_output = expand(OUTPUTDIR + "03_sims/data_for_sims_output.txt"),

#     output:
#         anndata_for_sims_output = expand(OUTPUTDIR + "03_sims/anndata_for_sims_output.txt"),
    
#     params:
#         matrix = config[""][""],
        