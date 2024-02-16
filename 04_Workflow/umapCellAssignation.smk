rule umapCellAssignation:
    input:
        sims_output = expand(OUTPUTDIR + "03_sims/output_sims.txt"),

    output:
        umapAssignation_output = expand(OUTPUTDIR + "03_sims/umapAssignation_output.txt"),

    params:
        sample_id = config["reference_sims"]["sample_id"],
        matrix = config["reference_sims"]["output_name_matrix"],
    
    conda :
        CONTAINER + "seurat.yaml"

    message:
        "Cells assignation on UMAP"

    script:
        SCRIPTDIR + "umapCellAssignation.R"
