#----------------------------------------
# Data preparation for SIMS
#----------------------------------------

rule data_for_sims:
    input:
        seurat_output = expand(OUTPUTDIR + "02_seurat/seurat_output.txt"),
        normalisation_output = expand(OUTPUTDIR + "02_seurat/normalisation_output.txt")
    
    output: 
        data_for_sims_output = expand(OUTPUTDIR + "03_sims/data_for_sims_output.txt"),

    params:
        # sims_rule = config["rules"]["sims_rule"],
        # General
        # r_script = config["reference_sims"]["r_script"],
        sample_id = config["reference_sims"]["sample_id"],
        output_name_ref_metadata = config["reference_sims"]["output_name_ref_metadata"],
        output_name_ref_matrix = config["reference_sims"]["output_name_ref_matrix"],
        output_name_matrix = config["reference_sims"]["output_name_matrix"],
        norm_method = config["seurat"]["norm_method"],
        norm_scale_factor = config["seurat"]["norm_scale_factor"],
        ### Reference Arlotta
        # arlotta_metadata = config["reference_sims"]["arlotta_metadata"],
        # arlotta_matrix = config["reference_sims"]["arlotta_matrix"],
        # arlotta_cells = config["reference_sims"]["arlotta_cells"],
        # arlotta_features = config["reference_sims"]["arlotta_features"],
        # developmental_time = config["reference_sims"]["developmental_time"],
        ### Reference Allen mouse whole cortex
        # allen_metadata = config["reference_sims"]["allen_metadata"],
        # allen_matrix = config["reference_sims"]["allen_matrix"],
        # allen_genes = config["reference_sims"]["allen_genes"],
        # allen_cells = config["reference_sims"]["allen_cells"],
        # Reference MouseGastrulation
        sample_list = expand("{sample_id.id}", sample_id = sample_id.itertuples()),
        mousegastrulation_samples = config["reference_sims"]["mousegastrulation_samples"].split(','),
        mousegastrulation_matrix = config["reference_sims"]["mousegastrulation_matrix"],
        data_matrix = config["reference_sims"]["data_matrix"],

    conda:
        CONTAINER + "preparation_sims.yaml"

    message:
        "Create matrix for sims (CSV format)"

    script:
        SCRIPTDIR + config["reference_sims"]["r_script"]

#----------------------------------------
# Use the csv create in the precedent
# rule to create an anndata file 
# for SIMS
#----------------------------------------

rule anndata_for_sims:
    input:
        data_for_sims_output = expand(OUTPUTDIR + "03_sims/data_for_sims_output.txt"),

    output:
        anndata_for_sims_output = expand(OUTPUTDIR + "03_sims/anndata_for_sims_output.txt"),
    
    params:
        sims_rule = config["rules"]["sims_rule"],
        sample_id = config["reference_sims"]["sample_id"],
        output_name_ref_matrix = config["reference_sims"]["output_name_ref_matrix"],
        output_name_matrix = config["reference_sims"]["output_name_matrix"],
        output_name_ref_metadata = config["reference_sims"]["output_name_ref_metadata"],
        cells_column = config["reference_sims"]["cells_column"],
        output_name_ref_matrix_metadata = config["reference_sims"]["output_name_ref_matrix_metadata"],

    conda:
        CONTAINER + "preparation_sims.yaml"

    message:
        "Create matrix for sims (H5AD format)"
    
    shell:
        """
        if [ {params.sims_rule} = "FALSE" ]
        then
            touch {output.anndata_for_sims_output}
        else
            python 03_Script/to_anndata_file.py {params.sample_id} {params.output_name_ref_matrix} {params.output_name_matrix} {params.output_name_ref_metadata} {params.cells_column} {params.output_name_ref_matrix_metadata}
            touch {output.anndata_for_sims_output}
        fi
        """