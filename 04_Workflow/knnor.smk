#----------------------------------------
# KNNOR : K-Nearest Neighbor
# OveRsampliing method : For imbalanced
# datasets
#----------------------------------------
rule knnor:
    input: 
        data_for_sims_output = expand(OUTPUTDIR + "03_sims/data_for_sims_output.txt"),
    
    output:
        knnor_output = expand(OUTPUTDIR +"03_sims/knnor_output.txt"),

    params:
        knnor_rule = config["rules"]["knnor_rule"],
        sample_id = config["reference_sims"]["sample_id"],
        matrix = config["reference_sims"]["output_name_ref_matrix"],
        metadata = config["reference_sims"]["output_name_ref_metadata"],
        class_label = config["knnor"]["class_label"],
        max_oversampling = config["knnor"]["max_oversampling"],
        cells_column = config["reference_sims"]["cells_column"],
    
    conda: 
        CONTAINER + "knnor.yaml"

    message:
        "Oversampling data with KNNOR"

    shell:
        """
        if [ {params.knnor_rule} = "FALSE" ]
        then
            touch {output.knnor_output}
        else
            
            python 03_Script/imbalanced_datasets.py {params.sample_id} {params.matrix} {params.metadata} {params.class_label} "{params.max_oversampling}" {params.cells_column}
            touch {output.knnor_output}
        fi
        """