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
    
    conda: 
        CONTAINER + "knnor.yaml"

    message:
        "Oversampling data with KNNOR"

    shell:
        """
            python 03_Script/imbalenced_datasets.py
            touch {output.knnor_output}
        """