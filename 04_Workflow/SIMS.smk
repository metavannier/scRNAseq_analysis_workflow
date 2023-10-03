#----------------------------------------
# SIMS : Scalable, Interpretable Models
# for Cell Annotation of large scale
# single-cell RNA-seq data
#----------------------------------------

rule sims:
    input:
        anndata_for_sims_output = expand(OUTPUTDIR + "03_sims/anndata_for_sims_output.txt"),

    output:
        sims_output = expand(OUTPUTDIR + "03_sims/output_sims.txt"),

    params:
        sims_rule = config["rules"]["sims_rule"],


        # A REECRiRE !
        # sims_h5ad_input = expand(REFERENCE + "SSp/sims_matrix10X_P30_WT_less_cells.h5ad"),
        # reference_0 = REFERENCE + "SSp/sims_matrix10X_P30_WT_less_cells.h5ad",
        # # reference_1 = REFERENCE + "SSp_SSs/P30_WT/sims_matrixSmart_seq_P30_WT.h5ad",
        # metadata_0 = REFERENCE + "SSp/ordered_metadata10X_P30_WT_less_cells.csv",
        # # metadata_1 = REFERENCE + "SSp_SSs/P30_WT/ordered_metadata_Smart-Seq_P30_WT.csv",
        # # matrix = REFERENCE + "SSp_SSs/P30_WT/sims_matrix_P30_WT.h5ad",
        # project_name = config["sims"]["project_name"],
        # class_label = config["sims"]["class_label"],
        # max_epoch = config["sims"]["max_epoch"],


    message:
        "Run SIMS for cell annotation"

    shell:
        """
        if [ {params.sims.rule} = "FALSE" ]
        then
            touch {output.sims.output}
        else
            export PATH="/home/lchabot/.local/bin:$PATH"
            export PATH="/scratch/lchabot/BIOINFO_PROJECT/workflow_scrnaseq_P5-30/.local/bin:$PATH"
            pip install --upgrade pip setuptools wheel
            pip install protobuf==3.20.*
            pip install --use-pep517 git+https://github.com/braingeneers/SIMS.git
            python 03_Script/sims.py {params.reference_0} {params.metadata_0} {params.project_name} {params.class_label} {params.max_epoch}
            touch {output.sims.output}
        fi
        """

# export Path => Chemin ou est localiser l'instalation python et pip