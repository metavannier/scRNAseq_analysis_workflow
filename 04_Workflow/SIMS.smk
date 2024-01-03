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
        sample_id = config["reference_sims"]["sample_id"],
        reference_matrix = config["reference_sims"]["output_name_ref_matrix"],
        project_name = config["sims"]["project_name"],
        class_label = config["sims"]["class_label"],
        num_workers = config["sims"]["num_workers"],
        monitor = config["sims"]["monitor"],
        patience = config["sims"]["patience"],
        max_epoch = config["sims"]["max_epoch"],
        matrix = config["reference_sims"]["output_name_matrix"],
        key = config["sims"]["key"],

    message:
        "Run SIMS for cell annotation"

    shell:
        """
        if [ {params.sims_rule} = "FALSE" ]
        then
            touch {output.sims_output}
        else
            pip install --use-pep517 git+https://github.com/braingeneers/SIMS.git
            python 03_Script/sims.py {params.sample_id} {params.reference_matrix} {params.project_name} {params.class_label} \
            {params.num_workers} {params.monitor} {params.patience} {params.max_epoch} {params.matrix} {params.key}
            touch {output.sims_output}
        fi
        """



    # Si ne marche pas de ma façon : A REFAIRE COMME CA !!!
    # shell:
    #     """
    #     if [ {params.sims.rule} = "FALSE" ]
    #     then
    #         touch {output.sims.output}
    #     else
    #         export PATH="/home/lchabot/.local/bin:$PATH"
    #         export PATH="/scratch/lchabot/BIOINFO_PROJECT/workflow_scrnaseq_P5-30/.local/bin:$PATH"
    #         pip install --upgrade pip setuptools wheel
    #         pip install protobuf==3.20.*
    #         pip install --use-pep517 git+https://github.com/braingeneers/SIMS.git
    #         python 03_Script/sims.py {params.reference_0} {params.metadata_0} {params.project_name} {params.class_label} {params.max_epoch}
    #         touch {output.sims.output}
    #     fi
    #     """

    # export Path => Chemin ou est localiser l'instalation python et pip