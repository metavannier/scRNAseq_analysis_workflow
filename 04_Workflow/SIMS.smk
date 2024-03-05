#----------------------------------------
# SIMS : Scalable, Interpretable Models
# for Cell Annotation of large scale
# single-cell RNA-seq data
#----------------------------------------

rule sims_training:
    input:
        anndata_for_sims_output = expand(OUTPUTDIR + "03_sims/anndata_for_sims_output.txt"),

    output:
        sims_training_output = expand(OUTPUTDIR + "03_sims/output_sims_training.txt"),

    params:
        sims_rule = config["rules"]["sims_rule"],
        sample_id = config["reference_sims"]["sample_id"],
        reference_matrix = config["reference_sims"]["output_name_ref_matrix_metadata"],
        project_name = config["sims"]["project_name"],
        class_label = config["sims"]["class_label"],
        num_workers = config["sims"]["num_workers"],
        monitor = config["sims"]["monitor"],
        patience = config["sims"]["patience"],
        max_epoch = config["sims"]["max_epoch"],
        matrix = config["reference_sims"]["output_name_matrix"],
        key = config["sims"]["key"],
        epoch = config["sims"]["epoch"],

    message:
        "Run SIMS for cell annotation"

    shell:
        """
        if [ {params.sims_rule} = "FALSE" ]
        then
            touch {output.sims_training_output}
        else
            pip install --use-pep517 git+https://github.com/braingeneers/SIMS.git@c3cc547e9223e979fdc70f9af2ae932b729da88e
            python 03_Script/sims_training.py {params.sample_id} {params.reference_matrix} {params.project_name} {params.class_label} \
            {params.num_workers} {params.monitor} {params.patience} {params.max_epoch} {params.matrix} {params.key} {params.epoch}
            touch {output.sims_training_output}
        fi
        """

rule sims_prediction:
    input:
        sims_training_output = expand(OUTPUTDIR + "03_sims/output_sims_training.txt"),

    output:
        sims_prediction_output = expand(OUTPUTDIR + "03_sims/output_sims_prediction.txt"),

    params:
        sims_rule = config["rules"]["sims_rule"],
        sample_id = config["reference_sims"]["sample_id"],
        matrix = config["reference_sims"]["output_name_matrix"],
        epoch = config["sims"]["epoch"],

    message:
        "Run SIMS for cell annotation"

    shell:
        """
        if [ {params.sims_rule} = "FALSE" ]
        then
            touch {output.sims_prediction_output}
        else
            pip install --use-pep517 git+https://github.com/braingeneers/SIMS.git@c3cc547e9223e979fdc70f9af2ae932b729da88e
            python 03_Script/sims_prediction.py {params.sample_id} {params.matrix} {params.epoch}
            touch {output.sims_prediction_output}
        fi
        """
