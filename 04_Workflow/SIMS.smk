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
        "Run SIMS for traning"

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
        # sims_prediction_report = report(expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_data_matrix_prediction.csv", sample_id = SAMPLE_ID), caption = REPORT + "label.rst", category = "03 sims"),
        # sims_prediction_unknown_report = report(expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_data_matrix_prediction_filtered.csv", sample_id = SAMPLE_ID), caption = REPORT + "label.rst", category = "03 sims"),

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

rule unknown_prediction:
    input:
        sims_prediction_output = expand(OUTPUTDIR + "03_sims/output_sims_prediction.txt"),

    output:
        unknown_prediction_output = expand(OUTPUTDIR + "03_sims/unknown_prediction_output.txt"),

    params:
        sample_id = config["reference_sims"]["sample_id"],
        matrix = config["reference_sims"]["output_name_matrix"],
        threshold = config["sims"]["threshold"],
        pred_filtered = config["sims"]["pred_filtered"],

    conda:
        CONTAINER + "eval_prediction.yaml"

    message:
        "Run SIMS for cell unknown prediction"

    script:
        SCRIPTDIR + "unknown_prediction.R"

rule evaluate_prediction:
    input:
        unknown_prediction_output = expand(OUTPUTDIR + "03_sims/unknown_prediction_output.txt"),

    output:
        evaluate_prediction_output = expand(OUTPUTDIR + "03_sims/evaluate_prediction_output.txt"),

    params:
        mousegastrulation_samples = config["reference_sims"]["mousegastrulation_samples"].split(','),
        sample_id = config["reference_sims"]["sample_id"],
        matrix = config["reference_sims"]["output_name_matrix"],
        pred_filtered = config["sims"]["pred_filtered"],
        confusion_matrix = config["sims"]["confusion_matrix"],

    conda:
        CONTAINER + "eval_prediction.yaml"

    message:
        "Run SIMS for cell annotation evaluation"

    script:
        SCRIPTDIR + "evaluate_prediction.R"