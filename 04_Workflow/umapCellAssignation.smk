rule umapCellAssignation:
    input:
        sims_prediction_output = expand(OUTPUTDIR + "03_sims/output_sims_prediction.txt"),

    output:
        umapAssignation_output = expand(OUTPUTDIR + "03_sims/umapAssignation_output.txt"),
        umap_sims_report = report(expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_sims.pdf", sample_id = SAMPLE_ID), caption = REPORT + "umap_sims.rst", category = "03 sims"),
        umap_sims_threshold_report = report(expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_sims_threshold.pdf", sample_id = SAMPLE_ID), caption = REPORT + "umap_sims.rst", category = "03 sims"),
        umap_per_labels_report = report(expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_per_labels.pdf", sample_id = SAMPLE_ID), caption = REPORT + "umap_sims.rst", category = "03 sims"),

    params:
        sample_id = config["reference_sims"]["sample_id"],
        matrix = config["reference_sims"]["output_name_matrix"],
        threshold = config["sims"]["threshold"],
        rds_assignation = config["sims"]["rds_assignation"],
        pred_filtered = config["sims"]["pred_filtered"],

    conda:
        CONTAINER + "seurat.yaml"

    message:
        "Cells assignation on UMAP"

    script:
        SCRIPTDIR + "umapCellAssignation.R"
