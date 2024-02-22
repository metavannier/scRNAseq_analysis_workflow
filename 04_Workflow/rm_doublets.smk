#----------------------------------------
# Doublet detection
#----------------------------------------

rule doublets:
    input:
        normalisation_output = expand(OUTPUTDIR + "02_seurat/normalisation_output.txt"),

    output:
        doublets_output = expand(OUTPUTDIR + "rm_doublet/doublets_output.txt"),
        doublets_report = report(expand(OUTPUTDIR + "rm_doublet/{sample_id}/{sample_id}_rm_doublet_report.html", sample_id = SAMPLE_ID), caption = REPORT + "rm_doublet.rst", category = "rm_doublet"),

    conda:
        CONTAINER + "rm_doublets.yaml"

    params:
        # seurat_data = expand(OUTPUTDIR + "02_seurat/{sample_id}/{sample_id}_filtered" + config["rm_doublets"]["normalisation"] + "_seurat_object.rds",sample_id = SAMPLE_ID),
        sample_id = expand("{sample_id.id}", sample_id = sample_id.itertuples()),
        normalisation = config["rm_doublets"]["normalisation"]

    message:
        "Run Doublet calls for removing the doublets"

    script:
        SCRIPTDIR + "rm_doublets_reports_compilation.R"