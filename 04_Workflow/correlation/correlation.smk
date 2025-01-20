rule dotplot:
    input:
        umapAssignation_output = expand(OUTPUTDIR + "02_seurat/normalisation_output.txt"),

    output:
        merge_rds_output = expand(OUTPUTDIR + "05_correlation/dotplot_output.txt"),

    params:
        input_count_matrices = config["correlation"]["input_count_matrices"],
        gene_X = config["correlation"]["gene_X"],
        gene_Y = config["correlation"]["gene_Y"],
        time = config["correlation"]["time"],

    conda:
        CONTAINER + "dotplot.yaml"

    message:
        "Run the merging of rds files"

    script:
        SCRIPTDIR + "correlation/dotplot.py"
