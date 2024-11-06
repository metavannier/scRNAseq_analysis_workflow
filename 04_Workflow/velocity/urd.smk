rule merge_rds:
    input:
        umapAssignation_output = expand(OUTPUTDIR + "03_sims/umapAssignation_output.txt"),

    output:
        merge_rds_output = expand(OUTPUTDIR + "04_velocity/merge_rds_output.txt"),

    conda:
        CONTAINER + "merge.yaml"

    params:
        sample_id = expand("{sample_id.id}", sample_id = sample_id.itertuples()),

    message:
        "Run the merging of rds files"

    script:
        SCRIPTDIR + "merge/merge_rds.R"

rule urd:
    input:
        merge_rds_output = expand(OUTPUTDIR + "04_velocity/merge_rds_output.txt"),

    output:
        urd_output = expand(OUTPUTDIR + "04_velocity/urd_output.txt"),
    
    params:
        urd_object = expand(OUTPUTDIR + "04_velocity/merged_urd_object.rds"),
        seed = config["urd"]["seed"],

    conda:
        CONTAINER + "urd.yaml"

    message:
        "Run URD package to analyse the data set"

    script:
        SCRIPTDIR + "velocity/urd_reports_compilation.R"