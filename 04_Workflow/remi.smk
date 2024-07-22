rule remi:
    input:
        ref_cellranger_output = expand(OUTPUTDIR + "01_cellranger/ref_cellranger_output.txt"),

    output: 
        remi_output = expand(OUTPUTDIR + "remi_output.txt"),

    conda:
        CONTAINER + "preparation_sims.yaml"

    message:
        "TEST REMI"

    shell:
        """
        Rscript 03_Script/remi.R
        touch {output.remi_output}
        """