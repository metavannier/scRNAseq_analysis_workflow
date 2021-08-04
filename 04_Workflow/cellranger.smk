##########################
# Count with cellranger  #
##########################

rule ref_cellranger:
    output:
        ref_cellranger = directory(REF + config["reference"]["ref_cellranger"]),
        out_ref = REF + config["reference"]["ref_cellranger"] + "/fasta/genome.fa"
    params:
        link_ref_fasta = config["reference"]["link_ref_fasta"],
        link_ref_gff = config["reference"]["link_ref_gff"],
        path_ref = REF,
        ref_name = config["reference"]["ref_cellranger"],
        ref_version = config["reference"]["ref_version"]
    singularity:
        CONTAINER + "cellranger.sif"
    message: 
        "Building the reference transcriptome"
    shell:
        SCRIPTDIR + "RefForCellranger.sh {params.link_ref_fasta} {params.link_ref_gff} {output.ref_cellranger} {params.path_ref} {params.ref_name} {params.ref_version}"
        
rule cellranger:
    input:
        ref_cellranger = REF + config["reference"]["ref_cellranger"] + "/",
        fastqs = config["cellranger"]["fastqs"]
    output:
        out_cellranger = report(OUTPUTDIR + "00_cellranger/" + config["sample"]["sname"] + "/outs/web_summary.html", caption = ROOTDIR + REPORT + "cellranger_summary.rst", category="01 cell ranger"),
    singularity:
        CONTAINER + "cellranger.sif"
    message: 
        "Counting with cellranger"
    shell:
        """
        cellranger count \
        --id={config[sample][sname]} \
        --expect-cells={config[cellranger][cells]} \
        --fastqs={input.fastqs} \
        --sample={config[sample][sname]} \
        --transcriptome={input.ref_cellranger} \
        --jobmode=local \
        --localcores={config[cellranger][localcores]} \
        --localmem={config[cellranger][localmem]}
        mv {config[sample][sname]}/* 05_Output/00_cellranger/{config[sample][sname]}/
        """