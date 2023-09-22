
# AGRR = config["cellranger"]["sagrr"].split(',') => A voir 

###########################################################################################
# optimizing and assembling genome annotations for 3’ single-cell RNA-sequencing analysis #
###########################################################################################

###MUST BE IMPROVED (Need to do a mapping before and change the header from the fasta reference ensembl)

rule reference_enhancer:
    input:
        multiqc_output = expand(OUTPUTDIR + "00_clean/multiqc_output.txt"),

    output: 
        reference_enhancer_output = expand(OUTPUTDIR + "01_cellranger/reference_enhancer_output.txt"),

    params:
        ref_gtf = config["reference"]["ref_gtf"],
        reference_enhancer_rule = config["rules"]["reference_enhancer_rule"],

    conda:
        CONTAINER + "ReferenceEnhancer.yaml"

    shell:
        """
        if [ {params.reference_enhancer_rule} = "FALSE" ]
        then
            touch {output.reference_enhancer_output}
        else
            chmod +x {SCRIPTDIR}ReferenceEnhancer.R
            {SCRIPTDIR}ReferenceEnhancer.R {params.ref_gtf}
            touch {output.reference_enhancer_output}
        fi
        """    

########################################
# Build the reference with cellranger  #
########################################

rule ref_cellranger:
    input:
        multiqc_output = expand(OUTPUTDIR + "00_clean/multiqc_output.txt"),

    output:
        ref_cellranger_output = expand(OUTPUTDIR + "01_cellranger/ref_cellranger_output.txt"),

    params:
        ref_cellranger_rule = config["rules"]["ref_cellranger_rule"],
        ref_cellranger = directory(REF + config["reference"]["ref_cellranger"]),
        link_ref_fasta = config["reference"]["link_ref_fasta"],
        link_ref_gtf = config["reference"]["link_ref_gtf"],
        path_ref = REF,
        ref_name = config["reference"]["ref_cellranger"],
        ref_version = config["reference"]["ref_version"],

    singularity:
        CONTAINER + "cellranger.sif"

    message: 
        "Building the reference transcriptome"

    shell:
        """
        if [ {params.ref_cellranger_rule} = "FALSE" ]
        then
            touch {output.ref_cellranger_output}
        else
            {SCRIPTDIR}RefForCellranger.sh {params.link_ref_fasta} {params.link_ref_gtf} {params.ref_cellranger} {params.path_ref} {params.ref_name} {params.ref_version}
            cd ../
            touch {output.ref_cellranger_output}
        fi
        """

##########################
#    Cell ranger count   #
##########################

rule cellranger:
    input:
        ref_cellranger_output = expand(OUTPUTDIR + "01_cellranger/ref_cellranger_output.txt"),
    
    output:
        cellranger_output = expand(OUTPUTDIR + "01_cellranger/cellranger_output.txt"),
    
    params:
        sample_id = expand("{sample_id.id}", sample_id = sample_id.itertuples()),
        sample_name = expand("{sample_name.name}", sample_name = sample_name.itertuples()),
        fastqs = config["cellranger"]["fastqs"],
        ref_cellranger = REF + config["reference"]["ref_cellranger"] + "/",

    singularity: 
        CONTAINER + "cellranger.sif"

    message:
        "Cell Ranger Count"

    shell:
        """
        id=({params.sample_id})
        name=({params.sample_name})
        len=${{#id[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do cellranger count \
        --id=${{id[$i]}} \
        --fastqs={params.fastqs} \
        --transcriptome={params.ref_cellranger} \
        --sample=${{name[$i]}} \
        --expect-cells={config[cellranger][expect_cells]} \
        --localcores={config[cellranger][localcores]} \
        --localmem={config[cellranger][localmem]}
        mv ${{id[$i]}}/outs/* 05_Output/01_cellranger/${{id[$i]}}/outs/
        rm -r ${{id[$i]}}/
        mv 05_Output/01_cellranger/${{id[$i]}}/outs/web_summary.html 05_Output/01_cellranger/${{id[$i]}}/outs/${{id[$i]}}_web_summary.html
        touch {output.cellranger_output}
        done
        """