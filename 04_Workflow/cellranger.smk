# SAMPLES = config["sample"]["sname"].split(',')
# AGRR = config["cellranger"]["sagrr"].split(',')
# NPROJECT = config["sample"]["nproject"]


########################################
# Build the reference with cellranger  #
########################################

rule ref_cellranger:
    input:
        multiqc_output = expand(OUTPUTDIR + "00_clean/multiqc_output.txt"),

    output:
        ref_cellranger_output = expand(OUTPUTDIR + "01_cellranger/ref_cellranger_output.txt"),

    params:
        ref_cellranger = directory(REF + config["reference"]["ref_cellranger"]),
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
        """
        ref_cellranger_output=({output.ref_cellranger_output})
        {SCRIPTDIR}RefForCellranger.sh {params.link_ref_fasta} {params.link_ref_gff} {params.ref_cellranger} {params.path_ref} {params.ref_name} {params.ref_version}
		echo "ref_cellranger step is FINISH" > ${{ref_cellranger_output}}
        """

##########################
#     Demultiplexing     #
##########################

rule multi:
    input:
        ref_cellranger_output = expand(OUTPUTDIR + "01_cellranger/ref_cellranger_output.txt"),

    output:
        multiplexing_output = expand(OUTPUTDIR + "01_cellranger/multiplexing_output.txt"),

    singularity:
        CONTAINER + "cellranger.sif"

    message: 
        "Demultiplexing & Cell Ranger"

    shell:
        """
        multiplexing_output=({output.multiplexing_output})
        cellranger multi  \
        --id={config[multi][id]}  \
        --csv={config[multi][config]} 
	    echo "Demultiplexing step is FINISH" > ${{multiplexing_output}}
        """

# rule cellranger:
#     input:
#         ref_cellranger = REF + config["reference"]["ref_cellranger"] + "/",
#         fastqs = config["cellranger"]["fastqs"]
#     output:
#         out_cellranger = report(expand(OUTPUTDIR + "01_cellranger/{samples}/outs/{samples}_web_summary.html", samples=SAMPLES), caption = REPORT + "cellranger_summary.rst", category="01 cell ranger"),
#     params:
#         samples = config["sample"]["sname"].split(','),
#     singularity:
#         CONTAINER + "cellranger.sif"
#     message: 
#         "Counting with cellranger"
#     shell:
#         """
#         sample=({params.samples})
#         len=${{#sample[@]}}
#         for (( i=0; i<$len; i=i+1 ))
#         do cellranger count \
#         --id=${{sample[$i]}} \
#         --expect-cells={config[cellranger][cells]} \
#         --fastqs={input.fastqs} \
#         --sample=${{sample[$i]}} \
#         --transcriptome={input.ref_cellranger} \
#         --jobmode=local \
#         --localcores={config[cellranger][localcores]} \
#         --localmem={config[cellranger][localmem]}
#         mv ${{sample[$i]}}/outs/* 05_Output/01_cellranger/${{sample[$i]}}/outs/
#         rm -r ${{sample[$i]}}/
#         mv 05_Output/01_cellranger/${{sample[$i]}}/outs/web_summary.html 05_Output/01_cellranger/${{sample[$i]}}/outs/${{sample[$i]}}_web_summary.html
#         done
#         """

# rule cellrangeraggr:
#     input:
#         out_cellranger = expand(OUTPUTDIR + "01_cellranger/{samples}/outs/{samples}_web_summary.html", samples=SAMPLES),
#     output:
#         aggrcsv = ROOTDIR + "/aggregation.csv",
#         out_aggregate = report(expand(OUTPUTDIR + "01_cellranger/" + NPROJECT + "/outs/aggregate_web_summary.html"), caption = REPORT + "cellranger_summary.rst", category="01 cell ranger"),
#     params:
#         samples = config["cellranger"]["sagrr"].split(','),
#         batch = config["cellranger"]["batch"].split(','),
#         categorie = config["cellranger"]["categorie"].split(','),
#         outfolder = config["sample"]["nproject"],
#         normagrr = config["cellranger"]["normagrr"],
#         outaggr = OUTPUTDIR + "01_cellranger/"
#     singularity:
#         CONTAINER + "cellranger.sif"
#     message: 
#         "Aggregate with cellranger"
#     shell:
#         """
#         sample=({params.samples})
#         batch=({params.batch})
#         categorie=({params.categorie})
#         aggrcsv=({output.aggrcsv})
#         outfolder=({params.outfolder})
#         normagrr=({params.normagrr})
#         outaggr=({params.outaggr})
#         echo "sample_id,molecule_h5,batch,statut" >> ${{aggrcsv}}
#         len=${{#sample[@]}}
#         for (( i=0; i<$len; i=i+1 ))
#         do echo "${{sample[$i]}},05_Output/01_cellranger/${{sample[$i]}}/outs/molecule_info.h5,${{batch[$i]}},${{categorie[$i]}}" >> ${{aggrcsv}}
#         done
#         cellranger aggr \
#         --id=${{outfolder}} \
#         --csv=${{aggrcsv}} \
#         --normalize=${{normagrr}}
#         mv ${{outfolder}}/outs/count/* ${{outaggr}}${{outfolder}}/outs/
#         mv ${{outfolder}}/outs/web_summary.html ${{outaggr}}/${{outfolder}}/outs/aggregate_web_summary.html
#         rm -r ${{outfolder}}
#         """
