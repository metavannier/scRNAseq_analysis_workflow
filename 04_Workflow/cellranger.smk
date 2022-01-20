SAMPLES = config["sample"]["sname"].split(',')
AGRR = config["cellranger"]["sagrr"].split(',')

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
        out_cellranger = report(expand(OUTPUTDIR + "01_cellranger/{samples}/outs/{samples}_web_summary.html", samples=SAMPLES), caption = REPORT + "cellranger_summary.rst", category="01 cell ranger"),
    params:
        samples = config["sample"]["sname"].split(','),
    singularity:
        CONTAINER + "cellranger.sif"
    message: 
        "Counting with cellranger"
    shell:
        """
        sample=({params.samples})
        len=${{#sample[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do cellranger count \
        --id=${{sample[$i]}} \
        --expect-cells={config[cellranger][cells]} \
        --fastqs={input.fastqs} \
        --sample=${{sample[$i]}} \
        --transcriptome={input.ref_cellranger} \
        --jobmode=local \
        --localcores={config[cellranger][localcores]} \
        --localmem={config[cellranger][localmem]}
        mv ${{sample[$i]}}/* 05_Output/01_cellranger/${{sample[$i]}}/
        rm -r ${{sample[$i]}}/
        mv 05_Output/01_cellranger/${{sample[$i]}}/outs/web_summary.html 05_Output/01_cellranger/${{sample[$i]}}/outs/${{sample[$i]}}_web_summary.html/
        done
        """

rule cellrangeraggr:
    input:
        out_cellranger = expand(OUTPUTDIR + "01_cellranger/{agrr}/outs/{samples}_web_summary.html", agrr=AGRR, samples=SAMPLES),
    output:
        aggrcsv = ROOTDIR + "/aggregation.csv",
    params:
        samples = config["cellranger"]["sagrr"].split(','),
        batch = config["cellranger"]["batch"].split(','),
        categorie = config["cellranger"]["categorie"].split(','),
        outfolder = config["sample"]["nproject"],
        normagrr = config["cellranger"]["normagrr"],
        outaggr = OUTPUTDIR + "01_cellranger/"
    singularity:
        CONTAINER + "cellranger.sif"
    message: 
        "Aggregate with cellranger"
    shell:
        """
        sample=({params.samples})
        batch=({params.batch})
        categorie=({params.categorie})
        aggrcsv=({output.aggrcsv})
        outfolder=({params.outfolder})
        normagrr=({params.normagrr})
        outaggr=({params.outaggr})
        echo "sample_id,molecule_h5,batch,statut" >> ${{aggrcsv}}
        len=${{#sample[@]}}
        for (( i=0; i<$len; i=i+1 ))
        do echo "${{sample[$i]}},05_Output/01_cellranger/${{sample[$i]}}/outs/molecule_info.h5,${{batch[$i]}},${{categorie[$i]}}" >> ${{aggrcsv}}
        done
        cellranger aggr \
        --id=${{outfolder}} \
        --csv=${{aggrcsv}} \
        --normalize=${{normagrr}}
        mkdir ${{outaggr}}${{outfolder}}
        mkdir ${{outaggr}}${{outfolder}}/outs
        mv ${{outfolder}}/outs/count/* ${{outaggr}}${{outfolder}}/outs/
        mv ${{outfolder}}/outs/web_summary.html ${{outaggr}}/${{outfolder}}/outs/web_summary.html
        rm -r ${{outfolder}}
        """