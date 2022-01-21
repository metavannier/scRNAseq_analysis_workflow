rule bcf:
    input:
        ref = REF + config["reference"]["ref_cellranger"] + "/fasta/genome.fa",
        bam = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/possorted_genome_bam.bam"
    output:
        bcf = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.bcf"
    conda:
        CONTAINER + "demuxlet.yaml"
    message: 
        "Run bcftools for demuxlet"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} {input.bam} | bcftools call -mv -Ob -o {output.bcf}
        """

rule demuxlet:
    input:
        bam = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/possorted_genome_bam.bam",
        demuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.bcf"
    output:
        demuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines"
    conda:
        CONTAINER + "demuxlet.yaml"
    params:
        spooled = config["demuxlet"]["spooled"],
    message: 
        "Run demuxlet to deconvolute sample identity and identify multiplets"
    shell:
        """
        if [ {params.spooled} = "TRUE" ]
        then
            demuxlet --sam {input.bam} --vcf {input.demuxlet} --field GT --geno-error 0.01 --alpha 0 --alpha 0.5 --out {output.demuxlet}
        else
            echo "tmp" > {output.demuxlet}
        fi
        """

# rule assignement:
#     input:
#         tsne = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/analysis/tsne/2_components/projection.csv",
#         demuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.best"
#     output:
#         tabdemuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.tsv"
#     conda:
#         CONTAINER + "demuxlet.yaml"
#     message: 
#         "Plot sample assignments and doublet predictions"
#     script:
#         SCRIPTDIR + "assignment.R"
