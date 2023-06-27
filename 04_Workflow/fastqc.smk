# ----------------------------------------------
# FastQC : Check quality of reads files
# ----------------------------------------------

rule fastqc:
    input:
        expand(RAWDATA + "{rawsample}_{library}_{types}.fastq.gz", rawsample=RAWSAMPLE, library=LIBRARY, types=TYPES),

    output:
        fastqc_output = expand(OUTPUTDIR + "00_fastqc/fastqc_output.txt"),

    message:
        "Quality with fastqc"

    conda:
        CONTAINER + "fastqc.yaml"

    shell:
        """
        export PATH="/$PWD/.local/bin:$PATH"
		    fastqc_output=({output.fastqc_output})
        fastqc {input} --outdir {OUTPUTDIR}00_fastqc/
		    echo "FastQC step is FINISH" > ${{fastqc_output}}
        """

# ----------------------------------------------
# MultiQC to check the reads trimmed quality
# ----------------------------------------------

rule multiqc:
  input:
    fastqc_zip = expand(OUTPUTDIR + "00_fastqc/{rawsample}_{library}_{types}_fastqc.zip", rawsample=RAWSAMPLE, library=LIBRARY, types=TYPES)
  output:
    raw_multi_html = report(OUTPUTDIR + "00_fastqc/raw_multiqc.html", caption = REPORT + "multiqc.rst", category="00 quality report"), 
  params:
    multiqc_output_raw = OUTPUTDIR + "00_fastqc/raw_multiqc_data"
  conda:
    CONTAINER + "multiqc.yaml"
  benchmark:
    BENCHMARK + "multiqc.benchmark.txt"
  log:
    LOG + "multiqc.log"
  message: 
    "Synthetize quality with multiqc"
  shell: 
    """
    multiqc -n {output.raw_multi_html} {input.raw_qc} #run multiqc
    rm -rf {params.multiqc_output_raw} #clean-up
    """ 
