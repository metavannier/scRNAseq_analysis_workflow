# ----------------------------------------------
# FastQC to check the row reads quality
# ----------------------------------------------

rule fastqc:
  input:
    RAWDATA + "{fastq}.fastq"
  output:
    html = OUTPUTDIR + "00_fastqc/{fastq}_fastqc.html",
    zip = OUTPUTDIR + "00_fastqc/{fastq}_fastqc.zip"
  message: 
    "Quality with fastqc"
  wrapper:
    "0.35.2/bio/fastqc"

# ----------------------------------------------
# MultiQC to check the reads trimmed quality
# ----------------------------------------------

rule multiqc:
  input:
    raw_qc = expand("{outputdir}00_fastqc/{fastq}_fastqc.zip", outputdir=OUTPUTDIR, fastq=FASTQ),
  output:
    raw_multi_html = report(OUTPUTDIR + "00_fastqc/raw_multiqc.html", caption = REPORT + "multiqc.rst", category="00 quality report"), 
  params:
    multiqc_output_raw = OUTPUTDIR + "00_fastqc/raw_multiqc_data"
  conda:
    CONTAINER + "multiqc.yaml"
  message: 
    "Synthetize quality with multiqc"
  shell: 
    """
    multiqc -n {output.raw_multi_html} {input.raw_qc} #run multiqc
    rm -rf {params.multiqc_output_raw} #clean-up
    """ 
