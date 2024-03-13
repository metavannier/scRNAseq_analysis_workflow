# ----------------------------------------------
# FastQC : Check quality of reads files
# ----------------------------------------------

rule fastqc:
	input:
		expand(RAWDATA + "{rawsample}_{library}_{types}.fastq.gz", rawsample=RAWSAMPLE, library=LIBRARY, types=TYPES)

	output:
		fastqc_output = expand(OUTPUTDIR + "00_clean/fastqc_output.txt")

	message:
		"Quality with fastqc"

	conda:
		CONTAINER + "fastqc.yaml"

	shell:
		"""
		export PATH="/$PWD/.local/bin:$PATH"
		fastqc_output=({output.fastqc_output})
		fastqc {input} --outdir {OUTPUTDIR}00_clean/
		touch ${{fastqc_output}}
		"""

# ----------------------------------------------
# Summarize analysis results in a single report
# ----------------------------------------------

rule multiqc:
	input:
		fastqc_output = expand(OUTPUTDIR + "00_clean/fastqc_output.txt")

	output:
		multiqc_output = expand(OUTPUTDIR + "00_clean/multiqc_output.txt"),
		raw_multiqc_html = report(expand(OUTPUTDIR + "00_clean/raw_multiqc.html", sample_id = SAMPLE_ID), caption = REPORT + "multiqc.rst", category = "00 quality"),

	params:
		raw_qc = expand(OUTPUTDIR + "00_clean/{rawsample}_{library}_{types}_fastqc.zip", rawsample=RAWSAMPLE, library=LIBRARY, types=TYPES),

	message:
		"Synthetize quality with multiqc"

	conda:
		CONTAINER + "multiqc.yaml"

	shell:
		"""
		multiqc -n {params.raw_multiqc_html} {params.raw_qc}
		touch {output.multiqc_output}
		"""

# ----------------------------------------------
# Trimmomatic: trimming reads and removing adapter sequences
# ----------------------------------------------

# rule trimmomatic:
# 	input:
# 		multiqc_output = expand(OUTPUTDIR + "00_clean/multiqc_output.txt"),

# 	output:
# 		trimmomatic_output = expand(OUTPUTDIR + "00_clean/trimmomatic_output.txt"),

# 	params:
# 		raw_qc = expand(OUTPUTDIR + "00_clean/{rawsample}_{library}_{types}_fastqc.zip", rawsample=RAWSAMPLE, library=LIBRARY, types=TYPES),
# 		raw_trimmed=expand( "05_Output/02_trimmomatic/{samples}_{run}.trimmed.fastq", samples=SAMPLES, run=RUN),
# 		raw_untrimmed=expand( "05_Output/02_trimmomatic/{samples}_{run}un.trimmed.fastq", samples=SAMPLES, run=RUN)
# 		adaptaters=config["trimmomatic"]["adaptaters"],
# 		leading=config["trimmomatic"]["leading"],
# 		trailing=config["trimmomatic"]["trailing"],
# 		minlen=config["trimmomatic"]["minlen"]

# 	conda: 
# 		CONTAINER + "trimmomatic.yaml"

# 	shell:
# 		"""
# 		raw_qc=({params.raw_qc})
# 		raw_trimmed=({params.raw_trimmed})
# 		raw_untrimmed=({params.raw_untrimmed})
# 		len=${{#raw_qc[@]}}
# 		for (( i=0; i<$len; i=i+2 ))
# 		do trimmomatic PE -threads 4 ${{raw_qc[$i]}} ${{raw_qc[$i+1]}} ${{raw_trimmed[$i]}} ${{raw_untrimmed[$i]}} ${{raw_trimmed[$i+1]}} ${{raw_untrimmed[$i+1]}} LEADING:20 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36
# 		done
# 		echo "Trimmomatic step is FINISH" > ${{multiqc_output}}
# 		"""
