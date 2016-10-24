#
# Helpers
#

bsub = scripts/bsub -K

raw-files = $(shell ls raw/sperm-small-??.cutadapt.fastq.gz)
fastqc-files = $(addprefix data/qc/,$(notdir ${raw-files:.fastq.gz=_fastqc.html}))

#
# Helper variables and definitions
#

directories = data data/qc

#
# QC
#

data/qc/%_fastqc.html: raw/%.fastq.gz | data/qc
	${bsub} -M4000 -R'select[mem>4000] rusage[mem=4000]' \
		"fastqc --outdir data/qc '$<'"
	@rm ${@:%.html=%.zip}

data/qc/multiqc_report.html: ${fastqc-files}
	multiqc --force --outdir data/qc $(sort $(dir $+))

#
# Directories
#

data/qc: data

${directories}:
	mkdir '$@'
