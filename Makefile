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

rm_engine = crossmatch
rm_threads = 32
rm_memlimit = 64000

#
# QC
#

data/qc/%_fastqc.html: raw/%.fastq.gz | data/qc
	${bsub} -M4000 -R'select[mem>4000] rusage[mem=4000]' \
		"fastqc --outdir data/qc '$<'"
	@rm ${@:%.html=%.zip}

## Generate the quality control report
data/qc/multiqc_report.html: ${fastqc-files}
	multiqc --force --outdir data/qc $(sort $(dir $+))

#
# Directories
#

data/qc: data

${directories}:
	mkdir '$@'

.DELETE_ON_ERROR:

.DEFAULT_GOAL := show-help
# See <https://gist.github.com/klmr/575726c7e05d8780505a> for explanation.
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)";echo;sed -ne"/^## /{h;s/.*//;:d" -e"H;n;s/^## //;td" -e"s/:.*//;G;s/\\n## /---/;s/\\n/ /g;p;}" ${MAKEFILE_LIST}|LC_ALL='C' sort -f|awk -F --- -v n=$$(tput cols) -v i=19 -v a="$$(tput setaf 6)" -v z="$$(tput sgr0)" '{printf"%s%*s%s ",a,-i,$$1,z;m=split($$2,w," ");l=n-i;for(j=1;j<=m;j++){l-=length(w[j])+1;if(l<= 0){l=n-i-length(w[j])-1;printf"\n%*s ",-i," ";}printf"%s ",w[j];}printf"\n";}'|more $(shell test $(shell uname) == Darwin && echo '-Xr')
