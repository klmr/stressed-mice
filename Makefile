SHELL := $(shell which bash)

#
# Helpers
#

bsub = scripts/bsub -K

raw-files = $(shell ls raw/sperm-small-??.cutadapt.fastq.gz)
fastqc-files = $(addprefix data/qc/,$(notdir ${raw-files:.fastq.gz=_fastqc.html}))

#
# Helper variables and definitions
#

directories = data/reference data/annotation data data/qc data/index data/repeat-quant

rm_engine = crossmatch
rm_threads = 32
rm_memlimit = 64000

long-raw-files = $(shell ls raw/*long*.cutadapt.fastq.gz)

genome-reference = data/reference/Mus_musculus.GRCm38.dna.primary_assembly.fa
repeat-reference = data/reference/Mus_musculus.GRCm38.75.repeats.fa
repeat-annotation = data/annotation/Mus_musculus.GRCm38.75.repeats.gtf
flanking-repeat-reference = data/reference/Mus_musculus.GRCm38.75.repeats-flanking.fa
short-repeat-index = data/index/Mus_musculus.GRCm38.75.repeats-short
long-repeat-index = data/index/Mus_musculus.GRCm38.75.repeats-long

repeat-quant = $(addprefix data/repeat-quant/,$(subst .cutadapt.fastq.gz,/quant.sf,$(notdir ${long-raw-files})))

#
# Annotation and reference
#

${genome-reference}: | data/reference
	curl -o $@.gz 'ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'
	gunzip '$@.gz'

${repeat-reference}: ${repeat-annotation} ${genome-reference}
	./scripts/gtf-to-fasta '$<' '$(lastword $+)' '$@'

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
# Expression quantification
#

.PHONY: repeat-indices
## Generate Salmon indices for repeats
repeat-indices: ${short-repeat-index} ${long-repeat-index}

.PRECIOUS: ${short-repeat-index}
${short-repeat-index}: ${repeat-reference} | data/index
	${bsub} -n2 -M8000 -R'span[hosts=1] select[mem>8000] rusage[mem=8000]' \
		"salmon index --type quasi --kmerLen 25 \
		--transcripts '$<' --index '$@'"

.PRECIOUS: ${long-repeat-index}
${long-repeat-index}: ${repeat-reference} | data/index
	${bsub} -n2 -M8000 -R'span[hosts=1] select[mem>8000] rusage[mem=8000]' \
		"salmon index --type quasi --kmerLen 31 \
		--transcripts '$<' --index '$@'"

.PHONY: repeat-quant
## Quantify repeat expression
repeat-quant: ${repeat-quant}

.PRECIOUS: ${repeat-quant}
data/repeat-quant/%/quant.sf: raw/%.cutadapt.fastq.gz ${long-repeat-index} | data/repeat-quant
	${bsub} -n8 -R'span[hosts=1]' -M12000 -R'select[mem>12000] rusage[mem=12000]' \
		"${SHELL} -c 'salmon quant --index $(lastword $^) --libType U \
		-r <(gunzip -c $<) -o ${@:%/quant.sf=%}'"

data/repeat-quant/samples.tsv: ${repeat-quant}
	./scripts/collect-samples $+ > '$@'

.PHONY: repeat-de
## Perform differential expression analysis on the repeat elements
repeat-de: data/repeat-quant/genes-sperm-vs-zygote.tsv

data/repeat-quant/genes-sperm-vs-zygote.tsv: data/repeat-quant/samples.tsv
	${bsub} -M1000 -R'select[mem>1000] rusage[mem=1000]' \
		"./scripts/differential-expression --prefix '$(dir $@)' sperm/zygote '$<'"

#
# Directories
#

data/reference: data
data/annotation: data
data/qc: data
data/index: data
data/repeat-quant: data

${directories}:
	mkdir '$@'

.DELETE_ON_ERROR:

.DEFAULT_GOAL := show-help
# See <https://gist.github.com/klmr/575726c7e05d8780505a> for explanation.
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)";echo;sed -ne"/^## /{h;s/.*//;:d" -e"H;n;s/^## //;td" -e"s/:.*//;G;s/\\n## /---/;s/\\n/ /g;p;}" ${MAKEFILE_LIST}|LC_ALL='C' sort -f|awk -F --- -v n=$$(tput cols) -v i=19 -v a="$$(tput setaf 6)" -v z="$$(tput sgr0)" '{printf"%s%*s%s ",a,-i,$$1,z;m=split($$2,w," ");l=n-i;for(j=1;j<=m;j++){l-=length(w[j])+1;if(l<= 0){l=n-i-length(w[j])-1;printf"\n%*s ",-i," ";}printf"%s ",w[j];}printf"\n";}'|more $(shell test $(shell uname) == Darwin && echo '-Xr')
