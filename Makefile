SHELL := $(shell which bash)

#
# Helpers
#

bsub = scripts/bsub -K
memreq = -M$1 -R'select[mem>$1] rusage[mem=$1]'

#
# Helper variables and definitions
#

rm-threads = 32
rm-mem = 24000

samples = $(shell tail -n+2 ./supporting/samples.tsv | tr $$'\t' '@')

get-sample-field = $(word $2,$(subst @, ,$1))
get-sample-file = $(call get-sample-field,$1,3)
get-sample-id = $(call get-sample-field,$1,1)

$(foreach i,${samples},$(eval sample_id_$(call get-sample-file,$i)=$(call get-sample-id,$i)))

$(foreach i,${samples},$(eval sample_file_$(call get-sample-id,$i)=$(call get-sample-file,$i)))

long-raw-files = $(foreach i,${samples},$(call get-sample-field,$i,3))

genome-reference = data/reference/Mus_musculus.GRCm38.dna.primary_assembly.fa
repeat-reference = data/reference/Mus_musculus.GRCm38.79.repeats.fa
repeat-annotation = data/annotation/Mus_musculus.GRCm38.79.repeats.gtf
rm-repeat-annotation = data/RepeatMasker/$(notdir ${genome-reference}).out
flanking-repeat-reference = data/reference/Mus_musculus.GRCm38.79.repeats-flanking.fa
long-repeat-index = data/index/Mus_musculus.GRCm38.79.repeats-long
protein-coding-annotation = data/annotation/Mus_musculus.GRCm38.79.gtf

repeat-quant = $(addprefix data/repeat-quant/,$(addsuffix /quant.sf,$(foreach i,${long-raw-files},${sample_id_$i})))

#
# Annotation and reference
#

.PHONY: RepeatMasker
## Mask/annotate repeats
RepeatMasker: ${rm-repeat-annotation}

${rm-repeat-annotation}: ${genome-reference} | data/RepeatMasker
	${bsub} -n${rm-threads} -R'span[hosts=1]' $(call memreq,${rm-mem}) \
		"RepeatMasker -pa ${rm-threads} -nolow -species mouse -dir ${@D} $<"

${genome-reference}: | data/reference
	curl -o $@.gz 'ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'
	gunzip '$@.gz'

${repeat-reference}: ${repeat-annotation} ${genome-reference}
	${bsub} $(call memreq,4000) "./scripts/gtf-to-fasta '$<' '$(lastword $+)' '$@'"

${protein-coding-annotation}: | data/annotation
	curl -o $@.gz 'ftp://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/Mus_musculus.GRCm38.79.gtf.gz'
	gunzip '$@.gz'

#
# QC
#

data/qc/%_fastqc.html: raw/%.fastq.gz | data/qc
	${bsub} $(call memreq,4000) "fastqc --outdir data/qc '$<'"
	@rm ${@:%.html=%.zip}

#
# Expression quantification
#

.PHONY: repeat-indices
## Generate Salmon indices for repeats
repeat-indices: ${long-repeat-index}/header.json

.PRECIOUS: ${long-repeat-index}/header.json
${long-repeat-index}/header.json: ${repeat-reference} | data/index
	${bsub} -n2 -R'span[hosts=1]' $(call memreq,24000) \
		"salmon index --type quasi --kmerLen 31 --perfectHash \
		--transcripts '$<' --index '$(dir $@)'"

.PHONY: repeat-quant
## Quantify repeat expression
repeat-quant: ${repeat-quant}

.SECONDEXPANSION:

.PRECIOUS: ${repeat-quant}
data/repeat-quant/%/quant.sf: $${sample_file_$$*} ${long-repeat-index}/header.json | data/repeat-quant
	${bsub} -n8 -R'span[hosts=1]' $(call memreq,64000) \
		"${SHELL} -c 'salmon quant --index $(dir $(lastword $^)) --libType U \
		-r <(gunzip -c $<) -o ${@:%/quant.sf=%}'"

data/repeat-quant/samples.tsv: ./supporting/samples.tsv | data/repeat-quant
	cut -f1-2 $< | paste - <(echo File; tr ' ' '\n' <<< '${repeat-quant}') > $@

#
# Differential expression and downstream analysis
#

de-repeat-genes = data/repeat-quant/genes-ms-vs-co.tsv

.PHONY: repeat-de
## Perform differential expression analysis on the repeat elements
repeat-de: ${de-repeat-genes}

${de-repeat-genes}: data/repeat-quant/samples.tsv ${repeat-quant}
	${bsub} $(call memreq,4000) \
		"./scripts/differential-expression --prefix '$(dir $@)' ms/co '$<'"

#
# Directories
#

directories = data data/reference data/annotation data/RepeatMasker data/qc data/index data/repeat-quant

${directories}: | $$(dir $$@)
	mkdir $@

./:
	@# Nothing to be done.

.DELETE_ON_ERROR:

.DEFAULT_GOAL := show-help
# See <https://gist.github.com/klmr/575726c7e05d8780505a> for explanation.
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)";echo;sed -ne"/^## /{h;s/.*//;:d" -e"H;n;s/^## //;td" -e"s/:.*//;G;s/\\n## /---/;s/\\n/ /g;p;}" ${MAKEFILE_LIST}|LC_ALL='C' sort -f|awk -F --- -v n=$$(tput cols) -v i=19 -v a="$$(tput setaf 6)" -v z="$$(tput sgr0)" '{printf"%s%*s%s ",a,-i,$$1,z;m=split($$2,w," ");l=n-i;for(j=1;j<=m;j++){l-=length(w[j])+1;if(l<= 0){l=n-i-length(w[j])-1;printf"\n%*s ",-i," ";}printf"%s ",w[j];}printf"\n";}'|more $(shell test $(shell uname) == Darwin && echo '-Xr')
