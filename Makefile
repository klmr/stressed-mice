include bsub.make
include help.make
include dirs.make

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
de-gene-list = raw/KG_GV_protein_lncRNA_DESeq2_name.txt

repeat-quant = $(addprefix data/repeat-quant/,$(addsuffix /quant.sf,$(foreach i,${long-raw-files},${sample_id_$i})))

$(call dirs,data data/reference data/annotation data/RepeatMasker data/qc data/index data/repeat-quant)

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

de-repeat-genes = data/repeat-quant/genes-MSUS-vs-Control.tsv

.PHONY: repeat-de
## Perform differential expression analysis on the repeat elements
repeat-de: ${de-repeat-genes}

${de-repeat-genes}: data/repeat-quant/samples.tsv ${repeat-quant}
	${bsub} $(call memreq,4000) \
		"./scripts/differential-expression --prefix '$(dir $@)' MSUS/Control '$<'"

.PHONY: te-changes
## Plot the differential expression of transposable elements
te-changes: data/repeat-quant/te-changes.pdf

data/repeat-quant/te-changes.pdf: ${de-repeat-genes}
	${bsub} $(call memreq,4000) \
		"./scripts/plot-de-genes \
		--gene-expression $(subst genes-,vsd-,$<) \
		--genes $< \
		--samples data/repeat-quant/samples.tsv \
		--annotation ${repeat-annotation} \
		$@"

te-co-expression = data/repeat-quant/co-expression-ks-ecdf-MSUS-vs-Control.pdf

.PHONY: te-co-expression
## Test for coexpression of upregulated TEs with protein-coding genes and lncRNAs
te-co-expression: ${te-co-expression}

${te-co-expression}: ${de-repeat-genes} ${te-annotation} ${protein-coding-annotation}
	${bsub} $(call memreq,8000) \
		"./scripts/co-expression-test \
		--te-list $< \
		--te-annotation ${repeat-annotation} \
		--de-genes ${de-gene-list} \
		--p-annotation ${protein-coding-annotation} \
		$@"


.DELETE_ON_ERROR:
