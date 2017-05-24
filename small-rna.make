include bsub.make
include help.make
include dirs.make

$(call dirs,data data/fastq data/trna-quant data/fastq/cca data/trna-quant/full)

small-rna-raw-files = $(shell ls raw/sperm-small-??.cutadapt.fastq.gz)
trna-fastq = $(addprefix data/fastq/,$(notdir ${small-rna-raw-files:cutadapt.fastq.gz=trna.fastq.gz}))
cca-trim-fastq = $(addprefix data/fastq/cca/,$(notdir ${small-rna-raw-files:cutadapt.fastq.gz=cca-trim.fastq.gz}))
trna-quant = $(addprefix data/trna-quant/,$(notdir ${trna-fastq:.trna.fastq.gz=}))
trna-full-quant = $(addprefix data/trna-quant/full/,$(notdir ${trna-fastq:.trna.fastq.gz=}))

trna-reference = data/reference/mm10-tRNAs.fasta
trna-index = data/index/mm10-tRNAs

${trna-reference}: | data/reference
	curl -o $@ 'http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.fa'

.PRECIOUS: ${trna-index}/header.json
${trna-index}/header.json: ${trna-reference} | data/index
	salmon index --type quasi --kmerLen 15 --perfectHash --transcripts $< --index $(dir $@)

.PHONY: trna-fastq
## Create the smallRNA-seq FASTQ files for tRNA-derived fragments.
trna-fastq: ${trna-fastq}

data/fastq/%.trna.fastq.gz: raw/%.cutadapt.fastq.gz | data/fastq
	${bsub} $(call memreq,500) \
		"gunzip -c $< | scripts/retrieve-trfs | gzip -c > $@"

.PHONY: cca-trim-fastq
## Create the smallRNA-seq FASTQ files with tailing CCA trimmed off.
cca-trim-fastq: ${cca-trim-fastq}
data/fastq/cca/%.cca-trim.fastq.gz: raw/%.cutadapt.fastq.gz | data/fastq/cca
	${bsub} $(call memreq,500) \
		"gunzip -c $< | scripts/trim-cca | gzip -c > $@"

.PHONY: trna-quant
## Quantify tRNA-derived fragments from reads ending in CCA.
trna-quant: ${trna-quant}
data/trna-quant/%: data/fastq/%.trna.fastq.gz ${trna-index}/header.json | data/trna-quant
	${bsub} $(call memreq,1000) \
		"bash -c 'salmon quant --index ${trna-index} --libType U -r <(gunzip -c $<) -o $@'"

.PHONY: trna-de
## Perform DE analysis for tRNA-derived fragments.
trna-de: data/trna-quant/trf-de-results.tsv

data/trna-quant/trf-de-results.tsv: ${trna-quant}
	./scripts/trf-de

.PHONY: trna-full-quant
## Quantify tRNA-derived fragments from all reads.
trna-full-quant: ${trna-full-quant}
data/trna-quant/full/%: data/fastq/cca/%.cca-trim.fastq.gz ${trna-index}/header.json | data/trna-quant/full
	${bsub} $(call memreq,2000) \
		"bash -c 'salmon quant --index ${trna-index} --libType U -r <(gunzip -c $<) -o $@'"

.DELETE_ON_ERROR:

# vim: ft=make
