include bsub.make
include help.make
include dirs.make

$(call directories,data data/fastq data/trna-quant)
$(call directories,foo foo/bar)

small-rna-raw-files = $(shell ls raw/sperm-small-??.cutadapt.fastq.gz)
trna-fastq = $(addprefix data/fastq/,$(notdir ${small-rna-raw-files:cutadapt.fastq.gz=trna.fastq.gz}))
trna-quant = $(addprefix data/trna-quant/,$(notdir ${trna-fastq:.trna.fastq.gz=}))

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
		"gunzip -c $< | perl scripts/retrieve-trfs | gzip -c > $@"

.PHONY: trna-quant
## Quantify tRNA-derived fragments.
trna-quant: ${trna-quant}
data/trna-quant/%: data/fastq/%.trna.fastq.gz ${trna-index}/header.json | data/trna-quant
	${bsub} $(call memreq,1000) \
		"bash -c 'salmon quant --index ${trna-index} --libType U -r <(gunzip -c $<) -o $@'"

.DELETE_ON_ERROR:

# vim: ft=make
