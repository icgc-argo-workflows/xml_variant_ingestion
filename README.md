### Input data:

Mapping file: `test/data/sample_new.tsv`

XML file: `test/data/example.xml`

### Reference files:

`test/reference/`

`hg19.fa` *Not provided, download needed* https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

`hg19.fa.fai`

`hg38.fa.gz` *Not provided, download needed* https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

`hg38.fa.fai`

`hg38.dict`

`hg19ToHg38.over.chain.gz` https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

### Usage

1. Download `hg19.fa` and `hg38.fa.gz` to `test/reference/`
2. Retrive `api_token`
3. Workflow command
```bash
nextflow run main.nf -profile docker,test,rdpc_qa --api_token $api_token
```
