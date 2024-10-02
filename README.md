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

Minimal Test
```bash
nextflow run main.nf \
    -profile docker,test,rdpc_qa \
    --api_token NNNNNNN \
    --hg19_ref_fa /PATH/TO/hg19.fa \
    --hg19_ref_fai /PATH/TO/hg19.fa.fai \
    --hg19_to_hg38_chain /PATH/TO/hg19ToHg38.over.chain.gz \
    --hg38_ref_fa /PATH/TO/hg38.fa.gz \
    --hg38_ref_fai /PATH/TO/hg38.fa.fai \
    --hg38_ref_dict /PATH/TO/hg38.dict
```
Test in RDPC-QA
```bash
nextflow run main.nf \
    -profile docker,rdpc_qa,test_rdpc_qa \
    --api_token NNNNNNN \
    --xml /PATH/TO/XML.FILE \
    --experiment_info_tsv /PATH/TO/MAPPING.TSV \
    --hg19_ref_fa /PATH/TO/hg19.fa \
    --hg19_ref_fai /PATH/TO/hg19.fa.fai \
    --hg19_to_hg38_chain /PATH/TO/hg19ToHg38.over.chain.gz \
    --hg38_ref_fa /PATH/TO/hg38.fa.gz \
    --hg38_ref_fai /PATH/TO/hg38.fa.fai \
    --hg38_ref_dict /PATH/TO/hg38.dict
```

#### Test data (Original)
```bash
nextflow run main.nf \
    -profile docker,rdpc_qa,test_rdpc_qa \
    --api_token NNNNNNN \
    --xml ./test/data/example.xml \
    --experiment_info_tsv ./test/data/sample_new.tsv \
    --hg19_ref_fa /PATH/TO/hg19.fa \
    --hg19_ref_fai /PATH/TO/hg19.fa.fai \
    --hg19_to_hg38_chain /PATH/TO/hg19ToHg38.over.chain.gz \
    --hg38_ref_fa /PATH/TO/hg38.fa.gz \
    --hg38_ref_fai /PATH/TO/hg38.fa.fai \
    --hg38_ref_dict /PATH/TO/hg38.dict
```

#### Test data ( Test1 - ORD-0000001-01.xml - short-variants/copy-number-alterations	)
```bash
nextflow run main.nf \
    -profile docker,rdpc_qa,test_rdpc_qa \
    --api_token NNNNNNN \
    --xml ./test/data/ORD-0000001-01.xml \
    --experiment_info_tsv ./test/data/sample_1.tsv \
    --hg19_ref_fa /PATH/TO/hg19.fa \
    --hg19_ref_fai /PATH/TO/hg19.fa.fai \
    --hg19_to_hg38_chain /PATH/TO/hg19ToHg38.over.chain.gz \
    --hg38_ref_fa /PATH/TO/hg38.fa.gz \
    --hg38_ref_fai /PATH/TO/hg38.fa.fai \
    --hg38_ref_dict /PATH/TO/hg38.dict
```

#### Test data ( Test2 - ORD-0000002-01.xml - short-variants )
```bash
nextflow run main.nf \
    -profile docker,rdpc_qa,test_rdpc_qa \
    --api_token NNNNNNN \
    --xml ./test/data/ORD-0000002-01.xml \
    --experiment_info_tsv ./test/data/sample_2.tsv \
    --hg19_ref_fa /PATH/TO/hg19.fa \
    --hg19_ref_fai /PATH/TO/hg19.fa.fai \
    --hg19_to_hg38_chain /PATH/TO/hg19ToHg38.over.chain.gz \
    --hg38_ref_fa /PATH/TO/hg38.fa.gz \
    --hg38_ref_fai /PATH/TO/hg38.fa.fai \
    --hg38_ref_dict /PATH/TO/hg38.dict
```

#### Test data ( Test3 - ORD-0000003-01.xml - short-variants )
```bash
nextflow run main.nf \
    -profile docker,rdpc_qa,test_rdpc_qa \
    --api_token NNNNNNN \
    --xml ./test/data/ORD-0000003-01.xml \
    --experiment_info_tsv ./test/data/sample_3.tsv \
    --hg19_ref_fa /PATH/TO/hg19.fa \
    --hg19_ref_fai /PATH/TO/hg19.fa.fai \
    --hg19_to_hg38_chain /PATH/TO/hg19ToHg38.over.chain.gz \
    --hg38_ref_fa /PATH/TO/hg38.fa.gz \
    --hg38_ref_fai /PATH/TO/hg38.fa.fai \
    --hg38_ref_dict /PATH/TO/hg38.dict
```

#### Test data ( Test4 - ORD-0000004-01.xml - short-variants )
```bash
nextflow run main.nf \
    -profile docker,rdpc_qa,test_rdpc_qa \
    --api_token NNNNNNN \
    --xml ./test/data/ORD-0000004-01.xml \
    --experiment_info_tsv ./test/data/sample_4.tsv \
    --hg19_ref_fa /PATH/TO/hg19.fa \
    --hg19_ref_fai /PATH/TO/hg19.fa.fai \
    --hg19_to_hg38_chain /PATH/TO/hg19ToHg38.over.chain.gz \
    --hg38_ref_fa /PATH/TO/hg38.fa.gz \
    --hg38_ref_fai /PATH/TO/hg38.fa.fai \
    --hg38_ref_dict /PATH/TO/hg38.dict
```

#### Test data ( Test5 - ORD-0000005-01.xml - short-variants/copy-number-alterations/rearrangements )
```bash
nextflow run main.nf \
    -profile docker,rdpc_qa,test_rdpc_qa \
    --api_token NNNNNNN \
    --xml ./test/data/ORD-0000005-01.xml \
    --experiment_info_tsv ./test/data/sample_5.tsv \
    --hg19_ref_fa /PATH/TO/hg19.fa \
    --hg19_ref_fai /PATH/TO/hg19.fa.fai \
    --hg19_to_hg38_chain /PATH/TO/hg19ToHg38.over.chain.gz \
    --hg38_ref_fa /PATH/TO/hg38.fa.gz \
    --hg38_ref_fai /PATH/TO/hg38.fa.fai \
    --hg38_ref_dict /PATH/TO/hg38.dict
```
