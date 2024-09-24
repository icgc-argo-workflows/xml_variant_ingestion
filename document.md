# XML Variant Ingestion

The ARGO Data Platform accepts variant calling data in XML format (based on hg19 genome reference). The XML file will be converted to VCF file followed by lift over to GRCh38 reference genome. For details, please see the latest version of the [ARGO xml_variant_ingestion](https://github.com/icgc-argo-workflows/dna-seq-processing-wfs/releases) workflow.

## Inputs
* Submitted XML file(s)
* Mapping file
* GRCh38 as human reference genome
* Genome liftover chain file
* Genome reference used to call variant in the XML file

## Processing
* Submitted variant calling (XML) are converted into VCFs based on variant types (copy number alteration, rearrangement, short variant).
* [Picard:liftovervcf](https://gatk.broadinstitute.org/hc/en-us/articles/27007978536219-LiftoverVcf-Picard) is used to lift the variant calling to GRCh38 reference genome.

## Output
* [Raw SNV Calls](https://docs.icgc-argo.org/docs/data/variant-calls#raw-snv-calls) and [VCF Index](https://docs.icgc-argo.org/docs/data/variant-calls#vcf-index)
* [Raw SV Calls](https://docs.icgc-argo.org/docs/data/variant-calls#raw-sv-calls) and [VCF Index](https://docs.icgc-argo.org/docs/data/variant-calls#vcf-index)
* [Raw CNV Calls](https://docs.icgc-argo.org/docs/data/variant-calls#raw-cnv-calls) and [VCF Index](https://docs.icgc-argo.org/docs/data/variant-calls#vcf-index)

## Workflow Diagram