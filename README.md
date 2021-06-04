# RNASeq Analysis Worklow

This workflow can be run in with the whole genome, or with a specific list of genes of interest.

## Initiate Nextflow Workflows

### Required Tools

This pipeline uses [Nextflow](https://www.nextflow.io/docs/latest/index.html), a bioinformatics workflow tool and [Singularity](https://sylabs.io/docs/), a containerization tool.

Make sure both tools rae installed before running this pipeline. If running on a HPC cluster then load required modules.

```
module load nextflow/20.01.0 singularity/3.5.3
```

### RNA Design File

The design file must named design.txt and be in tab seperated format for the workflows. All RNA workflows can be run usin the same design file format.  You can run in single-end mode with blank cells in the FqR2 column.

| SampleID | SampleGroup | SubjectID | SampleName | FqR1 | FqR2 |
|---|---|---|---|---|---|
| SRR1551069 |	monocytes |	53 |	53_Monocytes |	SRR1551069_1.fastq.gz |	SRR1551069_2.fastq.gz |
| SRR1551068 | neutrophils |	53 |	53_Neutrophils |	SRR1551068_1.fastq.gz	|	SRR1551068_2.fastq.gz |
| SRR1551055 | monocytes |	21 |	21_Monocytes |	SRR1551055_1.fastq.gz |	SRR1551055_2.fastq.gz |
| SRR1551054 | neutrophils |	21 |	21_Neutrophils |	SRR1551054_1.fastq.gz	|	SRR1551054_2.fastq.gz |
| SRR1551048 |	monocytes |	20 |	20_Monocytes |	SRR1551048_1.fastq.gz	|	SRR1551048_2.fastq.gz |
| SRR1551047 |	neutrophils |	20 |	20_Neutrophils |	SRR1551047_1.fastq.gz	|	SRR1551047_2.fastq.gz |
| SRR1550987 |	monocytes |	44 |	44_Monocytes |	SRR1550987_1.fastq.gz	|	SRR1550987_2.fastq.gz |
| SRR1550986 |	neutrophils |	44 |	44_Neutrophils |	SRR1550986_1.fastq.gz	|	SRR1550986_2.fastq.gz |


### RNA Parameters
* **--input**
  * directory containing the design file and fastq files
  * default is set to *'${basedir}/fastq'*
  * eg: **--input '/project/shared/bicf_workflow_ref/workflow_testdata/rnaseq/fastq'**
* **--output**
  * directory for the analysis output
  * default is set to *'${basedir}/analysis'*
  * eg: **--output '${basedir}/output'**
* **--genome**
  * directory containing all reference files for the various tools. This includes the genome.fa, gencode.gtf, genenames.txt, ect.
  * default is set for use on UTSW BioHPC.
  * eg: **--genome '/project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref'**
* **--stranded**
  * option for -s flag in featurecount used in geneabundance calculations
  * default is set to *'0'*
  * eg: **--stranded '0'**
* **--pairs**
  * select either 'pe' (paired-end) or 'se' (single-end) based on read inputs. Select 'pe' when both R1 and R2 are present. If only R1, then select 'se'.
  * default is set to *'pe'*
  * eg: **--pairs 'pe'**
* **--align**
  * select the algorithm/tool for alignment from 'hisat' or 'star'
  * default is set to *'hisat'*
  * eg: **--align 'hisat'**
* **--markdups**
  * select either picard (Mark Duplicates) or null (do not Mark Duplicates)
  * default is set to *'picard'*
  * eg: **--align 'picard'**

### RNA Run Workflow Testing

Human PE

```
module load nextflow/20.01.0 singularity/3.5.3
base=$repoClonedDirectory
datadir='/project/shared/bicf_workflow_ref/workflow_testdata/rnaseq'

nextflow -C ${base}/nextflow.config run ${base}/workflow/main.nf --design  ${datadir}/design.rnaseq.txt --input ${datadir} --output analysis

```

Mouse SE

```
module load nextflow/20.01.0 singularity/3.5.3
base=$repoClonedDirectory
datadir='/project/shared/bicf_workflow_ref/workflow_testdata/rnaseq'

nextflow -C ${base}/nextflow.config run -with-dag flowchart.png -with-timeline mouse_timeline.html -with-report mouse_report.html ${base}/workflow/main.nf --design  ${datadir}/mouse_se.design.txt --input ${datadir} --pairs se --output analysis

```

