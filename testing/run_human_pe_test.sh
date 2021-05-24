#!/bin/bash

module load nextflow/20.01.0 singularity/3.5.3
base='/project/BICF/BICF_Core/s166458/rnaseq_astrocyte'
datadir='/project/shared/bicf_workflow_ref/workflow_testdata/rnaseq'

nextflow -C ${base}/nextflow.config run ${base}/workflow/main.nf --design  ${datadir}/design.rnaseq.txt --input ${datadir} --output analysis
