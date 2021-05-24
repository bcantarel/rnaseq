#!/bin/bash

module load nextflow/20.01.0 singularity/3.5.3
base='/project/BICF/BICF_Core/s166458/rnaseq_astrocyte'
datadir='/project/shared/bicf_workflow_ref/workflow_testdata/rnaseq'

nextflow -C ${base}/nextflow.config run -with-dag flowchart.png -with-timeline mouse_timeline.html -with-report mouse_report.html ${base}/workflow/main.nf --design  ${datadir}/mouse_se.design.txt --input ${datadir} --pairs se --output analysis
