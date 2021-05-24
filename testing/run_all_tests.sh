#!/bin/bash

baseDir="`dirname \"$0\"`"
mkdir ${baseDir}/mouse_se_test
cd ${baseDir}/mouse_se_test/
cp ${baseDir}/run_mouse_se_test.sh ${baseDir}/mouse_se_test/run_test.sh
ln -s /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-*.img .
sbatch -p 32GB,super run_test.sh

mkdir ${baseDir}/human_pe_test
cd ${baseDir}/human_pe_test/
cp ${baseDir}/run_human_pe_test.sh ${baseDir}/human_pe_test/run_test.sh
ln -s /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-*.img .
sbatch -p 32GB,super run_test.sh
