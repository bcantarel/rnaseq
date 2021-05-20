#!/bin/bash

baseDir="`dirname \"$0\"`"

cd ${baseDir}/mouse_se_test/
sbatch -p 32GB,super run_test.sh
cd ${baseDir}/human_pe_test/
sbatch -p 32GB,super run_test.sh
