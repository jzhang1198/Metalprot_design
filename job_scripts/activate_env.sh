#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd

# The sole purpose of this bash script is to set up the virtual environment to run your jobs in. 

conda activate Metalprot_design

$@ ${SGE_TASK_ID}
