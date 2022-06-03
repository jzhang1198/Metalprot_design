#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd

# The sole purpose of this bash script is to set up the virtual environment to run your jobs in. 

conda activate Metalprot_learning

$@ ${SGE_TASK_ID}
