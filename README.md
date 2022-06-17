# Metalprot_design

This repository contains code for designing metal binding sites into helical bundles. There are two main steps in the design algorithm implemented here. In the first, all potential metal binding sites in an input backbone are enumerated. In the second, the metal binding properties of each site are predicted. With these predictions in hand, the user can filter based on desired properties. A more detailed description of each step is provided below. 

## Overview

### Site Enumeration
Site enumeration begins by identifying all helical residues using DSSP. Then, all pairs of $C_\alpha$ s within a user-defined cutoff distance are identified. After constructing imaginary $C_\alpha$ - $C_\beta$ bond vectors, pairs of contacting residues are filtered based on the relative orientations of these bond vectors. This filtering step is implemented to remove pairs of residues that do not satisfy the basic geometric constraints required for the function of metal binding. With the filtered pairs in hand, all possible cliques of size 2, 3, and 4 within the user-defined distance cutoff are enumerated. These cliques are, in essence, the putative metal binding sites. 

### Prediction
After site enumeration, the metal binding properties of each putative binding core are predicted using the neural networks from [Metalprot_learning](https://github.com/lonelu/Metalprot_learning). In the current state of the code, the x, y, and z coordinates of the metal are predicted using the coordinate regressor, along with a confidence metric. In the future, I hope to incorporate predictions from the classifier.

## Installation
This repo is a pip installable module. To run any of the code contained in <code>job_scripts</code>, you will need to pip install this repo to a virtual environment. To do this, you will need to first clone the repo to your local device using the command below.
```
git clone https://github.com/jzhang1198/Metalprot_design.git
```

Create a new conda environment and install the dependencies necessary to run the code.
```
conda env create -f environment.yml
```

If you had problems installing [PyPivoter](https://github.com/rckormos/PyPivoter) to the designated environment, you may need to manually install it. To do this, run
```
pip install --target=<path2sitepackages> pypivoter==0.1.0
```

After installing all dependencies, install the module. Make sure you in the root directory and have your new environment activated.
```
pip install .
```

## Running Jobs
Note that this code can only be run on Windows, linux, and potentially Ubuntu systems due to some persistent platform-specific bugs. The code does not currently work on MacOSX.

## Running Site Enumeration 
All the code for running site enumeration can be found in <code>job_scripts/load.py</code>. You will need to open this code with a text editor and define the following variables: <code>PDB_DIR</code>, <code>DISTANCE_CUTOFF</code>, <code>HELIX_NO_CUTOFF</code>. <code>PDB_DIR</code> should be a path to the directory containing your input pdb files. <code>DISTANCE_CUTOFF</code> is a positive float that defines the maximum distance between $C_\alpha$ atoms of a given binding core. For reference, most metal binding sites found in nature have $C_\alpha$ distances no greater than 12 $\AA$. Therefore, this is a good default value. <code>HELIX_NO_CUTOFF</code> is a positive integer that defines the minimum number of contributing helices in a binding site. For example, a value of 3 will result in enumeration of all possible binding cores containing 3 or more helices. If you are not concerned with the number of contributing helices, set this variable to <code>None<code>.

With the above defined, you are ready to run the code. If you are submitting jobs on Wynton, you will need to first edit the shell script in <code>job_scripts/activate_env.sh</code>. The only thing you need to add is a command for activating the virtual envinment containing all dependencies. This ensures that the cores running your job are using the correct virtual environment. After doing so, run the following command:
```
./run_jobs.py job_name path2output ./job_scripts/load.py -d SGE -n no_jobs -t HH:MM:SS -m memory
``` 

<code>job_name</code> is, quite simply, the user-defined name of the job. <code>path2output</code> is the path where output files should be written to. The -d flag specifies what the job distributor is, which for Wynton, is SGE. The -n flag specifiecs how many cores to parallelize the job across. As an example, if you have a total of 12 structures, you may want to distribute each of these structures to one of 12 cores to speed things up. The number of jobs should be no more than the total number of structures you have. The -t flag specifies the run time. For reference, an 80 amino acid helical bundle takes around 20 seconds to run. If you are using a large cutoff distance and/or a large structure, you may want to allocate more time to prevent your jobs from being killed. Lastly, the -m flag specifies how much memory, in GB, to allocate for each task. You shouldn't need to play around with this unless you have a large structure (if then). 

If you are running the code locally, you will still need to use the wrapper script, but the command will look a little different.
```
./run_jobs.py job_name path2output ./job_scripts/load.py
``` 

The output of these jobs are pickled files containing the flattened backbone distance matrices, residue numbers and chain IDs, and the path to the corresponding pdb file for all enumerated cores. They are named with the following convention: pdb_id_site_df.pkl. Before moving onto prediction, you will need to merge these files using the 

## Running Prediction
Before running prediction on enumerated sites, you will need to define the 