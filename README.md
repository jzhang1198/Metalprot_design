# Metalprot_design

Code for running metal site prediction on input backbones. DeGrado lab 2022.

## Installation
This repo is a pip installable module. To run any of the code contained in <code>job_scripts</code>, you will need to pip install this repo to a virtual environment. To do this, you will need to first clone the repo to your local device using the command below.
```
git clone [link to Metalprot_design]
```

Create a new conda environment and install the dependencies necessary to run the code.
```
conda env create -f environment.yml
```

If you had problems installing [PyPivoter](https://github.com/rckormos/PyPivoter) to the designated environment, you may need to manually install it. To do this, run
```
pip install --target=<path2sitepackages> pypivoter
```

After installing dependencies, install the module. Make sure you in the root directory and have your new environment activated.
```
pip install .
```

## Running Jobs
All jobs can be run using the wrapper script <code>run_jobs.py</code> in the root directory. The general syntax for running jobs from the command line using <code>run_jobs.py</code> is as follows.
```
./run_jobs.py job_name path2output ./job_scripts/job_script.py [options]
``` 

<code>job_name</code> is simply the user-defined name of the job. <code>path2output</code> is the path where output files should be written to. Refer to the docstring for more information on the options.

Keep in mind that there are a few variables you will need to define within each script under <code>job_scripts</code>. These variables are in all caps. 