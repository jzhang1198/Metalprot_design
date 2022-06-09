# Metalprot_design

Code for running metal site prediction on input backbones. DeGrado lab 2022.

## Installation
This repo is a pip installable module. To run any of the code contained in <code>job_scripts</code>, you will need to pip install this repo to a virtual environment. To do this, you will need to execute the following steps:

1) Clone the repo if you haven't done so already. Use the following command: <code>git clone [link to Metalprot_design]</code>

2) In the root directory of the cloned repo, create a new conda environment entitled <code>Metalprot_design</code> using the following command: <code>conda env create -f environment.yml</code>. This environment contains all necessary dependencies. 
    * Note: I had some issues installing PyPivoter. You can manually install this dependency using pip --target

3) Pip install the repo by running <code>pip install .</code>.

## Running Jobs
The scripts in <code>job_scipts</code> execute enumeration of possible metal binding sites in input structures and prediction of metal binding for each site. 