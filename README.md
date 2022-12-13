# PyTemperatureDLS
This code is designed to provide an analysis tool for temperature dependent dynamic light scattering data collected ALV CGS3.

The file *...scp* is an example code to control the instruments.

Several samples can be analyzed. The software is reading in all data within one specified folder. Each subfolder is considered as one sample. The folder name is treated as sample name.

A refractive index of n=1.332 is assumed to determine q.

Different fit approaches can be easily tested and applied to the data. The temperature dependent diffusion coefficients are displayed for the different  
Up to now, the different scattering vectors are treated independently, no global fits (respecting the q and tau dependence) are performed.

xelatex has to be installed to obtain an automatized generated pdf containing an overview of the different fits for each sample.

The file *DLSanalysis.py* is the main script. *DLSLib.py* contains supporting functions. The files *ilt.py* and *ldp.py* are for the Contin-like analysis and are modified based on https://github.com/caizkun/pyilt.
