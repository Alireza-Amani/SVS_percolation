# SVS_percolation

## Overview
This repository, SVS_percolation, contains the source code and associated scripts required to reproduce the simulations presented in the manuscript titled "Process-based simulations of percolation from various landfill final covers in a cold climate".

## Content
The repository contains the following notable files and directories:

- `SVS_source_code`: This directory contains the source code for the SVS land surface model. The code must be compiled to obtain the exexcutable file for SVS.
- `data`: This directory contains data required to run the simulation and plotting/evaluation scripts.
- `output`: This directory contains output files or results from the project's simulations.
- `svspaper_package`: This is the package hosting the Python wrap for SVS and related utility functions.
- `evaluations.ipynb`: A Jupyter notebook used for calculating evaluation metrics related to the project.
- `plot_SimVsObs.ipynb`: A Jupyter notebook used for plotting simulation versus observation data.
- `provide_information.py`: A Python script that provides utility for the project. This is where the path to the newly compiled SVS executable and other necessary parameters and information for running the simulations should be provided.
- `run_all_simulations.ipynb`: A Jupyter notebook that runs all the simulations related to the project.

## Getting Started
1. Clone this repository to your local machine.
2. Navigate into the `SVS_source_code` directory and compile the source code.
3. Update the `provide_information.py` script with the path to the new SVS executable and any other required parameters or information.
4. Run the `run_all_simulations.ipynb` Jupyter notebook to start the simulations.
5. Use the `evaluations.ipynb` and `plot_SimVsObs.ipynb` notebooks as needed to calculate evaluation metrics and plot simulation results against observation data.

## Notes
For compiling the SVS source code the reader is advised to refer to the [official guide](https://wiki.usask.ca/pages/viewpage.action?pageId=1885438549). <br><br>
In our case, we compiled the code with **gfortran** version "GNU Fortran (GCC) 9.3.0" on macOS (12.6.5). If you have gfortran installed, compiling it with would be as simple as running the following command in your terminal: `make gfortran OUT=SVS_EXEC_CUSTOM_NAME`. 

## License
This project is licensed under the MIT License.

 
