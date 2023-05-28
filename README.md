# SVS_percolation Repository

Welcome to the SVS_percolation repository. This repository hosts the source code and associated scripts necessary to reproduce the simulations presented in the manuscript, "Process-based simulations of percolation from various landfill final covers in a cold climate".

## Contents of the Repository

- `SVS_source_code`: This directory contains the source code for the SVS land surface model. You need to compile this code to obtain the executable file for SVS.
- `data`: This directory stores the data required to run the simulations and the scripts for plotting and evaluation.
- `output`: This is the directory where the output files or results from the simulations will be stored.
- `svspaper_package`: This package hosts the Python wrapper for SVS along with some utility functions.
- `evaluations.ipynb`: A Jupyter notebook for calculating evaluation metrics associated with the project.
- `plot_SimVsObs.ipynb`: A Jupyter notebook for plotting simulation data versus observation data.
- `provide_information.py`: A Python script that is used to provide necessary information for running the simulations, including the path to the SVS executable and other required parameters.

## Getting Started

1. Clone this repository to your local machine:<br> ```git clone https://github.com/Alireza-Amani/SVS_percolation.git```
3. Navigate to the `SVS_source_code` directory and compile the source code.
4. Update the `provide_information.py` script with the path to the new SVS executable and any other required parameters or information.
5. Run the `run_all_simulations.ipynb` Jupyter notebook to start the simulations.
6. Utilize the `evaluations.ipynb` and `plot_SimVsObs.ipynb` notebooks as needed to calculate evaluation metrics and plot simulation results against observation data.

## Important Note

For compiling the SVS source code, please refer to the [official guide](https://wiki.usask.ca/pages/viewpage.action?pageId=1885438549). In this project, the code was compiled using gfortran version "GNU Fortran (GCC) 9.3.0" on macOS (12.6.5). If you have gfortran installed, you can compile the code by running the following command in your terminal: `make gfortran OUT=SVS_EXEC_CUSTOM_NAME`.

## License

This project is licensed under the MIT License, permitting use, modification, and distribution with certain conditions.

We appreciate your interest in this project and hope you find it useful for your research. For any further questions or assistance, please feel free to reach out.
