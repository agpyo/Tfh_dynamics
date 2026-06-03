# Reciprocal TFH-B Cell Dynamics Shape Selection in the Germinal Center

This repository contains the MATLAB code used for the simulations presented in:

Pyo et al., "Reciprocal TFH-B Cell Dynamics Shape Selection in the Germinal Center"

## System Requirements

### Software
- MATLAB R2022a or later
- Statistics and Machine Learning Toolbox

### Hardware
- No non-standard hardware required

## Installation

1. Download or clone this repository.
2. Open MATLAB.
3. Navigate to the repository directory.

Typical installation time: less than 1 minute.

## Running the Demo

Run the main simulation script:

run main_simulation.m

Model parameters can be modified at the beginning of the script.

Expected output includes:
- B cell population dynamics
- TFH population dynamics
- Affinity distributions
- Clonal lineage abundances
- Diversity statistics

Expected runtime for the demo is typically less than 1 minute on a standard desktop or laptop computer.

## Using the Code

Users can modify model parameters within `main_simulation.m` to run custom simulations, including:
- Number of affinity classes
- Number of TFH clones
- Mutation parameters
- Selection parameters
- Population turnover rates
- Simulation duration

After modifying parameters, rerun:

run main_simulation.m

## Model Description

The model simulates coupled germinal center B cell and T follicular helper (TFH) cell dynamics, including:
- Antigen acquisition
- TFH-mediated selection
- B cell proliferation and apoptosis
- Somatic hypermutation
- Clonal lineage evolution
- Reciprocal TFH-B cell feedback

Simulation outputs include population dynamics, affinity distributions, lineage abundances, and diversity statistics.

## Contact

Andrew G. T. Pyo
agpyo@stanford.edu
