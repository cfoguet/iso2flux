# Iso2Flux
Version: 0.2

## Short Description

Open source software for steady state 13C flux analysis

## Description

Iso2Flux is a Python-based flux analysis software that allows to perform 13C Metabolic Flux Analysis on a sub-network of a large scale model. Iso2flux uses constraint-based modelling to compute steady state fluxes across a large metabolic network and uses such flxues to predict 13C distribution across a subser of the newtork. Then, given a set of 13C data Iso2flux can iteratively find the steady state flux distributions that are most consistent with such fluxes. 

## Key features

- Constraint-Based Modelling
- 13C Metabolic Flux Analysis

##Functionality

- Post-processing
- Statistical Analysis
- Workflows 

##Approaches

- Isotopic Labelling Analysis / 13C

##Instrument Data Types

- MS

## Screenshots


## Tool Authors

- Carles Foguet (University of Barcelona)
- Pedro de Atauri (University of Barcelona)
- Vitaly Selinvanov (University of Barcelona)
- Marta Cascante (Univesity of Barcelona)

## Container Contributors

- [Pablo Moreno](https://github.com/pcm32) 
- [Carles Foguet](https://github.com/cfoguet) 


## Website

- https://github.com/cfoguet/iso2flux


## Git Repository

- https://github.com/cfoguet/iso2flux

## Usage Instructions
Iso2flux takes the following inputs:
experimental_data_file= Name of the processed Metabolights file containing isotopologues distribution.
tracing_model= Name of the file(s) describing the label propagation model. Should either be a xlsx file with 2 tabs (one for defining metabolites where label can be simulated and one describing reactions that can propagate label) or 2 CSV files containing the same information (the two files names should be in a string separated by a coma)
sbml_model= Name of the SBML file describing the constraint based model that will be used
time= Incubation time that will be evaluated
factor= Experimental factor or condition that will be evaluated
parameters= Name of the CSV file defining additional parameters for Iso2flux (Optional)
constraints= Name of the file containing additional constraints for the constraint model (Optional)
quick_analysis Disables the confidence interval analysis (Optional)
Iso2Flux generates the follwing outputs:
unconstrained_fluxes.csv:Flux variability analysis before tracer data is integrated
best_label.csv:The simulated label distribution that is most consistent with experimental measurements
best_fluxes.csv:The flux distribution that is most consistent with experimental measurments
confidence.csv:The lower and upper limits of the confidence intervals

## Publications
- Carles Foguet, Pedro de Atauri, Vitaly A. Selivanov, Marta Cascante. Iso2Flux: A new software for 13C fluxomics developed in the framework of PhenoMeNal
