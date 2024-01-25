# MiSelf
Python script for modeling evolutionary dynamics of selfish mitochondrial genomes.

This document provides guidance on formatting input files for the MiSelf script (to see example input files, see supplementary data files of Gitschlag et al. 2024).

STATISTICS INPUT:

The only statistics input is called as t_values_at_95.csv, a two-column t-value table with degrees of freedom (deg_free) and t-values (t_value).

INTRA-ORGANISMAL SELECTION DATA:

Intra-organismal selection data is called as suborganismal_selection_raw_data.csv, with pairs of adjacent columns per genotype: genotype_gen1 for parents and genotype_gen2 for progeny (one biological replicate parent-progeny lineage per row), with the left-most column indexing biological replicates. Data should be expressed as raw mutant (heteroplasmic) frequency measurements.

ORGANISMAL SELECTION DATA:

If log ratio data are used as input organismal selection data, LN(het/hom) where het and hom are heteorplasmic and homoplasmic (wildtype) fractions of the population, respectively, input file is called as genotype_organismal_selection_data.csv, with columns indexing biological replicates (a replicate competing population, see Gitschlag et al. 2024 Methods) and each row representing a generational time point.

If raw mutant frequency data are used as input organismal selection data, input file is called as organismal_selection_raw_data.csv, with columns indexing biological replicate competing populations followed by non-competing control populations, labeled genotype_c_replicate (competing) or genotype_nc_replicate (non-competing).

MUTANT FREQUENCY MEASUREMENTS:

Mutant (heteroplasmic) frequency measurements are called as heteroplasmy_dists.csv, with genotypes and replicates indexed by column and row, respectively.

PARAMETERS FOR SIMULATING FITNESS FUNCTIONS AND FREQUENCY DISTRIBUTIONS:

For simulating the fitness functions and mutant frequency distributions among heteroplasmic hosts, input file is called as parameters_for_fitness_and_distribution_functions.csv. Parameters are indexed by column, with the drift parameter labeled “n (drift)” and selection parameters labeled alpha, beta, gamma, delta, and epsilon. A drift parameter specified as “high_n” instructs the script to use the variable sim_n (variable bottleneck sizes up to 1000 in approximately 2x increments from 62) as drift parameter values for the selection parameters in the corresponding row.
