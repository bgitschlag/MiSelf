# MiSelf
Python script for modeling evolutionary dynamics of selfish mitochondrial genomes.

This document provides guidance on formatting input files for the MiSelf script (to see example input files, see supplementary data files of Gitschlag et al. 2024).

INTRA-ORGANISMAL SELECTION DATA:

Intra-organismal selection data is called as {genotype}_intra_organismal_source_data.csv, with pairs of adjacent columns per genotype: genotype_gen1 for parents and genotype_gen2 for progeny (one biological replicate parent-progeny lineage per row), with the left-most column indexing biological replicates. Data should be expressed as raw mutant (heteroplasmic) frequency measurements.

ORGANISMAL SELECTION DATA:

Input inter-organismal selection data are expressed as log-ratios, LN(het/hom), where het and hom are heteorplasmic and homoplasmic (wildtype) fractions of the population, respectively. Input file is called as {genotype}_organismal_source_data.csv, with columns indexing biological replicates (a replicate competing population, see Gitschlag et al. 2024 Methods) and each row representing a generational time point.

MUTANT FREQUENCY MEASUREMENTS:

Mutant (heteroplasmic) frequency measurements are called as mutant_frequency_samples_source_data.csv, with genotypes and replicates indexed by column and row, respectively.
