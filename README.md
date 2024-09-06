# MiSelf
This document provides guidance on modeling evolutionary dynamics of selfish mitochondrial genomes, using the MiSelf Python script.

The MiSelf script is prepared to be run from the terminal, from a working directory that contains the sub-directory /Source_data, which will require three input data files (example input data can be downloaded from https://github.com/bgitschlag/MiSelf/Gitschlag_et_al_2024_SOURCE_DATA.zip):

1. INTRA-ORGANISMAL SELECTION DATA. The files are called as /Source_data/{genotype}_intra_organismal_source_data.csv. Each data file should consist of a three-column table with column labels 'Replicate', 'Parent' (for parents), and 'Adult progeny' (for progeny). Each row consists of one biological replicate parent-progeny lineage. Data should be expressed as float values between 0 and 1, representing raw mutant (heteroplasmic) frequency measurements.

2. INTER-ORGANISMAL SELECTION DATA. The files are called as /Source_data/{genotype}_organismal_source_data.csv. Each data file should consist of a table with columns indexing biological replicates (one competed population per replicate, see Gitschlag et al. 2024 Methods) and labeled sequentially as 'Rep 1' to 'Rep n' where n is the number of replicate populations. Rows index generational time points which are numbered sequentially in the column 'Generation'. Inter-organismal selection data are expressed as log-ratios, LN(het/hom), where het and hom are heteorplasmic and homoplasmic (wildtype) fractions of the population, respectively. As part of the original experiment, non-competed control populations were run in parallel, and their population-wide mutant frequency measurements are used as estimates for the mean homoplasmic (wildtype) fraction of the competed populations (see Gitschlag et al. 2024 Methods).

3. MUTANT FREQUENCY MEASUREMENTS. The file is called as /Source_data/mutant_frequency_samples_source_data.csv. This file should consist of a single table with genotype names as column labels. Genotypes and biological replicates are indexed by column and row, respectively.

Customizable settings (such as sample sizes and model input parameters for simulated data) can be found in lines 20-55 of the MiSelf script (https://github.com/user-attachments/assets/b69e4943-3b2b-4d87-bde2-07f53ac189e5)
