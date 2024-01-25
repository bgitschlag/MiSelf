# This script will perform population-genetic modeling of a Wright-Fisher process to calculate the most
# evolutionarily stable population-wide frequency distribution of a heteroplasmic mutation, given
# intra-organismal & organismal fitness functions (defined by model parameters), a neutral drift parameter,
# and a sample of mutant frequencies among heteroplasmic individuals. Model parameters are estimated by
# maximum likelihood using experimental data and their uncertainty is estimated by parametric bootstrapping.

import os
import numpy as np
import pandas as pd
import math
from scipy.linalg import eig
from scipy.optimize import minimize
from scipy.special import betainc
from scipy.stats import norm
from scipy.special import expit


# Calculate organismal fitness & average mutant frequency per individual, per generation:
def organismal_fitness_calculator(genotype):

    # Generate appendable lists and replicate-counters:
    c_rep = 1
    nc_rep = 1
    comp_z = []
    noncomp = []
    w_org = []
    logs = []
    measurements = []

    # Extract mutant frequency data (non-competing lines) from source file and reformat as numpy array:
    while genotype + f'_nc_{nc_rep}' in organism_data:
        noncomp.append(organism_data[genotype + f'_nc_{nc_rep}'].dropna())
        nc_rep += 1
    noncomp_z = np.array(sum(noncomp).T / len(noncomp))

    # Reformat mutant frequency data (competing lines) as numpy array, calculate fitness & log ratio of genotypes:
    while genotype + f'_c_{c_rep}' in organism_data:
        noncomp_ave = sum(np.array(noncomp)) / len(noncomp)
        z_raw = organism_data[genotype + f'_c_{c_rep}'].dropna()
        comp_z.append(z_raw)
        mut_array = np.array(z_raw)
        while mut_array[-1] == 0:
            mut_array = mut_array[:-1]
            noncomp_ave = noncomp_ave[:-1]
        measurements.append(len(noncomp_ave))
        timepoints = np.arange(0, len(noncomp_ave), 1)
        log_r = np.log((mut_array/noncomp_ave) / (1 - (mut_array/noncomp_ave)))
        ave_log = sum(log_r.T) / len(noncomp_ave)
        mid_t = sum(timepoints) / len(noncomp_ave)
        cov_tl = sum((log_r * timepoints).T)/len(noncomp_ave) - (ave_log*mid_t)
        var_t = sum(timepoints**2) / len(noncomp_ave) - (mid_t**2)
        slope_l = cov_tl / var_t

        # Fitness calculated as the exponential function of the slope of the log-ratio over time:
        w_org.append(np.exp(slope_l))
        inl_log = ave_log - slope_l * mid_t
        while len(log_r) < len(noncomp_z):
            log_r = np.append(log_r, slope_l * len(log_r) + inl_log)
        logs.append(log_r)
        c_rep += 1

    return w_org, np.array(comp_z), noncomp_z, np.array(logs, dtype=object).T, np.array(measurements)


# Find mutant frequency in generation n+1 given frequency of gen n and intra-organismal selection:
def intra_selection(g, d, e, z):
    fitness = abs(d) * (1-expit(z*g + e))
    covar = (fitness-1) * (z * (1-z))
    exp_fit = z*fitness + (1-z)

    # If z is an array, assign a nonzero exp_w at z=1, to avoid divide-by-zero error in delta_z:
    if isinstance(z, np.ndarray):
        exp_fit[z == 1] = 1
    delta_z = covar / exp_fit

    return np.array(z + delta_z, dtype=float)


# Simulate neutral drift by generating an n-sized binomial distribution around input frequency z:
def drift(n, z):
    drift_matrix = []
    for k in range(0, n+1):
        binomial = np.array(math.comb(n, k) * ((1 - z) ** (n - k)) * (z ** k))
        binomial = np.nan_to_num(binomial)
        binomial = np.array(binomial, dtype=float)
        drift_matrix.append(binomial)
    return np.array(drift_matrix)


# Simulate organismal selection using the cumulative density function of the beta distribution:
def org_selection(z, a, b):
    fit_cost = betainc(a, b, z)
    fitness = np.array(1 - fit_cost, dtype=float)
    fitness = np.nan_to_num(fitness)
    return fitness


# Find relative organismal fitness (cumulative density function of beta distribution) for a given frequency:
def frequency_of_fitness(w, a, b, itl_z):
    def neg_log_dif(x):
        fit = 1 - betainc(a, b, x)
        return abs(w - fit)
    freq_of_fit = minimize(neg_log_dif, itl_z, method="Nelder-Mead", bounds=np.array([[0.0, 1.0]]))
    return float(freq_of_fit['x'])


# Using selection & drift parameters, find mean organismal fitness & most stable mutant frequency distribution:
def population_genetics_model(z, n, a, b, g, d, e):

    # Generate principal submatrix of 2D probability distribution generated by multilevel selection & drift:
    matrix = (drift(n, intra_selection(g, d, e, z))[1:-1].T[:-1] * org_selection(z, a, b)[:-1]).T
    matrix[np.where(matrix == -0.0)] += 0.0
    eigen = eig(matrix)
    val = list(np.real(eigen[0]))
    vec = np.real(eigen[1])
    fit = max(val)
    prb = vec[:, val.index(fit)] / sum(vec[:, val.index(fit)])

    # If prb is not a true eigenvector, or contains a negative value, simulate evolution to find eigenvector:
    if np.any(prb < 0.0) or max(abs(np.dot(matrix, prb) - (prb * fit))) > 1e-15:
        gen = 1
        prb = abs(prb)
        while max(abs(np.dot(matrix, prb) - (prb * fit))) > 1e-15:
            evolvec = np.dot(matrix, prb)
            fit = sum(evolvec)
            prb = evolvec / fit
            gen += 1
    prb = np.concatenate((prb, np.zeros(1)))

    # Generate a smooth probability density distribution by applying gaussian smoothing to the eigenvector:
    gaussian = np.zeros(1000)
    for bin in np.arange(0, n, 1):
        gaussian += (norm.pdf(np.linspace(0, 1.0, 1000), loc=z[bin], scale=1 / n) * prb[bin])
    gaussian[0] = 0.0

    return fit, gaussian, prb


# Draw a sample of mutant frequencies from a given input probability-density distribution:
def sample_density_dist(dens, steps, draws):
    cumulative = np.cumsum(dens)
    cumulative *= 1 / cumulative[-1]
    sample = np.interp(np.random.uniform(size=draws), cumulative, steps)
    return sample


# Calculate log-likelihood of actual data given the expected data predicted by specific parameter values:
def parameters_likelihood(expected, actual, sample_size):
    ssq = sum((expected - actual)**2)
    var = ssq / sample_size
    likelihood = (sample_size/2) * (-np.log(2*np.pi) - np.log(var) - 1)
    return likelihood


# Calculate the log-likelihood of the frequency distribution corresponding to a set of parameters:
def distribution_likelihood(sample, densities):
    likelihoods = np.interp(np.array(sample), np.linspace(0, 1.0, 1000), densities)
    likelihoods[likelihoods < 0] = 0
    likelihood = sum(list(np.log(likelihoods)))
    return likelihood


# Calculate total log-likelihood of intra-organismal & organismal selection, and mutant frequency distribution:
def likelihood(z, n, a, b, g, d, e, i, z1, z2, log_r, duration, timepts, z_obs):

    # Log-likelihood of progeny mutant frequency given their expected frequency:
    exp_z2 = intra_selection(g, d, e, z1)
    l_intra = parameters_likelihood(exp_z2, z2, len(z2))

    # Log-likelihood of organismal selection data given expected vs actual log ratio of genotypes:
    popgen = population_genetics_model(z, n, a, b, g, d, e)
    r = sum(popgen[2] * drift(n, intra_selection(g, d, e, z))[0])
    het_frxn = i
    l_org = 0
    exp_lr = []

    # For all time-points in the competition experiment, calculate the expected log ratio of genotypes:
    for generation in range(0, duration):
        hom_frxn = 1 - het_frxn
        exp_lr.append(np.log(het_frxn / hom_frxn))
        progeny_het = popgen[0] * het_frxn
        progeny_hom = r * het_frxn
        progeny_tot = progeny_het + progeny_hom
        het_frxn = progeny_het / (progeny_tot + hom_frxn)
    for lr in log_r.T:
        dur_lr = timepts[np.where(log_r.T == lr)[0][0]]
        l_org += parameters_likelihood(np.array(exp_lr), lr, dur_lr)

    # Log-likelihood of sampled mutant frequencies given the distribution predicted by model parameters:
    l_dist = distribution_likelihood(z_obs, popgen[1])

    # Total log-likelihood of all data given the predictions of the model parameters:
    likelihood = l_intra + l_org + l_dist

    return likelihood


# Find maximum-likelihood values of all continuous parameters given input data and assumed drift parameter (n):
def maximize_likelihood(z, n, z1, z2, log_r, duration, timepts, z_obs, itl):
    def nl(x):
        lkhd = likelihood(z, n, x[0], x[1], x[2], x[3], x[4], x[5], z1, z2, log_r, duration, timepts, z_obs)
        return -lkhd
    m = minimize(nl, itl, method="Nelder-Mead",
                 bounds=np.array([[1e-16, np.inf], [1e-16, np.inf], [-np.inf, np.inf], [-np.inf, np.inf],
                                  [-np.inf, np.inf], [0.01, 0.99]]), options={"adaptive": True})
    return m['x'][0], m['x'][1], m['x'][2], m['x'][3], m['x'][4], m['x'][5], -m['fun']


# Find max-likelihood values of continuous parameters given input data across for a grid of drift parameter values:
def max_likelihood_values(z1, z2, log_r, duration, timepts, z_obs, min_n, max_n, itl):
    ml_par = pd.DataFrame(columns=['N', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'intl het/hom', 'log(L)'])

    # Generate table (Pandas dataframe) of parameter values for input data and a given drift parameter value:
    def mlp_for_n(n):
        z = np.arange(1/n, 1 + 1/n, 1/n)
        while len(z) > n:
            z = z[:-1]
        ml = maximize_likelihood(z, n, z1, z2, log_r, duration, timepts, z_obs, itl)
        mlp = pd.DataFrame(zip([n], [ml[0]], [ml[1]], [ml[2]], [ml[3]], [ml[4]], [ml[5]], [ml[6]]),
                           columns=['N', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'intl het/hom', 'log(L)'])
        print(mlp)
        return mlp

    # Grid search across n (drift parameter) values, generating max-likelihood parameter estimates for each n:
    for n in np.arange(min_n, max_n+1, 1):
        ml_par = pd.concat([ml_par, mlp_for_n(n)])

    # Expand grid search to larger n (increments of 125*sqrt(2)^integer) until likelihood peaks or n=1000:
    if int(ml_par[ml_par['log(L)'] == max(ml_par['log(L)'])]['N']) > 95:
        x = 1
        xtnd_n = np.concatenate((np.array([125]), np.int_(np.around(np.array([125 * np.sqrt(2) ** x])))))
        for n in xtnd_n:
            ml_par = pd.concat([ml_par, mlp_for_n(n)])
        while int(ml_par[ml_par['log(L)'] == max(ml_par['log(L)'])]['N']) >= xtnd_n[-2]:
            x += 1
            xtnd_n = np.concatenate((xtnd_n, np.int_(np.around(np.array([125 * np.sqrt(2) ** x])))))
            ml_par = pd.concat([ml_par, mlp_for_n(xtnd_n[-1])])
            if x == 6:
                break
    return ml_par


# Generate simulated (bootstrap) data, given input parameters and input data:
def bootstrap(z, z_mu, n, a, b, g, d, e, i, ipool, iboot, dboot, log_r, replicates, duration):
    popgen = population_genetics_model(z, n, a, b, g, d, e)

    # From intra-organismal selection & drift parameters, simulate parent & progeny mutant frequencies:
    z1 = sample_density_dist(popgen[1], z_mu, iboot)
    z2 = []
    for z1_i in z1:
        intra = drift(n, intra_selection(g, d, e, z1_i))
        density = n * (intra[1:] / sum(intra[1:]))
        draws = sample_density_dist(np.interp(z_mu, z, density), z_mu, ipool)
        z2.append(sum(draws)/ipool)

    # From model parameters and empirical variance, simulate organismal selection data (log ratio of genotypes):
    r = sum(popgen[2] * drift(n, intra_selection(g, d, e, z))[0])
    het_frxn = i
    exp_lr = []

    # Calculate expected log ratios of genotypes across a given organismal selection experiment:
    for generation in range(0, duration):
        hom_frxn = 1 - het_frxn
        exp_lr.append(np.log(het_frxn / hom_frxn))
        progeny_het = popgen[0] * het_frxn
        progeny_hom = r * het_frxn
        progeny_tot = progeny_het + progeny_hom
        het_frxn = progeny_het / (progeny_tot + hom_frxn)

    # Simulate observed (bootstrap) log ratios of genotypes across the same organismal selection experiment:
    sd_org = 0
    for lr in log_r.T:
        sd_org += (sum((np.array(exp_lr)-lr)**2) / duration)**0.5
    sd_org = sd_org/len(log_r.T)
    org_sim = []
    for i in range(0, replicates):
        log = []
        for timepoint in range(0, duration):
            log.append(float(np.random.normal(np.array(exp_lr)[timepoint], scale=sd_org, size=1)))
        org_sim.append(np.array(log))

    # Simulate sample draws from the predicted mutant frequency distribution (final returned item):
    return z1, np.array(z2), np.array(org_sim).T, sample_density_dist(popgen[1], z_mu, dboot)


# Generate tables of simulated (bootstrap) data & of their corresponding max-likelihood parameter values:
def boots_ml(path, z, z_mu, n, a, b, g, d, e, i, ipool, iboot, dboot, log_r, reps, duration, tmpts, min_n, max_n):

    # Generate and define output files and data structures:
    folder = path + f'/{genotype}_simulated_data_for_bootstrapping_{iboot}intra_{reps}org/'
    if not os.path.exists(folder):
        os.makedirs(folder)
    distrib = folder + f'/{genotype}_dist_bootstraps.csv'
    parent = folder + f'/{genotype}_parent_bootstraps.csv'
    progeny = folder + f'/{genotype}_progeny_bootstraps.csv'
    if os.path.exists(distrib):
        dist_sims = list(np.array(pd.read_csv(distrib)).T)
        z1_sims = list(np.array(pd.read_csv(parent)).T)
        z2_sims = list(np.array(pd.read_csv(progeny)).T)
    else:
        z1_sims = []
        z2_sims = []
        dist_sims = []
    org_folder = folder + f'/{genotype}_simulated_organismal_data/'
    if not os.path.exists(org_folder):
        os.makedirs(org_folder)
    ml_par = pd.DataFrame(columns=['N', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'intl het/hom', 'log(L)'])

    # For each in a given number of replicates, simulate a data set & find its max-likelihood parameter values:
    for s in range(0, simulations):
        print(f'data simulated under N={n} alpha={a} beta={b} gamma={g} delta={d} epsilon={e}')

        # Define the mutant frequency values used for the intra-organismal parents & distribution sample:
        if dboot == iboot:
            sim = bootstrap(z, z_mu, n, a, b, g, d, e, i, ipool, iboot, dboot, log_r, reps, duration)
            z_sample = sim[0]
        else:
            sim = bootstrap(z, z_mu, n, a, b, g, d, e, i, ipool, iboot, dboot-iboot, log_r, reps, duration)
            z_sample = np.concatenate((sim[0], sim[3]))
        z1_sims.append(sim[0])
        z2_sims.append(sim[1])
        dist_sims.append(z_sample)

        # Export simulated organismal selection data (log ratio of genotypes) to csv output file:
        org_rep = 1
        while os.path.exists(org_folder + f'/{genotype}_org_rep_{org_rep}.csv'):
            org_rep += 1
        pd.DataFrame(sim[2]).to_csv(org_folder + f'/{genotype}_org_rep_{org_rep}.csv')

        # Build table of maximum-likelihood parameter values given simulated (bootstrap) data:
        ml = max_likelihood_values(sim[0], sim[1], sim[2], duration, tmpts, z_sample,
                                   min_n, max_n, np.array([a, b, g, d, e, i]))
        ml_i = ml[ml['log(L)'] == max(ml['log(L)'])]
        print(f"ml_i = {ml_i}")
        ml_par = pd.concat([ml_par, ml_i])

    # Save simulated (bootstrap) data to csv output files:
    pd.DataFrame(np.array(z1_sims).T).to_csv(folder + f'/{genotype}_parent_bootstraps.csv', index=False)
    pd.DataFrame(np.array(z2_sims).T).to_csv(folder + f'/{genotype}_progeny_bootstraps.csv', index=False)
    pd.DataFrame(np.array(dist_sims).T).to_csv(folder + f'/{genotype}_dist_bootstraps.csv', index=False)

    return ml_par


# From a given set of model parameters, derive the fitness functions and mutant frequency descriptive results:
def fxns_of_parameters(genotype, parameters, output, z, drift_columns):

    # Define lists and tables for appending output results:
    z_mu_sq = z ** 2
    zs = np.arange(0, 1.005, 0.005)
    intra = [zs, np.ones(201)]
    org = [zs, np.ones(201)]
    pdf = [np.concatenate((np.array([0]), zs))]
    mean_z_table = pd.DataFrame(columns=drift_columns)
    stdv_z_table = pd.DataFrame(columns=drift_columns)

    # Read each set of parameters from a table containing a list of different combinations of parameter values:
    for ind in parameters.index:
        a = float(parameters['alpha'][ind])
        b = float(parameters['beta'][ind])
        g = float(parameters['gamma'][ind])
        d = float(parameters['delta'][ind])
        estr = parameters['epsilon'][ind]
        if isinstance(estr, str):
            e = eval(estr)
        else:
            e = float(estr)
        intra = np.vstack([intra, np.interp(zs, z, abs(d) * (1-expit(z*g+e)))])
        org = np.vstack([org, np.interp(zs, z, org_selection(z, a, b))])

        # Define drift parameters used for calculating mutant frequency data (distributions, mean, std dev):
        if genotype == 'simulated_genotype':
            drift = list(map(int, drift_columns[1:]))
            mean_z = ['N/A']
            stdv_z = ['N/A']
        else:
            _n_ = int(parameters['N'][ind])
            log2_scale_n = _n_ * (2 ** np.arange(4))
            np.put(log2_scale_n, np.where(log2_scale_n > 1000), 1000)
            drift = np.array([int(x) for x in list(log2_scale_n)])
            mean_z = [_n_]
            stdv_z = [_n_]

        # For each value of n (drift parameter), generate distributions, mean, & std dev of mutant frequencies:
        for n in drift:
            z_tot = np.arange(1/n, 1 + 1/n, 1/n)
            while len(z_tot) > n:
                z_tot = z_tot[:-1]
            popgen = population_genetics_model(z_tot, n, a, b, g, d, e)
            pdf = np.vstack([pdf, np.concatenate((np.array([n]), np.interp(zs, z, popgen[1])))])
            prob_dist = popgen[1] / sum(popgen[1])
            exp_z = sum(z_mu * prob_dist)
            exp_z_sq = sum(z_mu_sq * prob_dist)
            var_z = exp_z_sq - (exp_z ** 2)
            mean_z.append(exp_z)
            stdv_z.append(np.sqrt(var_z))

        # Append mean & std dev of mutation frequencies to their respective tables:
        mean_z_table = pd.concat([mean_z_table, pd.DataFrame([mean_z], columns=drift_columns)])
        stdv_z_table = pd.concat([stdv_z_table, pd.DataFrame([stdv_z], columns=drift_columns)])

    # Save the fitness functions and mutant frequency results (distributions, mean, std dev) to output files:
    pd.DataFrame(np.array(org).T).to_csv(output + f'/{genotype}_organismal.csv')
    pd.DataFrame(np.array(intra).T).to_csv(output + f'/{genotype}_intra_organismal.csv')
    pd.DataFrame(np.array(pdf).T).to_csv(output + f'/{genotype}_distribution.csv')
    mean_z_table.to_csv(output + f'/{genotype}_mean_z.csv')
    stdv_z_table.to_csv(output + f'/{genotype}_stdv_z.csv')

    return "finished compiling fitness functions and frequency data for input parameters"


# # For a given set of parameter values, calculate the rate of loss and decline of heteroplasmy prevalence:
def invasion(z, n, g, d, e, fitness, probabilities, generations):
    r = sum(probabilities * drift(n, intra_selection(g, d, e, z))[0])
    het_frxn = np.array([1.0])
    for k in generations:
        progeny_het = fitness * het_frxn[-1]
        progeny_hom = r * het_frxn[-1]
        progeny_tot = progeny_het + progeny_hom
        het_frxn = np.concatenate((het_frxn, np.array([progeny_het / (progeny_tot + (1 - het_frxn[-1]))])))
    return het_frxn, r

if __name__ == '__main__':
    # ADJUSTABLE SETTINGS:
    # Exclude 'simulated_genotype' from genotypes list if no downstream analysis is desired for a simulated genotype:
    genotypes = ['uaDf5', 'mptDf2', 'mpt4', 'mpt2', 'mptDf3']
    # genotypes = ['uaDf5']

    # Simulated experimental and model parameters:
    org_sel_replicates = 8
    distribution_reps = 320

    # Bottleneck and alpha through gamma, plus initial population-wide heteroplasmic fraction:
    simulated_parameters = [50, 12, 6, 8, 2, -6, 0.5]

    # Single-organism bottleneck frequency (for population-genetic simulations):
    org_bttlnk_freq = 0.5

    # Variable bottleneck sizes (for simulating the effects of drift):
    sim_n = [62, 125, 250, 500, 1000]

    # Thousand-increment range of mutant frequencies from 0 to 1, for downstream computations:
    z_mu = np.linspace(0, 1.0, 1000)

    # Number of progeny pooled and averaged per parent, from the intra-organismal selection experiment:
    ipool = 3

    # Lower & upper bounds of drift parameter search space (see max_likelihood_values for extended search space):
    min_bottleneck = 4
    max_bottleneck = 100

    # Number of bootstrap replicates:
    simulations = 100


    # IMPORT DATA FOR ANALYSIS:
    # Select file directory (file_dir) containing the subsequent csv files for downstream analysis:
    file_dir = '/USERS/WORKING_DIRECTORY'
    intr_org_data = pd.read_csv(file_dir + '/suborganismal_selection_raw_data.csv')
    organism_data = pd.read_csv(file_dir + '/organismal_selection_raw_data.csv')
    distributions = pd.read_csv(file_dir + '/heteroplasmy_dists.csv')
    t_values = pd.read_csv(file_dir + '/t_values_at_95.csv')
    theor_par = pd.read_csv(file_dir + f'/parameters_for_fitness_and_distribution_functions.csv')


    # WORKFLOW AND OUTPUTS FOR DATA ANALYSIS:
    for genotype in genotypes:
        results = file_dir + f'/{genotype}_results'
        if not os.path.exists(results):
            os.makedirs(results)
        theor_dir = results + f'/{genotype}_variable_N_on_specified_selection_parameters'
        if not os.path.exists(theor_dir):
            os.makedirs(theor_dir)


        # EXTRACT & FORMAT DATA FOR DOWNSTREAM MODELING & ANALYSIS:
        # For simulated genotypes, some data (e.g. observed parent and progeny frequencies) is not applicable:
        if genotype == 'simulated_genotype':
            z_obsrvd = "N/A"
            z_parent = "N/A"
            z_progeny = "N/A"
            w_intra = "N/A"

            # Use uaDf5 (or other empirical genotype) for competition experiment parameters (duration, timepoints):
            comp_data = organismal_fitness_calculator('uaDf5')
            org_w = comp_data[0]
            replicates = org_sel_replicates
            dboot = distribution_reps
            iboot = distribution_reps
            fxns_of_parameters(genotype, theor_par, theor_dir, z_mu, ["N for row"] + list(map(str, sim_n)))

        # For empirical genotypes, use data from input files corresponding to each specific genotype:
        else:
            z_obsrvd = distributions[genotype].dropna()
            z_parent = np.array(intr_org_data[genotype + '_gen1'].dropna())
            z_progeny = np.array(intr_org_data[genotype + '_gen2'].dropna())
            w_intra = (z_progeny / z_parent) / ((1 - z_progeny) / (1 - z_parent))
            comp_data = organismal_fitness_calculator(genotype)
            org_w = comp_data[0]
            replicates = len(org_w)
            dboot = len(z_obsrvd)
            iboot = len(z_parent)
        noncomp_z = comp_data[2]
        timepoints = comp_data[4]
        duration = len(noncomp_z)
        timepts_sim = np.ones(replicates) * duration

        # Format org selection data as 2D array (& convert to log ratios if starting with raw mutant frequencies):
        org_log_ratio_data = file_dir + f'/{genotype}_organismal_selection_data.csv'
        if os.path.exists(org_log_ratio_data):
            log_ratio = np.array(pd.read_csv(org_log_ratio_data))
            log_ratio = log_ratio[:, 1:]
        else:
            log_ratio = comp_data[3]


        # ESTIMATE PARAMETER VALUES OF THE MODEL GIVEN RAW DATA (OR SELECT PARAMETER VALUES, IF APPLICABLE):
        if genotype == 'simulated_genotype':
            n = simulated_parameters[0]
            a = simulated_parameters[1]
            b = simulated_parameters[2]
            g = simulated_parameters[3]
            d = simulated_parameters[4]
            e = simulated_parameters[5]
            i = simulated_parameters[6]
        else:
            if not os.path.exists(results + f'/{genotype}_max_likelihood_parameters_for_empirical_data.csv'):
                ml = max_likelihood_values(z_parent, z_progeny, log_ratio, duration, timepoints, z_obsrvd,
                                           min_bottleneck, max_bottleneck,
                                           np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.5]))
                pd.DataFrame(ml).to_csv(results + f'/{genotype}_max_likelihood_parameters_for_empirical_data.csv')
            ml = pd.read_csv(results + f'/{genotype}_max_likelihood_parameters_for_empirical_data.csv')
            n = int(ml[ml['log(L)'] == max(ml['log(L)'])]['N'])
            a = float(ml[ml['N'] == n]['alpha'])
            b = float(ml[ml['N'] == n]['beta'])
            g = float(ml[ml['N'] == n]['gamma'])
            d = float(ml[ml['N'] == n]['delta'])
            e = float(ml[ml['N'] == n]['epsilon'])
            i = float(ml[ml['N'] == n]['intl het/hom'])
            z = np.arange(1 / n, 1 + 1 / n, 1 / n)
            while len(z) > n:
                z = z[:-1]

            # ESTIMATE CONFIDENCE IN MODEL PARAMETERS BY SIMULATING DATA FOR BOOTSTRAPPING (IF NOT ALREADY DONE):
            if not os.path.exists(results + f'/{genotype}_dist{dboot}_intra{iboot}_org{replicates}_sim_results.csv'):
                pd.DataFrame(columns=['N', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'intl het/hom', 'log(L)']).\
                    to_csv(results + f'/{genotype}_dist{dboot}_intra{iboot}_org{replicates}_sim_results.csv')
                bs = pd.read_csv(results + f'/{genotype}_dist{dboot}_intra{iboot}_org{replicates}_sim_results.csv')
                bml = boots_ml(results, z, z_mu, n, a, b, g, d, e, i, ipool, iboot,
                               dboot, log_ratio, replicates, duration, timepts_sim, min_bottleneck, max_bottleneck)
                bs = pd.concat([bs, bml])
                bs.to_csv(results + f'/{genotype}_dist{dboot}_intra{iboot}_org{replicates}_sim_results.csv',
                          index=False)
            else:
                bs = pd.read_csv(results + f'/{genotype}_dist{dboot}_intra{iboot}_org{replicates}_sim_results.csv')

            # GENERATE & SAVE FREQ & FITNESS FUNCTIONS OF BOOTSTRAP PARAMETERS:
            log_2_n_scale = ['N for row', 'N', '2N or 1k', '4N or 1k', '8N or 1k']
            fxns_of_parameters(genotype, ml[ml['log(L)'] == max(ml['log(L)'])], theor_dir, z_mu, log_2_n_scale)
            theory_boots_dir = results + f'/{genotype}_variable_N_on_bootstraps'
            if not os.path.exists(theory_boots_dir):
                os.makedirs(theory_boots_dir)
            fxns_of_parameters(genotype, bs, theory_boots_dir, z_mu, log_2_n_scale)


        # GENERATE HISTOGRAM OF EMPIRICAL MUTANT FREQUENCIES:
        z_hist = np.histogram(z_obsrvd, bins=np.arange(1/50, 1 + 1/25, 1/25))
        pd.DataFrame(zip(z_hist[1], z_hist[0]),
                     columns=['frequency', 'density']).to_csv(results + f'/{genotype}_z_histogram.csv')


        # SAVE POPULATION GENETIC STATS (FITNESS, FREQUENCY, DRIFT, MUTANT INVASION), RAW & SIMULATED DATA:
        print(f"{genotype} max L: N={n}, alpha={a}, beta={b}, gamma={g}, delta={d}, epsilon={e}")
        # Specified mutant frequency distribution following single-organism bottleneck (for simulation purposes):
        org_btlk_prb = np.zeros(n)
        org_btlk_prb[np.where(abs(z-org_bttlnk_freq) == min(abs(z-org_bttlnk_freq)))[0][0]] = 1

        popgen = population_genetics_model(z, n, a, b, g, d, e)
        sim = pd.read_csv(results + f'/{genotype}_dist{dboot}_intra{iboot}_org{replicates}_sim_results.csv')
        z_steps = np.arange(0, 1.005, 0.005)

        # Create arrays of max-likelihood organismal & intra-organismal fitness functions & frequency distributions:
        org_sims = [z_steps, np.ones(201), org_selection(z_steps, a, b)]
        int_sims = [z_steps, np.ones(201), abs(d) * (1 - expit(z_steps * g + e))]
        den_sims = [z_steps, np.interp(z_steps, z_mu, popgen[1])]

        # Create dataframe of population-genetic stats (fitness, drift parameter, loss rate, invasion dynamics):
        stat_gens = np.concatenate((np.arange(0, 100), np.array([125, 250, 500, 1000, 2500, 5000, 10000])))
        invastat = invasion(z, n, g, d, e, popgen[0], popgen[2], stat_gens[1:])

        genpts = []
        for generation in stat_gens:
            genpts.append(f"generation {generation}")
        stats = pd.DataFrame(columns=['line', 'empirical'])
        stats = pd.concat([stats, pd.DataFrame(zip(["fitness", "n (drift)", "loss rate"],
                                                   [popgen[0], n, invastat[1]]), columns=['line', 'empirical'])])
        stats = pd.concat([stats, pd.DataFrame(zip(np.array(genpts), invastat[0]), columns=['line', 'empirical'])])

        # Max-likelihood organismal and intra-organismal fitness:
        ml_org_w = org_selection(z_mu, a, b)
        ml_intra_w = np.array(abs(d) * (1 - expit(z_mu * g + e)))

        # Calculate the mutant frequency at which multilevel selection forces are balanced:
        ml_bal_sel_w = (ml_intra_w - (1 - ml_org_w)) - 1
        idx_mlvl_bal = np.abs(ml_bal_sel_w).argmin()
        z_at_mlvl_bal = [z_mu[idx_mlvl_bal]]

        # Calculate the mutant frequency at which intra-organismal balancing selection is reached:
        idx_intra_bal = np.abs(ml_intra_w - 1).argmin()
        z_at_intra_bal = [z_mu[idx_intra_bal]]

        # Calculate the frequency at which organismal fitness is half that of wildtype organisms:
        idx_hf_w = np.abs(ml_org_w - 0.5).argmin()
        z_at_half_w = [frequency_of_fitness(0.5, a, b, np.array([z_mu[idx_hf_w]]))]

        # Calculate probabilistic & statistical results on mutant frequency (mean, std dev):
        z_mu_sq = z_mu ** 2
        exp_z = sum(z_mu * (popgen[1] / sum(popgen[1])))
        exp_z_sq = sum(z_mu_sq * (popgen[1] / sum(popgen[1])))
        var_z = exp_z_sq - (exp_z ** 2)
        mean_z = [exp_z]
        stdv_z = [np.sqrt(var_z)]

        # Repeat above calculations for each bootstrap data set:
        for ind in sim.index:
            n_sim = int(sim['N'][ind])
            a_sim = float(sim['alpha'][ind])
            b_sim = float(sim['beta'][ind])
            g_sim = float(sim['gamma'][ind])
            d_sim = float(sim['delta'][ind])
            e_sim = float(sim['epsilon'][ind])
            z_sim = np.arange(1 / n_sim, 1 + 1 / n_sim, 1 / n_sim)
            while len(z_sim) > n_sim:
                z_sim = z_sim[:-1]
            popgen_s = population_genetics_model(z_sim, n_sim, a_sim, b_sim, g_sim, d_sim, e_sim)
            invasim = invasion(z_sim, n_sim, g_sim, d_sim, e_sim, popgen_s[0], popgen_s[2], stat_gens[1:])
            simstat = np.concatenate((np.array([popgen_s[0], n_sim, invasim[1]]), invasim[0]))
            stats[f'sim rep {ind + 1}'] = simstat

            # Append organismal and intra-organismal fitness, and frequency distribution, to respective arrays:
            org_sims = np.vstack([org_sims, org_selection(z_steps, a_sim, b_sim)])
            int_sims = np.vstack([int_sims, abs(d_sim) * (1 - expit(z_steps * g_sim + e_sim))])
            den_sims = np.vstack([den_sims, np.interp(z_steps, z_mu, popgen_s[1])])

            # Max-likelihood organismal and intra-oragnismal fitness:
            ml_org_w_sim = org_selection(z_mu, a_sim, b_sim)
            ml_intra_w_sim = np.array(abs(d_sim) * (1 - expit(z_mu * g_sim + e_sim)))

            # Calculate the mutant frequency at which multilevel selection forces are balanced:
            ml_bal_sel_w_sim = (ml_intra_w_sim - (1 - ml_org_w_sim)) - 1
            idx_mlvl_bal_sim = np.abs(ml_bal_sel_w_sim).argmin()
            z_at_mlvl_bal.append(z_mu[idx_mlvl_bal_sim])

            # Calculate the mutant frequency at which intra-organismal balancing selection is reached:
            idx_intra_bal_sim = np.abs(ml_intra_w_sim - 1).argmin()
            z_at_intra_bal.append(z_mu[idx_intra_bal_sim])

            # Calculate the frequency at which organismal fitness is half that of wildtype organisms:
            idx_hf_w_sim = np.abs(ml_org_w_sim - 0.5).argmin()
            z_at_half_w.append(frequency_of_fitness(0.5, a_sim, b_sim, np.array([z_mu[idx_hf_w_sim]])))

            # Calculate probabilistic & statistical results on mutant frequency (mean, std dev, prob distribution):
            prob_dist = popgen_s[1] / sum(popgen_s[1])
            exp_z_sim = sum(z_mu * prob_dist)
            exp_z_sq_sim = sum(z_mu_sq * prob_dist)
            var_z_sim = exp_z_sq_sim - (exp_z_sim ** 2)
            mean_z.append(exp_z_sim)
            stdv_z.append(np.sqrt(var_z_sim))

        # Save results generated above to output files:
        if genotype == 'simulated_genotype':
            stat_dir = results + f'/{genotype}_simulated_data_for_bootstrapping_{iboot}intra_{replicates}org'
        else:
            stat_dir = results
        if not os.path.exists(stat_dir):
            os.makedirs(stat_dir)
        pd.DataFrame(np.array(z_at_mlvl_bal)).to_csv(stat_dir + f'/{genotype}_freq_at_multilevel_balancing.csv')
        pd.DataFrame(np.array(z_at_intra_bal)).to_csv(stat_dir + f'/{genotype}_freq_at_intra_org_balancing.csv')
        pd.DataFrame(np.array(z_at_half_w)).to_csv(stat_dir + f'/{genotype}_freq_when_mutant_loses_half_fitness.csv')
        pd.DataFrame(np.array([popgen[1]]).T).to_csv(stat_dir + f'/{genotype}_max_l_dist.csv')
        pd.DataFrame(np.array(mean_z)).to_csv(stat_dir + f'/{genotype}_mean_z.csv')
        pd.DataFrame(np.array(stdv_z)).to_csv(stat_dir + f'/{genotype}_stdv_z.csv')
        pd.DataFrame(np.array(org_sims).T).to_csv(stat_dir + f'/{genotype}_org_fitness.csv')
        pd.DataFrame(np.array(int_sims).T).to_csv(stat_dir + f'/{genotype}_intra_fitness.csv')
        pd.DataFrame(np.array(den_sims).T).to_csv(stat_dir + f'/{genotype}_distributions.csv')
        stats.to_csv(stat_dir + f'/{genotype}_population_stats_on_maintenance_of_heteroplasmic_state.csv')