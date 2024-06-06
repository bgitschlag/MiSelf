# This script will perform population-genetic modeling of a Wright-Fisher process to calculate the most
# evolutionarily stable population-wide frequency distribution of a heteroplasmic mutation, given
# intra-organismal & organismal fitness functions (defined by model parameters), a neutral drift parameter,
# and a sample of mutant frequencies among heteroplasmic individuals. Model parameters are estimated by
# maximum likelihood using experimental data and their uncertainty is estimated by parametric bootstrapping.

import os
import numpy as np
import pandas as pd
import math
import time
from scipy import stats
from scipy.linalg import eig
from scipy.optimize import minimize
from scipy.special import betainc
from scipy.special import expit


# ***ATTENTION: THESE ARE THE ADJUSTABLE SETTINGS FOR RUNNING THIS SCRIPT***
# Do you want to model a theoretical genotype with pre-set parameters? Enter "False" if analyzing empirical data:
simulate_genotype = True
# Do you want to re-estimate the model parameters by simulating data?
run_bootstrapping = True
# Do you want to fit a distribution to the intra-organismal selection data that is constrained on [0,1]?
beta_distribution = False

# Genotype(s) for downstream analysis:
if simulate_genotype:
    genotypes = ['simulated_genotype']
else:
    genotypes = ['uaDf5', 'mpt4', 'mptDf2', 'mpt2', 'mptDf3']
# Empirical genotype to use for statistical distributions:
empr_gntp = 'uaDf5'
# Parameters for theoretical modeling [N, alpha, beta, gamma, delta, epsilon, initial heteroplasmic fraction]:
sim_params = [50, 12, 4, 8, 2, -6, 0.5]
# Simulated experiment parameters:
# Inter-organismal selection duration (number of generations of simulated data):
org_dur = 8
# Inter-organismal replicate populations:
org_reps = 8
# Number of simulated draws from the stationary mutant frequency distribution:
dist_reps = 80
# Number of parents for intra-organismal selection (also drawn from stationary mutant frequency distribution):
intr_reps = 30
# Number of bins for empirical mutant frequency histogram:
bins = 24
# Number of progeny pooled and averaged per parent, for the intra-organismal selection experiment:
ipool = 3
# Lower & upper bounds of drift parameter search space (see max_likelihood_search for extended search space):
min_n = 10
max_n = 96
# Organism effective population sizes to model, for modeling invasion dynamics:
pop_sizes = [10, 100, 1000, 10000]
# Number of bootstrap replicates:
bootstrap_reps = 10


# SPECIFY WORKING FILE DIRECTORY AND IMPORT DATA FOR ANALYSIS:
dir = '/WORKING_DIRECTORY'
distributions = pd.read_csv(dir + '/Source_data/mutant_frequency_samples_source_data.csv')


# Specify selection data & dimensions (duration, replicates):
def data_processing(genotype, empir, dir, org_dur, org_reps, dist_reps, intr_reps):

    # For simulated genotypes, some data (e.g. pertaining to empirical selection measurements) are not applicable:
    if genotype == 'simulated_genotype':
        log_r = np.array(pd.read_csv(dir + f'/Source_data/{empir}_organismal_source_data.csv',
                                     usecols=lambda x: x != 'Generation')[:org_dur], dtype=float)
        intra_sel_data = ["N/A", "N/A", "N/A", "N/A"]
        data_dimensions = [org_dur, org_reps, dist_reps, intr_reps]

    # For empirical genotypes, use data from input files corresponding to each specific genotype:
    else:
        log_r = pd.read_csv(dir + f'/Source_data/{genotype}_organismal_source_data.csv')
        del log_r['Generation']
        log_r = log_r.transform(pd.to_numeric, errors='coerce')
        log_r = np.array(log_r, dtype=float)
        intr_org_data = pd.read_csv(dir + f'/Source_data/{genotype}_intra_organismal_source_data.csv')
        z_obs = distributions[genotype].dropna()
        z_par = np.array(intr_org_data['Parent'])
        z_pro = np.array(intr_org_data['Adult progeny'])
        w_int = (z_pro / z_par) / ((1 - z_pro) / (1 - z_par))
        intra_sel_data = [z_obs, z_par, z_pro, w_int]
        data_dimensions = [len(log_r), log_r.shape[1], len(z_obs), len(z_par)]

    return intra_sel_data, data_dimensions, log_r


# Generate a set of discrete, evenly-spaced mutant frequencies between 0 & 1 for downstream modeling:
def z_steps(n):
    z = np.arange(1 / n, 1 + 1 / n, 1 / n)
    while len(z) > n:
        z = z[:-1]
    return z


# Find mutant frequency in generation n+1 given frequency in gen n and intra-organismal selection:
def intra_selection(g, d, e, z):
    weighted_fitness = z * (d * (1 - expit(z * g + e)))
    return weighted_fitness / (weighted_fitness + (1 - z))


# Model neutral drift by generating an n-sized binomial distribution centered on input frequency z:
def drift(n, z):
    binom = []
    for k in range(0, n + 1):
        binom.append(np.array(np.nan_to_num(math.comb(n, k) * ((1 - z) ** (n - k)) * (z ** k)), dtype=float))
    return np.array(binom)


# Model organismal selection using the cumulative density function of the beta distribution:
def org_selection(z, a, b):
    fitness = np.array(1 - betainc(a, b, z), dtype=float)
    return np.nan_to_num(fitness)


# Using selection & drift parameters, find mean organism fitness and most stable mutant frequency distribution:
def population_genetics_model(z, n, a, b, g, d, e):

    # Generate principal submatrix of 2D probability distribution generated by multilevel selection & drift:
    intra_organism_dynamics = drift(n, intra_selection(g, d, e, z))
    matrix = (intra_organism_dynamics[1:-1].T[:-1] * org_selection(z, a, b)[:-1]).T
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

    # Generate a smooth probability-density distribution from the eigenvector:
    gaussian = np.zeros(1000)
    for bin in np.arange(0, n, 1):
        gaussian += (stats.norm.pdf(np.linspace(0, 1.0, 1000), loc=z[bin], scale=1 / n) * prb[bin])

    return fit, gaussian, prb, intra_organism_dynamics


# Sample mutant frequencies from a given input probability-density distribution:
def sample_density_dist(dens, steps, draws):
    cumulative = np.cumsum(dens)
    cumulative *= 1 / cumulative[-1]
    return np.interp(np.random.uniform(size=draws), cumulative, steps)


# Find log-likelihood of actual data given expected data from the model:
def parameters_likelihood(expected, actual, sample_size):
    ssq = sum((expected - actual) ** 2)
    var = ssq / sample_size
    return (sample_size / 2) * (-np.log(2 * np.pi) - np.log(var) - 1)


# Find log-likelihood of a sampled for a given frequency distribution given by the model:
def distribution_likelihood(sample, densities):
    likelihoods = np.interp(np.array(sample), np.linspace(0, 1.0, 1000), densities)
    likelihoods[likelihoods < 0] = 0
    return sum(list(np.log(likelihoods)))


# Find log-likelihood of intra- and inter-organismal selection and mutant frequency distribution:
def likelihood(z, n, a, b, g, d, e, het_frxn, z1, z2, log_r, dur, z_obs):

    # Log-likelihood of progeny mutant frequencies given their expected frequencies:
    exp_z2 = intra_selection(g, d, e, z1)
    if beta_distribution:
        l_intra = 0
        for zi in range(len(z1)):
            m, v = exp_z2[zi], exp_z2[zi] * (1 - exp_z2[zi]) / (ipool * n)
            aa = m * ((m * (1-m) / v) - 1)
            bb = aa * ((1 - m) / m)
            l_intra += np.log(stats.beta.pdf(z2[zi], aa, bb))
    else:
        l_intra = parameters_likelihood(exp_z2, z2, len(z2))

    # Log-likelihood of inter-organismal selection data given expected vs actual log ratio of genotypes:
    popgen = population_genetics_model(z, n, a, b, g, d, e)
    r = sum(popgen[2] * popgen[3][0])
    l_org = 0
    exp_log_r = []
    for generation in range(0, dur):
        hom_frxn = 1 - het_frxn
        exp_log_r.append(np.log(het_frxn / hom_frxn))
        progeny_het = popgen[0] * het_frxn
        progeny_hom = r * het_frxn
        het_frxn = progeny_het / (progeny_het + progeny_hom + hom_frxn)
    for lr in log_r.T:
        lr = lr[~(np.isnan(lr))]
        l_org += parameters_likelihood(np.array(exp_log_r[:len(lr)]), lr, len(lr))

    # Log-likelihood of sampled mutant frequencies for a distribution given by the model:
    l_dist = distribution_likelihood(z_obs, popgen[1])

    return l_intra + l_org + l_dist


# Find maximum-likelihood values of continuous parameters given input data & a specified drift parameter (n):
def maximize_likelihood(z, n, z1, z2, log_r, dur, z_obs, itl):

    def nl(x):
        lkhd = likelihood(z, n, x[0], x[1], x[2], x[3], x[4], x[5], z1, z2, log_r, dur, z_obs)
        return -lkhd
    m = minimize(nl, itl, method="Nelder-Mead",
                 bounds=np.array([[1e-16, np.inf], [1e-16, np.inf], [-np.inf, np.inf], [0.0, np.inf],
                                  [-np.inf, np.inf], [0.01, 0.99]]), options={"adaptive": True})

    return m['x'][0], m['x'][1], m['x'][2], m['x'][3], m['x'][4], m['x'][5], -m['fun']


# Find max-likelihood values of continuous parameters given input data across for a grid of drift parameter values:
def max_likelihood_search(phase, z1, z2, log_r, dur, z_obs, min_n, max_n, itl):

    # Generate table of parameter values for input data and a given drift parameter value:
    def mlp_for_n(n, itl):
        start_time = time.time()
        ml = maximize_likelihood(z_steps(n), n, z1, z2, log_r, dur, z_obs, itl)
        mlp = pd.DataFrame(zip([n], [ml[0]], [ml[1]], [ml[2]], [ml[3]], [ml[4]], [ml[5]], [ml[6]]),
                           columns=['N', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'intl het/hom', 'log(L)'])
        print(f'N={n} alpha={round(ml[0], 4)} beta={round(ml[1], 4)} gamma={round(ml[2], 4)} '
              f'delta={round(ml[3], 4)} epsilon={round(ml[4], 4)} intl het/hom={round(ml[5], 4)} '
              f'log(L)={round(ml[6], 4)}, N={n} took {round(time.time() - start_time, 2)} seconds')
        return mlp, np.array([ml[0], ml[1], ml[2], ml[3], ml[4], ml[5]]), ml[6]

    # Grid search across n (drift parameter) values, generating max-likelihood parameter estimates for each n:
    ml_par = pd.DataFrame(columns=['N', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'intl het/hom', 'log(L)'])
    for n in np.arange(min_n, max_n + 1, 1):
        if phase == "bootstrapping":
            ml_par = merge(ml_par, mlp_for_n(n, itl)[0])
        else:
            par_adp, par_unb = mlp_for_n(n, itl), mlp_for_n(n, np.array([1.0, 1.0, 0.0, 2.0, 0.0, 0.5]))
            if par_adp[2] >= par_unb[2]:
                ml_par = merge(ml_par, par_adp[0])
                itl = par_adp[1]
            else:
                ml_par = merge(ml_par, par_unb[0])
                itl = np.array([1.0, 1.0, 0.0, 2.0, 0.0, 0.5])

    # Expand search space of n (to 100, then increments of 125*sqrt(2)^int) until max likelihood peaks or n=1000:
    if int(ml_par.loc[ml_par['log(L)'] == max(ml_par['log(L)'])].iloc[0]['N']) >= max_n - 1:
        for n in range(max_n + 1, max_n + 5):
            ml_par = pd.concat([ml_par, mlp_for_n(n, itl)[0]])
    if int(ml_par.loc[ml_par['log(L)'] == max(ml_par['log(L)'])].iloc[0]['N']) > max_n:
        x = 1
        xtnd_n = np.concatenate((np.array([125]), np.int_(np.around(np.array([125 * np.sqrt(2) ** x])))))
        for n in xtnd_n:
            ml_par = pd.concat([ml_par, mlp_for_n(n, itl)[0]])
        while int(ml_par.loc[ml_par['log(L)'] == max(ml_par['log(L)'])].iloc[0]['N']) >= xtnd_n[-2]:
            x += 1
            xtnd_n = np.concatenate((xtnd_n, np.int_(np.around(np.array([125 * np.sqrt(2) ** x])))))
            ml_par = pd.concat([ml_par, mlp_for_n(xtnd_n[-1], itl)[0]])
            if x == 6:
                break

    return ml_par


# Generate simulated data given input parameters and input data:
def bootstrap(popgen, z, z_mu, n, g, d, e, i, ipool, iboot, dboot, log_r, reps, dur):

    # From intra-organismal selection & drift parameters, simulate parent & progeny mutant frequencies:
    z1 = sample_density_dist(popgen[1], z_mu, iboot)
    z2 = []
    for z1_i in z1:
        intra = drift(n, intra_selection(g, d, e, z1_i))
        density = n * (intra[1:] / sum(intra[1:]))
        draws = sample_density_dist(np.interp(z_mu, z, density), z_mu, ipool)
        z2.append(sum(draws) / ipool)

    # From model parameters and empirical variance, simulate organismal selection data (log ratio of genotypes):
    r = sum(popgen[2] * drift(n, intra_selection(g, d, e, z))[0])
    het_frxn = i
    exp_lr = []
    for generation in range(0, dur):
        hom_frxn = 1 - het_frxn
        exp_lr.append(np.log(het_frxn / hom_frxn))
        progeny_tot = popgen[0] * het_frxn
        het_frxn = (progeny_tot * (1 - r)) / (progeny_tot + hom_frxn)
    s_org = 0
    for lr in log_r.T:
        lr = lr[np.logical_not(np.isnan(lr))]
        s_org += (sum((np.array(exp_lr[:len(lr)]) - lr) ** 2) / len(lr)) ** 0.5
    org_sim = []
    for i in range(0, reps):
        ln = []
        for timepoint in range(0, dur):
            ln.append(float(np.random.normal(np.array(exp_lr)[timepoint], scale=s_org / reps, size=1)[0]))
        org_sim.append(np.array(ln))

    # Simulate sample draws from the predicted mutant frequency distribution (final returned item):
    return z1, np.array(z2), np.array(org_sim).T, sample_density_dist(popgen[1], z_mu, dboot)


# Generate tables of simulated (bootstrap) data & of their corresponding max-likelihood parameter values:
def max_likelihood_estimates_of_bootstraps(genotype, dir, labels, popgen, z, z_mu, n, a, b, g, d, e, i,
                                           log_l, ipool, iboot, dboot, log_r, reps, dur, min_n, max_n):

    # Generate and define output files and data structures:
    if genotype == 'simulated_genotype':
        data_type, l = "ground_truth", "N/A"
        initialize = np.array([1.0, 1.0, 0.0, 2.0, 0.0, 0.5])
    else:
        data_type, l = "empirical", log_l
        initialize = np.array([a, b, g, d, e, i])
    folder = dir + f'/{genotype}_simulated_data_for_bootstrapping/'
    org_folder = folder + f'/{genotype}_simulated_organismal_data/'
    if not os.path.exists(org_folder):
        os.makedirs(org_folder)
    if os.path.exists(folder + f'/{genotype}_dist_bootstraps.csv'):
        dist_sims = list(np.array(pd.read_csv(folder + f'/{genotype}_dist_bootstraps.csv')).T)
        z1_sims = list(np.array(pd.read_csv(folder + f'/{genotype}_parent_bootstraps.csv')).T)
        z2_sims = list(np.array(pd.read_csv(folder + f'/{genotype}_progeny_bootstraps.csv')).T)
        ml_par = pd.DataFrame(columns=['data'] + labels + ['log(L)'])
    else:
        dist_sims, z1_sims, z2_sims = [], [], []
        ml_par = pd.DataFrame([[data_type, n, a, b, g, d, e, i, l]], columns=['data'] + labels + ['log(L)'])

    # For each in a number of simulated data-set replicates, simulate a data set & estimate its parameter values:
    for s in range(0, bootstrap_reps):
        sim = bootstrap(popgen, z, z_mu, n, g, d, e, i, ipool, iboot, dboot - iboot, log_r, reps, dur)
        z1_sims.append(sim[0])
        z2_sims.append(sim[1])
        dist_sims.append(np.concatenate((sim[0], sim[3])))
        boot_rep = 1
        while os.path.exists(org_folder + f'/{genotype}_org_rep_{boot_rep}.csv'):
            boot_rep += 1
        pd.DataFrame(sim[2]).to_csv(org_folder + f'/{genotype}_org_rep_{boot_rep}.csv')
        ml = max_likelihood_search("bootstrapping", sim[0], sim[1], sim[2], dur,
                                   np.concatenate((sim[0], sim[3])), min_n, max_n, initialize)
        ml_i = ml[ml['log(L)'] == max(ml['log(L)'])]
        ml_i.insert(loc=0, column='data', value=f"bootstrap")
        ml_par = merge(ml_par, ml_i)
        print(f"data simulated for N={n} alpha={a} beta={b} gamma={g} delta={d} epsilon={e}; max-L={ml_i}")

    # Save simulated (bootstrap) data to csv files:
    pd.DataFrame(np.array(z1_sims).T).to_csv(folder + f'/{genotype}_parent_bootstraps.csv', index=False)
    pd.DataFrame(np.array(z2_sims).T).to_csv(folder + f'/{genotype}_progeny_bootstraps.csv', index=False)
    pd.DataFrame(np.array(dist_sims).T).to_csv(folder + f'/{genotype}_dist_bootstraps.csv', index=False)

    return ml_par


# Define parameter values for bootstrap data, by estimating (first time) or retrieving (if estimated already):
def bootstrap_parameters(genotype, dir, labels, popgen, z, z_mu, n, a, b, g, d, e,
                         i, l, ipool, dboot, iboot, log_r, dur, reps, min_n, max_n):

    def bootstraps_output(output_name):
        sim_results = dir + f'/{genotype}_dist{dboot}_intra{iboot}_org{reps}_{output_name}.csv'
        if not os.path.exists(sim_results):
            pd.DataFrame(columns=['data'] + labels + ['log(L)']).to_csv(sim_results)
        return pd.read_csv(sim_results)

    bs = bootstraps_output('sim_results')
    bml = max_likelihood_estimates_of_bootstraps(genotype, dir, labels, popgen, z, z_mu, n, a, b, g, d, e,
                                                 i, l, ipool, iboot, dboot, log_r, reps, dur, min_n, max_n)
    bs = merge(bs, bml)
    bs.to_csv(dir + f'/{genotype}_dist{dboot}_intra{iboot}_org{reps}_sim_results.csv', index=False)

    return bs


# Generate range of values of n for each set of selection parameters, to model drift:
def vary_n(selection_parameters, par_labels, spec_n):
    par_table = pd.DataFrame(columns=par_labels)
    if 'N' in selection_parameters:
        selection_parameters = selection_parameters.drop(['N'], axis=1)
    n_set = spec_n * (2 ** np.arange(4))
    np.put(n_set, np.where(n_set > 1000), 1000)
    for p in range(len(selection_parameters)):
        for n in n_set:
            all_par = pd.DataFrame([[n] + list(selection_parameters.iloc[p])], columns=par_labels)
            par_table = merge(par_table, all_par)
    return par_table


# From a given set of model parameters, derive the fitness functions and mutant frequency statistical results:
def functions_of_parameters(genotype, parameters, labels, dir, z, z_mu, instruction):
    if 'log(L)' in parameters:
        labels = labels + ['log(L)']
    na = ["N/A"] * len(labels)
    i_w = pd.DataFrame(zip(labels + list(z), na + list(np.ones(201))), columns=['variable', 'neutral'])
    o_w = pd.DataFrame(zip(labels + list(z), na + list(np.ones(201))), columns=['variable', 'neutral'])
    pdf = pd.DataFrame(zip(labels + list(z)), columns=['variable'])
    stat_mmts = [[], []]
    z_mu_sq = z_mu ** 2
    # For a table of model parameter values, generate mutant frequency & fitness data:
    for p in range(len(parameters)):
        n = int(parameters.iloc[p]['N'])
        a = float(eval(str(parameters.iloc[p]['alpha'])))
        b = float(eval(str(parameters.iloc[p]['beta'])))
        g = float(eval(str(parameters.iloc[p]['gamma'])))
        d = float(eval(str(parameters.iloc[p]['delta'])))
        e = float(eval(str(parameters.iloc[p]['epsilon'])))
        i = float(parameters.iloc[p]['intl het/hom'])
        par = [n, a, b, g, d, e, i]
        if 'data' in parameters:
            par = [parameters.iloc[p]['data']] + par
        if 'log(L)' in parameters:
            par = par + [float(parameters.iloc[p]['log(L)'])]
        popgen = population_genetics_model(z_steps(n), n, a, b, g, d, e)
        i_w[f'parameter_set_{p + 1}'] = par + list(np.interp(z, z_mu, d * (1 - expit(z_mu * g + e))))
        o_w[f'parameter_set_{p + 1}'] = par + list(np.interp(z, z_mu, org_selection(z_mu, a, b)))
        pdf[f'parameter_set_{p + 1}'] = par + list(np.interp(z, z_mu, popgen[1]))
        exp_z = sum(z_mu * (popgen[1] / sum(popgen[1])))
        stat_mmts[0].append(exp_z)
        stat_mmts[1].append(np.sqrt(sum(z_mu_sq * (popgen[1] / sum(popgen[1]))) - (exp_z ** 2)))

    # Save fitness & mutant frequency results to output files:
    parameters['mean_z'] = stat_mmts[0]
    parameters['stdv_z'] = stat_mmts[1]
    parameters.to_csv(dir + f'/{genotype}_parameters_and_frequency_stats_{instruction}.csv')
    i_w.to_csv(dir + f'/{genotype}_intra_organismal_fitness_{instruction}.csv')
    o_w.to_csv(dir + f'/{genotype}_organismal_fitness_{instruction}.csv')
    pdf.to_csv(dir + f'/{genotype}_distributions_{instruction}.csv')

    return print(f"finished compiling fitness functions and frequency data for {genotype}_{instruction}")


# For a given set of parameter values, calculate the rate of loss and decline of heteroplasmy prevalence:
def invasion(z, n, g, d, e, fitness, probabilities, generations):
    r = sum(probabilities * drift(n, intra_selection(g, d, e, z))[0])
    het_frxn = np.array([1.0])
    for k in generations:
        progeny_het = fitness * het_frxn[-1]
        progeny_hom = r * het_frxn[-1]
        progeny_tot = progeny_het + progeny_hom
        het_frxn = np.concatenate((het_frxn, np.array([progeny_het / (progeny_tot + (1 - het_frxn[-1]))])))
    return het_frxn, r


# Median time until invasion of a heteroplasmic population by a de novo wildtype lineage due to loss of mutation:
def invasion_timing(fitness, loss_rate, population_sizes):
    s = (1 / fitness) - 1
    inv_labels, inv_times = [], []
    for n_e in population_sizes:
        inv_labels.append(f"med inv time, N={n_e}")
        inv_times.append(int(round(np.log(0.5) / (-2 * s * loss_rate * n_e))))
    return inv_labels, inv_times


# Generate a new data frame to fill with data, or append to existing data frame, as appropriate:
def merge(append_to, to_append):
    if append_to.empty:
        append_to = to_append
    else:
        append_to = pd.concat([append_to, to_append])
    return append_to


# Generate and save empirical mutant frequency histogram:
def histogram(res, z_obs, bins):
    hg = np.histogram(z_obs, bins=np.arange(1 / (2 * bins), 1 + 1 / bins, 1 / bins))
    pd.DataFrame(zip(hg[1], hg[0]), columns=['freq', 'density']).to_csv(res + f'/{genotype}_histogram.csv')
    return


# Find relative organism fitness (cumulative density function of beta distribution) for a given frequency:
def z_at_w(w, a, b, itl_z):
    def neg_log_dif(x):
        fit = 1 - betainc(a, b, x)
        return abs(w - fit)
    freq_of_fit = minimize(neg_log_dif, itl_z, method="Nelder-Mead", bounds=np.array([[0.0, 1.0]]))
    return float(freq_of_fit['x'][0])


# Generate statistical descriptions of population genetic values (fitness, frequency, invasion dynamics) from data:
def pop_gen_stats(popgen, z, z_bins, n, a, b, g, d, e, gens):
    org = org_selection(z, a, b)
    intra = np.array(d * (1 - expit(z * g + e)))
    i_ml_bal = np.abs((intra - (1 - org)) - 1).argmin()
    i_in_bal = np.abs(intra - 1).argmin()
    i_half_w = np.abs(org - 0.5).argmin()
    return invasion(z_bins, n, g, d, e, popgen[0], popgen[2], gens[1:]), i_in_bal, i_ml_bal, i_half_w, org[i_in_bal]


# Specify model parameter values given a set of data (or return pre-specified parameter values if applicable):
def specify_parameters(genotype, dir, sim_par, z_obs, z1, z2, lr, dur, min_n, max_n):
    if genotype == 'simulated_genotype':
        parameters = sim_par + ['N/A']
    else:
        if not os.path.exists(dir + f'/{genotype}_max_likelihood_parameters_for_empirical_data.csv'):
            m = max_likelihood_search("initial_data_analysis", z1, z2, lr, dur, z_obs, min_n,
                                      max_n, np.array([1.0, 1.0, 0.0, 2.0, 0.0, 0.5]))
            pd.DataFrame(m).to_csv(dir + f'/{genotype}_max_likelihood_parameters_for_empirical_data.csv')
        ml_tbl = pd.read_csv(dir + f'/{genotype}_max_likelihood_parameters_for_empirical_data.csv')
        ml_par = ml_tbl.loc[ml_tbl['log(L)'] == max(ml_tbl['log(L)'])]
        ml = ml_par.iloc[0]
        parameters = [int(ml['N']), float(ml['alpha']), float(ml['beta']), float(ml['gamma']),
                      float(ml['delta']), float(ml['epsilon']), float(ml['intl het/hom']), float(ml['log(L)'])]
        print(f"{genotype} max L: N={parameters[0]}, alpha={parameters[1]}, beta={parameters[2]}, "
              f"gamma={parameters[3]}, delta={parameters[4]}, epsilon={parameters[5]}")
    return parameters


# Save population genetic statistics results to output files:
def save_pop_gen_results(genotype, res, eq_z_i, eq_z_m, zhw, w_o_i_bal, st):
    eq_z_i.to_csv(res + f'/{genotype}_freq_at_intra_org_balancing.csv')
    eq_z_m.to_csv(res + f'/{genotype}_freq_at_multilevel_balancing.csv')
    zhw.to_csv(res + f'/{genotype}_freq_at_half_fitness.csv')
    st.to_csv(res + f'/{genotype}_pop_stats_on_maintenance_of_heteroplasmic_state.csv')
    w_o_i_bal.to_csv(res + f'/{genotype}_organismal_fitness_at_intra_org_balancing_freq.csv')
    return print(f"{genotype} population genetic results saved to output folder")


if __name__ == '__main__':
    # WORKFLOW AND OUTPUTS FOR DATA ANALYSIS:
    z_mu, z_200 = np.linspace(0, 1.0, 1000), np.linspace(0, 1.0, 201)
    for genotype in genotypes:
        res = dir + f'/Results/{genotype}_results'
        if not os.path.exists(res):
            os.makedirs(res)


        # EXTRACT & FORMAT DATA FOR DOWNSTREAM MODELING & ANALYSIS:
        dat_prc = data_processing(genotype, empr_gntp, dir, org_dur, org_reps, dist_reps, intr_reps)
        dboot, iboot = dat_prc[1][2], dat_prc[1][3]
        duration, replicates, log_ratio = dat_prc[1][0], dat_prc[1][1], dat_prc[2]
        if not genotype == 'simulated_genotype':
            histogram(res, dat_prc[0][0], bins)


        # ESTIMATE MODEL PARAMETER VALUES GIVEN RAW DATA (OR SELECT PRE-SPECIFIED PARAMETER VALUES, IF APPLICABLE):
        labels = ['N', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'intl het/hom']
        par = specify_parameters(genotype, res, sim_params, dat_prc[0][0], dat_prc[0][1],
                                 dat_prc[0][2], log_ratio, duration, min_n, max_n)
        n, a, b, g, d, e, i, l = par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]
        par_tbl = pd.DataFrame([[n, a, b, g, d, e, i]], columns=labels)


        # GENERATE & SAVE FUNCTIONS OF PARAMETERS (FITNESS & FREQUENCY DISTRIBUTIONS & DESCRIPTIVE STATISTICS):
        if genotype == 'simulated_genotype':
            dtp = 'ground_truth'
            labeling = f'a{round(a, 3)}_b{round(b, 3)}_g{round(g, 3)}_d{round(d, 3)}_e{round(e, 3)}'
            res = res + f'/{genotype}_results_for_{labeling}'
            if not os.path.exists(res):
                os.makedirs(res)
            b_res = res + f'/{genotype}_bootstrap_results_for_{labeling}/' \
                          f'{genotype}_data_for_{iboot}_intra_{replicates}_org'
            if not os.path.exists(b_res):
                os.makedirs(b_res)
            functions_of_parameters(genotype, par_tbl, labels, res, z_200, z_mu, "theoretical_parameters")
        else:
            dtp = 'empirical'
            b_res = res
            par_var_n = vary_n(par_tbl, labels, n)
            functions_of_parameters(genotype, par_var_n, labels, res, z_200, z_mu, "log2_scale_var_N_to_1k")
            functions_of_parameters(genotype, par_tbl, labels, res, z_200, z_mu, "function_of_max_l_parameters")


        # GENERATE POPULATION-GENETIC STATISTICS ON SELECTION & INVASION DYNAMICS ON EMPIRICAL DATA:
        gens = np.concatenate((np.arange(0, 100), np.array([125, 250, 500, 1000, 2500, 5000, 10000])))
        pgm = population_genetics_model(z_steps(n), n, a, b, g, d, e)
        pgs = pop_gen_stats(pgm, z_mu, z_steps(n), n, a, b, g, d, e, gens)
        c_zin, c_zml = ['data', 'intra_eqlbrm_freq'], ['data', 'mlvl_eqlbrm_freq']
        c_zhw, w_ibl = ['data', 'freq_half_w'], ['data', 'w_org_at_intra_bal']
        w, loss = pgm[0], pgs[0][1]
        eq_z_i = pd.DataFrame([[dtp, z_mu[pgs[1]]]], columns=c_zin)
        eq_z_m = pd.DataFrame([[dtp, z_mu[pgs[2]]]], columns=c_zml)
        w_o_i_bal = pd.DataFrame([[dtp, pgs[4]]], columns=w_ibl)
        zhw = pd.DataFrame([[dtp, z_at_w(0.5, a, b, np.array([z_mu[pgs[3]]]))]], columns=c_zhw)
        gnrtns = []
        for generation in gens:
            gnrtns.append(f"generation {generation}")
        inv = invasion_timing(w, loss, pop_sizes)
        st = pd.DataFrame(zip(["fitness", "N", "loss rate"] + inv[0], [w, n, loss] + inv[1]), columns=['line', dtp])
        st = pd.concat([st, pd.DataFrame(zip(np.array(gnrtns), pgs[0][0]), columns=['line', dtp])])
        save_pop_gen_results(genotype, res, eq_z_i, eq_z_m, zhw, w_o_i_bal, st)


        # ESTIMATE MODEL PARAMETERS VALUES & GENERATE POPULATION GENETIC STATISTICS ON SIMULATED (BOOTSTRAP) DATA:
        if run_bootstrapping:
            print("bootstrapping has begun")
            bsp = bootstrap_parameters(genotype, b_res, labels, pgm, z_steps(n), z_mu, n, a, b, g, d, e, i, l,
                                       ipool, dboot, iboot, log_ratio, duration, replicates, min_n, max_n)
        sim_res_csv = (b_res + f'/{genotype}_dist{dboot}_intra{iboot}_org{replicates}_sim_results.csv')
        if os.path.exists(sim_res_csv):
            sim_res = pd.read_csv(sim_res_csv)
            eq_z_m = pd.read_csv(res + f'/{genotype}_freq_at_multilevel_balancing.csv')
            eq_z_i = pd.read_csv(res + f'/{genotype}_freq_at_intra_org_balancing.csv')
            zhw = pd.read_csv(res + f'/{genotype}_freq_at_half_fitness.csv')
            w_o_i_bal = pd.read_csv(res + f'/{genotype}_organismal_fitness_at_intra_org_balancing_freq.csv')
            st = pd.read_csv(res + f'/{genotype}_pop_stats_on_maintenance_of_heteroplasmic_state.csv')
            for p in range(len(sim_res)):
                sim = sim_res.iloc[p]
                if p >= 1:
                    n, a, b = int(sim['N']), float(sim['alpha']), float(sim['beta'])
                    g, d, e = float(sim['gamma']), float(sim['delta']), float(sim['epsilon'])
                    popgen = population_genetics_model(z_steps(n), n, a, b, g, d, e)
                    pgs = pop_gen_stats(popgen, z_mu, z_steps(n), n, a, b, g, d, e, gens)
                    inv = invasion_timing(popgen[0], pgs[0][1], pop_sizes)
                    w, loss = popgen[0], pgs[0][1]
                    z_half_w = z_at_w(0.5, a, b, np.array([z_mu[pgs[3]]]))
                    eq_z_i = merge(eq_z_i, pd.DataFrame([[f'bootstrap {p}', z_mu[pgs[1]]]], columns=c_zin))
                    eq_z_m = merge(eq_z_m, pd.DataFrame([[f'bootstrap {p}', z_mu[pgs[2]]]], columns=c_zml))
                    zhw = merge(zhw, pd.DataFrame([[f'bootstrap {p}', z_half_w]], columns=c_zhw))
                    w_o_i_bal = merge(w_o_i_bal, pd.DataFrame([[f'bootstrap {p}', pgs[4]]], columns=w_ibl))
                    st[f'bootstrap {p}'] = np.concatenate((np.array([popgen[0], n, pgs[0][1]] + inv[1]), pgs[0][0]))
            functions_of_parameters(genotype, sim_res, ['data'] + labels, res, z_200, z_mu, "bootstrap")


        # SAVE POPULATION-GENETIC STATISTICAL RESULTS TO OUTPUT FILES:
        save_pop_gen_results(genotype, res, eq_z_i, eq_z_m, zhw, w_o_i_bal, st)
