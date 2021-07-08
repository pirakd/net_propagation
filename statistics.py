import numpy as np
import scipy.stats
import scipy
from scipy.stats import norm


class StatRes:
    def __init__(self, p_value=None, directionality=None, z_score=None, name=None):
        self.p_value = p_value
        self.directionality = directionality
        self.z_score = z_score
        self.name = name


def get_stat_test_func(name):
    if name == 'empirical_mean_diff':
        return empirical_mean_diff
    elif name == 'two_sample_z_test':
        return two_sample_z_test
    elif name == 'man_whit_U_test':
        return man_whit_U_test


def empirical_mean_diff(experiment_scores, elements_scores, n_draws=5000) -> StatRes:
    empirical_mean_diff.name = 'Empirical mean diff'

    if not isinstance(elements_scores, np.ndarray):
        elements_scores = np.array(elements_scores)

    mean_scores = np.mean(experiment_scores)
    # mean_stds = [np.std(x, ddof=1) for x in experiment_scores]
    random_scores =\
        np.array([np.mean(elements_scores[np.random.randint(0, len(elements_scores), len(experiment_scores))]) for x in range(n_draws)])
    random_scores = np.array(random_scores).transpose()
    p_val, direction = get_sample_p_values(mean_scores, random_scores, two_tailed=True)
    return StatRes(p_value=p_val[0], directionality=direction[0], name=empirical_mean_diff.name)


def man_whit_U_test(experiment_scores, elements_scores) -> StatRes:
    man_whit_U_test.name = 'Mann–Whitney U test'
    p_vals = scipy.stats.mannwhitneyu(experiment_scores, elements_scores).pvalue
    direction = np.mean(experiment_scores) > np.mean(elements_scores)
    return StatRes(p_value=p_vals, directionality=direction, name=man_whit_U_test.name)


def two_sample_z_test(sample_1, sample_2, mu_diff=0) -> StatRes:

    two_sample_z_test.name = 'two_sample_z_test'
    mu_1, mu_2 = np.mean(sample_1), np.mean(sample_2)
    sd_1, sd_2 = np.std(sample_1, ddof=1), np.std(sample_2, ddof=1)
    n_1, n_2 = len(sample_1), len(sample_2)
    pooledSE = np.sqrt(sd_1**2/n_1 + sd_2**2/n_2)
    z = ((mu_1 - mu_2) - mu_diff)/pooledSE
    directionality = z>0
    pval = 2*(1 - norm.cdf(abs(z)))
    return StatRes(p_value=pval, directionality=directionality , z_score=z, name=two_sample_z_test.name)


def bh_correction(p_values):
    p_vals_rank = scipy.stats.rankdata(p_values, 'max') - 1
    p_vals_rank_ord = scipy.stats.rankdata(p_values, 'ordinal') - 1

    p_values_sorted = np.zeros_like(p_vals_rank)
    p_values_sorted[p_vals_rank_ord] = np.arange(len(p_vals_rank_ord))

    p_vals = p_values * (len(p_values) /(p_vals_rank+1))
    adj_p_vals_by_rank = p_vals[p_values_sorted]

    p_vals_ordered = np.minimum(adj_p_vals_by_rank, np.minimum.accumulate(adj_p_vals_by_rank[::-1])[::-1])
    adj_p_vals = p_vals_ordered[p_vals_rank]
    return adj_p_vals


def get_sample_p_values(original_network_scores, random_networks_scores, two_tailed=False):
    if len(random_networks_scores.shape) != 2:
        random_networks_scores = random_networks_scores[:, np.newaxis]
    if not (isinstance(original_network_scores,list) or  (isinstance(original_network_scores, np.ndarray))):
        original_network_scores = [original_network_scores]
    n_experiments = random_networks_scores.shape[0]
    sorted_scores = np.sort(random_networks_scores, axis=0)
    gene_score_rank = np.array(
        [np.searchsorted(sorted_scores[:, i], original_network_scores[i], side='left')
         for i in range(random_networks_scores.shape[1])])
    step = 1 / (n_experiments+1)
    p_values_stepes = np.arange(1, n_experiments + 2, 1) * step
    if not two_tailed:
        p_values = p_values_stepes[gene_score_rank]
        return p_values
    else:
        higher_than_half = gene_score_rank >= (n_experiments / 2)
        p_values = np.zeros(gene_score_rank.shape)
        p_values[np.logical_not(higher_than_half)] = p_values_stepes[gene_score_rank[np.logical_not(higher_than_half)]]
        p_values[higher_than_half] =\
            1-(p_values_stepes[gene_score_rank[higher_than_half] - 1])
        p_values = p_values * 2
        return p_values, higher_than_half
