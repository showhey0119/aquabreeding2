'''
A module for phenotype
'''

import numpy as np
import aquabreeding2 as aq


def snp_index(n_snp, gblup):
    '''
    Randomly select neutral SNPs and SNPs with effects

    Args:
        n_snp (int): No. causal SNPs
        gblup (int): No. neutral SNPs

    Returns:
        - ndarray: Index of neutral SNPs
        - ndarray: Index of SNPs with effects
    '''
    if gblup is None:
        return None, np.arange(n_snp, dtype=np.int32)
    tmp_index = np.arange(n_snp + gblup, dtype=np.int32)
    np.random.shuffle(tmp_index)
    index_neu = tmp_index[:gblup]
    index_eff = tmp_index[gblup:]
    index_neu.sort()
    index_eff.sort()
    return index_neu, index_eff
# snp_index


def calculate_p_1p(gen_mat, n_hap):
    '''
    Calculate E[p(1 - p)] from a genotype matrix

    Args:
        gen_mat (numpy.ndarray): Genotype matrix (no. samples x np. loci)
        n_hap (int): No. haplotypes

    Returns:
        float: E[p(1 - p)]
    '''
    def cal_freq(arr1d, n_hap):
        '''
        Calculate p(1 - p)

        Args:
            arr1d (numpy.ndarray): List of genotypes
            n_hap (int): No. haplotypes

        Return:
            float: p(1 - p)
        '''
        freq_1 = np.sum(arr1d) / n_hap
        return freq_1 * (1.0 - freq_1)
    return np.mean(np.apply_along_axis(cal_freq, 0, gen_mat, n_hap))
# calculate_p_1p


def effect_size_variance(v_g, gen_eff):
    '''
    Calculate the variance of effect size

    Args:
        v_g (float): Additive genetic variance
        gen_eff (numpy.ndarray): Genotype matrix of the causal SNPs

    Returns:
        float: Variance of effect size
    '''
    n_gen, c_gen = np.shape(gen_eff)
    # E[p(1-p)]
    ep_1p = calculate_p_1p(gen_eff, 2*n_gen)
    # Variance of effect size
    return v_g/c_gen/ep_1p/2.0
# effect_size_variance


def variance_effect(var_p, h2_p, pro_gen):
    '''
    Calculate the variance of effect size

    Args:
        var_p (float): Variance of phenotype
        h2_p (float): Heritability
        pro_gen (numpy.ndarray): Genotype matrix of causal SNPs

    Returns:
        - float: Additive genetic variance
        - float: Residual variance
        - float: Variance of effect size
    '''
    # Additive/residual variance
    v_g = h2_p * var_p
    v_e = (1.0 - h2_p) * var_p
    # Variance of effect size
    v_s = effect_size_variance(v_g, pro_gen)
    return v_g, v_e, v_s
# variance_effect


def genmat_effect_all(pro_inf, index_eff):
    '''
    Combine genotype matrix with effects of all populations

    Args:
        pro_inf (list): List of progeny populations
        index_eff (numpy.ndarray): Index of SNPs with effects

    Returns:
        numpy.ndarray: combined genotype matrix with effects
    '''
    res_mat = None
    for p_inf in pro_inf:
        if res_mat is None:
            res_mat = p_inf.gen_mat.copy()
        else:
            res_mat = np.concatenate([res_mat, p_inf.gen_mat], axis=0)
    return res_mat[:, index_eff]
# genmat_effect_all


class PhenotypeInfo:
    '''
    Class for phenotype information

    Phenotypic values,  true/estimated breeding values,
    variance components are stored

    Args:
        n_population (int): No. breeding populations
        mean_phenotype (float): Mean phenotype
        var_phenotype (float): Variance of phenotype
        heritatility (float): Heritability

    Attributes:
        n_population (int): No. breeding populations
        mean_pv (float): Mean phenotype
        var_p (float): Variance of phenotype
        h2_p (float): Heritability
        v_g (float): True additive genetic variance
        v_e (float): True residual variance
        v_s (float): Variance of effect_size
        effect_size (numpy.ndarray): Effect size
        pheno_v (numpy.ndarray): Phenotypic values
        true_bv (numpy.ndarray): True breeding value
        hat_bv (numpy.ndarray): Estimated breeding value
        hat_beta (numpy.ndarray): Estimated fixed effects
        hat_vg (float): Estimated additvie genetic variance
        hat_ve (float): Estimated residual variance
        index_neu (numpy.ndarray): Index of neutral SNPs
        index_eff (numpy.ndarray): Index of causal SNPs
        _first_gen (bool): Check if the fist generation or not
    '''
    def __init__(self, n_population, mean_phenotype, var_phenotype,
                 heritability):
        '''
        Constructor
        '''
        self.n_population = n_population
        self.mean_pv = [mean_phenotype, mean_phenotype]
        self.var_p = var_phenotype
        self.h2_p = heritability
        self.v_g = None
        self.v_e = None
        self.v_s = None
        self.effect_size = None
        self.pheno_v = [None] * self.n_population
        self.true_bv = [None] * self.n_population
        self.hat_bv = [None] * self.n_population
        self.hat_beta = [None] * self.n_population
        self.hat_vg = [None] * self.n_population
        self.hat_ve = [None] * self.n_population
        # Index of netural and causal SNPs
        self.index_neu = None
        self.index_eff = None
        # to normalize true breeding value as zero
        self._first_gen = [True] * self.n_population
    # __init__

    def calculate_bv(self, method, par_inf, pro_inf, n_snp, gblup, x_i,
                     train_gen, train_phe):
        '''
        Calculate phenotype and breeding value

        Args:
            method (str): 'BLUP', 'GBLUP', or 'no'
            par_inf (PopulationInfo): Founder population
            pro_inf (list): List of PpulationInfo for progenies
            n_snp (int): No. causal SNPs
            gblup (int): No. neutral SNPs
            x_i (int): index of progeny populations
            train_gen (numpy.ndarray): genotype matrix of training population
            train_phe (numpy.ndarray): phenotypes of training population

        Note:
            The first half of self.pheno_v contains phenotypes of females, and
            seconf half of them are phenotypes of males
        '''
        # Set index of neutral and causal SNPs
        if self.index_eff is None:
            self.index_neu, self.index_eff = snp_index(n_snp, gblup)
        # Set effect size
        if self.effect_size is None:
            gen_eff_all = genmat_effect_all(pro_inf, self.index_eff)
            self.v_g, self.v_e, self.v_s = variance_effect(self.var_p,
                                                           self.h2_p,
                                                           gen_eff_all)
            self.effect_size = np.random.normal(loc=0.0,
                                                scale=np.sqrt(self.v_s),
                                                size=n_snp)
        # Genotypte matrix of causal SNPs
        gen_eff = pro_inf[x_i].gen_mat[:, self.index_eff]
        n_progeny = pro_inf[x_i].n_f + pro_inf[x_i].n_m
        # true bv
        self.true_bv[x_i] = gen_eff @ self.effect_size.T
        if self._first_gen[x_i]:
            self.mean_pv[x_i] -= np.mean(self.true_bv[x_i])
            self._first_gen[x_i] = False
        # Phenotype
        rand_norm = np.random.normal(loc=0.0, scale=np.sqrt(self.v_e),
                                     size=n_progeny)
        self.pheno_v[x_i] = self.mean_pv[x_i] + self.true_bv[x_i] + rand_norm
        # Numerator relationship matrix
        aq.nrm_cpp(par_inf, pro_inf[x_i])
        # BLUP
        if method == 'BLUP':
            # Breeding value estimation
            aq.bv_estimation(self, pro_inf[x_i].a_mat, x_i)
        # GBLUP
        elif method == 'GBLUP':
            # Genomic relationship matrix
            aq.convert_gmatrix(pro_inf[x_i],
                               pro_inf[x_i].gen_mat[:, self.index_neu],
                               0)
            # Breeding value estimation
            aq.bv_estimation(self, pro_inf[x_i].g_mat, x_i)
        # genomic prediction by GBLUP
        elif method == 'GP':
            # Genomic relationship matrix
            aq.convert_gmatrix(pro_inf[x_i],
                               pro_inf[x_i].gen_mat[:, self.index_neu],
                               0)
            # G matrix for genomic prediction
            gen_all = np.concatenate([train_gen[:, self.index_neu],
                                     pro_inf[x_i].gen_mat[:, self.index_neu]],
                                     axis=0)
            aq.convert_gmatrix(pro_inf[x_i], gen_all, 1)
            aq.bv_estimation2(self, pro_inf[x_i].g_mat2, x_i, train_phe)
    # calclaate_phenotype
# PhenotypeInfo


def main():
    '''
    main
    '''
    print('A module for phenotype')
# main


if __name__ == '__main__':
    main()
