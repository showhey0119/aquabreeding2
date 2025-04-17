'''
Module for simulating aquaculture breeding

Aquabreeding2 allows simulating aquaculture breeding.

A founder population is generated together with coalesent simulation,
implemented in msprime.

Progenies are produced following some implemented mating schemes or
user's own setting.

The individuals with large phenotypic/breeding values are selected.
Mass selection, within-family selection, and family selection are implemented.
The selected ones are used as parentsin the next generation.

Phenotype, true/estimated breeding value, inbreeding coefficient,
and variance components are output.

Note:
    This is a beta version.

Todo:
    * Dominance and epistasis
    * Multiple phenotypes
    * Founders from a population with exponential growth or bottleneck
'''

import sys
import numpy as np
import aquabreeding2 as aq


class AquaBreeding:
    '''
    Class for aquaculture breeding

    Args:
        founder_size (tuple): (Nos. females and  males) in each population
        n_population (int): No. breeding population
        chrom (tuple): Chrom num, chrom len (bp), female cM/Mb, male cM/Mb
        mean_phenotype (float): Mean phenotype
        var_phenotype (float): Variance of phenotype
        heritability (float): Heritability
        n_wild (int): No. individuals in a extra natural population,
                     default None

    Attributes:
        chrom (tuple): Chrom num, chrom len (bp), female cM/Mb, male cM/Mb
        par_inf (PopulationInfo): Founder/parental population
        pro_inf (PopulationInfo): Progeny population
        n_snp (int): No. causal SNPs
        phe_inf (PhenotypeInfo): Phenotype information
        cross_inf (numpy.ndarray): Index pairs of female and male parents
        gblup (int): No. neutral SNPs
    '''
    def __init__(self, founder_size, n_population, chrom, mean_phenotype,
                 var_phenotype, heritability, n_wild=None):
        '''
        constructor
        '''
        # check argument
        aq.check_founder(founder_size, n_population)
        aq.check_tuple(chrom, 'chrom', 4)
        self.founder_size = founder_size
        self.n_population = n_population
        # chromosome info
        self.chrom = chrom
        # parental population
        self.par_inf = []
        for i in range(self.n_population):
            self.par_inf.append(aq.PopulationInfo(self.founder_size[i],
                                                  self.chrom, n_wild))
            self.par_inf[i].init_founder()
        # progeny population
        self.pro_inf = [None] * self.n_population
        # Phenotype info
        self.phe_inf = aq.PhenotypeInfo(self.n_population, mean_phenotype,
                                        var_phenotype, heritability)
        # No. causal SNPs
        self.n_snp = None
        # No. neutral SNPs
        self.gblup = None
        # Mating design
        self.cross_inf = [None] * self.n_population
    # __init__

    def snp(self, model, n_snp, gblup=None, n_pop=None, fst=None,
            n_female=None, n_male=None):
        '''
        Generate SNPs

        Args:
            model (str): 'WF' or 'SP'
            n_snp (int): No. causal SNPs
            gblup (int): No. neutral SNPs
            n_pop (int): No. populations
            fst (float): Average Fst value among populations
            n_female (tuple): Nos. females in populations
            n_male (tuple): Nos. males in populations
        '''
        # check args
        if n_pop is not None:
            sys.exit('Under construction')
            check_tuple(n_female, 'n_female', n_pop)
            check_tuple(n_male, 'n_male', n_pop)
            if sum(n_female) != self.par_inf.n_f:
                sys.exit(f'Sum of n_female is not {self.par_inf.n_f}')
            if sum(n_male) != self.par_inf.n_m:
                sys.exit(f'Sum of n_male is not {self.par_inf.n_m}')
        aq.check_model(model)
        self.n_snp = n_snp
        self.gblup = gblup
        # coalescent simulation
        aq.coalescent_simulation(model, self.par_inf, self.n_snp, self.gblup,
                                 n_pop, fst, n_female, n_male)
    # snp

    def mating_design(self, design, select_size=None):
        '''
        Set mating design

        Args:
            design (unknown): If str, design should be '1x(int)' for partial
                              factorial cross or 'full' factorial mating.
                              If numpy.ndarray, design contains two columns:
                              index of female parents and index of male
                              parents.
            select_size (tuple): Set nos. selected females and males
        '''
        if select_size is None:
            select_size = [None] * self.n_population
        aq.check_tuple(design, 'design', self.n_population)
        aq.check_tuple(select_size, 'select_size', self.n_population)
        for i in range(self.n_population):
            if select_size[i] is None:
                select_size[i] = (self.par_inf[i].n_f, self.par_inf[i].n_m)
            else:
                aq.check_tuple(select_size[i], 'select_size2', 2)
            self.cross_inf[i] = aq.set_mating_design(design[i], select_size[i])
    # mating_design

    def mating(self, progeny_size):
        '''
        Mate founder/parental individuals to produce progenies

        Args:
            progeny_size (tuple): Nos. female and male progenies
        '''
        aq.check_tuple(progeny_size, 'progeny_size', self.n_population)
        for i in range(self.n_population):
            aq.check_tuple(progeny_size[i], 'progeny_size2', 2)
            if self.pro_inf[i] is None:
                self.pro_inf[i] = aq.PopulationInfo(progeny_size[i], self.chrom)
                self.pro_inf[i].init_progeny()
            else:
                self.pro_inf[i].change_size(progeny_size[i])
            aq.start_mating(self.cross_inf[i], self.par_inf[i], self.pro_inf[i])
    # mating

    def breeding_value(self, method):
        '''
        Calculate phenotype and breeding value

        Args:
            method (str): If 'BLUP', numerator relationship matrix is used
                          to estimate breeding values.  If 'GBLUP', genomic
                          relationship matrix is used.  If 'no', breeding
                          values are not estimated. Default 'BLUP'
        '''
        aq.check_tuple(method, 'method', self.n_population)
        # genotype matrix
        for i in range(self.n_population):
            self.pro_inf[i].genotype_matrix(self.n_snp, self.gblup)
        # Calculate phenotype and breeding value
        for i in range(self.n_population):
            aq.check_method(method[i])
            self.phe_inf.calculate_bv(method[i], self.par_inf[i], self.pro_inf,
                                      self.n_snp, self.gblup, i)
    # breeding_value

    def selection(self, target, method, top_prop=None, n_family=None,
                  select_size=None, max_r=None):
        '''
        Select parents of next generation

        Args:
            target (str): Selection based on breeding value ('bv'),
                          phenotypic value ('phenotype'), or random
                          ('random')
            method (str): How to select from progenies such as mass
                          selection ('mass'), within-family selection
                          ('within-family'),  family selection ('family'),
                          or selection based on A matrix ('RvalueA') or
                          G matrix ('RvalueG')
            top_prop (float): Select progenies with top X% breeding values
                              in within-family selection. Set 0.0 < top_prop
                              <= 1.0.
            n_family (int): No. families to be selected
            select_size (tulple): Number of selected founders, default: None
            max_r (float): R' values among selected individuals are less than
                           max_r

        Returns:
            int: if 0, excuted correctly. if 1, terminated irregurally
        '''
        if top_prop is None:
            top_prop = [1.0] * self.n_population
        if n_family is None:
            n_family = [None] * self.n_population
        if max_r is None:
            max_r = [None] * self.n_population
        if select_size is None:
            select_size = [None] * self.n_population
        aq.check_tuple(select_size, 'select_size', self.n_population)
        for i in range(self.n_population):
            if select_size[i] is None:
                select_size[i] = (self.par_inf[i].n_f, self.par_inf[i].n_m)
            else:
                aq.check_tuple(select_size[i], 'select_size2', 2)
            self.par_inf[i].change_size(select_size[i])
            check_r = aq.start_selection(self.par_inf[i], self.pro_inf[i],
                                         self.phe_inf, target[i], method[i],
                                         self.cross_inf[i], top_prop[i],
                                         n_family[i], select_size[i], max_r[i], i)
        return check_r
    # selection

    def get_ibd(self):
        '''
        Output inbreeding coefficient

        Returns:
            numpy.ndarray: Inbreeding coefficient
        '''
        return np.diag(self.pro_inf.a_mat) - 1.0
    # get_ibd

    def get_phenotype(self):
        '''
        Output phenotype

        Returns:
            numpy.ndarray: Phenotype
        '''
        return self.phe_inf.pheno_v
    # get_phenotype

    def get_true_bv(self):
        '''
        Output true breeding value

        Returns:
            numpy.ndarray: True breeding value
        '''
        return self.phe_inf.true_bv
    # get_true_bv

    def get_ebv(self):
        '''
        Output estimated breeding value

        Returns:
            numpy.ndarray: Estimated breeding value
        '''
        return self.phe_inf.hat_bv
    # get_ebv

    def variance_component(self):
        '''
        Outut estimated additive and residual variance

        Returns:
            tuple: Additive and residual variance
        '''
        return self.phe_inf.hat_vg, self.phe_inf.hat_ve
    # variance_component

    def wild_parent(self, num):
        '''
        Replance parents with wild individuals before mating

        Args:
            num (tuple): No. famale/male parents relpaced by wild individuals
        '''
        check_tuple(num, 'num', 2)
        self.par_inf.replace_wild(num)
    # wild_parent

    def get_mean_phenotype(self):
        '''
        Output mean phenotypes of all breeding populations
        
        Returns:
            list: mean phenotypes
        '''
        res_p = []
        for i in range(self.n_population):
            res_p.append(np.mean(self.phe_inf.pheno_v[i]))
        return res_p
    # get_mean_phenotype

    def get_mean_ibd(self):
        '''
        Output meaninbreeding coefficient of all breeding populations

        Returns:
            list: mean ibd
        '''
        res_i = []
        for i in range(self.n_population):
            res_i.append(np.mean(np.diag(self.pro_inf[i].a_mat) - 1.0))
        return res_i
    # get_mean_ibd

    def merge_population(self, target):
        '''
        Merge two breeding populations and remove second one
        Keep selected females from a population and males from
        another population

        Args:
            target (list): Index (starts with 1) of two breeding populations
        '''
        aq.merge_pop(par_inf, target)
        self.n_population -= 1
# AquaBreeding


def main():
    '''
    main
    '''
    print('A module for simulating aquacuture breeding')
# main


if __name__ == '__main__':
    main()
