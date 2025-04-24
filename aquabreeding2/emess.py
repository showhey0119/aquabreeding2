'''
A module for error messages
'''

import sys


def check_model(model):
    '''
    Check the input demographic model for coalescent simulaitn
    '''
    if model not in ['WF', 'SP', 'buri']:
        sys.exit('model should be \'WF\' \'buri\' or \'SP\'')
# check_model


def check_method(method):
    '''
    Check a method for estimating breeding value

    Args:
        methos (str): method
    '''
    if method not in ['BLUP', 'GBLUP', 'GP', 'no']:
        sys.exit('method should be BLUP, GBLUP, GP or no')
# check_method


def check_tuple(obj, tag, l_tuple):
    '''
    Check the arguments of AquaBreeding class

    Check if an argument is tuple and
    if the length of the tuple is correct

    Args:
        obj (unknown): A focal variable
        tag (str): Name of the variable
        l_tuple (int): Correct tuple length
    '''
    if not isinstance(obj, (tuple, list)):
        sys.exit(f'{tag} should be tuple/list')
    if len(obj) != l_tuple:
        if tag in ('founder_size',  'progeny_size'):
            e_mess = '(No. males, No. females)'
        if tag == 'chrom':
            e_mess = '(Chrom no., chrom len (bp), male cM/Mb, female cM/Mb)'
        if tag in ('n_female', 'n_male'):
            e_mess = 'Length should be equal to be n_pop'
        sys.exit(f'Length of {tag} is incorrect\n{e_mess}')
# check_tuple


def check_founder(founder_size, n_population):
    '''
    Check the founder_size arugument
    '''
    l_founder = len(founder_size)
    if l_founder != n_population:
        sys.exit('Length of founder_size should be equal to n_population')
    for f_s in founder_size:
        check_tuple(f_s, 'founder_size', 2)
# check_founder


def main():
    '''
    main
    '''
    print('A module of aquabreeding')
# main


if __name__ == '__main__':
    main()
