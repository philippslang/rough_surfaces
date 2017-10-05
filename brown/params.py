import collections as co


SAPARAMS = ['hrms', 'dimensions', 'hurst', 'lambda_L_over_lambda_0', 'lambda_L_over_lambda_1']
SADEFAUL = [0.01, (1.0, 1.0), 0.8, None, None]


'''Minimalistic parameter set for self affine surfaces.'''
SelfAffineParameters = co.namedtuple('SelfAffineParameters', SAPARAMS)


def self_affine_default_parameters():
    '''
    Provides some (arbitrary) defaults for a self affine surfaces.

    >>> saparams = self_affine_default_parameters()
    >>> saparams.hrms == SADEFAUL[0]
    True
    >>> saparams.dimensions
    (1.0, 1.0)
    >>> saparams.hurst
    0.8
    '''
    return SelfAffineParameters(*SADEFAUL)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
