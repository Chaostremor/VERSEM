import numpy as np

def stf_delta_function(t):
    """ This function takes in the time stepping vector and outputs
    a zero vector and the time step that is 1/100th of all time steps.

    :param t: 1D ``numpy`` array

    :rtype: delta source_time function vector

    """

    # Length of t
    nt = len(t)

    # initialize new array
    s = np.zeros(nt)

    # Find index where to :math:`t_i = t(n_t/100)
    s[int(nt/100)] = 1 

    return s