def stf_ricker_int(t,f):
    """Modified function from seismo-live.org. Computes the
    Integral of the Ricker Wavelet
    """

    # delta 
    dt = t[1] - t[0]

    # Origin shift index
    t0 = 1/f / dt

    # Ricker factor
    a_ricker = 4 * f

    s = -2 * a_ricker * (t-t[int(t0)]) * np.exp(-(a_ricker * ((t-t[int(t0)]))) ** 2) 

    return s