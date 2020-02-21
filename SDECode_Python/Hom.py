import numpy as np


def compute_homogeneous(sigma, alpha, period, p):
    # Define the discretization
    n_integral = 10000
    discr = np.linspace(0.0, period, n_integral+1)
    h = period / n_integral

    # Integrate the fluctuating potential
    p_eval = p(discr)
    integrandZHat = np.exp(p_eval / sigma)
    integrandZ = np.exp(-p_eval / sigma)
    ZHat = np.trapz(integrandZHat, dx = h)
    Z = np.trapz(integrandZ, dx = h)
    k = period * period / (Z * ZHat)

    # Compute the homogenized coefficients
    hom_param = np.zeros(np.size(alpha)+1)
    hom_param[0:-1] = k * alpha
    hom_param[-1] = k * sigma
    return hom_param
