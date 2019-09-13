import numpy as np


def compute_homogeneous(param, period, grad_p):
    n_integral = 10000
    discr = np.linspace(0.0, period, n_integral+1)
    h = period / n_integral
    sigma = param[-1]

    z = np.zeros(2)
    grad_p_eval = np.zeros(2)
    for i in range(0, n_integral):
        grad_p_eval[0] = grad_p(discr[i])
        grad_p_eval[1] = grad_p(discr[i+1])
        z[0] += h * (np.exp(grad_p_eval[0] / sigma) + np.exp(grad_p_eval[1] / sigma)) / 2.0
        z[1] += h * (np.exp(-grad_p_eval[0] / sigma) + np.exp(-grad_p_eval[1] / sigma)) / 2.0

    hom_param = np.zeros(2)
    k = period * period / (z[0] * z[1])
    hom_param[0] = k * param[0]
    hom_param[1] = k * param[1]

    return hom_param
