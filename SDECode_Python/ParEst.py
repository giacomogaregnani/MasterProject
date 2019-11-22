import numpy as np
import scipy.stats as stats
from scipy.signal import convolve
from scipy.special import gamma

def moving_average(a, n):
    ret = np.zeros(np.size(a))
    for i in range(0, n):
        ret[i] = ret[i-1] + (a[i] - a[0]) / n
    for i in range(n, np.size(ret)):
        ret[i] = ret[i-1] + (a[i] - a[i-n]) / float(n)
    return ret


def filter_trajectory(a, delta, T, beta=1):
    n = np.size(a)
    h = T / (n-1)
    t_vec = np.arange(n) * h
    if beta == 1:
        Cb = 1.0
    else:
        Cb = 1.0 / gamma((beta+1.0)/beta)
    k = Cb * np.exp(-1.0 * (t_vec ** beta) / delta) / (delta ** (1.0 / beta))
    ret = convolve(a, k, mode='full') * h
    ret=ret[0:n]
    return ret, k[0]


# This is just a naive implementation used for checking filter_trajectory
def kernel(t, s, delta):
    return 1.0 / delta * np.exp(-1.0 * (t - s) / delta)


def filter_trajectory_2(a, delta, T,beta=1):
    n = np.size(a)
    h = T / (n-1)
    t_vec = np.arange(n) * h
    ret = np.zeros(n)
    for i in range(0, n):
        ret[i] = h * np.sum(np.multiply(a[0:i], kernel(t_vec[i], t_vec[0:i], delta)))
    return ret


def get_density(x):
    kde = stats.gaussian_kde(x)
    x_sorted = np.sort(x)
    density = kde.evaluate(x_sorted)
    return x_sorted, density


class ParEst:
    def __init__(self, x, grad_v, h, is_vect=False):
        self.x = x
        self.grad_v = grad_v
        self.h = h
        self.x_diff = np.diff(x)
        self.n = np.size(x)
        self.is_vect = is_vect

    def phiPrime(self, y):
        return self.period / self.zHat * np.exp(self.p(y) / self.sigma) - 1.0

    def phi(self, y):
        nGrid = 10
        grid = np.linspace(0.0, y, nGrid+1)
        dxIntegral = y / nGrid
        func = np.exp(self.p(grid) / self.sigma)
        return self.period / self.zHat * np.trapz(func, dx=dxIntegral) - y

    def dtk(self, delta, t, s):
        return -1.0 / (delta ** 2) * np.exp(-1.0 * (t - s) / delta)

    def formula_verification(self, sigma, alpha, p, eps, delta, data, dw):
        self.p = p
        self.sigma = sigma
        self.period = 2.0*np.pi

        # Compute ZHat for the definition of Phi and its derivative
        nZIntegral = 10000
        zIntX = np.linspace(0.0, self.period, nZIntegral+1)
        zIntdx = self.period / nZIntegral
        p_eval = p(zIntX)
        integrandZHat = np.exp(p_eval / sigma)
        self.zHat = np.trapz(integrandZHat, dx=zIntdx)

        # Compute the stochastic integral
        nChunks = (self.n-1) / delta
        integralV = np.zeros(nChunks)
        firstIntegral = np.zeros(nChunks)
        secondIntegral = np.zeros(nChunks)
        thirdIntegral = np.zeros(nChunks)
        idx = 0
        for i in range(0, nChunks):
            for j in range(0, delta):
                evalPhiPrime = self.phiPrime(data[idx+j]/eps)
                firstIntegral[i] += self.grad_v(data[idx+j]) * (1.0 + evalPhiPrime) * self.h
                secondIntegral[i] += (1.0 + evalPhiPrime) * dw[idx+j]
                integralV[i] += self.grad_v(self.x[idx+j]) * self.h
            firstIntegral[i] *= -alpha
            secondIntegral[i] *= np.sqrt(2.0 * sigma)
            thirdIntegral[i] = -eps * (self.phi(data[idx+delta]/eps) - self.phi(data[idx]/eps))
            idx = idx + delta

        num1 = -1.0 / (delta * self.h) * np.sum(np.multiply(firstIntegral, integralV))
        num2 = -1.0 / (delta * self.h) * np.sum(np.multiply(secondIntegral, integralV))
        num3 = -1.0 / (delta * self.h) * np.sum(np.multiply(thirdIntegral, integralV))
        term1 = num1 / (self.h * self.den)
        term2 = num2 / (self.h * self.den)
        term3 = num3 / (self.h * self.den)
        print(term1, term2, term3, term1 + term2 + term3)

        return term1

    def drift(self, strat=False, Sigma=1, lapl=None):
        grad_v_eval = np.zeros(self.n-1)
        if self.is_vect:
            grad_v_eval = self.grad_v(self.x[0:-1])
        else:
            for i in range(0, self.n - 1):
                grad_v_eval[i] = self.grad_v(self.x[i])
        num_summand = np.multiply(self.x_diff, grad_v_eval)
        num = np.sum(num_summand)
        # Convert from Stratonovich to Ito
        if strat:
            correction = Sigma * self.h * np.sum(lapl(self.x[0:-1]))
            num -= correction
        self.den = np.sum(np.square(grad_v_eval))
        return -1.0 * num / (self.h * self.den)

    def drift_alternative(self, lapl, delta=1.0, Sigma=0):
        if Sigma==0:
            Sigma = self.diffusion() / delta
        grad_v_eval = self.grad_v(self.x[0:-1])
        den = np.sum(np.square(grad_v_eval)) * self.h
        num = Sigma * np.sum(lapl(self.x[0:-1])) * self.h
        return num / den, Sigma

    def diffusion(self):
        return np.sum(np.square(self.x_diff)) / (2.0 * self.h * self.n)

    def driftBayesian(self, a_pr, var_pr):
        grad_v_eval = np.zeros(self.n - 1)
        if self.is_vect:
            grad_v_eval = self.grad_v(self.x[0:-1])
        else:
            for i in range(0, self.n - 1):
                grad_v_eval[i] = self.grad_v(self.x[i])
        int_num = np.sum(np.multiply(self.x_diff, grad_v_eval))
        int_denom = self.h * np.sum(np.square(grad_v_eval))
        den = 1.0 + var_pr**2.0 * int_denom
        num = a_pr - var_pr**2.0 * int_num
        a_post = num / den
        var_post = 1.0 / (var_pr**(-2.0) + int_denom)
        return a_post, var_post
