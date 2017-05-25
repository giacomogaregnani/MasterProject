%% Example with lognormal

hMax = 0.1;
P = @(h, p) 1 - logncdf(hMax, log(h) - 0.5 * log(1 + h.^(2 * p - 2)), 0.5 * log(1 + h.^(2 * p - 2)));

