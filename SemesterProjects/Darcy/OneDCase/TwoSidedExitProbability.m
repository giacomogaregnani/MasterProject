function p = TwoSidedExitProbability(x0,xh,sigma,h,Bounds)
% Computes the EXIT PROBABILITY between two consecutive timesteps

p = max(exp(-2 * (Bounds(2) - xh) * (Bounds(2) - x0) / (sigma^2 * h)), ...
    exp(-2 * (xh - Bounds(1)) * (x0 - Bounds(1))/ (sigma^2 * h)));

