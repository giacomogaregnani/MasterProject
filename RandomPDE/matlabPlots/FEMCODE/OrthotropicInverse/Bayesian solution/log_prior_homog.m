function p = log_prior_homog(param)

if param(2) < 0 || param(2) > 0.5
    p = -inf;
else
    p = -0.5 * param(1)^2;
end

return