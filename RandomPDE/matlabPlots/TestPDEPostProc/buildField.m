function kappa = buildField(theta)

kappa = @(x) (x >= 0) .* (x <= 0.1);
for i = 1 : 9
    coeff = theta(i);
    kappa = @(x) kappa(x) + coeff .* (x > i/10) .* (x <= (i+1)/10);   
end

return
