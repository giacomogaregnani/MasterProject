function kappa = buildField_Andrea(theta)

% kappa = @(x) 1 * (x >= 0) .* (x <= 0.1);
% for i = 2 : 10
%     coeff = exp(theta(i - 1));
%     kappa = @(x) kappa(x) + coeff .* (x > (i-1)/10) .* (x <= i/10);   
% end

if length(theta) == 1
    kappa = @(x) exp(theta) * x.^0;
end

if length(theta) == 2
    kappa = @(x) exp(theta(1)) * (x < 0.5) + exp(theta(2)) * (x >= 0.5);
end

if length(theta) == 4
   kappa = @(x) 0;
   for i = 1 : 4
       kappa = @(x) kappa(x) + exp(theta(i)) .* (x >= (i-1)/4) .* (x < i / 4);
   end
end

if length(theta) == 9
    kappa = @(x) exp(theta(1)) .* (x == 0);
    for i = 2 : 10
       kappa = @(x) kappa(x) + exp(theta(i-1)) * (x <= i / 10) .* (x > (i-1) / 10); 
    end
end

return
