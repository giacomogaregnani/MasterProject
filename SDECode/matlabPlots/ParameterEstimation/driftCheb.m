function eval = driftCheb(x, alpha, nParam)

pol = cell(nParam+1, 1);
pol{1} = @(x) 1 * x.^0;
pol{2} = @(x) x;
for i = 3 : nParam+1
   pol{i} = @(x) 2 * x .* pol{i-1}(x) - pol{i-2}(x);
end

eval = zeros(size(x));
for i = 1 : nParam
    eval = eval + alpha(i) * pol{i+1}(x);
end