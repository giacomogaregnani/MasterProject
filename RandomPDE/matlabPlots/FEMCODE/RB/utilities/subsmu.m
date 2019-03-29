function expr = subsmu(expr,mu)
for p=1:size(mu,2)
    expr = subs(expr,['mu' num2str(p)], mu(p));
end
expr = double(expr);
end