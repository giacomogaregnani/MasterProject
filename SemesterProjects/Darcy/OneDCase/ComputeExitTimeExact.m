function tau = ComputeExitTimeExact(X0, V, g, Bounds, BoundCond)
% Computes the EXACT EXIT TIME in the 1D case in case f = -V'. 

l = Bounds(1);
r = Bounds(2);

psi = @(x) 2 * ( V(l) - V(x) ) ./ (g(x).^2); 

FirstTermArg = @(x,y) exp(-psi(x)) .* exp(psi(y)) ./ (g(y).^2);
LimSupInt = @(x) x;
FirstTerm = integral2(FirstTermArg,l,X0,l,LimSupInt,'AbsTol',1e-8,'RelTol',1e-8);

SecondTermArg = @(y) exp(-psi(y));
SecondTerm = integral(SecondTermArg,l,X0,'AbsTol',1e-8,'RelTol',1e-8);

if BoundCond(2) == 0
    c1 = 2 * integral2(FirstTermArg,l,r,l,LimSupInt) / integral(SecondTermArg,l,r);
elseif BoundCond(2) == 1
    c1Arg = @(z) exp(psi(z)) ./ (g(z).^2);
    c1 = 2 * integral(c1Arg,l,r);
end

tau = -2 * FirstTerm + c1 * SecondTerm;

end