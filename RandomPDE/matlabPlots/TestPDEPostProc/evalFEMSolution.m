function y = evalFEMSolution(x, u, xEval)

y = interp1(x, u, xEval);
y(isnan(y)) = 0;