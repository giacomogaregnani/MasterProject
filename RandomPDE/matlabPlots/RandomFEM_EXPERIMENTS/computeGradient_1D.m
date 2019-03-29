function g = computeGradient_1D(x, u)

g = (u(2:end) - u(1:end-1)) ./ (x(2:end) - x(1:end-1));

end