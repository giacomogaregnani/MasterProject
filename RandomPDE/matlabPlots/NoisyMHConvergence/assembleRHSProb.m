function F = assembleRHSProb(f, x)

F = 0.5 * f(x(2:end-1)) .* (x(3:end) - x(1:end-2))';
