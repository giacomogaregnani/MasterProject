function y = explicitEulerProb(f, N, y, h)

for i = 1 : N
    H = h + h^1.5 * (2 * rand(1) - 1);
    y = y + H * f(H, y);
end

end

