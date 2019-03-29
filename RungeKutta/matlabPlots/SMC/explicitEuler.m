function y = explicitEuler(f, N, y, h)

for i = 1 : N
    y = y + h * f(h, y);
end

end