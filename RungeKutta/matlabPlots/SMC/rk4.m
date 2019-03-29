function y = rk4(f, N, y, h)

for i = 1 : N
    k1 = f(h, y);
    k2 = f(h, y + h/2 * k1);
    k3 = f(h, y + h/2 * k2);
    k4 = f(h, y + h * k3);
    y = y + h/6 * k1 + h/3 * k2 + h/3 * k3 + h/6 * k4;
end

end