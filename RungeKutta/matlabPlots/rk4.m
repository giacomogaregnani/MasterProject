function y = rk4(fun, h, t, y)

k1 = fun(t, y);
k2 = fun(t + h/2, y + h/2 * k1);
k3 = fun(t + h/2, y + h/2 * k2);
k4 = fun(t + h, y + h * k3);
y = y + h/6 * k1 + h/3 * k2 + h/3 * k3 + h/6 * k4;

end