h = 0.1;
x = [0, h, 2*h];
phi = @(xx) interp1(x, [0, 1, 0], xx);

figure
hold on
plot(xx, phi(xx))

for i = 1 : 1000
X = x + 0.5 * h * (rand(1, 3)-0.5);

Phi = @(xx) interp1(X, [0, 1, 0], xx);

xx = linspace(min(x(1), X(1)), max(x(3), X(3)), 1000);
plot(xx, phi(xx) - Phi(xx), 'color', 0.6 * ones(1, 3)) 
end