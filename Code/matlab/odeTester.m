roberF = @(t, y) [-0.04 * y(1) + 10e4 * y(2) * y(3); ...
    0.04 * y(1) - 10e4 * y(2) * y(3) - 3e7 * y(2)^2; ...
    3e7 * y(2)^2];
options = odeset('InitialStep', 1e-10);
[t, y] = ode23s(roberF, [0, 10e11], [1; 0;0], options);


for i = 1 : 3
    figure
    semilogx(t, y(:, i))
end

figure
semilogy(t(2 : end) - t(1 : end - 1))
