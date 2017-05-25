clc; clear; close all;

% Obtain data
suppString = {'q15', 'q2', 'q25', 'q3'};
q = [1.5 2 2.5 3];
p = 2;

mark = 'o+<>';

figure

for i = 1 : length(q)
    results = dlmread(['error' suppString{i} '.txt']);
    h = results(:, 1);
    err = results(:, 2);
    thOrder = min(q(i) - 0.5, p);
    loglog(h, err, 'marker', mark(i));
    hold on
    leg{2 * i - 1} = ['q = ', num2str(q(i))];
    if q(i) - 0.5 <= p
        loglog(h, h.^thOrder, 'k', 'marker', mark(i));
        leg{2 * i} = ['slope = ', num2str(thOrder)];
    end
end

legend(leg, 'Location', 'SE')