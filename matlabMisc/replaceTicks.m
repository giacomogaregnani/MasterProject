function replaceTicks(linOrLogX, linOrLogY)

% Generate figure and remove ticklabels
set(gca,'yticklabel',[], 'xticklabel', []) %Remove tick labels

% Reset the ytick labels in desired font
yTicks = get(gca, 'ytick');
xTicks = get(gca, 'xtick');
ax = axis; 

if linOrLogX == 'lin'
    HorizontalOffset = 0.01 * (xTicks(end) - xTicks(1));
else
    HorizontalOffset = 5 * min(xTicks) / 10;
end

for i = 1 : length(yTicks)
    if linOrLogY == 'lin'
        text(ax(1) - HorizontalOffset, yTicks(i), ['$' num2str(yTicks(i)) '$'],...
            'HorizontalAlignment','Right','interpreter', 'latex');
    else
        text(ax(1) - HorizontalOffset, yTicks(i), ['$' num2str(yTicks(i), '%10.0e\n') '$'],...
            'HorizontalAlignment','Right','interpreter', 'latex');
    end
end

% Reset ylabel
yLab = get(gca, 'ylabel');
set(gca, 'ylabel', []);

if linOrLogY == 'lin'
    yPosYLab = 0.5 * (yTicks(end) - yTicks(1));
else
    yPosYLab = 10 ^ (0.5 * (log10(yTicks(1)) - log10(yTicks(end))));
end

if linOrLogX == 'log'
    text(ax(1) - 7.6 / 5 * HorizontalOffset, yPosYLab, ['$', yLab.String, '$'], ...
        'HorizontalAlignment', 'Center', 'interpreter', 'latex', 'FontSize', 10);
else
    text(ax(1) - 12 * HorizontalOffset, yPosYLab, ['$', yLab.String, '$'], ...
        'HorizontalAlignment', 'Center', 'interpreter', 'latex', 'FontSize', 10);
end

% Reset the xtick labels in desired font
minY = min(yTicks);

if linOrLogY == 'lin'
    verticalOffset = 0.01 * (yTicks(end) - yTicks(1));
else
    verticalOffset = 5 * minY / 10;
end

for xx = 1 : length(xTicks)
    if linOrLogX == 'lin'
        text(xTicks(xx), minY - verticalOffset, ['$' num2str(xTicks(xx)) '$'],...
            'HorizontalAlignment', 'Center', 'interpreter', 'latex', 'FontSize', 10);
    else
        text(xTicks(xx), minY - verticalOffset, ['$' num2str(xTicks(xx), '%10.0e\n'), '$'],...
            'HorizontalAlignment', 'Center', 'interpreter', 'latex', 'FontSize', 10);
    end
end

% Reset xlabel
xLab = get(gca, 'xlabel');
set(gca, 'xlabel', []);

if linOrLogX == 'lin'
    xPosXLab = 0.5 * (xTicks(end) - xTicks(1));
else
    xPosXLab = 10 ^ (0.5 * (log10(xTicks(end)) - log10(xTicks(1))));
end

if linOrLogY == 'log'
    text(xPosXLab, minY - 8 / 5 * verticalOffset, ['$' xLab.String '$'], ...
        'HorizontalAlignment', 'Center', 'interpreter', 'latex', 'FontSize', 10);
else
    text(xPosXLab, minY - 2 * verticalOffset, ['$' xLab.String '$'], ...
        'HorizontalAlignment', 'Center', 'interpreter', 'latex', 'FontSize', 10);
end


