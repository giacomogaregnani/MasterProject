function placeLegendExtreme(fig, margin, loc)
%% placeLegendExtreme(fig, margin, loc)
% fig must be a handle to an existing figure
% margin is the desired distance in pixels from the axes boundary
% loc is the desired location of the figure.
% valid values for loc = {'NW', 'NE', 'SW', 'SE', 'N', 'S', 'W', 'E'}

ax = get(fig, 'currentaxes');
leg = get(ax, 'legend');
ax.Units = 'pixels';
leg.Units = 'pixels';
legendPosition = leg.Position;

switch loc
    case 'NW'
        set(leg, 'Position', [ax.Position(1)+margin, ax.Position(2)+ax.Position(4)-legendPosition(4)-margin, ...
            legendPosition(3), legendPosition(4)]);
    case 'NE'
        set(leg, 'Position', [ax.Position(1)+ax.Position(3)-legendPosition(3)-margin, ax.Position(2)+ax.Position(4)-legendPosition(4)-margin, ...
            legendPosition(3), legendPosition(4)]);
    case 'SW'
        set(leg, 'Position', [ax.Position(1)+margin, ax.Position(2)+margin, ...
            legendPosition(3), legendPosition(4)]);
    case 'SE'
        set(leg, 'Position', [ax.Position(1)+ax.Position(3)-legendPosition(3)-margin, ax.Position(2)+margin, ...
            legendPosition(3), legendPosition(4)]);
    case 'N'
        set(leg, 'Position', [ax.Position(1)+ax.Position(3)/2-legendPosition(3)/2-margin, ax.Position(2)+ax.Position(4)-legendPosition(4)-margin, ...
            legendPosition(3), legendPosition(4)]);
    case 'S'
        set(leg, 'Position', [ax.Position(1)+ax.Position(3)/2-legendPosition(3)/2-margin, ax.Position(2)+margin, ...
            legendPosition(3), legendPosition(4)]);
    case 'E'
        set(leg, 'Position', [ax.Position(1)+ax.Position(3)-legendPosition(3)-margin, ax.Position(2)+ax.Position(4)/2-legendPosition(4)/2, ...
            legendPosition(3), legendPosition(4)]);
    case 'W'
        set(leg, 'Position', [ax.Position(1)+margin, ax.Position(2)+ax.Position(4)/2-legendPosition(4)/2, ...
            legendPosition(3), legendPosition(4)]);
    otherwise
        fprintf('valid values for loc = {NW, NE, SW, SE, N, S, W, E}\n')
end
return