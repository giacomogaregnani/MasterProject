function [W, lambda] = plotDensity(filename, tol, zoom)

[W, lambda, x] = computeDensity(filename, tol);

figure
if size(x, 2) == 4
    if nargin == 3
        hold on
        xMin = zoom(1); xMax = zoom(2);
        yMin = zoom(3); yMax = zoom(4);
        for i = 1 : size(x, 1)
            if x(i, 1) < xMax && x(i, 1) >= xMin ...
                    && x(i, 2) < yMax && x(i, 2) >= yMin
                X = x(i, 1:2);
                R = x(i, 3:4);
                fill([X(1); X(1)+R(1); X(1)+R(1); X(1)     ], ...
                    [X(2); X(2);      X(2)+R(2); X(2)+R(2)], ...
                    W(i), 'edgecolor', 'none');
            end
        end    
        axis([zoom]);
        box on
        colorbar
    else
        hold on
        for i = 1 : size(x, 1)
            X = x(i, 1:2);
            R = x(i, 3:4);
            fill([X(1); X(1)+R(1); X(1)+R(1); X(1)     ], ...
                [X(2); X(2);      X(2)+R(2); X(2)+R(2)], ...
                W(i), 'edgecolor', 'none');
        end
    end
else
    ver = [1 1 0;
        0 1 0;
        0 1 1;
        1 1 1;
        0 0 1;
        1 0 1;
        1 0 0;
        0 0 0];
    fac = [1 2 3 4;
        4 3 5 6;
        6 7 8 5;
        1 2 8 7;
        6 7 1 4;
        2 3 5 8];
    axis equal
    for i = 1 : size(x, 1)
        cube = [ver(:,1)*x(i, 4)+x(i, 1),ver(:,2)*x(i, 5)+x(i, 2),ver(:,3)*x(i, 6)+x(i, 3)];
        p = patch('Faces',fac,'Vertices',cube,'FaceVertexCData',W(i),'EdgeColor', 'none');
        p.FaceColor = 'flat';
    end
    view(57, 15)
end

colorbar

return