function plotBoxCollectionColor(x, Col, transp)

figure
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

linCol = 0.15 + 0.85 * (Col - min(Col)) / (max(Col) - min(Col));
for i = 1 : size(x, 1)
    cube = [ver(:,1)*x(i, 4)+x(i, 1),ver(:,2)*x(i, 5)+x(i, 2),ver(:,3)*x(i, 6)+x(i, 3)];
    if transp
        p = patch('Faces',fac,'Vertices',cube,'FaceVertexCData',Col(i), ...
            'EdgeColor', 'none', 'FaceAlpha', linCol(i));
    else
        p = patch('Faces',fac,'Vertices',cube,'FaceVertexCData',Col(i), ...
            'EdgeColor', 'none');
    end
    p.FaceColor = 'flat';
end
view(57, 15)