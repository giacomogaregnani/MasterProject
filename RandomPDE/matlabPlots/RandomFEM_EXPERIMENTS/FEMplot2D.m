function FEMplot2D(MESH, U, cont, lev, nContF, nContF2)

% figure

if nargin == 1
    pdeplot(MESH.vertices,[], MESH.elements, 'mesh','on');
    view(2)
    box off
    
elseif nargin == 2
    pdeplot(MESH.vertices,[], MESH.elements,'xydata',U(1:MESH.numVertices),'xystyle','interp',...
        'Contour','on','colorbar','on','mesh','off');
    lighting phong
    view(2)
    colormap parula
    
elseif nargin == 3
    pdeplot(MESH.vertices,[], MESH.elements,'xydata',U(1:MESH.numVertices),'xystyle','interp',...
        'Contour', cont, 'colorbar','on', 'mesh', 'off');
    lighting phong
    view(2)
    colormap parula
    
elseif nargin == 4
    pdeplot(MESH.vertices, [], MESH.elements, 'xydata', U(1:MESH.numVertices), 'xystyle', 'flat',...
        'Contour', cont, 'Levels', lev, 'colorbar', 'on', 'mesh', 'off');
    lighting phong
    view(2)
    colormap parula
    
elseif nargin >= 5 && nargin <= 6
    uInt = pdeInterpolant(MESH.vertices, MESH.elements, U);
    x = linspace(0, 1, nContF);
    [XX, YY] = meshgrid(x, x);
    UU = uInt.evaluate([XX(:)'; YY(:)']);
    UU = reshape(UU, nContF, nContF);
    contourf(XX, YY, UU, nContF2);
    figure
    surf(XX, YY, UU, 'edgecolor', 'none')

end

xlim([0,1])
ylim([0,1])