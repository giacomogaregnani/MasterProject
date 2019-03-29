function FEMplot2D(MESH, U, cont)

figure


if nargin == 1
    pdeplot(MESH.vertices,[], MESH.elements, 'mesh','on');
    view(2)
    
elseif nargin == 2
    pdeplot(MESH.vertices,[], MESH.elements,'xydata',U(1:MESH.numVertices),'xystyle','interp',...
        'Contour','on','colorbar','on','mesh','off');
    lighting phong
    view(2)
    colormap parula
    
elseif nargin == 3
   pdeplot(MESH.vertices,[], MESH.elements,'xydata',U(1:MESH.numVertices),'xystyle','interp',...
        'Contour', cont, 'colorbar','on','mesh','off');
    lighting phong
    view(2)
    colormap parula
end

xlim([0,1])
ylim([0,1])


end