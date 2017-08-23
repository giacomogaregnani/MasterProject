function FEMplot2D(MESH, U)

figure


if nargin == 1
    pdeplot(MESH.vertices,[], MESH.elements(1:3,:), 'mesh','on');
    view(2)
    
elseif nargin == 2
pdeplot(MESH.vertices,[], MESH.elements(1:3,:),'xydata',U(1:MESH.numVertices),'xystyle','interp',...
    'zdata',U(1:MESH.numVertices),'zstyle','continuous',...
    'colorbar','on','mesh','off');
lighting phong
view(2)
colormap parula
end


end