function mesh = renumber(mesh)
% changes node numbering to cover 1:N and removes unused nodes
[i, ~, mesh.elem(:)] = unique(mesh.elem(:),'legacy');
mesh.node = mesh.node(i,:);
end

