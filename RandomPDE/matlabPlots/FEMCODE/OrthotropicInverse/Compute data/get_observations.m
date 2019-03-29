function values = get_observations(mesh, sol, locations, femspace)
values = zeros(3,length(locations),2);
edge_north = @(x) (x(:,2) > mesh.box(4) - mesh.mesh_size*1e-3);
edge_west = @(x) (x(:,1) < mesh.box(1) + mesh.mesh_size*1e-3);
edge_south = @(x) (x(:,2) < mesh.box(3) + mesh.mesh_size*1e-3);
edge_east = @(x) (x(:,1) > mesh.box(2) - mesh.mesh_size*1e-3);

dofs_north = find(edge_north(mesh.node(femspace.bd_dof,:))==1);
[a, b] = sort(mesh.node(femspace.bd_dof(dofs_north),1));
sol_north = sol(femspace.bd_dof(dofs_north(b)),:);
sol_north_locations = interp1(a, sol_north, locations);
dofs_west = find(edge_west(mesh.node(femspace.bd_dof,:))==1);
[a, b] = sort(mesh.node(femspace.bd_dof(dofs_west),2));
sol_west = sol(femspace.bd_dof(dofs_west(b)),:);
sol_west_locations = interp1(a, sol_west, locations);
dofs_east = find(edge_east(mesh.node(femspace.bd_dof,:))==1);
[a, b] = sort(mesh.node(femspace.bd_dof(dofs_east),2));
sol_east = sol(femspace.bd_dof(dofs_east(b)),:);
sol_east_locations = interp1(a, sol_east, locations);
values(1,:,:) = sol_west_locations;
values(2,:,:) = sol_north_locations;
values(3,:,:) = sol_east_locations;
end