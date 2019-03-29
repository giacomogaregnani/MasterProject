function loc = get_loc(mesh, femspace)
%GET_LOC gets coordinates of all nodes of a femspace
loc = zeros(femspace.ndof, 2);
for i=1:femspace.ldof
  loc(get_dof(mesh,'all',femspace,i), :) = ...
    get_rc(mesh, [], 'all', [], femspace.nodelambda(i,:));
end
end

