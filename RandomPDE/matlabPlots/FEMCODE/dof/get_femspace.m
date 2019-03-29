function femspace= get_femspace(mesh, elemtype)
%DOF Summary of this function goes here
%   Detailed explanation goes here

h = str2func([elemtype,'_dof']);
femspace = h(mesh);
end

