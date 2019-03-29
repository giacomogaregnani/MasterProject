function dof = get_dof(mesh, whe, femspace, i)
%GET_DOF Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) -1;
NT = size(mesh.elem, 1);
if ~isfield(femspace,'elem2dof'),
    BD = 0;
else
    BD = size(femspace.elem2dof,2);
end

if strcmp(femspace.elemtype,'p0d') || ...
    strcmp(femspace.elemtype,'p1d') || ...
    strcmp(femspace.elemtype,'p2d') || ...
    strcmp(femspace.elemtype,'p3d')
  dof = (i: femspace.ldof : NT*femspace.ldof).';

	if ~strcmp(whe,'all'), dof = dof(whe); end
else
	if (i <= dim+1) %% NODE DOF
		if strcmp(whe,'all')
			dof = mesh.elem(:,i);
		else
			dof = mesh.elem(whe,i);
		end
    elseif (i<= dim+1 + BD) %% BOUNDARY DOF
        if strcmp(whe,'all'),
			dof = femspace.elem2dof(:,i-dim-1);
		else
			dof = femspace.elem2dof(whe,i-dim-1);
        end
    else %% INTERIOR DOF
        j = i-(dim+1)-BD;
        if strcmp(whe,'all')
            dof = femspace.lastbdnode + (NT*(j-1)+1 : NT*j)';
        else
            dof = femspace.lastbdnode + whe;
        end
	end
end
end