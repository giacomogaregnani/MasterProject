function simpplot_qsol(mesh,sol)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
%[mesh, ~] = deperiodize(mesh);
bary = get_rc(mesh);
d=size(mesh.node,2);
switch d
	case 1
		error('not yet implemented');
	case 2
		quiver(bary(:,1), bary(:,2), sol(:,1), sol(:,2));
		axis equal; axis off;
	case 3
		quiver3(bary(:,1), bary(:,2), bary(:,3), sol(:,1), sol(:,2), sol(:,3));
	otherwise
		error('Unimplemented dimension.');
end
