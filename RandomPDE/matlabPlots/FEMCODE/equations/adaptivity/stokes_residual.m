function res = stokes_residual(mesh, vp, ufemspace, usol, pfemspace, psol)
%STOKES_RESIDUAL computes the Stokes residual

if vp.a ~= 1
    error('Stokes residual was not yet updated for general viscosity');
end

%% initialization
dim = size(mesh.elem,2) - 1;
if strcmp(vp.elemtype,'p1b')
  ufemspace = get_femspace(mesh,'p1');
	usol = usol(1:size(mesh.node,1),:,:);
end

if ~isfield(mesh,'volume')
  mesh.volume = simplex_volume(mesh);
end

%% COMPUTE RESIDUALS
jumpResidual = stokes_jump_residual;
elemResidual = stokes_elem_residual;
divResidual = stokes_div_residual;	

%% ASSIGN JUMP RESIDUALS TO ELEMENTS
elem2subsim = auxstructure(mesh,'elem2subsim');
res = 0;
for ii=1:dim+1
	res = res + jumpResidual(elem2subsim(:,ii),:,:); 
end

%% COMPUTE RESIDUALS
diam = simplex_diameter(mesh);
res = 1/2 * bsxfun(@times, res, diam) + ...
	bsxfun(@times,elemResidual,diam.^2) + ...
	divResidual;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function eres = stokes_elem_residual
	%STOKES_ELEM_RESIDUAL computes interior residuals of a Stokes FEM sol.
	%
	% res = stokes_elem_residual
	% computes the interior residuals of the discrete solution over the
	% elements. For every k=1:NT we compute
	%
	%  || f - \nabpa p^h + \Delta v^h ||_{L^2(element num. k)}

	deg = max([ufemspace.deg - 2, pfemspace.deg - 1, 1]);	
	[lambda, weight] = quadpts(dim, deg);	
	eres = 0;
	for j=1:numel(weight)
		fx = evalf(mesh, 'all', [], lambda(j,:), vp.f);
		for i=1:dim % i-th coordinate
			ares = fx(:,i,:) - evalf(mesh, 'all', pfemspace, ...
				lambda(j,:), psol, i);
			for k=1:dim
				der = zeros(dim,1); der(k) = 2; % second derivative in x_k
				ares = ares + evalf(mesh, 'all', ufemspace, ...
					lambda(j,:), usol(:,i,:), der);
			end
			eres = eres + weight(j) * ares.^2;
		end
	end
	eres = bsxfun(@times, eres, mesh.volume);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function dres = stokes_div_residual
	%DIV_RESIDUAL Summary of this function goes here
	%   Detailed explanation goes here
	
	deg = ufemspace.deg - 1;	
	[lambda, weight] = quadpts(dim, deg);
	dres = 0;
	for i = 1:numel(weight)
		ares = 0;
		for j = 1:dim
			ares = ares + evalf(mesh, 'all', ufemspace, lambda(i,:), usol(:,j,:), j);
		end
		dres = dres + weight(i) * ares.^2;
	end
	dres = bsxfun(@times, dres,mesh.volume);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function jres = stokes_jump_residual
	%STOKES_JUMP_RESIDUAL computes jump residuals of a Stokes FEM solution
	%
	% res = stokes_jump_residual
	% computes the jump residuals of the discrete solution over the edges
	% (faces). For every k=1:NE, where NE is the number of edges (faces) we
	% compute
	%
	%      || [ dv^h / dn ]_{jump} ||_{L^2(e)}^2
	%
	% where e is the k-th edge (given by nodes subsim(k)).
	
	[normals, subsim, subsim2elem] = get_normals(mesh);	
	NE  = size(subsim,1);
	deg = max([ufemspace.deg - 1, 1]);
	
	%% Arrays lam1 and lam2 - computation
	% these arrays are of size NE x dim+1 with values in {1, ..., dim+1}.
	% They provide a connection between the node numbers in subsim(k,:) and the
	% node numbers of mesh.elem(subsim2elem(k,1),:) resp.
	% mesh.elem(subsim2elem(k,2),:)
	%
	% lam1(k,:) gives the local indices of mesh.elem(subsim2elem(k,1),:) in
	% subsim(k,:), while the one additional node has local index dim+1.
	%
	% EXAMPLE:
	% IF
	%   subsim(k,:) = [4, 8, 9]
	%   mesh.elem(subsim2elem(k,1),:) = [9, 2, 4, 8]
	%   mesh.elem(subsim2elem(k,2),:) = [8, 4, 9, 7]
	% THEN
	%   lam1(k,:) = [3, 4, 1, 2]
	%   lam2(k,:) = [2, 1, 3, 4]
	%
	% COMMENT: the current implementation is vectorized but not the fastest
	% possible
	lam1 = repmat(dim+1,[NE,dim+1]);
	lam2 = repmat(dim+1,[NE,dim+1]);
	for i=1:dim+1
		for j=1:dim
			ix1 = subsim(:,j) == mesh.elem(subsim2elem(:,1),i);
			ix2 = subsim(:,j) == mesh.elem(subsim2elem(:,2),i);
			lam1(ix1,i) = j;
			lam2(ix2,i) = j;
		end
	end
	clear ix1 ix2;
	
	% quadrature formula on subsims (last weight is set tozero)
	[lambda, weight] = quadpts(dim-1, deg);
	lambda(:,dim+1) = 0;
	
	%% residual computation
	jres = 0;
	for j=1:numel(weight)
		lambdaj = lambda(j,:); % defines integration point x_{lambdaj}
		lambdaj1 = lambdaj(lam1);
		lambdaj2 = lambdaj(lam2);
		for i=1:dim
			ares = 0;
			for k=1:dim
				% here we add:       [ dv^h_i(x_{lambdaj})/dx_k * n_k ]_jump
				ares = ares + bsxfun(@times, ...
					evalf(mesh, subsim2elem(:,1), ufemspace, ...
					lambdaj1, usol(:,i,:), k) - ...
					evalf(mesh, subsim2elem(:,2), ufemspace, ...
					lambdaj2, usol(:,i,:), k), ...
					normals(:,k));
			end
			jres = jres + ares.^2 * weight(j);
		end
	end
	jres = bsxfun(@times, jres,subsimplex_volume(mesh,[],[],[],subsim));
	end
end