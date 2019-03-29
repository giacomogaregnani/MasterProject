function res = poisson_residual(mesh, vp, femspace, sol)
%POISSON_RESIDUAL computes the Poisson residual 

%% COMPUTE RESIDUALS

dim = size(mesh.elem, 2) - 1;
alambda = vp.alambda;
aweight = vp.aweight;
aNQ = numel(aweight);

%% ASSEMBLE FLUX AT QUADRATURE POINTS
g = []; % INIT   NT * -- * -- * NQ
for k=1:aNQ
	f = evalf(mesh, 'all', femspace, alambda(k,:), vp.f);
	gg = 0; % INIT
	for i=1:dim
		ggg = []; % INIT
		for j=1:dim
			if vp.fully_discrete % fully discrete tensor
				vpakij = vp.a(:,i,j,k);
			elseif isnumeric(vp.a) 
				vpakij = vp.a(i,j);
			else
				vpakij = vp.a(get_rc(mesh,[],'all',[],alambda(k,:)),i,j);
			end
			%  a \nabla u    OR     a (\nabla u - f)
			gggadd = vpakij .* ...
				(evalf(mesh, 'all', femspace, alambda(k,:), sol, i) ...
				- vp.derivatives * f(:,i+1,:));
			ggg = cat(2,ggg,gggadd);
		end
		gg = gg + ggg;
	end
	g = cat(4, g, gg);
end

%% RESTORE FLUX TO DISCOTINUPUS FE
[rfemspace, restored] = restore(mesh, g, vp.elemtype);
clear g gg ggg;

if ~isfield(mesh,'volume')
	mesh.volume = simplex_volume(mesh);
end

jumpResidual = poisson_jump_residual;
elemResidual = poisson_elem_residual;

%% ASSIGN JUMP RESIDUALS TO ELEMENTS
elem2subsim = auxstructure(mesh,'elem2subsim');
res = 0;
for i=1:dim+1
	res = res + jumpResidual(elem2subsim(:,i),:,:); 
end

%% COMPUTE RESIDUALS
diam = simplex_diameter(mesh);
res = 1/2*bsxfun(@times, res, diam)+ bsxfun(@times, elemResidual, diam.^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function eres = poisson_elem_residual
	eres = 0;
	for kk=1:aNQ
		f = evalf(mesh, 'all', femspace, alambda(kk,:), vp.f);
		addres = f(:,1,:);
		for ii=1:dim
			addres = addres + ...
				evalf(mesh, 'all', rfemspace, alambda(kk,:), restored(:,ii,:), ii);
		end
		eres = eres + addres.^2 * aweight(kk);
	end
	eres = bsxfun(@times, eres, mesh.volume);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function jres = poisson_jump_residual
	% quadrature formula on subsims (last weight is set tozero)
	[normals, subsim, subsim2elem] = get_normals(mesh);
	NE  = size(subsim,1);
	[slambda, sweight] = quadpts(dim-1, rfemspace.deg);
	slambda(:,dim+1) = 0;
	sNQ = numel(sweight);
	
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
	for ii=1:dim+1
		for jj=1:dim
			ix1 = subsim(:,jj) == mesh.elem(subsim2elem(:,1),ii);
			ix2 = subsim(:,jj) == mesh.elem(subsim2elem(:,2),ii);
			lam1(ix1,ii) = jj;
			lam2(ix2,ii) = jj;
		end
	end
	clear ix1 ix2;
	
	isbdsubsim = auxstructure(mesh,'isbdsubsim');
	
	%% residual computation
	jres = 0; % INIT
	for jj=1:sNQ
		lambdaj = slambda(jj,:); % defines integration point x_{lambdaj}
		lambdaj1 = lambdaj(lam1);
		lambdaj2 = lambdaj(lam2);
		ares = 0; % INIT
		for ii=1:dim
			ares = ares + bsxfun(@times, ...
				evalf(mesh, subsim2elem(:,1), rfemspace, ...
				lambdaj1, restored(:,ii,:)) - ...
				evalf(mesh, subsim2elem(:,2), rfemspace, ...
				lambdaj2, restored(:,ii,:)), ...
				normals(:,ii));
% TODO: Different behaviour for boundary edges
%			if strcmp(vp.bc,'natural')
%				ares(isbdsubsim,:,:,:) = ares(isbdsubsim,:,:,:) + ...
%                    bsxfun(@times, ...
%					evalf(mesh, subsim2elem(isbdsubsim,1), rfemspace, ...
%					lambdaj1(isbdsubsim,:), restored(:,ii,:)), ...
%					normals(isbdsubsim,ii));
%			end
		end
		jres = jres + ares.^2 * sweight(jj);
	end
	jres = bsxfun(@times, jres, subsimplex_volume(mesh,[],[],[],subsim));
	end
end

