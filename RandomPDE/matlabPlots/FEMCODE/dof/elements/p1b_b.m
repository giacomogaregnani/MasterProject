function res = p1b_b( mesh, whe, lambda, bf, der )
%P1B_B2 Summary of this function goes here
%   Detailed explanation goes here

% CAN BE OPTIMIZED TO LESS OPERATIONS WHEN some fields order are > 1

if strcmp(whe,'all')
	whe = (1:size(mesh.elem,1))';
end
NT = numel(whe);
NL = size(lambda,1);
dim = size(mesh.elem, 2) -1;
CD = (dim+1)^(dim+1);

nder = sum(der);
d= zeros(1,sum(der));
for i=1:dim
	d(sum(der(1:i-1))+1 : sum(der(1:i))) = i;
end

if (bf <= dim + 1) && (nder == 0)
	res = repmat(lambda(:,bf), [NT/NL, 1] );
elseif (bf <= dim + 1) && (nder == 1)
	res = mesh.dlambda(whe,d,bf);
elseif (bf == dim + 2) && (nder <= dim + 1)
	wdlambda = get_subsimplices(1:(dim+1), nder-1);
	wlambda  = flipdim(get_subsimplices(1:(dim+1),dim - nder),1);
	
	res = zeros(NT,1);
	for k=1:size(wdlambda,1)
		resa = CD;
		for j=1:size(wlambda, 2)
			resa = resa .* lambda(:, wlambda(k,j));
		end
				
		per = perms(wdlambda(k,:));
		resb = zeros(NT,1);
		for l = 1:size(per,1)
			resc = ones(NT,1);
			for j = 1:size(per,2)
				resc = resc .* mesh.dlambda(whe, d(j), per(l,j));
			end
			resb = resb + resc;
		end
		
		res = res + resa.*resb;
	end
else
	res = zeros(NT,1);
end
 



