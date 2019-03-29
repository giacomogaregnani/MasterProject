function res = p2_b(mesh, whe, lambda, bf, der)
%P1_B2 Summary of this function goes here
%   Detailed explanation goes here

if strcmp(whe,'all')
	NT = size(mesh.elem,1);
else
	NT = numel(whe);
end
NL = size(lambda,1);
dim = size(mesh.elem, 2) - 1;
id=find(der);

tr = get_subsimplices(1:dim+1, 1);
if (bf > dim + 1)
	i1 = tr(bf-(dim+1),1);
	i2 = tr(bf-(dim+1),2);
end

if (sum(der)==0)
	if (bf <= dim +1) 
		res = repmat(lambda(:,bf).*(2*lambda(:,bf)-1), [NT/NL, 1] );
	else 
		res = repmat(4*lambda(:,i1).*lambda(:,i2), [NT/NL, 1] );
	end
elseif (sum(der)==1)
	if (bf <= dim+1)
		if strcmp(whe,'all')
			res = (4*lambda(:,bf)-1) .* mesh.dlambda(:,id,bf);
		else
			res = (4*lambda(:,bf)-1) .* mesh.dlambda(whe,id,bf);
		end
	else
		if strcmp(whe,'all')
			res = 4*lambda(:,i1) .* mesh.dlambda(:,id,i2) ...
				+ 4*lambda(:,i2) .* mesh.dlambda(:,id,i1);
		else
			res = 4*lambda(:,i1) .* mesh.dlambda(whe,id,i2) ...
				+ 4*lambda(:,i2) .* mesh.dlambda(whe,id,i1);
		end
	end
elseif (sum(der)==2) && (numel(id) == 1)
	if (bf <= dim+1)
		if strcmp(whe,'all')
			res = 4 * mesh.dlambda(:,id,bf).^2;
		else
			res = 4 * mesh.dlambda(whe,id,bf).^2;
		end
	else
		if strcmp(whe,'all')
			res = 8 * mesh.dlambda(:,id,i1) .* mesh.dlambda(:,id,i2);
		else
			res = 8 * mesh.dlambda(whe,id,i1) .* mesh.dlambda(whe,id,i2);
		end
	end
elseif (sum(der)==2) && (numel(id) == 2)
	if (bf <= dim+1)
		if strcmp(whe,'all')
			res = 4 * mesh.dlambda(:,id(1),bf) .* ...
				      mesh.dlambda(:,id(2),bf);
		else
			res = 4 * mesh.dlambda(whe,id(1),bf) .* ...
				      mesh.dlambda(whe,id(2),bf);
		end
	else
		if strcmp(whe,'all')
			res = 4 * mesh.dlambda(:,id(1),i1) .* ...
				      mesh.dlambda(:,id(2),i2) ...
				+ 4 * mesh.dlambda(:,id(2),i1) .* ...
				      mesh.dlambda(:,id(1),i2);
		else
			res = 4 * mesh.dlambda(whe,id(1),i1) .* ...
				      mesh.dlambda(whe,id(2),i2) ...
				+ 4 * mesh.dlambda(whe,id(2),i1) .* ...
				      mesh.dlambda(whe,id(1),i2);
		end
	end
else
	res = zeros(NT,1);
end


end
