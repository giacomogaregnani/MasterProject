function res = p3_b(mesh, whe, lambda, bf, der)
%P1_B2 Summary of this function goes here
%   Detailed explanation goes here

if strcmp(whe,'all')
	NT = size(mesh.elem,1);
	whe = (1:NT)';
else
	NT = numel(whe);
end
dim = size(mesh.elem, 2) - 1;
id=find(der);

L1 = (dim+1);

LNE = (dim+1)*dim/2;
L2 = (dim+1) + 2*LNE;

if (bf <= L1)
	l1 = lambda(:,bf);
	if (sum(der)>0)
		m1 = mesh.dlambda(whe,id,bf);
	end
end
if (bf>L1) && (bf <= L2)
	tr = get_subsimplices(1:dim+1, 1);
	tr = [tr; flipdim(tr,2)];
	tr = tr(bf-L1,:);
    
	l1 = lambda(:,tr(1));
	l2 = lambda(:,tr(2));
	if (sum(der)>0)
		m1 = mesh.dlambda(whe,id,tr(1));
		m2 = mesh.dlambda(whe,id,tr(2));
	end
end
if (bf>L2)
    tr = get_subsimplices(1:dim+1, 2);
	tr = tr(bf-L2,:);
	l1 = lambda(:,tr(1));
	l2 = lambda(:,tr(2));
    l3 = lambda(:,tr(3));
	if (sum(der)>0)
		m1 = mesh.dlambda(whe,id,tr(1));
		m2 = mesh.dlambda(whe,id,tr(2));
		m3 = mesh.dlambda(whe,id,tr(3));
	end
end

NL = size(l1,1);
if (sum(der)==0)
	if (bf <= L1) 
		res = repmat((1/2) * l1 .*  (3*l1-1) .* (3*l1-2), [NT/NL, 1] );
    elseif (bf <= L2) 
		res = repmat((9/2) * ...
            l1 .*  (3*l1 - 1) .* l2, [NT/NL, 1] );
    else
        res = repmat(27 * l1 .* l2 .* l3, [NT/NL, 1]);
	end
elseif (sum(der)==1)
    if (bf <= L1)
		res = (1 + l1 .* (27/2* l1 - 9)) .* m1;
    elseif (bf <= L2)
		res = ((9/2) * l1 .* (3*l1 - 1)) .* m2 ...
			+ (27*l1 -9/2) .* l2 .* m1;
    else
		res = 27 * ( ...
			l1 .* l2 .* m3 + ...
			l1 .* l3 .* m2 + ...
			l2 .* l3 .* m1);
    end
elseif (sum(der)==2) && (numel(id) == 1)
	if (bf <= L1)
        res = (27* l1 - 9) .* m1.^2;
	elseif (bf <= L2)
		res = 27 * l2 .* m1.^2 + (54*l1-9) .* m1 .* m2;
	else
		res = 2*27 * ( ...
			l1 .* m2 .* m3 + ...
			l2 .* m1 .* m3 + ...
			l3 .* m1 .* m2);
	end
elseif (sum(der)==2) && (numel(id) == 2)
    error('not yet implemented');
elseif (sum(der)==3)
    error('not yet implemented');
else
	res = zeros(NT,1);
end

end
