function res = p4_b(mesh, whe, lambda, bf, der)
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

L1 = (dim+1);
NE = (dim+1)*dim/2;

L2 = (dim+1) + 3*NE;
NF = (dim+1)*dim*(dim-1)/6;

L3 = (dim+1) + 3*NE + 3*NF;

if (bf>L1) && (bf <= L2)
    tr = get_subsimplices(1:dim+1, 1);
    i = mod(bf-L1-1,NE)+1;
	i1 = tr(i,1);
	i2 = tr(i,2);
    switch ceil((bf-L1)/NE)
        case 1
            f=@(x,y)((16/3) * x.*(4*x-1).*(2*x-1).*y);
            g=@(x,y,dx,dy)((16/3)*(  (x.*(24*x -12)+1).*dx.*y  + x.*(4*x-1).*(2*x-1).*dy  ));
        case 2
            f=@(x,y)(4*x.*(4*x-1).*y.*(4*y-1));
            g=@(x,y,dx,dy)(4*(  (8*x-1).*dx.*y.*(4*y-1) + (8*y-1).*dy.*x.*(4*x-1)   ));
        case 3
            f=@(y,x)((16/3) * x.*(4*x-1).*(2*x-1).*y);
            g=@(y,x,dy,dx)((16/3)*(  (x.*(24*x -12)+1).*dx.*y  + x.*(4*x-1).*(2*x-1).*dy  ));
    end      
end
if (bf>L2) && (bf<=L3) && (dim>1)
    tr = get_subsimplices(1:dim+1, 2);
    i = mod(bf-L2-1,NF)+1;
	i1 = tr(i,1);
	i2 = tr(i,2);
    i3 = tr(i,3);
    f=@(x,y,z)(32 * x.*(4*x-1).*y.*z);
    g=@(x,y,z,dx,dy,dz)(32*(  (8*x-1).*dx.*y.*z + x.*(4*x-1).*dy.*z + x.*(4*x-1).*y.*dz    ));
    switch ceil((bf-L2)/NF)
        case 2
            f=@(x,y,z)(f(y,x,z));
            g=@(x,y,z,dx,dy,dz)(g(y,x,z,dy,dx,dz));
        case 3
            f=@(x,y,z)(f(z,x,y));
            g=@(x,y,z,dx,dy,dz)(g(z,x,y,dz,dx,dy));
    end         
end
if (bf>L3) && (dim>2)
    tr = get_subsimplices(1:dim+1, 3);
    i=bf-L3;
    i1 = tr(i,1);
	i2 = tr(i,2);
    i3 = tr(i,3);
    i4 = tr(i,4);
end

if (sum(der)==0)
	if (bf <= L1) 
		res = repmat((1/3) * ...
            lambda(:,bf) .* ...
            (4*lambda(:,bf)-1) .* ...
            (2*lambda(:,bf)-1) .* ...
            (4*lambda(:,bf)-3), [NT/NL, 1] );
    elseif (bf <= L2) 
		res = repmat(f(lambda(:,i1), lambda(:,i2)), [NT/NL, 1]);
    elseif (bf <= L3)
        res = repmat(f(lambda(:,i1), lambda(:,i2), lambda(:,i3)), ...
            [NT/NL, 1]);
    else
        res =  repmat(4^4 * lambda(:,i1) .* lambda(:,i2) .* ...
            lambda(:,i3) .* lambda(:,i4), [NT/NL, 1]);
	end
elseif (sum(der)==1)
    if (bf <= L1)
		if strcmp(whe,'all')
			res = (-1 + lambda(:,bf).* (44/3 + lambda(:,bf) .* ...
                (-48 + (128 *lambda(:,bf))/3))) .* ...
                mesh.dlambda(:,id,bf);
		else
			res = (-1 + lambda(:,bf).* (44/3 + lambda(:,bf) .* ...
                (-48 + (128 *lambda(:,bf))/3))) .* ...
                mesh.dlambda(whe,id,bf);
        end
    elseif (bf <= L2)
		if strcmp(whe,'all')
			res = g(lambda(:,i1), lambda(:,i2), ...
                mesh.dlambda(:,id,i1), mesh.dlambda(:,id,i2));
		else
			res = g(lambda(:,i1), lambda(:,i2), ...
                mesh.dlambda(whe,id,i1), mesh.dlambda(whe,id,i2));
        end
    elseif (bf<=L3)
        if strcmp(whe,'all')
			res = g(lambda(:,i1), lambda(:,i2), lambda(:,i3), ...
                mesh.dlambda(:,id,i1), mesh.dlambda(:,id,i2), ...
                mesh.dlambda(:,id,i3));
        else
			res = g(lambda(:,i1), lambda(:,i2), lambda(:,i3), ...
                mesh.dlambda(whe,id,i1), mesh.dlambda(whe,id,i2), ...
                mesh.dlambda(whe,id,i3));
        end
    else
        if strcmp(whe,'all')
			res = 4^4 * ( ...
                lambda(:,i1) .* lambda(:,i2) .* lambda(:,i3) .* mesh.dlambda(:,id,i4) + ...
                lambda(:,i1) .* lambda(:,i2) .* lambda(:,i4) .* mesh.dlambda(:,id,i3) + ...
                lambda(:,i1) .* lambda(:,i3) .* lambda(:,i4) .* mesh.dlambda(:,id,i2) + ...
                lambda(:,i2) .* lambda(:,i3) .* lambda(:,i4) .* mesh.dlambda(:,id,i1));
        else
			res = 4^4 * ( ...
                lambda(:,i1) .* lambda(:,i2) .* lambda(:,i3) .* mesh.dlambda(whe,id,i4) + ...
                lambda(:,i1) .* lambda(:,i2) .* lambda(:,i4) .* mesh.dlambda(whe,id,i3) + ...
                lambda(:,i1) .* lambda(:,i3) .* lambda(:,i4) .* mesh.dlambda(whe,id,i2) + ...
                lambda(:,i2) .* lambda(:,i3) .* lambda(:,i4) .* mesh.dlambda(whe,id,i1));
        end
    end
elseif (sum(der)==2) && (numel(id) == 1)
	if (bf <= dim+1)
        error('not yet implemented');
	end
elseif (sum(der)==2) && (numel(id) == 2)
	if (bf <= dim+1)
        error('not yet implemented');
	end
elseif (sum(der)==3)
    error('not yet implemented');
elseif (sum(der)==4)
    error('not yet implemented');
else
	res = zeros(NT,1);
end


end
