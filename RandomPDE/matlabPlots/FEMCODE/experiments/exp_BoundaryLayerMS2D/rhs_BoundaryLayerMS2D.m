function out = rhs_BoundaryLayerMS2D(x, der)
%RHS_POISSONMS_2D Summary of this function goes here
%   Detailed explanation goes here

if (nargin > 1) && (sum(der) > 0)
	error('Derivatives not implemented');
else

    scal = 1e-3;
    
    a = 1-exp(-1/scal);
    b = @(x,y) exp(-(1-x).*(1-y)./scal);
    f_add = @(x,y) (x+y).*(1-b(x,y))/a + (2*y.*(1-y) + 2*x.*(1-x)).*b(x,y)./a + ...
        1/scal*(x.*y.*(1-y).^2 + x.*y.*(1-x).^2 - x.*y.*(1-x) - x.*y.*(1-y)).*b(x,y)./a;
    out = f_add(x(:,1),x(:,2));
	
end

end

