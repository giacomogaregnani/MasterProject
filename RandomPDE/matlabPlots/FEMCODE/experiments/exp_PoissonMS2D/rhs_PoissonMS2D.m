function out = rhs_PoissonMS2D(x, der)
%RHS_POISSONMS_2D Summary of this function goes here
%   Detailed explanation goes here

if (nargin > 1) && (sum(der) > 0)
	error('Derivatives not implemented');
else

%    % periodic
   out = 0.*x(:,1) + 0.*x(:,2);
    
%     % locally periodic
%     scal = 1;
% 
% 	out = scal*(1./(x(:,1).^3+5)+1./(x(:,2).^2+0.05+2*(x(:,1).*x(:,2)+1)))*4*pi^2.*sin(2*pi*x(:,1)).*sin(2*pi*x(:,2))...
%     + (0+...%1./(x(:,1).^3+5)
%     + scal*3*x(:,1).^2./(x(:,1).^3+5).^2)*2*pi.*cos(2*pi*x(:,1)).*sin(2*pi*x(:,2))...
%     + (0+...%1./(x(:,2).^2+0.05+2*(x(:,1).*x(:,2)+1))+
%     scal*(2*(x(:,2)+x(:,1)))./(x(:,2).^2+0.05+2*(x(:,1).*x(:,2)+1)).^2)*2*pi.*sin(2*pi*x(:,1)).*cos(2*pi*x(:,2));

	
end

end

