function r = x2par_expe( x )
%X_TO_PARAM Summary of this function goes here
%   Detailed explanation goes here

% THIS INFO IS ALSO IN ANOTHER SCRIPT
r = -x(:,1).^2/8 + (3-x(:,2))/3;
% 1/2 + (sin(x(:,1))+1)/2 .* (cos(2*pi*x(:,2)/3)+1)/2;
r = r - floor(r);
end

