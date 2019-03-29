function r = x2par_expd( x )
%X_TO_PARAM Summary of this function goes here
%   Detailed explanation goes here

r = (sin(x(:,1))+1)/2 .* (sin(2*pi*x(:,2)/3)+1)/2;
end

