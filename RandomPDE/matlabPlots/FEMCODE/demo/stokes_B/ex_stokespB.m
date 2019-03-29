function xx = ex_stokespB(x)
%EX_STOKESPB Summary of this function goes here
%   Detailed explanation goes here
y=x(:,2);
x=x(:,1);

xx = sin(x.^2 + y.^2) - 0.56129039832190474467; % normalization

end

