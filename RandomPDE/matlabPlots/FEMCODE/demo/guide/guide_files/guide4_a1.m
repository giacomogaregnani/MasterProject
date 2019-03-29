function tensor = guide4_a1(x, k ,l)
if k~=l
%    tensor = zeros(size(x,1),1);
   % tensor = 0.15*(cos(x(:,1)*2*pi) + sin(x(:,2)*2*pi));
   tensor = zeros(size(x,1),1);
else
  %  tensor = 1 + 0.3*(sin(x(:,k)*2*pi) + sin(x(:,2)*2*pi));
  tensor = 64/9/sqrt(17) *  (sin(2*pi*x(:,1)) +9/8) .*  (cos(2*pi*x(:,2)) +9/8);
end
end