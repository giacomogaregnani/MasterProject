function [Theta, ThetaF, AX, Fq] = reduce_X(Theta, ThetaF, AX, Fq)
%REDUCE_X Summary of this function goes here
%   Detailed explanation goes here

NX = size(AX{1},1);

prf = [1 0 0; 1 0 0; 0 0 -2; 1 2 0; 1 0 -2; 0 2 0; 1 0 0; 1 0 0;
  0 0 2; 1 -2 0; 1 0 2; 0 -2 0; 1 0 0; 1 2 0; 1 0 -2; 1 -2 0; 1 0 2];
A1 = sparse(NX,NX);
Amu1 = sparse(NX,NX);
Amu2 = sparse(NX,NX);

s1 = sym('1');
s2 = sym('mu1');
s3 = sym('mu2');

for i=13:numel(Theta)
  if ~(prf(i-12,1)*s1 + prf(i-12,2)*s2 + prf(i-12,3)*s3 - Theta{i} == 0),
    error('bad precomputed values');
  end
  A1   = A1   + prf(i-12,1) * AX{i};
  Amu1 = Amu1 + prf(i-12,2) * AX{i};
  Amu2 = Amu2 + prf(i-12,3) * AX{i};
end
Theta(16:end) = [];
Theta{13} = s1;
Theta{14} = s2;
Theta{15} = s3;
AX(16:end) = [];
AX{13} = A1;
AX{14} = Amu1;
AX{15} = Amu2;
end

