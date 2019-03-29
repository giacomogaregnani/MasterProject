function [ThetaNew, ThetaFNew, AXNew, FqNew] = reduce_Y(Theta, ThetaF, AX, Fq)
%REDUCE_X Summary of this function goes here
%   Detailed explanation goes here
Q = numel(Theta);
QF = numel(ThetaF);

indices  = cell(Q,1);
indicesF = cell(QF,1);

ind = [];
ind2 = [];

for i=1:Q
  if ~isempty(strfind(char(Theta{i}),'/'))
    ind = [ind, i];
  else
    ind2 = [ind2, i];
  end
end

mapQ = zeros(Q,1);

n=0;
for i=ind
  for j=ind
    if j>=i, 
      n=n+1;
      ThetaNew{n} = Theta{i};
      indices{i} = [n,1];
      break;
    end
    s = simplify(Theta{j}/Theta{i});
    if isempty(strfind(char(s),'mu1')) && isempty(strfind(char(s),'mu2')) 
      %fprintf('%d:%d : %s\n',j,i, char(s));
      mapQ(i) = j;
      indices{i} = [indices{j}(1),1];
      break;
    end
  end
end

mu1 = sym('mu1');
mu2 = sym('mu2');
ThetaNew{n+1} = sym(1);
ThetaNew{n+2} = mu1;
ThetaNew{n+3} = mu2;
ThetaNew{n+4} = mu1*mu2;

for i=ind2
  k=zeros(4,1);
  k(1) = subsmu(Theta{i},[0,0]);
  k(2) = subsmu(diff(Theta{i},mu1),[0,0]);
  k(3) = subsmu(diff(Theta{i},mu2),[0,0]);
  k(4) = subsmu(diff(diff(Theta{i},mu1),mu2),[0,0]);
  indices{i} = [n+(1:4)',k];
end

for i=1:Q
  s = Theta{i};
  for j=1:size(indices{i},1)
    s=s - indices{i}(j,2)*ThetaNew{indices{i}(j,1)};
  end
  s = simplify(s);
  if s~=0
    error('simplification is not correct');
  end
end

ThetaFNew{1} = sym(1);
ThetaFNew{2} = mu1;
ThetaFNew{3} = mu2;
ThetaFNew{4} = mu1*mu2;

for i=1:QF
  k=zeros(4,1);
  k(1) = subsmu(ThetaF{i},[0,0]);
  k(2) = subsmu(diff(ThetaF{i},mu1),[0,0]);
  k(3) = subsmu(diff(ThetaF{i},mu2),[0,0]);
  k(4) = subsmu(diff(diff(ThetaF{i},mu1),mu2),[0,0]);
  indicesF{i} = [(1:4)',k];
end

for i=1:QF
  s = ThetaF{i};
  for j=1:size(indicesF{i},1)
    s=s - indicesF{i}(j,2)*ThetaFNew{indicesF{i}(j,1)};
  end
  s = simplify(s);
  if s~=0
    error('simplification is not correct');
  end
end


Qnew = numel(ThetaNew);
QFNew = numel(ThetaFNew);

for i=1:Qnew
  AXNew{i} = sparse(size(AX{1},1),size(AX{1},2));
end
for i=1:QFNew
  for k=1:size(Fq,2)
    FqNew{i,k} = zeros(size(Fq{1}));
  end
end

for i=1:Q
  for j=1:size(indices{i},1)
    AXNew{indices{i}(j,1)} = AXNew{indices{i}(j,1)} + indices{i}(j,2)*AX{i}; 
  end
end

for i=1:QF
  for j=1:size(indicesF{i},1)
    for k=1:size(Fq,2)
      FqNew{indicesF{i}(j,1),k} = FqNew{indicesF{i}(j,1),k} + indicesF{i}(j,2)*Fq{i,k};
    end
  end
end
end
