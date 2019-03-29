function bdgmsh = bdmesh2bdgmsh( bdmesh )
%BDMESH2BDGMSH Summary of this function goes here
%   Detailed explanation goes here

% assumptions: 
% 1. not periodic
% 2. 3D

d = size(bdmesh.node,2);
N = size(bdmesh.node,1);
if d ~= 3, error('unsupported dimension'); end

bdgmsh = struct;
bdgmsh.node = bdmesh.node;


[edge, elem2edge] = auxstructure(bdmesh, 'edge','elem2edge');
NE = size(edge,1);
bdgmsh.line = edge;

NT = size(bdmesh.elem,1);
for i=1:NT
  bdgmsh.lineloop{i} = [...
    double(elem2edge(i,1)) * (-1)^(bdmesh.elem(i,2) > bdmesh.elem(i,1)), ...
    double(elem2edge(i,2)) * (-1)^(bdmesh.elem(i,1) > bdmesh.elem(i,3)), ...
    double(elem2edge(i,3)) * (-1)^(bdmesh.elem(i,3) > bdmesh.elem(i,2)) ...
    ];
end

for i=1:3
  nnum = find(bdmesh.node(:,i) == 1);
  NN = numel(nnum);
  a = zeros(NE,1);
  for j=1:NN
    a = a + sum(edge == nnum(j), 2);
  end
  ned = find(a==2);
  NC = numel(ned);
  used = false(NC,1);
  bdgmsh.lineloop{NT+i} = ned(1); used(1) = true; con = edge(ned(1),2);
  for j=2:NC
    for k=1:NC
      if used(k), continue; end
      if edge(ned(k), 1) == con
        bdgmsh.lineloop{NT+i} = [bdgmsh.lineloop{NT+i}; ned(k)];
        used(k) = true;
        con = edge(ned(k), 2);
      elseif edge(ned(k), 2) == con
        bdgmsh.lineloop{NT+i} = [bdgmsh.lineloop{NT+i}; -ned(k)];
        used(k) = true;
        con = edge(ned(k), 1);        
      end
    end
  end  
end
for i= 1:NT+3
  bdgmsh.planesurface{i} = i;
end
bary = get_rc(bdmesh);
bdgmsh.surfaceloop{1} = [find(all(bary<=1,2)); NT+1; NT+2; NT+3];
bdgmsh.surfaceloop{2} = [find(bary(:,1)>=1); NT+1];
bdgmsh.surfaceloop{3} = [find(bary(:,2)>=1); NT+2];
bdgmsh.surfaceloop{4} = [find(bary(:,3)>=1); NT+3];

for i=1:4
  bdgmsh.volume{i} = i;
end




end

