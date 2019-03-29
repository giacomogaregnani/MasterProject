N = size(mesh.node,1);
bmesh = struct('node', zeros(0,3), 'elem', zeros(0,4));
for k=0:1
  for j=0:1
    for i=0:1
      m = (-1).^[i,j,k];
      n = i+2*j+4*k;
      bmesh.node = [bmesh.node; bsxfun(@times,m,mesh.node)];
      bmesh.elem = [bmesh.elem; mesh.elem + n*N];
    end
  end
end

bmesh = remove_duplicate_nodes(bmesh);