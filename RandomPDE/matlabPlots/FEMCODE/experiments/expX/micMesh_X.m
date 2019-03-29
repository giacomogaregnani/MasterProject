function mesh = micMesh_X(x)
par = x2par_X(x);
mesh = struct;
mesh.elem = [2,6,1; 5,6,2; 4,5,2; 4,2,3];
mesh.node =  [-1/2,0; par(1),par(2); 0,-1/2; 1/2,-1/2; 1/2,1/2; -1/2,1/2];
mesh = bisect(mesh, [2,3]);
mesh = uniformrefine(mesh);
mesh = periodize(mesh,[-1/2,1/2,-1/2,1/2]);
mesh = label(mesh);
end

