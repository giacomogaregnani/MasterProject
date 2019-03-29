%% generate mesh
mesh = structured_mesh([0,1,0,1,0,1],10);

%% display mesh
figure;
simpplot(mesh);


%% display cut mesh
figure;
clear options;
options.expr = 'mesh.bary(:,1)+mesh.bary(:,2)-mesh.bary(:,3)>=0.001';
simpplot(mesh, options);