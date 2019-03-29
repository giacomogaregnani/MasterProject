box = [-0.5,0.5,-0.5,0.5,-0.5,0.5];
mesh = structured_mesh(box, 9, struct('centre',true));
abary = abs(get_rc(mesh));
in = ((abary(:,1) < 1/6) & (abary(:,2) < 1/6)) | ((abary(:,1) < 1/6) & (abary(:,3) < 1/6)) | ((abary(:,2) < 1/6) & (abary(:,3) < 1/6));
mesh.elem(~in,:) = [];
mesh = renumber(mesh);
mesh = periodize(mesh, box);