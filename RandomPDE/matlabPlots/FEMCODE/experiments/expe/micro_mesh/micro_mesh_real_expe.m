function mesh = micro_mesh_real_expe( p, vp )
%MICRO_MESH_REAL_EXPE Summary of this function goes here
%   Detailed explanation goes here

eps = vp.epsilon;
k = floor(p/eps + 1/2);
shifts = [-0.5,-0.5; -0.5,0.5; 0.5,-0.5; 0.5,0.5];
num = 0;
bdgmsh = cell(4,1);
for i=1:4
	x = eps * (k+shifts(i,:));
	if ((x(1) >= 0) && (x(1)<=2) && (x(2)>=0) && (x(2) <= 3)) && ...
			~((x(1) >=1) && (x(2) > 1) && (x(2) < 2))
		num = num+1;
	else
		continue;
	end			
	[~,bdgmsh{num}] = pormat_expe(x,0);
end

bdgmsh_final = bdgmsh{1};
for i=5:8
	max_dist = Inf;
	for j=1:num
		new_dist = sum((eps*(k+shifts(j,:)+bdgmsh{j}.node(i,:)) - p).^2);
		if new_dist < max_dist
			bdgmsh_final.node(i,:) = bdgmsh{j}.node(i,:);
			max_dist = new_dist;
		end
	end
end

options.lc = 0.3;
options.data_loc = '~/repos/experiments/auxdata/';
options.gmsh = '~/gmsh/bin/./gmsh';

%% MESH GENERATION
mq = 0;
while mq < 0.3
	mesh = gmsh(bdgmsh_final, options);
	mq = min(simplex_quality(mesh));
	options.lc = options.lc / 1.2;
end	

end

