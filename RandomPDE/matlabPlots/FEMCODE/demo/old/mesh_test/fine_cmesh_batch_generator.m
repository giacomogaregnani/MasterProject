for i=0:5
	j=2^i;
	mesh = fine_cmesh_generator(j);
	save(['computed_data/cmesh_' num2str(j) '.mat'],'mesh');
end