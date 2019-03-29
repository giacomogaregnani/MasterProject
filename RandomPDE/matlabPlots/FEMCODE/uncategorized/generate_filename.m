function myfile = generate_filename(x, direc, nd)
%GENERATE_FILENAME generates filename based on coordinate and two numbers

dim = numel(x);

if dim == 2
	myfile =sprintf('dir%d_qp%d_x%.12d_y%.12d.mat', direc, nd, ...
		round(x(1)*3*2^30),...
		round(x(2)*3*2^30));
elseif dim == 3
	myfile =sprintf('%di%d_x%.12d_y%.12d_z%.12d.mat', direc, nd, ...
		round(x(1)*3*2^30),...
		round(x(2)*3*2^30),...
		round(x(3)*3*2^30));
else
    error('Unsupported dimension');
end
end

