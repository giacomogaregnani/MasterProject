function subsimplices = get_subsimplices(elem, dim)
%GET_SUBSIMPLICES gets all subsets of a prescribed size

d = size(elem, 2);
NT = size(elem,1);

% opposite edges
if ~isnumeric(dim) && strcmp(dim,'opposite')
	comb = flipdim(fastcombnk(d,d-1),1);
	dim = d-2;
elseif isnumeric(dim)
	comb = fastcombnk(d, dim+1);
else
	error('Wrong dim parameter');
end

if (NT == 1)
	subsimplices = reshape(elem(comb), size(comb));
else
	ncomb = size(comb,1);
	subsimplices = zeros(ncomb*NT, dim+1);
	for i=1:ncomb
		subsimplices((i-1)*NT+1: i*NT,:) = elem(:,comb(i,:));
	end
end

	function res = fastcombnk(A,B)
	switch A
		case 1
			switch B
				case 1
					res = 1;
				otherwise
					res = zeros(1,0);
			end
		case 2
			switch B
				case 1
					res = [1;2];
				case 2
					res= [1,2];
				otherwise
					res = zeros(1,0);
			end
		case 3
			switch B
				case 1
					res = [1;2;3];
				case 2
					res= [1,2; 1,3; 2,3];
				case 3
					res= [1,2,3];
				otherwise
					res = zeros(1,0);
			end
		case 4
			switch B
				case 1
					res = [1;2;3;4];
				case 2
					res= [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
				case 3
					res= [1,2,3; 1,2,4; 1,3,4; 2,3,4];
				case 4
					res= [1,2,3,4];
				otherwise
					res = zeros(1,0);
			end
		otherwise
			res = sortrows(combnk(1:A, B));
	end
	end
end