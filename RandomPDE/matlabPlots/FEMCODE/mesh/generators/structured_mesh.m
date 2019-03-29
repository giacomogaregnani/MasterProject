function  mesh  = structured_mesh( box, N, options )
%STRUCTURED_MESH generates structured *mesh* in a rectangle/cube 
%   This function generates 2D and 3D structured simplicial and
%   rectangular/parallelpipedal meshes in arbitrary rectangle or
%   parallelpiped aligned with coordinate axes and arbitrary number of
%   nodes in each dimension. Several types of discretizations and of
%   boundary conditions are available. The creation process is vectorized
%   and hence suitable for very fast generation even for large number of
%   degrees of freedom.
%
%   box - rectangle or parallelpiped description as 
%   [x1, x2, y1, y2(, z1, z2)],   x1 < x2, y1 < y2, z1 < z2
%
%   N - vector [N(1), N(2) (, N(3))]. We divide the box into 
%   N(1) x N(2) (x N(3)) smaller rectangles of parallelpipedes - let us
%   call them SUBBOXES. Then we divide each subbox into simplicial elements 
%   (if needed).
%
%   options - describes boundary conditions (options.bc), shape of elements
%   (options.shape) and the type of the mesh (options.type).    
%       options.bc = 'dirichlet', 'neumann', 'periodic', 'periodic-index'
%       options.centre =    0 (default - without subbox centers), 
%                           1 (with subbox center)
%       options.shape = 's' (default - simplical - centre = 0/1), 
%                       'q' (quadrilaterals / parallelpipeds - centre = 0)
%                       'p' (pyramids - dim = 3, centre = 1)
%       options.type = 1 (default - simple), 2, 3, ....


%% INITIALIZATION AND PARAMETER CHECKING
dim = round(numel(box)/2); % dimension and N, checking
if (dim<1) || (dim>3) || (2*dim ~= numel(box))
    error('Invalid dimensions of parameters box');
end
if (numel(N) == 1), N = repmat(N,[1,dim]); end
if (numel(N) ~= dim), error('Invalid N'); end

if sum(box(1:2:end) >= box(2:2:end)) % check condition x1<x2, y1<y2, z1<z2
    error('box is not given as lower-left and upper-top point/singular');
end
if nargin < 3, options = struct; end

if isfield(options,'periodic') && options.periodic
	options.bc='periodic';
end

if ~isfield(options,'bc')
	options.bc = 'dirichlet';  % use default settings
elseif strcmp(options.bc,'periodic') && sum(N<3)
    error('Too small N for periodic boundary conditions');    
end

if ~isfield(options,'type'),   options.type = 1;   end
if ~isfield(options,'centre'), options.centre = 0; end

if ~isfield(options,'shape')
    options.shape = 's';
elseif ~strcmp(options.shape,'s') && ... 
    ~(strcmp(options.shape,'q') && ~options.centre) &&...
    ~(strcmp(options.shape,'p') && (dim == 3) && options.centre)
        error('Unknown shape of elements');   
end

if (dim==1)
	if strcmp(options.bc,'periodic')
		mesh.node = linspace(box(1),box(2) - (box(2)-box(1))/N, N)';
		mesh.elem = [1:N; 2:N, 1]';
		mesh.periodic = true;
		mesh.box = box;
	else
		mesh.node = linspace(box(1),box(2),N+1)';
		mesh.elem = [1:N; 2:N+1]';
		mesh.periodic = false;
	end
	return;
end

%% Definition of auxiliary variables
x1 = box(1); x2 = box(2); xx = N(1);
y1 = box(3); y2 = box(4); yy = N(2);
if (dim == 3)
    z1 = box(5); z2 = box(6); zz = N(3);
end

if options.centre
    xx1 = x1 + (x2-x1) / (2*xx);   xx2 = x2 - (x2-x1) / (2*xx); 
    yy1 = y1 + (y2-y1) / (2*yy);   yy2 = y2 - (y2-y1) / (2*yy);
    if (dim==3),	
		zz1 = z1 + (z2-z1) / (2*zz);  zz2 = z2 - (z2-z1) / (2*zz); 
    end
end

mesh.periodic = strcmp(options.bc, 'periodic') == 1;
if mesh.periodic, mesh.box = box; end

if strcmp(options.bc,'periodic') 
	x=xx;   x2 = x2 - 1/x*(x2-x1);    
	y=yy;   y2 = y2 - 1/y*(y2-y1);    
    if (dim==3),
		z=zz;    z2 = z2 - 1/z*(z2-z1);   
    end
else
	x=xx+1; y=yy+1; 
    if (dim == 3), z=zz+1; end
end

%% AUXILIARY FUNCTIONS 
    function elem_correct(ind, direction, cor)
        num_elem = size(L, 1);
        for i= 1 : num_elem
            ind1 = ind - (num_elem -i);
            ind2 = find(loc(L(i,:),direction) == 1);
            elem(ind1, ind2) = elem(ind1, ind2) + cor;
        end
    end            

%% NODE GENERATION
if (dim == 2)
	mesh.node(:,2) = repval(linspace(y1, y2, y)', [x,1]);
	mesh.node(:,1) = repmat(linspace(x1, x2, x)', [y,1]);
elseif (dim == 3)
	mesh.node(:,3) = repval(linspace(z1, z2, z)', [x*y,1]);
	mesh.node(:,2) = repmat(repval(linspace(y1, y2, y)', [x, 1]), [z,1]);
	mesh.node(:,1) = repmat(linspace(x1, x2, x)', [y*z,1]);
end

%% CENTERS OF SUBBOXES GENERATION
if options.centre 
    if (dim == 2)        
        midnode(:,1) = repmat(linspace(xx1, xx2, xx)', [yy,1]);
        midnode(:,2) = repval(linspace(yy1, yy2, yy)', [xx,1]);
        mesh.node = [mesh.node; midnode];
    elseif (dim == 3)
        midnode(:,1) = repmat(linspace(xx1, xx2, xx)', [yy*zz,1]);
        midnode(:,2) = repmat(repval(linspace(yy1, yy2, yy)', ...
			[xx,1]), [zz,1]);
        midnode(:,3) = repval(linspace(zz1, zz2, zz)', [xx*yy,1]);
        mesh.node = [mesh.node; midnode];
    end
    clear midnode;
end

%% FIRST SUBBOX DEFINITION
V = [1, 2, x+2, x+1];
if (dim == 2) 
    V = [V ,x*y+1];
    loc = [0,0; 1,0; 1,1; 0,1; 
           0.5, 0.5];
else % (dim == 3) 
    V = [V, y*x+1, y*x+2, (y+1)*x+2, (y+1)*x+1, x*y*z+1]; 
    loc = [0,0,0; 1,0,0; 1,1,0; 0,1,0; 
           0,0,1; 1,0,1; 1,1,1; 0,1,1; 
           0.5,0.5,0.5];
end

%% ELEMENTS OF FIRST SUBBOX DEFINITION
%       4---------3         Y
%       |         |         |
%       |    5    |         |
%       |         |         |
%       1---------2         O---------X
if (dim == 2) 
    if strcmp(options.shape,'q')
        L = 1:4;
    elseif ~options.centre
        switch options.type 
            case 2
                L = [4 2 3; 2 4 1]; % [3 4 2; 1 2 4]
            otherwise % default
                L = [3 1 2; 1 3 4]; %[2 3 1; 4 1 3];
        end
    else 
        switch options.type 
            otherwise
                spec = 3; norm = 1:2;
                L = [1 2 5; 2 3 5; 3 4 5; 4 1 5];
        end
    end    
%          8---------7
%         /|        /|
%        / |       / |                Z  Y
%       5---------6  |                | /
%       |  4------|--3   AND 9 in     |/
%       | /       | /    the centre   O------X
%       |/        |/
%       1---------2	
else % (dim == 3) 
    if strcmp(options.shape,'q')
        L = 1:8;
    elseif strcmp(options.shape,'p')
        L = [1 2 3 4 9; 5 8 7 6 9;
             2 6 7 3 9; 1 4 8 5 9;
             2 1 5 6 9; 3 7 8 4 9];
        spec = 5; norm = 1:4;
    elseif ~options.centre
        switch options.type 
            otherwise % default
                L = [1 7 2 3; 1 7 3 4; 1 7 4 8;
                     1 7 8 5; 1 7 5 6; 1 7 6 2];
        end
    else 
        switch options.type 
            otherwise % default
                L = [1 3 4 9; 3 1 2 9; % pyramid 12349 and its opposite:		 
                     5 7 6 9; 7 5 8 9;  % pyramid 58769
                     2 7 3 9; 7 2 6 9;  % pyramid 26739 and its opposite:		 
                     1 8 5 9; 8 1 4 9;  % pyramid 14859		 
                     2 5 6 9; 5 2 1 9;  % pyramid 21569 and its opposite:		 
                     3 8 4 9; 8 3 7 9]; % pyramid 37849
                spec = 4; norm = 1:3;
        end
    end
end

%% FAST GENERATION OF OTHER ELEMENTS
elem = V(L);
if (dim == 2)
    if ~options.centre && ~strcmp(options.bc,'periodic')
        elem = repmat(elem, [x-1,1]) + repval((0:x-2)', size(elem));
        elem = repmat(elem, [y-1,1]) + repval((0:x:(y-2)*x)', size(elem));        
    elseif options.centre && ~strcmp(options.bc,'periodic')
        elem = repmat(elem, [x-1,1]) + repval((0:1:x-2)', size(elem));
        siz = size(elem,1);
		elem = repmat(elem,[y-1,1]);
        elem(:,spec) = elem(:,spec) + repval((0:x-1:(y-2)*(x-1))', ...
			[siz,numel(spec)]);
        elem(:,norm) = elem(:,norm) + repval((0:x:(y-2)*x)',       ...
			[siz,numel(norm)]);
    else % periodic
        elem = repmat(elem,[x*y,1]) + repval((0:x*y-1)', size(elem));
        N = size(L,1);  % Correct periodicity
        ind = N*x :N*x: N*x*y;             elem_correct(ind, 1, -x );        
        ind = N*(x*(y-1)+1) : N : N*x*y;   elem_correct(ind, 2, -x*y); 
    end
else % dim == 3
    if ~options.centre && ~strcmp(options.bc,'periodic')
        elem = repmat(elem, [x-1,1]) + repval((0:x-2)', size(elem));
        elem = repmat(elem, [y-1,1]) + repval((0:x:(y-2)*x)', size(elem));
        elem = repmat(elem, [z-1,1]) + ...
			repval((0:x*y:(z-2)*x*y)', size(elem));        
    elseif options.centre && ~strcmp(options.bc,'periodic')
        elem = repmat(elem,[x-1,1]) + repval((0:x-2)', size(elem));
        siz = size(elem,1);
        elem = repmat(elem, [y-1,1]);
        elem(:,norm) = elem(:,norm) + repval((0:x:(y-2)*x)',       ...
			[siz, numel(norm)]);
        elem(:,spec) = elem(:,spec) + repval((0:x-1:(y-2)*(x-1))', ...
			[siz, numel(spec)]);
        siz = size(elem,1);
        elem = repmat(elem, [z-1,1]); 
        elem(:,norm) = elem(:,norm) + ...
			repval((0:x*y:(z-2)*x*y)',                 [siz, numel(norm)]);
        elem(:,spec) = elem(:,spec) + ...
			repval((0:(x-1)*(y-1):(z-2)*(x-1)*(y-1))', [siz, numel(spec)]);        
    else % periodic
        elem = repmat(elem,[x*y*z,1]) + repval((0:x*y*z-1)', size(elem));
        N = size(L,1);      % Correct periodicity
        ind = N*x  :N*x: N*x*y*z;                   
		elem_correct(ind, 1, -x );        
        ind = repmat(N*(x*(y-1)+1):N:N*x*y,[1,z]) + ...
            repval(0:N*x*y:N*x*y*(z-1), [1,x]);     
		elem_correct(ind, 2, -x*y);        
    	ind = N*(x*y*(z-1)+1) : N : N*x*y*z;        
		elem_correct(ind, 3, -x*y*z);	
    end
end

if strcmp(options.shape,'s'),     mesh.elem  = elem;
elseif strcmp(options.shape,'q'), mesh.elemq = elem;
else                              mesh.elemp = elem;
end

end
