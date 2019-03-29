function output = detn(mat,a,b)
%DETN computes generalized determinant of an array of matrices
% 
% output = detn(mat,a,b) computes the generalized determinant of an array
% of (relatively small) matrices using the standard sum-product formula.
% Parameters a and b specify the subsets of rows and columns that should be
% considered.
%
% if numel(a) =  numel(b), then output(i) = determinant of mat(i,a,b).
%
% if numel(a) ~= numel(b), then output(i) = Euclidean norm of a vector
% consisting of the determinants of maximal square submatrices of 
% mat(i,a,b).

%% INITIALIZATION
% Number of Determinants
ND = size(mat,1);

% Implicit parameter settings
if nargin < 2, a = 1:size(mat,2); end
if nargin < 3, b = 1:size(mat,3); end

% matrix dimensions
s1 = numel(a);
s2 = numel(b);

%% DIMENSIONS 1,2,3 - OPTIMIZED (WITHOUT RECURSIVE CALLING)
if (s1 == s2) && (s1 <= 3)
	if (s1 == 1)
			output = reshape(mat(:,a,b), [ND, 1]);
	elseif (s1 == 2)
		output = mat(:,a(1),b(1)).*mat(:,a(2),b(2)) ...
			   - mat(:,a(2),b(1)).*mat(:,a(1),b(2));
	elseif (s1 == 3)
		output = mat(:,a(1),b(1)).*(mat(:,a(2),b(2)).*mat(:,a(3),b(3))  ...
			                      - mat(:,a(2),b(3)).*mat(:,a(3),b(2))) ...
			   + mat(:,a(1),b(2)).*(mat(:,a(2),b(3)).*mat(:,a(3),b(1))  ...
			                      - mat(:,a(2),b(1)).*mat(:,a(3),b(3))) ...
			   + mat(:,a(1),b(3)).*(mat(:,a(2),b(1)).*mat(:,a(3),b(2))  ... 
			                      - mat(:,a(2),b(2)).*mat(:,a(3),b(1)));
	end
end

%% GENERAL CASES (below optimized for low dimensions)
if (s1 == s2) && (s1 > 3)
	output = zeros(ND,1);
	for k=1:s1
		output = output + (-1)^(k+1) * mat(:,k,1) .* ...
		detn(mat, a([1:k-1, k+1:s1]), b(2:s1));
	end
	return;
elseif (s1 ~= s2)
	if s1 < s2
		B = combnk(b, s1);
		NC = size(B,1);
		A = repmat(a, [NC, 1]);
	else
		A = combnk(a,s2);
		NC = size(A, 1);
		B = repmat(b, [NC, 1]);
	end
	output = zeros(ND,1);
	for i=1:NC
		output = output + detn(mat, a(A(i,:)), b(B(i,:))).^2;
	end
	output = sqrt(output);
	return;
end

end