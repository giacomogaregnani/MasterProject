function A = subs_array(A, from, to)
%SUBS_ARRAY substitutes values in an array A using rule from -> to
%
% subs_array(A, from, to)
%   makes a subsitution in an array A by using the rule:
%      from(i) --->  to(i)
%   for every i.
%
% INPUT:
%   A      n-dimensional numerical array
%   from   a numerical array with NS elements of the same type as A
%          we assume that its elements are distinct
%   to     a numerical array with NS elements of the same type as A
% 
% EXAMPLES:
%   subs_array([1,2,1; 2,1,3], [1,2], [2,5])
%   
%   ans = [2,5,2; 5,2,3]
%
% ALGORITHM:
%   if number of substitutions is small (less than 10), we substitute 
%   values one by one in a vectorized fashion.
%
%   if number os substitutions is large, we sort the arrays and do a
%   non-vectorized sweep of the arrays.

%% argument check
if nargin==1, 
  return; 
elseif nargin<3
  error('not enough parameters');
end
if numel(from) ~= numel(to)
  error('substitution not correctly submitted');
end

%% initialization
NS = numel(from);
N  = numel(A);

if NS<10   % vectorized one-by-one substitution
  notUpdated = true(size(A));
  for i=1:NS
    update = notUpdated & (A==from(i));
    A(update) = to(i);
    notUpdated = notUpdated & ~update;
  end
  
else       % non-vectorized way processing all subsitutions at once 
  % sort A
  sizeA = size(A);
  [A, indA] = sort(A(:));
  
  % sort vfrom and vto
  [from, iv] = sort(from);
  to = to(iv);
  
  % sweep the arrays
  i=1; j=1;
  while (i <= N) && (j <= NS)
    if A(i) == from(j)
      A(i) = to(j);
      i = i+1;
    elseif A(i) < from(j)
      i = i+1;
    else
      j = j+1;
    end
  end
  
  % return to the original order and shape
  A(indA) = A; 
  A = reshape(A, sizeA);
end
end