function A = stack_blkdiag(A,dim)
if dim == 1
  return;
elseif dim == 2,
  A = blkdiag(A,A);
elseif dim == 3,
  A = blkdiag(A,A,A);
else
  error('stack_blkdiag only implemented for dim = 1:3');
end
end