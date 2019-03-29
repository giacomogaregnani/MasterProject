function [u, uTilde] = computeDecomposition(uSum, x, X)
% Computes the unique decomposition of a function in the FE sum space into the two
% components in the two FE summands. 

hatUp   = @(xx, xL, xR) (xx - xL) / (xR - xL);
hatDown = @(xx, xL, xR) (xR - xx) / (xR - xL);

N = length(x) - 1;
A = sparse(2*(N-1), 2*(N-1));

uTilde = 0; 
u = uTilde;

if x(2) < X(2)
    A(1, 1) = 1;
    A(1, 2) = hatUp(x(2), X(1), X(2));
    A(2, 1) = hatDown(X(2), x(2), x(3));
    A(2, 2) = 1;
    A(2, 3) = hatUp(X(2), x(2), x(3));
else
    A(1, 1) = hatUp(X(2), x(1), x(2));
    A(1, 2) = 1;
    A(2, 1) = 1; 
    A(2, 2) = hatDown(x(2), X(2), X(3));
    A(2, 4) = hatUp(x(2), X(2), X(3));
end

for i = 3 : N-1
    
    idx = 2*i-3;
    
    if x(i) < X(i)
        A(idx, idx)   = 1;
        A(idx, idx+1) = hatUp(x(i), X(i-1), X(i));
        A(idx, idx-1) = hatDown(x(i), X(i-1), X(i));
        A(idx+1, idx+1) = 1; 
        A(idx+1, idx+2) =  hatUp(X(i), x(i), x(i+1));
        A(idx+1, idx) =  hatDown(X(i), x(i), x(i+1));
    else
        A(idx, idx) = hatUp(X(i), x(i-1), x(i));
        A(idx, idx+1) = 1;
        A(idx, idx-2) = hatDown(X(i), x(i-1), x(i));
        A(idx+1, idx) = 1;
        A(idx+1, idx+1) = hatDown(x(i), X(i), X(i+1));
        A(idx+1, idx+3) = hatUp(x(i), X(i), X(i+1));
    end    
    
end

if x(end-1) < X(end-1)
    A(end-1, end-2) = hatDown(x(end-1), X(end-2), X(end-1));
    A(end-1, end-1) = 1;
    A(end-1, end) = hatUp(x(end-1), X(end-2), X(end-1));
    A(end, end-1) = hatDown(X(end-1), x(end-1), x(end)); 
    A(end, end) = 1;
else
    A(end-1, end-3) = hatDown(X(end-1), x(end-2), x(end-1));
    A(end-1, end-1) = hatUp(X(end-1), x(end-2), x(end-1));
    A(end-1, end) = 1;
    A(end, end-1) = 1;
    A(end, end) = hatDown(x(end-1), X(end-1), X(end)); 
end

RHS = uSum(2:end-1);

sol = A \ RHS;

u = [0; sol(1:2:end-1); 0];
uTilde = [0; sol(2:2:end); 0];

end

