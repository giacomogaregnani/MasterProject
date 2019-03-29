function A = assembleMatrix(k, x)

N = length(x) - 1;

uDiag = zeros(1, N-2);
lDiag = zeros(1, N-2);
mDiag = zeros(1, N-1);

% First variable
phiUp = @(xx) -k(xx) / (x(3) - x(2))^2;
uDiag(1) = integral(phiUp, x(2), x(3));
phiUp = @(xx) k(xx) / (x(2) - x(1))^2; 
phiDown = @(xx) k(xx) / (x(3) - x(2))^2;
mDiag(1) = integral(phiUp, x(1), x(2)) + integral(phiDown, x(2), x(3));

for i = 2 : N-2
    phiUp = @(xx) -k(xx) / (x(i+2) - x(i+1))^2;
    uDiag(i) = integral(phiUp, x(i+1), x(i+2));
    phiDown = @(xx) -k(xx) / (x(i+1) - x(i))^2;
    lDiag(i-1) = integral(phiDown, x(i), x(i+1));
    phiUp = @(xx) k(xx) / (x(i+1) - x(i))^2;
    phiDown = @(xx) k(xx) / (x(i+2) - x(i+1))^2;
    mDiag(i) = integral(phiUp, x(i), x(i+1)) + integral(phiDown, x(i+1), x(i+2));   
end
% Last variable
phiDown = @(xx) -k(xx) / (x(end-1) - x(end-2))^2;
lDiag(end) = integral(phiDown, x(end-2), x(end-1));
phiUp = @(xx) k(xx) / (x(end-1) - x(end-2))^2; 
phiDown = @(xx) k(xx) / (x(end) - x(end-1))^2;
mDiag(end) = integral(phiUp, x(end-2), x(end-1)) + integral(phiDown, x(end-1), x(end));

uDiag = [0, uDiag];
lDiag = [lDiag, 0];

A = spdiags([lDiag', mDiag', uDiag'], [-1, 0, 1], N-1, N-1);

end