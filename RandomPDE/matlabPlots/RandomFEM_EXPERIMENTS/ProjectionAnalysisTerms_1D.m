clc; clear; close all
%%

f = @(x) sin(2*pi*x);
hVec = 2.^[-1:-1:-7];

M = 1e2;
p = 3;

for i = 1 : length(hVec)
    
    h = hVec(i);
    x = [0.2, 0.2+h, 0.2+2*h];
    xFine = x(1)-h : 0.01*h^p : x(3)+h;
    
    errlocH1 = zeros(M, 1);
    errlocL2 = zeros(M, 1);
    
    F = f(x(2));
    phi = (xFine - x(1)) / h .* (xFine >= x(1)) .* (xFine <= x(2)) ...
         +(x(3) - xFine) / h .* (xFine > x(2))  .*  (xFine <= x(3));   
    gradPhi = 1 / h .* (xFine >= x(1)) .* (xFine <= x(2)) ...
             -1 / h .* (xFine > x(2))  .* (xFine <= x(3));
    
    parfor j = 1 : M
        X = x + h^p * (rand - 0.5);
        
        FPr = f(X(2));
        phiPr = (xFine - X(1)) / h .* (xFine >= X(1)) .* (xFine <= X(2)) ...
               +(X(3) - xFine) / h .* (xFine > X(2))  .*  (xFine <= X(3));    
        gradPhiPr = 1 / (X(2) - X(1)) .* (xFine >= X(1)) .* (xFine <= X(2)) ...
                   -1 / (X(3) - X(2)) .* (xFine > X(2))  .* (xFine <= X(3));
               
        errlocH1(j) = trapz(xFine, (F * gradPhi - FPr * gradPhiPr).^2);
        errlocL2(j) = trapz(xFine, (F * phi - FPr * phiPr).^2);
    end
    
    errH1(i) = mean(errlocH1);
    errL2(i) = mean(errlocL2);

end

loglog(hVec, errH1, 'o-');
hold on
loglog(hVec, errL2, 'o-');
loglog(hVec, hVec.^(p-2), 'k--');
loglog(hVec, hVec.^(p+1), 'k--');