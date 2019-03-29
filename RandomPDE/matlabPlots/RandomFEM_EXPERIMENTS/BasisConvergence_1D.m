clc; clear; close all
%%

hVec = 2.^[-1:-1:-7];

M = 1e1;

p = 2;

for i = 1 : length(hVec)
    
    h = hVec(i);
    x = [0, h, 2*h];
    xFine = x(1)-h : 0.1*h^p : x(3)+h;
    
    errlocH1 = zeros(M, length(xFine));
    errlocL2 = zeros(M, length(xFine));
    
    phi = (xFine - x(1)) / h .* (xFine >= x(1)) .* (xFine <= x(2)) ...
         +(x(3) - xFine) / h .* (xFine > x(2))  .* (xFine <= x(3));   
    gradPhi = 1 / h * (xFine >= x(1)) .* (xFine <= x(2)) ...
             -1 / h * (xFine > x(2))  .*  (xFine <= x(3));
    
    parfor j = 1 : M
        X = x + 0.5 * h^p * (rand - 0.5);
        phiPr = (xFine - X(1)) / h .* (xFine >= X(1)) .* (xFine <= X(2)) ...
               +(X(3) - xFine) / h .* (xFine > X(2))  .* (xFine <= X(3));    
        gradPhiPr = 1 / (X(2) - X(1)) .* (xFine >= X(1)) .* (xFine <= X(2)) ...
                   -1 / (X(3) - X(2)) .* (xFine > X(2))  .* (xFine <= X(3));
        errlocH1(j, :) = (gradPhi - gradPhiPr).^2;
        errlocL2(j, :) = (phi - phiPr).^2;
    end
    
    errH1(i) = trapz(xFine, mean(errlocH1));
    errL2(i) = trapz(xFine, mean(errlocL2));

end

loglog(hVec, errH1, 'o-');
hold on
loglog(hVec, errL2, 'o-');
loglog(hVec, hVec.^(p-2), 'k--');
loglog(hVec, 1e-1*hVec.^(2*p-1), 'k--');