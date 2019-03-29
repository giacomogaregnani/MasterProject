clc; clear; close all
%%

hVec = 2.^[0:-1:-5];

M = 1e3;
p = 1;

errloc = zeros(1, M);

for i = 1 : length(hVec)
    
    h = hVec(i);
    
    xFine = -h : h/1000 : 2*h+h;
    gradPhiPrMean = zeros(size(xFine));

    x = [0, h, 2*h];

    gradPhi = 1 / h .* (xFine >= x(1)) .* (xFine <= x(2)) ...
             -1 / h .* (xFine > x(2)) .* (xFine <= x(3));
    
    for j = 1 : M
        X = x + h^p * (rand - 0.5);
        gradPhiPr = 1 / (X(2) - X(1)) .* (xFine >= X(1)) .* (xFine <= X(2)) ...
                - 1 / (X(3) - X(2)) .* (xFine > X(2)) .* (xFine <= X(3));
        gradPhiPrMean = gradPhiPrMean + 1/M * gradPhiPr;
    end
   
    figure
    plot(xFine, gradPhiPrMean);
    hold on
    plot(xFine, gradPhi);
    
    errFunc = zeros(size(xFine));
    for j = 1 : length(xFine)
        if gradPhi(j) == 0
            errFunc(j) = abs(gradPhiPrMean(j));
        else
            errFunc(j) = abs(gradPhiPrMean(j) - gradPhi(j)) / gradPhi(j);
        end
    end
    
    err(i) = trapz(xFine, errFunc);
    
        
end

figure
loglog(hVec, err, 'o-');
hold on
loglog(hVec, hVec.^(p-1), 'k--');