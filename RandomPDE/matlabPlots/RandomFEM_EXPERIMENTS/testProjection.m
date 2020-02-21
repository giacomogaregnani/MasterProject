clc; clear; close all
%%

f = @(x) sin(2*pi*x);
NVec = 2.^(3:12);
hVec = 1 ./ NVec;

M = 1;
p = 1;

errH1MC = zeros(length(NVec), M);
errL2MC = zeros(length(NVec), M);
errH1 = zeros(size(NVec));
errL2 = zeros(size(NVec));

for i = 1 : length(hVec)
    
    x = linspace(0, 1, NVec(i)+1);
    v = f(x);
    gradv = computeGradient_1D(x, v);
    
    for j = 1 : M
        xTilde = [0, x(2:end-1) + hVec(i)^p * (rand(1, NVec(i)-1) - 0.5), 1];
        xTot = sort([xTilde, x(2:end-1)]);
        vTilde = interp1(x, v, xTilde);
        gradvTilde = computeGradient_1D(xTilde, vTilde);        
        for kk = 1 : length(xTot) - 1
            xm = (xTot(kk) + xTot(kk+1)) / 2;
            errH1MC(i, j) = errH1MC(i, j) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradv, xm) - gradFunction(xTilde, gradvTilde, xm))^2;
        end
    end
    
    errH1(i) = mean(sqrt(errH1MC(i, :)));
end

loglog(hVec, errH1, 'ko-');
hold on
loglog(hVec, 10 * hVec.^((p+1)/2), 'k--');