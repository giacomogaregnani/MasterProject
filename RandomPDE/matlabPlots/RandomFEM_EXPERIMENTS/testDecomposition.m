clc; clear; close all
%%

uEx = [];
isFuncEx = true;
field = @(x) x.^0;

% f = @(x) x.^0;
% uEx = @(x) -1 / 2 * x .* (x - 1);
% graduEx = @(x)  0.5 - x;
a = 2;
f = @(x) sin(a * pi * x);
uEx = @(x) 1 / (a^2 * pi^2) * sin(a * pi * x);
graduEx = @(x) 1 / (a * pi) * cos(a * pi * x);
% uEx = @(x) -sin(12*pi*x) .* exp(-(x-0.5).^2 / 0.01);
% graduEx = @(x) exp(-100*(x - 1/2).^2).*sin(12*pi*x).*(200*x - 100) - 12*pi*exp(-100*(x - 1/2).^2).*cos(12*pi*x);
% f = @(x) exp(-100*(x - 1/2).^2).*sin(12*pi*x).*(200*x - 100).^2 - 200*exp(-100*(x - 1/2).^2).*sin(12*pi*x) - ...
%     144*pi^2*exp(-100*(x - 1/2).^2).*sin(12*pi*x) - 24*pi*exp(-100*(x - 1/2).^2).*cos(12*pi*x).*(200*x - 100);
% f = @(x) x.^2 .* (x >= 0.5);
% field = @(x) 1  + (x>=0.2) .* (x<=0.4) + (x>=0.6) .* (x <= 0.8);

if isempty(uEx)
    
    isFuncEx = false;
    
    disp('Computing reference solution...');
    
    xEx = linspace(0, 1, 5e3+1);
    F = assembleRHS(f, xEx);
    A = assembleMatrix(field, xEx);
    
    uExVec = [0; A \ F; 0];
    uEx = @(xx) interp1(xEx, uExVec, xx);
    gradVec = computeGradient_1D(xEx, uExVec);
    graduEx = @(xx) gradFunction(xEx, gradVec, xx);
    
    disp('Computed reference solution');
    
end

%%

p = 1;

NVec = round(2.^(2:9));
hVec = 1 ./ NVec;

M = 10;

errLInfMC = zeros(length(NVec), M);
errLInf = zeros(size(hVec));
errLInfWeak = zeros(size(hVec));
errH1MC = zeros(length(NVec), M);
errH1 = zeros(size(hVec));

xFine = linspace(0, 1, 10000);

for i = 1 : length(hVec)
    
    disp(['iteration ' num2str(i) '/' num2str(length(hVec))])    
    N = NVec(i); h = hVec(i);
    x = linspace(0, 1, N+1);
    meanw = zeros(size(xFine))';
    meanwTilde = zeros(size(xFine))';
    
    for jj = 1 : M
        
        % Build perturbed grid
        xTilde = x;
        for kk = 2 : length(x)-1
            hBar = min(x(kk)-x(kk-1), x(kk+1)-x(kk));
            xTilde(kk) = x(kk) + 0.5 * hBar^p * (rand - 0.5);
        end
        xTot = sort([xTilde, x(2:end-1)]);
        
        % FEM solution on sum space
        FSum = assembleRHS(f, xTot);
        ASum = assembleMatrix(field, xTot);
        uSum = [0; ASum\FSum; 0];
        uSumEval = interp1(xTot, uSum, xFine)';
        gradUSum = computeGradient_1D(xTot, uSum');
        
        % Compute decomposition
        [w, wTilde] = computeDecomposition(uSum, x, xTilde);
        gradw = computeGradient_1D(x, w');
        gradwTilde = computeGradient_1D(xTilde, wTilde');
        wEval = interp1(x, w, xFine)';
        wTildeEval = interp1(xTilde, wTilde, xFine)';
        meanw = meanw + 1/M * wEval;
        meanwTilde = meanwTilde + 1/M * wTildeEval;
        
        % Compute Linf error and H1 error

        errLInfMC(i, jj) = max(abs(wEval - wTildeEval));
        for kk = 1 : length(xFine) - 1
            xm = (xFine(kk) + xFine(kk+1)) / 2;
            errH1MC(i, jj) = errH1MC(i, jj) + (xFine(kk+1) - xFine(kk)) * (gradFunction(x, gradw, xm) - gradFunction(xTilde, gradwTilde, xm))^2;
        end
    end
    
    errLInf(i) = mean(errLInfMC(i, :)); 
    errLInfWeak(i) = max(abs(meanw - meanwTilde));
    errH1(i) = mean(errH1MC(i, :).^(1/2));
    
    figure
    plot(xFine, meanw, xFine, meanwTilde);
    hold on
    plot(xFine, uSumEval, xFine, meanwTilde + meanw);
    legend({'w', '$\tilde w$', '$u^+$', '$w+\tilde w$'}, 'interpreter', 'latex')
end

%% Plot error

figure
loglog(hVec, errLInf, 'o-k')
hold on
loglog(hVec, errLInfWeak, 'x-k')
loglog(hVec, hVec.^((p+1)/2), 'k--')

figure
loglog(hVec, errH1, 'o-k')
hold on
loglog(hVec, hVec.^p, 'k--')

