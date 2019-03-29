clc; clear; close all
%%

uEx = [];
isFuncEx = true;
field = @(x) x.^0;

% f = @(x) x.^0;
% uEx = @(x) -1 / 2 * x .* (x - 1);
% graduEx = @(x)  0.5 - x;
a = 6;
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

errH1Ex = zeros(1, length(NVec));
errH1SumEx = zeros(1, length(NVec));
errH1Loc = zeros(M, length(hVec));
errH1ProjLoc = zeros(M, length(hVec));
errH1SumLoc = zeros(M, length(hVec));
errH1SumTildeLoc = zeros(M, length(hVec));
errH1DecompLoc = zeros(M, length(hVec));
errH1DecompProdLoc = zeros(M, length(hVec));
errH1DecompProdLoc2 = zeros(M, length(hVec));
errH1DecompProjLoc = zeros(M, length(hVec));
errH1DecompProjTildeLoc = zeros(M, length(hVec));
errH1DecompProjProjLoc = zeros(M, length(hVec));
errH1SumExLoc = zeros(M, 1);
errTest = zeros(M, length(hVec));
errTest2 = zeros(M, length(hVec));
errTest3 = zeros(M, length(hVec));

for i = 1 : length(hVec)
    
    N = NVec(i); h = hVec(i);
    x = linspace(0, 1, N+1);
    
    F = assembleRHS(f, x);
    A = assembleMatrix(field, x);
    u = [0; A \ F; 0];
    gradu = computeGradient_1D(x, u');
    
    meanw = zeros(size(x));
    meanwTilde = zeros(size(x));
    meanz = zeros(size(x));
    
    parfor jj = 1 : M
        
        % Build perturbed grid
        xTilde = x;
        
        for kk = 2 : length(x)-1
            hBar = min(x(kk)-x(kk-1), x(kk+1)-x(kk));
            xTilde(kk) = x(kk) + 0.2 * hBar^p * (rand - 0.5);
        end
        xTot = sort([xTilde, x(2:end-1)]);
        
        % FEM solution on perturbed grid
        FTilde = assembleRHS(f, xTilde);
        ATilde = assembleMatrix(field, xTilde);
        uTilde = [0; ATilde\FTilde; 0];
        graduTilde = computeGradient_1D(xTilde, uTilde');
        PiU = interp1(xTilde, uTilde, x);
        gradPiU = computeGradient_1D(x, PiU);
        PiTildeu = interp1(x, u, xTilde);
        gradPiTildeu = computeGradient_1D(xTilde, PiTildeu);
        
        % FEM solution on sum space
        FSum = assembleRHS(f, xTot);
        ASum = assembleMatrix(field, xTot);
        uSum = [0; ASum\FSum; 0];
        gradUSum = computeGradient_1D(xTot, uSum');
        
        % Temporary test
        PiUSum = interp1(xTot, uSum, x)';
        gradPiUSum = computeGradient_1D(x, PiUSum');
        PiTildeUSum = interp1(xTot, uSum, xTilde)';
        gradPiTildeUSum = computeGradient_1D(xTilde, PiTildeUSum');
        
        % Compute decomposition
        [w, wTilde] = computeDecomposition(uSum, x, xTilde);     
        meanw = meanw + 1/M * w;
        meanwTilde = meanwTilde + 1/M * wTilde;        
        gradw = computeGradient_1D(x, w');
        gradwTilde = computeGradient_1D(xTilde, wTilde');
        wTot = interp1(x, w, xTot)';
        wTildeTot = interp1(xTilde, wTilde, xTot)';
        PiwTilde = interp1(xTilde, wTilde, x);
        gradPiwTilde = computeGradient_1D(x, PiwTilde);
        PiTildew = interp1(x, w, xTilde);
        gradPiTildew = computeGradient_1D(xTilde, PiTildew);
        
        ratio = wTot ./ uSum;
        ratio = ratio(2:end-1);
%         gamma = mean(ratio(abs(ratio) < 10));
        gamma = 1/2;

        % Compute difference between stuff and other stuff
        z = w - gamma * u;
        meanz = meanz + 1/M * z;
        gradz = computeGradient_1D(x, z');
        zTilde = wTilde - (1 - gamma) * uTilde;
        gradzTilde = computeGradient_1D(xTilde, zTilde');
        PiTildez = interp1(x, z, xTilde);
        gradPiTildez = computeGradient_1D(xTilde, PiTildez);
        PizTilde = interp1(xTilde, zTilde, x);
        gradPizTilde = computeGradient_1D(x, PizTilde);
%         display(num2str(norm(zTilde) / norm(U)));
%         display(num2str(norm(z) / norm(u)));
        
        % Exact error computations
        eL2I = @(xx) interp1(x, z, xx)^2;
        for kk = 1 : length(xTot) - 1
            xm = (xTot(kk) + xTot(kk+1)) / 2;
            errH1Loc(jj, i) = errH1Loc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradu, xm) - gradFunction(xTilde, graduTilde, xm))^2;
            errH1ProjLoc(jj, i) = errH1ProjLoc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradu, xm) - gradFunction(xTilde, gradPiTildeu, xm))^2;
            errH1SumLoc(jj, i) = errH1SumLoc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradu, xm) - gradFunction(xTot, gradUSum, xm))^2;
            errH1SumTildeLoc(jj, i) = errH1SumTildeLoc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(xTilde, graduTilde, xm) - gradFunction(xTot, gradUSum, xm))^2;
            
            errH1DecompLoc(jj, i) = errH1DecompLoc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradw, xm) - gradFunction(xTilde, gradwTilde, xm))^2;
            errH1DecompProjLoc(jj, i) = errH1DecompProjLoc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradw, xm) - gradFunction(xTilde, gradPiTildew, xm))^2;
            errH1DecompProjTildeLoc(jj, i) = errH1DecompProjTildeLoc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(xTilde, gradwTilde, xm) - gradFunction(x, gradPiwTilde, xm))^2;
            errH1DecompProjProjLoc(jj, i) = errH1DecompProjProjLoc(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(xTilde, gradPiTildew, xm) - gradFunction(x, gradPiwTilde, xm))^2;
            errH1DecompProdLoc(jj, i) = errH1DecompProdLoc(jj, i) + ...
                (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradu, xm) - gradFunction(xTilde, graduTilde, xm)) * (gradFunction(x, gradw, xm) - gradFunction(xTilde, gradwTilde, xm));
            errH1DecompProdLoc2(jj, i) = errH1DecompProdLoc2(jj, i) + ...
                (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradu, xm) - gradFunction(xTilde, graduTilde, xm)) * (gradFunction(x, gradz, xm) - gradFunction(xTilde, gradzTilde, xm));    
            
            errTest(jj, i) = errTest(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(x, gradz, xm) - gradFunction(xTilde, gradPiTildez, xm))^2;
            errTest2(jj, i) = errTest2(jj, i) + (xTot(kk+1) - xTot(kk)) * (gradFunction(xTilde, gradzTilde, xm) - gradFunction(x, gradPizTilde, xm))^2;
            errTest3(jj, i) = errTest3(jj, i) + (xTot(kk+1) - xTot(kk)) * gradFunction(x, gradz, xm)^2;            
%             errTest3(jj, i) = errTest3(jj, i) + (xTot(kk+1) - xTot(kk))/6 * (eL2I(xTot(kk)) + 4*eL2I(xm) + eL2I(xTot(kk+1)));
        end
        
        errH1SumExLoc(jj) = integral(@(xx) (graduEx(xx) - gradFunction(xTot, gradUSum, xx)).^2, 0, 1);
        
    end
    
%     figure
%     plot(x, meanw, x, meanwTilde);
%     figure
%     plot(x, meanz);
    
    errH1SumEx(i) = sqrt(mean(errH1SumExLoc));
    
    % Compute exact errors
    if isFuncEx
        errH1Ex(i)  = sqrt(integral(@(xx) (graduEx(xx) - gradFunction(x, gradu, xx)).^2, 0, 1));
    else
        for kk = 1 : length(xEx) - 1
            xm = mean([xEx(kk), xEx(kk+1)]);
            errH1Ex(i) = errH1Ex(i) + (xEx(kk+1) - xEx(kk)) * (gradFunction(x, gradu, xm) - graduEx(xm))^2;
        end
        errH1Ex(i) = sqrt(errH1Ex(i));
    end
    
    disp(['iteration ' num2str(i) '/' num2str(length(hVec))])
    
end

if M > 1
    errH1 = sqrt(mean(errH1Loc));
    errH1Proj = sqrt(mean(errH1ProjLoc));
    errH1Sum = sqrt(mean(errH1SumLoc));
    errH1SumTilde = sqrt(mean(errH1SumTildeLoc));
    errH1Decomp = sqrt(mean(errH1DecompLoc));
    errH1DecompProd = mean(errH1DecompProdLoc);
    errH1DecompProd2 = mean(errH1DecompProdLoc2);
    errH1DecompProj = sqrt(mean(errH1DecompProjLoc));
    errH1DecompProjTilde = sqrt(mean(errH1DecompProjTildeLoc));
    errH1DecompProjProj = sqrt(mean(errH1DecompProjProjLoc));
    errTest = sqrt(mean(errTest));
    errTest2 = sqrt(mean(errTest2));
    errTest3 = sqrt(mean(errTest3));
else
    errH1 = sqrt(errH1Loc);
    errH1Proj = sqrt(errH1ProjLoc);
    errH1Sum = sqrt(errH1SumLoc);
    errH1SumTilde = sqrt(errH1SumTildeLoc);
    errH1Decomp = sqrt(errH1DecompLoc);
    errH1DecompProd = errH1DecompProdLoc;
    errH1DecompProd2 = errH1DecompProdLoc2;
    errH1DecompProj = sqrt(errH1DecompProjLoc);
    errH1DecompProjTilde = sqrt(errH1DecompProjTildeLoc);
    errH1DecompProjProj = sqrt(errH1DecompProjProjLoc);
    errTest = sqrt(errTest);
    errTest2 = sqrt(errTest2);
    errTest3 = sqrt(errTest3);
end

%%

enhanced = 2;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', enhanced);

W = 6.7; H = 6.7;

% H1 error
fig = createFigure(W, H, 'enhanced',enhanced);
loglog(hVec, errH1Ex, 'ks-', 'markersize', 4*enhanced);
hold on
loglog(hVec, errH1SumEx, 'kd-', 'markersize', 4*enhanced);
loglog(hVec, errH1, 'k<-', 'markersize', 4*enhanced);
loglog(hVec, errH1Sum, 'ko-', 'markersize', 4*enhanced);
loglog(hVec, errH1SumTilde, 'ko--', 'markersize', 4*enhanced);
loglog(hVec, errH1Decomp, 'k+-', 'markersize', 4*enhanced);
loglog(hVec, errH1DecompProj, 'ks--', 'markersize', 4*enhanced);
loglog(hVec, errH1DecompProjTilde, 'kd--', 'markersize', 4*enhanced);
loglog(hVec, errH1DecompProjProj, 'k+--', 'markersize', 4*enhanced);
loglog(hVec, errTest, 'bo-')
loglog(hVec, errTest2, 'ro-')
loglog(hVec, errTest3, 'co-')
legend({'$|u_h - u|_{H^1}$', '$E |u_h^+ - u|_{H^1}$', '$E|u_h - \tilde u_h|_{H^1}$', ...
        '$E |u_h^+ - u_h|_{H^1}$', '$E |u_h^+ - \tilde u_h|_{H^1}$', '$E |w_h - \tilde w_h|_{H^1}$', ...
        '$E |w_h - \tilde \Pi_h w_h|_{H^1}$', '$E |\Pi_h \tilde w_h - \tilde w_h|_{H^1}$', '$E |\Pi_h \tilde w_h - \tilde \Pi_hw_h|_{H^1}$'}, ...
       'location', 'se', 'interpreter', 'latex', 'fontsize', fontsizeTICK)
xlabel('$h$', 'interpreter', 'latex')
ylabel('$H^1$ error', 'interpreter', 'latex')

% Another figure
fig = createFigure(W, H, 'enhanced',enhanced);
loglog(hVec, errH1Sum.^2 + errH1SumTilde.^2, 'ks-')
hold on
loglog(hVec, errH1.^2, 'kd--')
loglog(hVec, errH1.^2 + errTest2.^2 + errTest.^2, 'ko--')
loglog(hVec, errTest.^2 + errTest2.^2, 'k<--')
legend('LHS', 'H1', 'bound', 'interp', 'location', 'se')


% A posteriori constants
fig = createFigure(W, H, 'enhanced',enhanced);
loglog(hVec, (errH1Sum.^2 + errH1SumTilde.^2) ./ errH1.^2, 'ks-', 'markersize', 4*enhanced);
hold on
loglog(hVec, (errH1Sum.^2 + errH1SumTilde.^2) ./ errH1DecompProd, 'kd-', 'markersize', 4*enhanced);
loglog(hVec, (errH1Sum.^2 + errH1SumTilde.^2) ./ (0.5 * errH1.^2 + errH1DecompProd2), 'k+-', 'markersize', 4*enhanced);
% loglog(hVec, errH1Sum.^2 ./ bound, 'ko-');
xlabel('$h$', 'interpreter', 'latex')
% legend({'$\beta$', '$\gamma$'}, 'location', 'nw', 'interpreter', 'latex', 'fontsize', fontsizeTICK)