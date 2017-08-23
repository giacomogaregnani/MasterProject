clc; clear; close all
% Just a test to see if everything is all right

% Problem data
kappa = @(x) 1;
f = @(x) sin(2 * pi * x)';

% Exact solution
uEx = @(x) 1 / (2 * pi)^2 * sin(2 * pi * x);% + x;
uExDer = @(x) 1 / (2 * pi) * cos(2 * pi * x);% + 1;

N = 3 * 2.^[1 : 8];
mesh = struct('N', [], 'xMin', 0, 'xMax', 1);

L2ErrVec = [];
L2ErrVecProb = [];
L2ErrVecDetProb = [];
VecErr = [];
H1ErrVec = [];
H1ErrVecProb = [];

rBC = 0;
nMC = 100;

k = 1;
L2ErrVecProb = cell(3, 1);

% Figure settings
W = 6;
H = 6;
enhanced = 1;
fontsizelb = 'normalsize';
fontsizetk = 'small';
fontsizeLAB = getLatexTextSize(fontsizelb, 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize(fontsizetk, 'enhanced', enhanced);

for p = [1 1.5 2]
    
    q = 1;
    for n = N
          
        display(n)
        
        mesh.N = n;
        mesh = buildMesh(mesh);
        
        if p == 1 && (n == N(1) || n == N(2) || n == N(3) || n == N(4))
            
            fig = createFigure(W, H, 'enhanced',enhanced);
            hold on
            for j = 1 : nMC
                meshProb = mesh;
                meshProb.xInt = meshProb.xInt + ...
                    meshProb.h^p * (-0.5 + rand(size(meshProb.xInt)));
                uProb = solveFwdProblemProb(meshProb, kappa, f, rBC);
                exMesh = uEx(mesh.x)';
                plot(mesh.x, uProb, 'color', [0.7, 0.7, 0.7]);
            end
            
            u = solveFwdProblem(mesh, kappa, f, rBC);
            plot(mesh.x, u, 'k', 'LineWidth', 1)
            
            set(gca, 'ytick', [])
            set(gca, 'xtick', [])
            
            axis square
            box on
            
            filename =  ['../../../../Reports/RandMesh/VERSION2/solution', num2str(q), '.eps'];
            print('-depsc2', filename)
            
            q = q + 1;
       else
            meshProb = mesh;
            meshProb.xInt = meshProb.xInt + ...
                meshProb.h^p * (-0.5 + rand(size(meshProb.xInt)));
            uProb = solveFwdProblemProb(meshProb, kappa, f, rBC);
            exMesh = uEx(mesh.x)';
        end
        
        
        quadPoints = linspace(mesh.xMin, mesh.xMax, 5*n);
        L2ErrProb = sqrt(trapz(quadPoints, ...
            (interp1(mesh.x, exMesh, quadPoints) - interp1(mesh.x, uProb, quadPoints)).^2));
        L2ErrVecProb{k} = [L2ErrVecProb{k} L2ErrProb];
        
    end
    
    k = k + 1;
end

%%
W = 12; H = 8;
fig = createFigure(W, H, 'enhanced',enhanced);

loglog(1./N, L2ErrVecProb{1}, 'ko-')
hold on
loglog(1./N, L2ErrVecProb{2}, 'k+-')
loglog(1./N, L2ErrVecProb{3}, 'k>-')

loglog(1./N, 0.5 * 1e-1 * 1./N, 'k--.')
loglog(1./N, 0.5 * 1e-1 * 1./N.^1.5, 'k-.')
loglog(1./N, 1e-1 * 1./N.^2, 'k--')
legend('p = 1', 'p = 1.5', 'p = 2', 'slope 1', 'slope 1.5', 'slope 2', 'Location', 'best')

axis tight

xlabel('$h$', 'interpreter', 'LaTex');
ylabel('$L^2$ error', 'interpreter', 'LaTeX');

set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);

filename =  ['../../../../Reports/RandMesh/VERSION2/L2Convergence.eps'];
print('-depsc2', filename)

% figure
% loglog(1./N, H1ErrVec, 'o-')
% hold on
% loglog(1./N, H1ErrVecProb, 'o-')
% loglog(1./N, 1./N, 'k--')
% loglog(1./N, 1./N.^(p-1), 'k')
% legend('det', 'prob', 'h', 'h^{p-1}')

