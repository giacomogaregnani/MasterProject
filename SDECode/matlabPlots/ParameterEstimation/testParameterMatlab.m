clc; clear; close all
%%

solutions = dlmread('testParamSol.txt');
x = dlmread('testParam.txt');
hom = x(1, 2:3);

T = 100;
N = size(solutions, 2)-1;
h = T / N;
t = linspace(0, T, N+1);

nEps = size(solutions, 1);

Sigma = hom(2);

mu = hom(1);
sigmaPrior = 0.1;

pX = @(x1, x2, a) -1/(2*Sigma*h) * (x2 - (1-a*h)*x1)^2;
prior = @(a) -1/(2*sigmaPrior^2) * (a-mu).^2;

aVec = linspace(-1, 3, 1000);
nA = length(aVec);

for i = 3
    
    plot(t, solutions(i, :))
    
    X = solutions(i, :);
    
    l = zeros(1, nA);
    post = zeros(1, nA);
    for k = 1 : nA
        for j = 1 : N
            l(k) = l(k) + pX(X(j), X(j+1), aVec(k));
            post(k) = post(k) + pX(X(j), X(j+1), aVec(k));
        end
        post(k) = post(k) + prior(aVec(k));
    end
    
    figure
    hold on
    plot(aVec, l, aVec, post)
    
    [mL, IL] = max(l); 
    [minL, ~] = min(post);
    [mP, IP] = max(post);     
    
    aVec(IL)
    x(i+1, 2)
    
    aVec(IP)
    x(i+1, 4)
    
%     plot([x(i, 2), x(i, 2)], [minL, mL], '--k') 
%     plot([x(i, 4), x(i, 4)], [minL, mP], '-k')
    
end



