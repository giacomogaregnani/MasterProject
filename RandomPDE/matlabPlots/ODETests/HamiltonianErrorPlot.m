clc; clear; close all
addpath('resultsHamiltonian2')

%%
p0 = 1.5; q0 = -pi;
H0 = p0^2/2 - cos(q0);

fileBase = 'newPlot';
fileNames = { [fileBase, '200.0.txt'], ...
    [fileBase, '100.0.txt'], ...
    [fileBase, '50.0.txt'], ...
    [fileBase, '25.0.txt']};

fileNamesDet = {[fileBase, '200.0det.txt'], ...
    [fileBase, '100.0det.txt'], ...
    [fileBase, '50.0det.txt'], ...
    [fileBase, '25.0det.txt']};

h = [0.2, 0.1, 0.05, 0.025];

for i = 1 : length(fileNames)
    data = dlmread(fileNames{i});
    T{i} = data(:, 1);
    err{i} = data(:, 2);
    %     err{i} = mean(abs(data(:, 2:end) - H0), 2);
    data = dlmread(fileNamesDet{i});
    TDet{i} = data(:, 1);
    %     errDet{i} = abs(data(:, 2) - H0);
    errDet{i} = data(:, 2);
end

%%
q = 2;
p = 2;
%
enhanced = 1;

fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('footnotesize', 'enhanced', enhanced);
W = 12; H = 5;

linestyles = {'-', '--', '-.', ':'};

% figure
createFigure(W, H, 'enhanced',enhanced);

for i = 1 : length(fileNames)
      
    % Compute constants (only for first iteration, they are constants...)
    if i == 1
        M = max(errDet{i}); % Max error fixed h
        
%         errCopy = err{i};
%         tCopy = T{i};
%         found = false;
%         while ~found
%             idx = find(errCopy > M, 1);
%             tAttempt = tCopy(idx);
%             if isempty(tAttempt)
%                 tBar(i) = T{i}(end);
%                 found = true;
%             else
%                 found = true;
%                 for j = idx : idx+100
%                     found = found * (errCopy(j) > M);
%                 end
%                 if found
%                     tBar(i) = tCopy(idx);
%                 end
%                 errCopy = errCopy(idx+100:end);
%                 tCopy = tCopy(idx+100:end);
%             end
%         end
        c1 = M / h(i)^q; % Compute C1
%         c2 = c1 / (h(i)^p * tBar(i)^0.5) % Compute C2
%         c3 = c1 / (h(i)^(2*p-1) * tBar(i))
    end
    
    % Plot error
    loglog(T{i}(1:2*i:end), err{i}(1:2*i:end));
    hold on
    %     loglog(T{i}(1:2*i:end), errDet{i}(1:2*i:end));
    xlim([1, max(TDet{i})])
    %     text(TDet{end}(3)/1.5, 0.6*max(errDet{i}), ['$h = ', num2str(h(i)), '$'], 'interpreter', 'laTeX')

    % Plot error estimate
    loglog(T{end}, c1*h(i)^q*ones(size(T{end})), 'k', 'linestyle', linestyles{i});
%         loglog(T{end}, c1*h(i)^q + c2 * T{end}.^(0.5) * h(i)^(p+q) + c3 * T{end} * h(i)^(2*p+q-1), 'k', 'linestyle', linestyles{i});

    % Plot reference slope
    tPlot = [2e2, 1e6];
    loglog(tPlot, 3e-3*tPlot.^(1/2), 'k')
    Text05 = text(1.4e3, 0.2, 'slope $1/2$', 'interpreter', 'latex', 'fontsize', fontsizeTICK);
    set(Text05,'Rotation', 13);
end

a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', fontsizeTICK)
set(gca, 'XTick', [1e1 1e2 1e3 1e4 1e5 1e6])
set(gca, 'XTickLabel', {'10^{1}', '10^{2}', '10^3', '10^4', '10^5', '10^6'})
a = get(gca,'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', fontsizeTICK)

xlabel('$t_n$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)
ylabel('$E|Q(Y_n) - Q(y_0)|$', 'interpreter', 'laTeX', 'fontsize', fontsizeLAB)

print -depsc2 ../../../Reports/PaperRTSRK_18/VERSION15/HamiltonianError2.eps



