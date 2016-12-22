function plotMultiRV(X, trueVal, xRef)
% Plot the distribution of a multi-valued RV
% Input the RV as a table with each row a realization
% Distribution on diagonals, scatter plots out of diagonals

figure

XFull = X;
X = unique(X, 'rows');

if nargin == 3
    xRefFull = xRef;
    xRef = unique(xRef, 'rows');
    [NRef, ~] = size(xRef);
    IndexRef = 1 : 2 : NRef;
end

[N, dim] = size(X);

AxLim = zeros(dim, 2);
for i = 1 : dim
    AxLim(i, :) = [min(X(:, i))-0.1*(abs(min(X(:, i)))), max(X(:, i))+0.1*(abs(max(X(:, i))))];
end

Index = 1 : 2 : N;


[ha, ~] = tight_subplot(dim, dim, [.04 .02], [.1 .05], [.1 .05]);
if (nargin == 1)
    for i = 1 : dim
        for j = 1 : dim
            axes(ha(dim * (i - 1) + j));
            if j == i
                [f, xi] = ksdensity(XFull(:, j));
                plot(xi, f);
                if i ~= 1
                    set(gca, 'YTickLabel', '')
                end
                if i ~= dim
                    set(gca, 'XTickLabel', '')
                end
                xlim(AxLim(i, :));
                
            else
                plot(X(Index, j), X(Index, i), '.', 'MarkerSize', 3)
                hold on
                if j ~= 1
                    set(gca, 'YTickLabel', '')
                end
                if i ~= dim
                    set(gca, 'XTickLabel', '')
                end
                xlim(AxLim(j, :))
                ylim(AxLim(i, :))
                
            end
        end
    end
end

if (nargin == 2)
    for i = 1 : dim
        for j = 1 : dim
            axes(ha(dim * (i - 1) + j));
            if j == i
                [f, xi] = ksdensity(XFull(:, j));
                plot(xi, f);
                if i ~= 1
                    set(gca, 'YTickLabel', '')
                end
                if i ~= dim
                    set(gca, 'XTickLabel', '')
                end
                xlim(AxLim(i, :));
                
            else
                plot(X(Index, j), X(Index, i), '.', 'MarkerSize', 3)
                hold on
                if j ~= 1
                    set(gca, 'YTickLabel', '')
                end
                if i ~= dim
                    set(gca, 'XTickLabel', '')
                end
                plot(trueVal(j), trueVal(i), 'go', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
                xlim(AxLim(j, :))
                ylim(AxLim(i, :))
                
            end
        end
    end
end


if (nargin == 3 && isempty(trueVal) == 0)
    for i = 1 : dim
        for j = 1 : dim
            axes(ha(dim * (i - 1) + j));
            if j == i
                [f, xi] = ksdensity(XFull(:, j));
                plot(xi, f);
                hold on
                [f, xi] = ksdensity(xRefFull(:, j));
                plot(xi, f);
                if i ~= 1
                    set(gca, 'YTickLabel', '')
                end
                if i ~= dim
                    set(gca, 'XTickLabel', '')
                end
                xlim(AxLim(i, :));
                set(gca, 'FontSize', 10)
            else
                plot(X(Index, j), X(Index, i), '.', 'MarkerSize', 1)
                hold on
                plot(xRef(IndexRef, j), xRef(IndexRef, i), '.', 'MarkerSize', 1)
                if j ~= 1
                    set(gca, 'YTickLabel', '')
                end
                if i ~= dim
                    set(gca, 'XTickLabel', '')
                end
                plot(trueVal(j), trueVal(i), 'go', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 4)
                xlim(AxLim(j, :))
                ylim(AxLim(i, :))
                set(gca, 'FontSize', 10)
            end
        end
    end
end

if (nargin == 3 && isempty(trueVal))
    for i = 1 : dim
        for j = 1 : dim
            axes(ha(dim * (i - 1) + j));
            if j == i
                [f, xi] = ksdensity(XFull(:, j));
                plot(xi, f);
                hold on
                [f, xi] = ksdensity(xRefFull(:, j));
                plot(xi, f);
                if i ~= 1
                    set(gca, 'YTickLabel', '')
                end
                if i ~= dim
                    set(gca, 'XTickLabel', '')
                end
                xlim(AxLim(i, :));
                
            else
                plot(X(Index, j), X(Index, i), '.', 'MarkerSize', 3)
                hold on
                plot(xRef(IndexRef, j), xRef(IndexRef, i), '.', 'MarkerSize', 3)
                if j ~= 1
                    set(gca, 'YTickLabel', '')
                end
                if i ~= dim
                    set(gca, 'XTickLabel', '')
                end
                xlim(AxLim(j, :))
                ylim(AxLim(i, :))
                
            end
        end
    end
end