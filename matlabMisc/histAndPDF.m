function [h, l] = histAndPDF(data, varargin)

hold on
argsHist = cell(0);
I = find(strcmp(varargin, 'NumBins'), 1);
if ~isempty(I)
   argsHist = [argsHist, varargin{I:I+1}];
end
I = find(strcmp(varargin, 'FaceColor'), 1);
if ~isempty(I)
   argsHist = [argsHist, varargin{I:I+1}];
end
I = find(strcmp(varargin, 'FaceAlpha'), 1);
if ~isempty(I)
   argsHist = [argsHist, varargin{I:I+1}];
end
h = histogram(data, 'Normalization', 'pdf', 'EdgeColor', 'none', argsHist{:});

argsLine = cell(0);
I = find(strcmp(varargin, 'LineStyle'), 1);
if ~isempty(I)
   argsLine = [argsLine, varargin{I+1}];
end
I = find(strcmp(varargin, 'LineColor'), 1);
if ~isempty(I)
   argsLine = [argsLine, 'color', varargin{I+1}];
end
I = find(strcmp(varargin, 'LineWidth'), 1);
if ~isempty(I)
   argsLine = [argsLine, varargin{I:I+1}];
end

[dens, xi] = ksdensity(data);
l = plot(xi, dens, argsLine{:});
maxDens = max(dens);