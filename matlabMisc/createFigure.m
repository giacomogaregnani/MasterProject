function fig = createFigure(W, H, varargin)

set(0,'DefaultTextFontname', 'CMUSerif')
set(0,'DefaultAxesFontName', 'CMUSerif')

% W, H in centimeter

p = inputParser;
addParamValue(p, 'figvarargin', {});
addParamValue(p, 'enhanced', 1);
parse(p,varargin{:});

figvarargin = p.Results.figvarargin;
enhanced = p.Results.enhanced;


fig = figure('units', 'centimeter', 'position', ...
    [0 0 enhanced*W enhanced*H], 'color', 'w', 'PaperPositionMode', 'auto', ...
    figvarargin{:});
sc  = get(0, 'screensize');
set(fig, 'units','pixels');
fp = get(fig, 'position');
set(fig, 'position', [(sc(3)/2-fp(3))/2, (sc(4)-fp(4))/2, fp(3:4)]);

end
