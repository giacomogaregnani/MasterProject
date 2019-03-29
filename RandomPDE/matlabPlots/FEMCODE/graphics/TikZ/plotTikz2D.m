function plotTikz2D(mesh, sol, opt)
%PLOTTIKZ2D Summary of this function goes here
if nargin < 3, opt = struct; end
if ~isfield(opt,'colormap'), opt.colormap = CMRmap(256); end
if ~isfield(opt,'border'), opt.border = 0; end
if ~isfield(opt,'color'), opt.color = 'white'; end
if ~isfield(opt,'figfile'), error('no figure file provided'); end
if ~isfield(opt,'tikzfile'), error('no tikz file provided'); end
if ~isfield(opt,'height'), opt.height = 2.5; end

% DEPERIODIZE
[mesh, updates] = deperiodize(mesh);
if isnumeric(sol)
	sol = sol(updates);
else
	sol = sol(mesh.node);
end

% PICTURE BOUNDS
xmin = min(mesh.node(:,1)) - opt.border;
xmax = max(mesh.node(:,1)) + opt.border;
xspan = xmax - xmin;
ymin = min(mesh.node(:,2)) - opt.border;
ymax = max(mesh.node(:,2)) + opt.border;
yspan = ymax - ymin;
maxspan = max(xspan, yspan);
xrelspan = xspan / maxspan;
yrelspan = yspan / maxspan;

% PLOT
simpplot_sol(mesh, sol);

% SET COLORMAP
colormap(opt.colormap);

% SET LIMITS OF THE AXES
xlim([xmin,xmax]);
ylim([ymin,ymax]);

% SET FIGURE - GENERAL PROPERTIES
set(gcf,'color',opt.color);

% SET SIZES
set(gcf,'pos',[100,100,250,250]);
apos = [0.2,0.2,xrelspan*0.6, yrelspan*0.6];
set(gca, 'pos', apos);

% MINIMUM AND MAXIMUM VALUES OF THE SOLUTION
valMin = min(sol); 
valMax = max(sol);
valSpan = valMax-valMin;
if (abs(valMin) < 1e-6 * valSpan), valMin = 0; end
if (abs(valMax) < 1e-6 * valSpan), valMax = 0; end

% COMPUTE WHERE TO PLACE TICKS
steps = reshape(((10.^(5:-1:-5))'*[5,2,1])',1,[]);
for s=1:numel(steps) % tck - ticks
  step = steps(s);
  minTick = ceil(valMin/step);
  maxTick = floor(valMax/step);
  if (maxTick - minTick + 1 < 3)
    continue
  end
  tck = (minTick:maxTick)*step;
  break;
end

% compute where to place contours
for s=1:10
  contStep = step / s;
  minCont = ceil(valMin / contStep);
  maxCont = floor(valMax / contStep); 
  if (maxCont - minCont + 1 < 10)
    continue
  end
  cont = (minCont:maxCont)*contStep;
  break;
end

% plot contours
tricontour2(mesh, sol, cont, 0.5, valMax+0.1*valSpan);
tricontour2(mesh, sol, tck,  1,   valMax+0.1*valSpan);

% DRAW COLORBAR
hc = colorbar;
cpos = [apos(1) + 1.05*apos(3), apos(2) + 0.05*apos(4), 0.05*apos(3), 0.9*apos(4)];
set(hc, ...
  'units', 'normalized', ...
  'location','eastoutside',...
  'position', cpos, ...
  'TickLength', [0.05*apos(3) 0.05*apos(3)], ...
  'ytick', cont, ...
  'yticklabel',repmat('',[numel(cont),1]))

% EXPORT FIGURE
export_fig([opt.folder, opt.figfile], '-m6', '-a3', '-zbuffer');

% EXPORT TIKZ FILE
hf = fopen([opt.folder, opt.tikzfile],'w');
fprintf(hf,'\\begin{scope}\n');
fprintf(hf, ['\\node[anchor=north west,inner sep=0]' ...
  ' at (0, 0) {\\includegraphics[height=%fcm]{' ...
  'figures/' opt.figfile '}};\n'], opt.height);
for s=1:numel(tck)
   tickheight = (tck(s) - valMin)/(valMax-valMin)*(0.9*opt.height) - 0.95*opt.height;
   fprintf(hf, ['\\node[anchor=west, inner sep = 0] ' ...
     'at (%f,%f) {\\tiny $' num2str(tck(s)) '$};\n'], ...
     opt.height*1.13, tickheight);
end
fprintf(hf,'\\end{scope}\n');
fclose(hf);
end

