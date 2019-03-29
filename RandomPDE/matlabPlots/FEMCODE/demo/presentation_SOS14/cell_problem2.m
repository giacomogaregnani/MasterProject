function cell_problem2

%% Define coarse problem
vp = struct('elemtype', 'p2', 'pelemtype', 'p1', 'bc', [0,0], 'a', 1, 'f', @(x)(repmat([0,-1], size(x,1),1)));
point = [0,2]; mesh = micro_mesh_expc(point); mesh.bdflag = 'dirichlet';


mag = '-m8'; aa = '-a2'; bgcolor = [0.1,0.3,1];
scale = 3;

mesh = uniformrefine(mesh);

  %% solve coarse problem
  [usol, ~, ~, ~, mesh, vp] = stokes(mesh, vp);
  
  
  [mesh2,update] = deperiodize(mesh);
  usol2 = usol(update,:);
  
  init_figure;
  quiver(mesh2.node(:,1), mesh2.node(:,2), scale*usol2(:,1), scale*usol2(:,2), 'color','white','autoscale','off');
  hold on; simpplot(mesh, struct('elevation',-1,'facecolor',bgcolor));
  axisset;
  export_fig(['v3.png'],'-transparent',mag,aa);
  close  
  
  function init_figure
    figure; set(gcf,'color','black','position',[100,100,270,180]);
  end

  function axisset
    xlim([-0.6,0.6]); ylim([-0.6,0.6]); axis equal; axis off;
  end
end
