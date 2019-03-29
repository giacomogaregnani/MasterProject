for it = 0:5
  Ndiv=2^it;
  saveFile = ['~/repos/experiments/comp/e/peps/ad_eps' num2str(Ndiv) '.mat'];
  saveFileRes = ['~/repos/experiments/comp/e/peps/ad_eps' num2str(Ndiv) 'res.mat'];
  
  try
    load(saveFile,'mesh','pormesh');
  catch
    lc = 0.1;
    pormesh = composite_mesh(@pormat_expe,'C',Ndiv,struct('periodic',true,'lc',lc));
    
    pormesh.bdflag = 'dirichlet';
    vp = struct('a', 1, 'f', @(x)(repmat([0,-1],[size(x,1),1])), 'elemtype', 'p2', 'pelemtype','p1',...
      'solver','iuzawapcg','iuzawapcgopt',struct('verbose',2,'adaptive',true),...
      'adapt',struct('maxit',20,'maxdof', 20*size(pormesh.elem,1)));
    
    % adaptive Stokes solver
    [usol, ufemspace, psol, pfemspace, pormesh, vp] = ...
      stokes_ad(pormesh, vp, struct('verbose',true));
    
    mesh = pormesh;
    mesh = cleanfields(mesh);
    N =  size(mesh.node,1);
    
    bdcoor = pormesh.node(auxstructure(pormesh,'isbdnode'), :);
    bdind = auxstructure(pormesh,'bdnode');
    myeps = 1e-5;
    
    for i=1:2*Ndiv
      for j=1:3*Ndiv
        ind = (Ndiv*bdcoor(:,1) > i-1 + myeps) & (Ndiv*bdcoor(:,1) < i - myeps) ...
          & (Ndiv*bdcoor(:,2) > j-1 + myeps) & (Ndiv*bdcoor(:,2) < j - myeps);
        thiscoor = bdcoor(ind,:);
        thisind  = bdind(ind,:);
        if (isempty(thisind))
          continue
        end
        
        ang = angle((thiscoor(:,1)-(i-1/2)/Ndiv)+1j*(thiscoor(:,2)-(j-1/2)/Ndiv));
        [~,order] = sort(ang);
        thiscoor = thiscoor(order,:);
        thisind  = thisind(order,:);
        bdgmsh = struct;
        bdgmsh.node = thiscoor;
        M = size(bdgmsh.node,1);
        bdgmsh.line = [(1:M)', [2:M,1]'];
        bdgmsh.lineloop = { 1:M };
        bdgmsh.planesurface = { 1 };
        opt = struct('lc',lc/Ndiv,'periodic',false);
        hole = gmsh(bdgmsh, opt);
        
        numNewElem = size(hole.elem,1);
        numCommonNode = numel(thisind);
        numNewNode = size(hole.node,1) - numCommonNode;
        hole.elem = subs_array(hole.elem, 1:size(hole.node,1), [thisind', N+1:N+numNewNode]);
        mesh.node = [mesh.node; hole.node(numCommonNode+1:end,:)];
        mesh.elem = [mesh.elem; hole.elem];
        N = N + numNewNode;
      end
    end
    save(saveFile,'mesh','pormesh','Ndiv');
  end
  
  vphom = struct('a',@ahom_e,'f',@(x)(repmat([0,0,-1],[size(x,1),1])),'elemtype','p1','solver','agmg');
  mesh.bdflag = 'neumann';
  [sol, femspace, mesh] = poisson(mesh,vphom);
  
  %% COMPUTE L2 error over the fluid part of the domain  
  dim = 2; deg = 2; NT = size(pormesh.elem,1);
  [lambda, weight] = quadpts(dim, deg);
  if ~isfield(pormesh,'volume'), mesh.volume = simplex_volume(pormesh); end
  res = 0;
  for i=1:numel(weight)
    res = res + weight(i) * ( ...
      evalf(mesh, 1:NT, femspace, lambda(i,:), sol) - ...
      evalf(pormesh, 'all', pfemspace, lambda(i,:), psol)).^2;
  end
  res = bsxfun(@times, sum(res,2), pormesh.volume);
  L2normF = sqrt(sum(res, 1));
  
  %% Compute L2 error over the solid part of the domain
  porbary = get_rc(pormesh);
  bary = get_rc(mesh);
  bary(1:NT,:) = [];
  L2normS  = 0;
  for i=1:2*Ndiv
    for j=1:3*Ndiv
      ind = (Ndiv*porbary(:,1) > i-1) & (Ndiv*porbary(:,1) < i) & ...
        (Ndiv*porbary(:,2) > j-1) & (Ndiv*porbary(:,2) < j);
      ind = find(ind);
      if isempty(ind)
        continue;
      end
      % average in pormesh
      [lambda, weight] = quadpts(dim, pfemspace.deg);
      res = 0;
      for k=1:numel(weight)
        res = res + weight(k) * ...
          evalf(pormesh, ind, pfemspace, lambda(k,:), psol);
      end
      aver = sum(bsxfun(@times, res, pormesh.volume(ind)), 1);
      aver = Ndiv^2 * aver /  (1 - 0.6*0.4); % eplicit area
      
      % L2 norm over the pore
      ind = (Ndiv*bary(:,1) > i-1) & (Ndiv*bary(:,1) < i) & ...
        (Ndiv*bary(:,2) > j-1) & (Ndiv*bary(:,2) < j);
      ind = find(ind) + NT;
      [lambda, weight] = quadpts(dim, deg);
      res = 0;
      for zz=1:numel(weight)
        res = res + weight(zz) * (evalf(mesh, ind, femspace, lambda(zz,:), sol) - aver).^2;
      end
      res = bsxfun(@times, sum(res,2), mesh.volume(ind));
      L2normS = L2normS + sum(res, 1);
    end
  end
  L2normS = sqrt(L2normS);
  fprintf('\n\n1 / EPSILON = %d\n',Ndiv);
  fprintf('L2 ERROR IN FLUID PART: %f, \nL2 ERROR IN SOLID PART: %f\n\n',L2normF,L2normS);
  save(saveFileRes,'-v7.3');
end