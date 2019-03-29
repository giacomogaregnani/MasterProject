function [mesh, father] = get_fine_ad4_X(RB)
%GET_FINE_AD3_X Summary of this function goes here
%   Detailed explanation goes here
maxdof = 7070; 
nref = 4;
dim=2;
filename = ['~/repos/anmc_fem/experiments/expX/mesh_fine_ad' num2str(nref) '.mat'];
try
  load(filename,'mesh','father');
catch
  mesh = RB.rmesh;
  [mesh, father] = bisect(mesh, [2,3], RB.R);
  [mesh, father] = uniformrefine(mesh, father);
  mesh = periodize(mesh,[-1/2,1/2,-1/2,1/2]);

  [paramsX, paramsY] = meshgrid(...
    linspace(RB.param.min(1),RB.param.max(1),5), ...
    linspace(RB.param.min(2),RB.param.max(2),5));
  params = [paramsX(:), paramsY(:)];
  NP = size(params,1);
  [meshes, res] = deal(cell(NP,1));
  
  % adaptive
  iter=1;
  while true
    for k=1:NP
      % GET TRANSFORMED MESH
      meshes{k} = mesh;
      par = params(k,:);
      R = numel(RB.C);
       done = [];
      for i=1:R
        C = subsmu(RB.C{i}, par);
        G = subsmu(RB.G{i}, par);
        for j=1:dim+1
          el = find(father==i);
          nodes = meshes{k}.elem(el, j);
          [newnodes, ind] = setdiff(nodes,done);
          done = union(done, newnodes);
          nodeval = get_rc(meshes{k}, [], el(ind), j, 1, struct('adjust_to','bary'));
          meshes{k}.node(newnodes,:) = bsxfun(@plus, C', nodeval * G');
        end
      end
      
      meshes{k}.bdflag = 'dirichlet';
      vp = struct('f',@fstokesmicro,'elemtype','p2','pelemtype','p1','a',1);
    
      % SOLVE PROBLEMS
      [usol, ufemspace, psol, pfemspace, meshes{k}, vp] = stokes(meshes{k}, vp);
      ndof = 2*ufemspace.ndof + pfemspace.ndof;
      fprintf('ITER: %d, DOF: %d/%d, PROBLEM: %d/%d\n', iter, ndof,maxdof,k, NP);
      
      % ESTIMATE
      res{k} = stokes_residual(meshes{k}, vp, ufemspace, usol, pfemspace, psol);
    end
    
    %% STOP?
    if (ndof >= maxdof) 
      break;
    end
    
    maxres = max(res{1},[],3);
    for k=2:NP
      maxres = max(maxres, max(res{k},[],3));
    end
    
    %% MARK
    mElements = markelem(maxres, struct('method','L2','unify',true,'theta',0.25));
    
    %% REFINE
    [mesh, father] = bisect(mesh, mElements, father);
    iter=iter+1;
  end
  mesh=cleanfields(mesh);
  save(filename,'mesh','father');
end
end

