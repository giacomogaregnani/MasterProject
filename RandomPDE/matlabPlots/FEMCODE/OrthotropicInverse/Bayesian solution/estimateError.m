function [mErr, sErr] = estimateError(mesh, locations, nErrEst, nFine, samples)

box = [0 1 0 1];
meshFine = structured_mesh(box, [nFine, nFine], struct('centre',false));
meshFine.box = box;
meshFine.mesh_size = box(4) / nFine;
meshFine = set_macro_bdflag(meshFine);
pde = define_pde();

err = [];

if nargin == 4
    samples = zeros(2, nErrEst);
    samples(1, :) = randn(1, nErrEst);
    samples(2, :) = 0.5 * rand(1, nErrEst);
end

for i = 1 : nErrEst
    display(['error estimation iteration: ' , num2str(i)])

    param(1) = samples(1, i); param(2) = samples(2, i);
    vp = get_variational_problem(pde, [exp(param(1)), param(2)]);
    
    [solCoarse, femspace, ~] = my_elasticity_2d(mesh, vp);
    resultsCoarse = get_observations(mesh, solCoarse, locations, femspace);
    
    [solFine, femspace, ~] = my_elasticity_2d(meshFine, vp);
    resultsFine = get_observations(meshFine, solFine, locations, femspace);
    
    err = [err, resultsFine(:) - resultsCoarse(:)];    
end

mErr = mean(err, 2);
sErr = cov(err');


return