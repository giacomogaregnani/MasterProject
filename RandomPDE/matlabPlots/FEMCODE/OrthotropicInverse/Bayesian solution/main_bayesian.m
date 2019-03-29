% Load data
clc; clear;
load('observations_numObs10.mat');

% Noise data
Observations = reshape(Observations, size(Observations, 2), size(Observations, 3), size(Observations, 4));
std_noise = 0.0005;
noise = normrnd(0,std_noise, [nr_edges, numObservations_perEdge, 2]);
noisy_observations = Observations + noise;

% Truth
truth = [7, 0.3];

% Input parameters
N = 5;
pde = define_pde();

% Macro mesh and femspace
dim = 2; box = [0 1 0 1];
mesh = structured_mesh(box,[N, N],struct('centre',false));
mesh.box = box;
mesh.mesh_size = box(4)/N;
mesh = set_macro_bdflag(mesh);
vp = get_variational_problem(pde, []);
femspace = get_femspace(mesh, 'p1');

% Prior and likelihood definitions
log_prior = @(p) log_prior_homog(p);

% Bayesian solution
Nsample = 1000;
% errEst = [true, true];
errEst = true;

for pr = errEst
    sample = run_bayesian(Nsample, noisy_observations, std_noise, locations, log_prior, mesh, false, pr, 10, 100, 10);
    figure
    plot(truth(1), truth(2), 'xk')
    hold on
    [dens, eval] = ksdensity([exp(sample(1, :)); sample(2, :)]');
    evalX = reshape(eval(:, 1), 30, 30);
    evalY = reshape(eval(:, 2), 30, 30);
    dens = reshape(dens, 30, 30);
    contour(evalX, evalY, dens, 'b');
    if pr
        save(['results_errEst_N_', num2str(N)])
    else
        save(['results_det_N_', num2str(N)])
    end
end