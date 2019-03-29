% DATAFILE for CMCS FEM 1D solver

data.flag_dirichlet = [1];
data.flag_neumann   = [2];
data.flag_robin     = [];


data.bcDir          = @(x,t,param)(1.*x.^0);
data.bcNeu          = @(x,t,param)(0.*x);
data.bcRob_alpha    = @(x,t,param)(0.*x);
data.bcRob_gamma    = @(x,t,param)(0.*x);
data.bcRob_fun      = @(x,t,param)(0.*x);

data.force          = @(x,t,param)(0 + 0.*x);

data.diffusion      = @(x,t,param)(4.5 + 0.*x);
data.transport      = @(x,t,param)(32 .* x.^3 - 16 .* x + 1);
data.reaction       = @(x,t,param)(0 + 0.*x);


data.uexact         = @(x,t,param)(0 .* x);
data.uxexact        = @(x,t,param)(0 .* x);


data.u0             = @(x,t,param)(0 .* x);


data.time.t0        = 0;
data.time.tf        = 0.5;
data.time.dt        = param(1);
data.time.theta     = param(2);