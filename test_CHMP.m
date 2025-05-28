%% Higher-Order Rotational Synchronization
%% ----------------------------------------
%% Test CHMP for rotational synchronization on a Erdos-Renyi random hypergraph.
%% ----------------------------------------
%% Author: Addie Duncan
%% ----------------------------------------
%% Input Parameters:
% n: Number of random rotations to generate.
% p: Edge probability for Erdos-Renyi model.
% q: Corruption probability.
% noise_model: Noise model, choose 'gaussian' or 'uniform'.
% sigma: Noise level.
% sampling: Number of cycles to sample from cycle set, default is 1 which considers all cycles.
% beta_init: Reweighting parameter initial value.
% beta_max: Reweighting parameter max value.
% iter: Number of iterations (determines rate of beta_t increase between beta_init and beta_max
%% ----------------------------------------
%% Outputs:
% CHMP average relative error
% Rotation recovery procrustes error (two methods, MST and GCW)


%% define parameters

% data parameters
parameters.n = 30; 
parameters.p = 1; 
parameters.q = 0.2; 

% noise parameters
noise_model = 'gaussian';
parameters.sigma = 0; 

% CHMP parameters
parameters.sampling = 1; % to use all possible cycles set parameter = 1
parameters.beta_init = 1;
parameters.beta_max = 40;
parameters.iter = 20;

parameters.rate = nthroot(parameters.beta_max/parameters.beta_init,parameters.iter);


%% generate data
disp('Generating random hypergraph')
tic;

% generate random rotations
rotations = generateRotations(parameters);

% generate random hypergraph (Erdos Renyi Model)
hypergraph = generateHypergraph(parameters,rotations,noise_model);

data_generation_time=toc; % save time for data generation


%% Recover hyperedge corruption errors using CHMP

% run CHMP
CHMP_out = CHMP(hypergraph,parameters);


%% Calculate Results

% calculate error for hyperedge corruption
CHMP_SVec = CHMP_out.SVec;
CHMP_ErrVec = hypergraph.ErrVec;
CHMP_num_edges = hypergraph.num_hyperedges;

CHMP_REL_ERR = sum(abs(CHMP_ErrVec - CHMP_SVec))/CHMP_num_edges;


%% Recover Rotations

% MST
[CHMP_R_est_MST,CHMP_MST_recovery_time] = MST_hypergraph(hypergraph,parameters,CHMP_out);
CHMP_rot_recovery_ERR_MST = procrustes_error_SO3(CHMP_R_est_MST, rotations.R_orig, parameters);

% GCW
[CHMP_R_est_GCW,CHMP_GCW_recovery_time] = GCW_hypergraph(hypergraph,parameters,CHMP_out);
CHMP_rot_recovery_ERR_GCW = procrustes_error_SO3(CHMP_R_est_GCW, rotations.R_orig, parameters);


%% Print results

disp("CHMP Error =")
disp(CHMP_REL_ERR)

disp("Procrustes Error (MST) =")
disp(CHMP_rot_recovery_ERR_MST)

disp("rocrustes Error (GCW) =")
disp(CHMP_rot_recovery_ERR_GCW)

