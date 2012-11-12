function run_shrimp_feed
% Wrapper to run a single simulation of predation on a shrimp (same
% conditions as V0.9.5 simulation.

% Root path for saving and loading simulation data
sim_root = '/Users/mmchenry/Dropbox/Projects/Holzmann/sims';

% Default parameter values
[sim,pred,prey] = default_params;

% Flow field is loaded, if previously defined
if 0 %~isempty(dir([sim_root filesep 'flow_data.mat']))
    disp(' ');disp('Loading flow data . . .');disp(' ');
    load([sim_root filesep 'flow_data.mat']);
    
% Otherwise, feeding_splines generates the 'flow_data' file
else
    makeFeedField(sim,pred,sim_root);
end

% Find 