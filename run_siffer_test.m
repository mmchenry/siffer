function run_siffer_test
% Wrapper to validate the siffer solver with simplified morphology and flow
% conditions.


%% Define paths

% Root directory for saving and loading simulation data
sim_path = '/Users/mmchenry/Dropbox/Projects/Holzmann/sims/validate';
%sim_path = '/Volumes/Flow HD/Dropbox/Projects/Holzmann/sims/validate';


%% Set parameter values

% Default parameter values
%[sim,pred,prey] = params_validate;
[sim,pred,prey] = params_danio;

run_const_accel = 1;
run_const_grad  = 0;

% Visualize flow
vis_flow = 0;

% Vizualize simulation results
vis_sim = 1;


%% Test code on a spatially uniform flow field with constant acceleration
% Used to text predictions of drag 

if run_const_accel

% Create spatially uniform flow field
%fl = flow_uniform_ramp(sim,pred,vis_flow);
%save([sim_path filesep 'uni_flow'],'fl')
load([sim_path filesep 'uni_flow.mat'])    

% Run simulation
tic
r = siffer(sim,prey,fl);

% Save results
save([sim_path filesep 'sim_results.mat'],'r')

% Report time
t_lapse = toc;
disp(['Sim complete in ' num2str(t_lapse) 's'])


% Visualize results
if vis_sim 
    figure;
    
    subplot(2,1,1)
    plot(r.t,1000.*r.pos(1,:),'k-')
    xlabel('time (s)');ylabel('position (mm)');
    title('Constant acceleration');
    
    subplot(2,1,2)
    plot(r.t,r.D(1,:),'r-',r.t,r.PF(1,:),'b-',r.t,r.AR(1,:),'g-')
    legend('D','PF','AR')
    xlabel('time (s)');ylabel('Force (N)');
end

end


%% Test code on a spatially uniform flow field with fixed spatial gradient

if run_const_grad

% Create spatially uniform flow field
%fl = flow_constant_gradient(sim,pred,vis_flow);
%save([sim_path filesep 'grad_flow'],'fl')
load([sim_path filesep 'grad_flow.mat'])    

% Run simulation
tic
r = siffer(sim,prey,fl);

% Save feeding flow field data
save([sim_path filesep 'sim_results.mat'],'r')

% Report time
t_lapse = toc;
disp(['Sim complete in ' num2str(t_lapse) 's'])


% Visualize results
if vis_sim 
    figure;
    plot(r.t,r.pos(:,1),'r--')
    xlabel('time');ylabel('position');title('Fixed gradient');
end


end







end

