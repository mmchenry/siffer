function fl = flow_noflow(sim)
% Creates flow field with zero flow.  Used for testing escape response
% simulations.


%% Define time and spatial domain

% Time vector
t = linspace(0,sim.dur,sim.num_time)';

% Spatial vectors
xs = linspace(sim.flow_lim(1),sim.flow_lim(2),sim.num_x)-sim.flow_lim(1);
ys = sim.flow_lim(3):mean(diff(xs)):sim.flow_lim(4);

% Mesh position for flow in global FOR
[Xg,Yg] = meshgrid(xs,ys);


%% Predator gape and position

% Store values
fl.t        = t;
fl.pos      = [0.*t 0.*t];  
fl.gape_spd = 0.*t;
fl.X        = Xg;
fl.Y        = Yg;
fl.U        = zeros(size(Xg,1),size(Xg,2),length(t));
fl.V        = fl.U;
fl.dUdx     = fl.U;
fl.dUdy     = fl.U;
fl.dUdt     = fl.U;
fl.dVdx     = fl.U;
fl.dVdy     = fl.U;
fl.dVdt     = fl.U;

end


