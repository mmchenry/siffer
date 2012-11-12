function [sim,pred,prey] = default_params
% Generates strctures with default parameter values for shrimp feeding
% simulations.


%% SIMULATION PARAMETERS

% Water density (kg m^-3)
sim.rho_water = 1000;

% Kinematic viscosity (Pa s)
sim.mu_water = 0.001;

% Maximum simulation duration (s)
sim.dur = .5; 

% Relative tolerance for the solver (dimensionless)
sim.reltol = 1e-1;

% Number of time values to define predator flow
sim.num_time = 100;

% Number of values along the x-axis to define pred flow
sim.num_x = 100;

% Tolerance for smoothing spline for defining pred flow
sim.sp_tol = 1e-3;

% Boundaries in x and y directions for defining pred flow
sim.flow_lim = [-.05 .1 -.1 .1];


%% PREY PARAMETERS 

% Prey diameter (m)
prey.diam = 2e-3; 

% Prey length (m)
prey.len = 20e-3;

% Prey density (kg m^-3)
prey.rho = 1000;

% Number of body segments defining morphology
prey.num_segs = 20;

% Find area and volume of body segments
[prey.s,prey.area,prey.vol,prey.mass,prey.wet_area] = getMorph(prey);

% Position of COM
prey.sCOM = 0.25 * prey.len;

% Added mass coefficient for a cylinder(?) (dimensionless)
prey.added_mass = 0.6; 

% Initial prey position
prey.pos0= [9.2e-3 0];

% Initial prey speed
prey.spd0 = [0 0];

% Vector of escape force (N)
prey.esc = [8e-3,10e-3,0];

% Sensitivity threshold
prey.thresh = 0.006;


%% PREDATOR PARAMETERS 

% Max approach speed (m/s)
pred.spd.max = 1;

% Time of max speed (s)
pred.spd.t_max = 30e-3;

% Shape factor for speed
pred.spd.alpha = 2;

% Max gape (m)
pred.gape.max = 20e-3;

% Time of max gape (s)
pred.gape.t_max = 30e-3;

% Shape factor
pred.gape.alpha = 2;

% Max dist (m)
pred.dist.max = 9e-3;

% Time of max dist (s)
pred.dist.t_max = 60e-3;

% Distance shape factor 
pred.dist.alpha = 2;

% Initial distance (m)
pred.dist.init = 9.2e-3;




end


function [sD,seg_area,seg_vol,mass,wet_area] = getMorph(prey)
% Body morphology of prey

% Normalized body position 
s = linspace(0,1,prey.num_segs)';

% Dimensional body position
sD = linspace(0,prey.len,prey.num_segs)';

% Peripheral shape of the prey body
width = flipud(-(s.^6).*226.67 + (s.^5).*671.83 - (s.^4).*742.67 + ...
                        (s.^3).*371.74 - (s.^2).*81.765 + s.*7.8304);

%width = flipud(width);
% x-Sectional area of segments
seg_area = ((width.*(prey.diam/2)).^2).*pi;     

% Wetted area of the segments
wet_area = trapz(sD,2*pi.*(width.*(prey.diam/2)));

% Volume of segments
seg_vol =  trapz(sD,seg_area);

% Mass
mass = sum(seg_vol) .* prey.rho;
end



