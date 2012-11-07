function r = siff(p)
% Runs a SIFF simulation using matlab's ODE solver functions

%global p


%% Parameter values

% Use default parameter values, if none are provided
if nargin<1
    
    % Load flow velocity data ('f')
    load('flow_field')

    % Function (below) for specifying parameter values
    p = default_params(f);
    
    clear f
end

% Scale parameters (have yet to do)

% % Define scaling constants
% sL = L;
% sT = 10^-3;
% sM = mass_body*10^6;
% sF = sM .* sL ./ sT^2;


%% Run simulation

% Solver options
%options = odeset('Events',@evnts,'RelTol',1e-1);
options = odeset('OutputFcn',@status_check,'RelTol',1e-1);
%options = odeset('RelTol',1e-1);

% Function that runs the numerical solver
[t,X] = solver(p,options);


% Store results
r.t          = t;
r.prey.pos   = X(:,1:2);
r.prey.vel   = X(:,3:4);
r.prey.accel = [diff(r.prey.vel(:,1))./diff(t) diff(r.prey.vel(:,2))./diff(t)];

r.pred.pos      = getPredPos(t,p);
r.pred.gape_spd = getGapeSpeed(t,p);
r.pred.gape     = getGape(t,p);

r.prey.fl_vel = getFlowVel(t,p,r.prey.pos,r.pred.pos,r.pred.gape);
r.prey.PF     = 0.*getPressureForce(r.t,p,r.prey.pos,r.pred.pos,r.prey.accel,...
                                r.pred.gape,r.prey.fl_vel);
r.prey.drag = getDrag(p,r.prey.fl_vel,r.prey.vel);

clear t X options


%% Plot results
figure

subplot(3,1,1)

ax1 = gca;
h1 = line(r.t,1e3*r.pred.pos(:,1),'Color','b','parent',ax1,'linestyle','--');
set(ax1,'XColor','b','YColor','b')
xlabel('time(s)');
ylabel('Displacement (mm)');
ax2 = axes('Position',get(ax1,'Position'),...
          'XAxisLocation','top',...
          'YAxisLocation','right',...
          'Color','none',...
          'YColor','black');
%set(ax2,'YColor','black')
ylabel('Speed(m/s)');

h2 = line(r.t,r.pred.gape_spd(:,1),'Color','black','parent',ax2);
h3 = line(r.t,1e3*r.pred.gape,'Color','b');
h = [h1(1) h2 h3];
legend(h,'x pos (mm)','gape(mm)','x spd (m/s)');
title('Kinematic characteristics');


subplot(3,1,2)
plot(r.t,r.prey.drag(:,1)+r.prey.PF(:,1),'k-')
xlabel('time (s)')
ylabel('Total x force (N)')

subplot(3,1,3)
plot(1000.*r.pred.pos(:,1),1000.*r.pred.gape./2,'r-',...
    1000.*r.pred.pos(:,1),-1000.*r.pred.gape./2,'r-')
hold on
plot(1000.*r.prey.pos(:,1),1000.*r.prey.pos(:,2),'b.')
hold off
xlabel('x - coord (mm)')
ylabel('y - coord (mm)')
axis equal



end


function [t,X] = solver(p,options)

global pred_position

% Acceleration in previous timestep
A_prev = [0 0];


% Solve governing equation
[t,X] = ode45(@gov_eqn,[0 p.dur],[p.X0; p.U0],options);

    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % Prey body position (inertial FOR)
        pos = X(1:2)';
        
        % Prey body velocity (inertial FOR)
        V = X(3:4)';
        
        % Position of the predator's gape (inertial FOR)
        pred_pos = getPredPos(t,p);
        
        pred_position = pred_pos;
        
        % Instantaneous gape diameter
        gape = getGape(t,p);
        
        % Flow velocity (intertial FOR)
        U = getFlowVel(t,p,pos,pred_pos,gape);
        
        % Thrust 
        %T = getThrust(t,p);
        
        % Pressure force
        PF = 0.*getPressureForce(t,p,pos,pred_pos,A_prev,gape,U);
        %PF = [-5e-3*gape./p.gape.max 0];
        
        % Drag
        D = getDrag(p,U,V);
        
        % Body acceleration
        accel = (D + PF)./...
                 (p.prey.mass + p.rho_water*sum(p.prey.vol)*p.added_mass);
   
        % Define output: speed
        dX(1,1) = V(1);
        dX(2,1) = V(2);
        
        % Define output: acceleration
        dX(3,1) = accel(1);
        dX(4,1) = accel(2);
        
        % Update previous body acceleration
        A_prev = accel;
        
    end
    


end

%function [value,isterminal,direction] = evnts(t,X)
function status = status_check(t,y,flag,varargin)

global pred_position

if ~isempty(t)
    pos = y(1);
    if (pos <= pred_position(1))
        status = 1;
    else
        status = 0;
    end
else
    status = 1;
end
    


%isterminal = 1;
%direction = 0;

%rel_pos = pos - pred_position(1);

% if (pos <= pred_position(1))
%     status = 1;
%     
% else
%     status = 0;;
% end
%value

% % Check that linkage geometry is possible
% if pos < pred_position(1)
%     value = 0;
%     direction = 0;
%     disp('Capture successful!')
% else
%     value = 1;
%     direction = 1;
% end
end


function [sD,seg_area,seg_vol,mass,wet_area] = getMorph(p)
% Body morphology of prey

% Normalized body position 
s = linspace(0,1,p.prey.num_segs)';

% Dimensional body position
sD = linspace(0,p.prey.len,p.prey.num_segs)';

% Peripheral shape of the prey body
width = flipud(-(s.^6).*226.67 + (s.^5).*671.83 - (s.^4).*742.67 + ...
                        (s.^3).*371.74 - (s.^2).*81.765 + s.*7.8304);

%width = flipud(width);
% x-Sectional area of segments
seg_area = ((width.*(p.prey.diam/2)).^2).*pi;     

% Wetted area of the segments
wet_area = trapz(sD,2*pi.*(width.*(p.prey.diam/2)));

% Volume of segments
seg_vol =  trapz(sD,seg_area);

% Mass
mass = sum(seg_vol) .* p.prey.rho;

end


function spd = getGapeSpeed(t,p)
% Speed of flow (inertial FOR) at mouth
    
spd = p.spd.max * ((t./p.spd.t_max).*...
                          (exp(1-(t./p.spd.t_max)))).^p.spd.alpha;

end


function gape = getGape(t,p)
% Gape diameter
gape = p.gape.max.*((t./p.gape.t_max).*...
                  (exp(1-(t./p.gape.t_max)))).^p.gape.alpha;
% Avoid zero values 
%gape = gape + realmin;
end


function pos = getPredPos(t,p)
    dist = 0*p.dist.init + p.dist.max.*((t./p.dist.t_max).*...
                         (exp(1-(t./p.dist.t_max)))).^p.dist.alpha;
    pos = [dist dist.*0];
end


function esc = getThrust(t,p)
    esc = p.esc;
end


function D = getDrag(p,U,V)
% Drag on prey body

for i = 1:size(U,1)
    if (sum(U(i,:))==0) 
        D(i,:) = [0 0];
    else
        % Flow speed in the prey FOR
        spd = sqrt( (U(i,2)-V(i,2)).^2 + (U(i,1)-V(i,1)).^2);
        
        % Reynolds number
        Re = spd * p.prey.diam / p.eta;
        
        % Drag coefficient [Kils 1979 (Re 300-15,000)
        % & Roi's mesurements (10,000-85,000)]
        Cd = (Re.^(-0.1703)).*0.0708;
        
        % Drag
        D(i,:) = 0.5 * Cd * sum(p.prey.wet_area) * ...
                  p.rho_water * spd .* (U(i,:)-V(i,:));
        
        clear Re Cd spd
    end
end

%D = [-1e-5 0];
end


function [P,C] = getPressureForce(t,p,pos,pred_pos,A_prev,gape,U)
% Calc pressure force

for i = 1:length(gape)
    if gape(i)==0
        P(i,:) = [0 0];
    else
        
        % Prey position relative to predator
        rel_pos = (pos(i,:) - pred_pos(i,:))./gape(i);
        
        % Unit vector for direction of prey WRT predator
        b_orient = rel_pos ./ abs(rel_pos);
        
        b_orient(isnan(b_orient)) = zeros(1,sum(isnan(b_orient)));
        
        % Speed at mouth
        spd_mouth = getGapeSpeed(t(i),p);
        
        % Position values for body
        Xs = repmat(rel_pos,p.prey.num_segs,1) + ...
                    [p.prey.s./gape(i) p.prey.s.*0] .* ...
                    repmat(b_orient,p.prey.num_segs,1);
        
%         % Scale position data for flow field
%         x = p.flow.x * gape(i);
%         y = p.flow.y * gape(i);
        
        % Flow speed in the prey FOR
        spd = sqrt( U(i,1).^2 + U(i,2).^2 );
        
        
        
        %TODO: modify this to do 2D
        
        % If the flow field doesn't encompass the prey . . .
        if max(max(p.flow.x(:))>Xs(:,1))==0
            dUdx = 0.*Xs;
        
        elseif (max(p.flow.x(:))>min(Xs(:,1))) && ...
               (max(p.flow.x(:))<max(Xs(:,1)))
            % min(max(p.flow.y(:))>Xs(:,2))==0 || ...
            % min(min(p.flow.y(:))<Xs(:,2))==0
         
            % . . . use the speed gradient at edge of flow field
%             dU = spd_mouth .* (p.flow.u(1,end)-p.flow.u(1,end-5));
%             dx = gape(i).*(p.flow.x(1,end)-p.flow.x(1,end-5));
            dx = gape(i).*(max(Xs(:,1)) - max(p.flow.x(:)));
            dU = - min(spd_mouth .* p.flow.u(:,end));
            dUdx = [repmat(dU/dx,length(Xs),1) zeros(length(Xs),1)];
        %elseif min(Xs(:,1))<=0
            %[gape(i) dU dx dU/dx]
        else
            % Otherwise, calcualte the speed at each segment
            Us(:,1) = spd_mouth .* interp2(p.flow.x,p.flow.y,p.flow.u,...
                Xs(:,1),Xs(:,2),'linear',0);
            Us(:,2) = spd_mouth .* interp2(p.flow.x,p.flow.y,p.flow.v,...
                Xs(:,1),Xs(:,2),'linear',0);
            
            % Velocity gradient along the body
            dUdx = [0 0; diff(Us(:,1))./diff(p.prey.s) diff(Us(:,2))./diff(p.prey.s)];
        end
        
        % Pressure gradient along the body
        dPdx(:,1) = p.rho_water * (1+p.added_mass) .* (A_prev(1) + spd.*dUdx(:,1));
        dPdx(:,2) = p.rho_water * (1+p.added_mass) .* (A_prev(2) + spd.*dUdx(:,2));
        
        % Pressure force, integrated long the length
        P(i,1) = -trapz(p.prey.s,dPdx(:,1).*p.prey.area);
        P(i,2) = -trapz(p.prey.s,dPdx(:,2).*p.prey.area);
        
        C = spd.*dUdx(:,1);
        
        clear x y Us dUdx dPdx rel_pos b_orient spd_mouth Xs dU dx
    end
end
%P = 0.*pos;
end



function U = getFlowVel(t,p,pos,pred_pos,gape)
% Flow in at the position of the prey in the inertial FOR

% Loop through gape values
for i = 1:length(gape)
    if gape(i)==0
        U(i,:) = [0 0];
    else
        % Prey position relative to predator
        rel_pos(i,:) = (pos(i,:) - pred_pos(i,:))./gape(i);
        
        % Speed at mouth
        spd_mouth = getGapeSpeed(t(i),p);
        
        % Flow speed in the inertial FOR
        U(i,1) = spd_mouth .* interp2(p.flow.x,p.flow.y,p.flow.u,...
                                    rel_pos(i,1),rel_pos(i,2),'linear',0);
        U(i,2) = spd_mouth .* interp2(p.flow.x,p.flow.y,p.flow.v,...
                                    rel_pos(i,1),rel_pos(i,2),'linear',0);
    end
end
end


function p = default_params(f)


%% Parameters

% Prey diameter (m)
p.prey.diam = 2e-3; 

% Prey length (m)
p.prey.len = 20e-3;

% Prey density (kg m^-3)
p.prey.rho = 1000;

% Number of body segments defining morphology
p.prey.num_segs = 20;

% Find area and volume of body segments
[p.prey.s,p.prey.area,p.prey.vol,p.prey.mass,p.prey.wet_area] = getMorph(p);

% Duration (s)
p.dur = .5; 

% Max approach speed (m/s)
p.spd.max = 1;

% Time of max speed (s)
p.spd.t_max = 30e-3;

% Shape factor for speed
p.spd.alpha = 2;

% Max gape (m)
p.gape.max = 20e-3;

% Time of max gape (s)
p.gape.t_max = 30e-3;

% Shape factor
p.gape.alpha = 2;

% Max dist (m)
p.dist.max = 9e-3;

% Time of max dist (s)
p.dist.t_max = 60e-3;

% Distance shape factor 
p.dist.alpha = 2;

% Initial distance (m)
p.dist.init = 9.2e-3;

% Vector of escape force (N)
p.esc = [8e-3,10e-3,0];

% Sensitivity threshold
p.th = 0.006;

% Water density (kg m^-3)
p.rho_water = 1000;

% Kinematic viscosity (m^2 s^-1)
p.eta = 1e-6;

% Added mass coefficient of a sphere (dimensionless)
p.added_mass = 0.6; 

% Initial prey position
p.X0 = [9.2e-3 0];

% Initial prey speed
p.U0 = [0 0];

% Flow data (dimensionless)
p.flow.x = f.x;
p.flow.y = f.y;
p.flow.u = f.u;
p.flow.v = f.v;

clear f
end