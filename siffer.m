function r = siffer(sim,prey,fl)
% Runs a SIFF simulation using matlab's ODE solver functions. The three
% input structions (sim,prey,fl) respectively provide parameter values for
% the simulation, prey and flow field.


%% Scale parameters
% Keeps parameter values at reasonable orders of magnitude.  Also passes
% only required parameters into the solver.

% Prey parameters
py.s        = prey.s            ./sim.sL;
py.r        = prey.r            ./sim.sL;
py.sCOM     = prey.sCOM         ./sim.sL;
py.mass     = prey.mass         ./sim.sM;
py.wet_area = prey.wet_area     ./sim.sL^2;
py.x_area   = prey.x_area       ./sim.sL^2;
py.vol      = prey.vol          ./sim.sL^3;

py.pos0     = [prey.pos0(1:2)./sim.sL prey.pos0(3)];
py.vel0     = [prey.vel0(1:2)./sim.sL.*sim.sT prey.pos0(3).*sim.sT];
py.add_mass = prey.add_mass;
py.Cd       = prey.Cd;

% Flow field
f.t         = fl.t              ./sim.sT;
f.pos       = fl.pos            ./sim.sL;
f.gape_spd  = fl.gape_spd       ./sim.sL .*sim.sT;
f.X         = fl.X              ./sim.sL;
f.Y         = fl.Y              ./sim.sL;
f.U         = fl.U              ./sim.sL .*sim.sT;
f.V         = fl.V              ./sim.sL .*sim.sT;
f.dUdx      = fl.dUdx           .*sim.sT;
f.dVdy      = fl.dVdy           .*sim.sT;
f.dUdt      = fl.dUdt           ./sim.sL .*sim.sT^2;
f.dVdt      = fl.dVdt           ./sim.sL .*sim.sT^2;

% Simulation parameters
s.rho_water = sim.rho_water     ./sim.sM .*sim.sL^3;
s.dur       = sim.dur           ./sim.sT;
s.reltol    = sim.reltol;

% Scaling parameters
s.sL = sim.sL;
s.sT = sim.sT;
s.sM = sim.sM;
s.sF = sim.sF;

% Overwrite previous structures, clear temporaries
prey = py; clear py
fl   = f;  clear f
sim  = s;  clear s

% Define useful global parameter
global pred_position


%% Run simulation

% Solver options
%options = odeset('Events',@evnts,'RelTol',1e-1);
options = odeset('OutputFcn',@status_check,'RelTol',sim.reltol);
%options = odeset('OutputFcn',@odeplot,'RelTol',sim.reltol);

% Prey midline coorindates in prey FOR
mid_x = prey.s - prey.sCOM;
mid_y = mid_x.*0;

% Body coorindates in prey FOR
marg_x = [prey.s; prey.s(end:-1:1)]-prey.sCOM;
marg_y = [prey.r; -prey.r(end:-1:1)];

% Solve governing equation
sol = ode45(@(t,X)gov_eqn(t,X),[0 sim.dur],[prey.pos0 prey.vel0],options);
%sol = ode15s(@(t,X)gov_eqn(t,X,A_prev),[0 sim.dur],[prey.pos0 prey.vel0],options);

 % Predator's position (used by status_check)
 pred_pos = interp1(fl.t,fl.pos(:,1),sol.x);

% Determine indicies with & without capture
if max(isnan(sol.y(1,:)))
    warning(['Nan values appeared in your solution -- the prey is likely ' ...
            'falling outside of your fluid domain'])
    capture = nan;
    idx = 1:find(isnan(sol.y(1,:)),1,'first');
elseif max(sol.y(1,:) < pred_pos) 
    capture = 1;
    idx = 1:find(sol.y(1,:)<pred_pos,1,'first');
else
    capture = 0;
    idx = 1:length(sol.y(1,:));
end

% Time to evaluate results
t = linspace(0,max(sol.x(idx)),length(sol.x(idx))*1.5);

% Find solution at equal intervals
[X,dX] = deval(sol,t);

% Evaluate governing equation for other variables 
% (Note: comment this out when not used)
[dX,D,PF,AR,m_x,m_y] = gov_eqn(t,X);

 % Predator's position (used by status_check)
 pred_pos = interp1(fl.t,fl.pos(:,1),t);

% Store results in re-scaled units
r.cap      = capture;
r.t        = t .* sim.sT;
r.pos      = [X(1:2,:).*sim.sL; X(3,:)];
r.marg_x   = [m_x.*sim.sL];
r.marg_y   = [m_y.*sim.sL];
r.vel      = [X(4:5,:).*sim.sL./sim.sT; X(6,:)./sim.sT];
r.acc      = [dX(4:5,:).*sim.sL./sim.sT^2; dX(6,:)./sim.sT^2];
r.D        = D .*sim.sF;
r.PF       = PF .*sim.sF;
r.AR       = AR .*sim.sF;
r.pred_pos = [pred_pos.*sim.sL; pred_pos.*0];


%% Governing equation 

    function [dX,D,PF,AR,m_x,m_y] = gov_eqn(t,X)
        % ODE of the dynamics of the system
        % Matricies arranged with time along columns, if evaluated over
        % more than a single instant
        
%         % Check inputs
%         if nargin<3 && length(t)>1
%             error(['You need to provide the body acceleration for the ' ...
%                 ' prey in order to evaluate the results over time']);
%         end
        
%         % Prey COM acceleration (if analyzing simulation results)
%         if length(t) > 1
%             accCOM_x   = bod_acc(1,:);
%             accCOM_y   = bod_acc(2,:);
%             accCOM_ang = bod_acc(3,:);
%         end
        
        % Prey COM position
        posCOM_x   = X(1,:);
        posCOM_y   = X(2,:);
        posCOM_ang = X(3,:);
        
        % Prey COM velocity
        velCOM_x   = X(4,:);
        velCOM_y   = X(5,:);
        velCOM_ang = X(6,:);
        
        % Loop through time for prey position
        for i = 1:length(t)
            % Rotation matrix
            M = [cos(posCOM_ang(i)) -sin(posCOM_ang(i)); ...
                 sin(posCOM_ang(i))  cos(posCOM_ang(i))];
            
            % Rotate
            py_mid  = M*[mid_x';  mid_y'];
            py_marg = M*[marg_x'; marg_y'];
            
            % Translate position into inertial FOR
            s_x(:,i) = py_mid(1,:)' + posCOM_x(i);
            s_y(:,i) = py_mid(2,:)' + posCOM_y(i);
            m_x(:,i) = py_marg(1,:)' + posCOM_x(i);
            m_y(:,i) = py_marg(2,:)' + posCOM_y(i);
            
            clear M py_mid py_marg
        end
        
        clear i
        
        %         % Velocity of prey's segments (rotation + translation)
        %         s_dxdt = ang_velCOM.*(prey.s-prey.sCOM) .* sin(ang_velCOM) + x_velCOM;
        %         s_dydt = ang_velCOM.*(prey.s-prey.sCOM) .* cos(ang_velCOM) + y_velCOM;
        
        % Velocity of prey's segments (just translation, n segments x m time vals))
        s_dxdt = s_x.*0 + repmat(velCOM_x,length(prey.s),1);
        s_dydt = s_y.*0 + repmat(velCOM_y,length(prey.s),1);
        
        % Interpolate flow conditions at segments
        for i = 1:length(t)
            
            % Bracket time around present time
            tmp = abs(t(i)-fl.t);
            tIdx = find(tmp==min(tmp),1,'first');
            tIdx = [max([1 tIdx-2]) min([length(fl.t) tIdx+2])];
            
            % Bracket region-of-interest around prey in x
            xOffset = range(m_x(:,i))/4;
            xMin    = min(m_x(:,i))-xOffset;
            xMax    = max(m_x(:,i))+xOffset;
            xIdx    = (fl.X(1,:)>xMin) & (fl.X(1,:)<xMax);
            
            % Bracket region-of-interest around prey in y
            yOffset = range(m_y(:,i))/4;
            yMin    = max([min(m_y(:,i))-yOffset min(fl.Y(:,1))]);
            yMax    = min([max(m_y(:,i))+yOffset max(fl.Y(:,1))]);
            yIdx    = (fl.Y(:,1)>yMin) & (fl.Y(:,1)<yMax); 
            
            % Coordinates for ROI
            X_t    = repmat(fl.X(yIdx,xIdx),[1,1,length(tIdx)]);
            Y_t    = repmat(fl.Y(yIdx,xIdx),[1,1,length(tIdx)]);
            
            % Check if prey outside of flow field
            if ~max(xIdx) || ~max(yIdx)
                
                % Alert
                warning('prey outside of flow field');
                
                % Flow velocity
                s_U(:,i) = 0.*s_x(:,i);
                s_V(:,i) = 0.*s_x(:,i);
                
                % Spatial gradient in velocity at segments
                s_dUdx(:,i) = 0.*s_x(:,i);
                s_dVdy(:,i) = 0.*s_x(:,i);
                
                % Flow acceleration at segments
                s_dUdt(:,i) = 0.*s_x(:,i);
                s_dVdt(:,i) = 0.*s_x(:,i);
                s_dUdt_prev(:,i) = 0.*s_x(:,i);
                s_dVdt_prev(:,i) = 0.*s_x(:,i);
                
            % If inside, continue with interpolation
            else  
                % Time step
                dt = mean(diff(fl.t(tIdx)));
                
                % Time values for current bracket
                for j = 1:length(tIdx)
                    T_t(:,:,j)  = fl.X(yIdx,xIdx).*0+fl.t(tIdx(j));
                end
                
                % Reduce flwo domain to ROI and time bracket
                U_t    = fl.U(yIdx,xIdx,tIdx);
                V_t    = fl.V(yIdx,xIdx,tIdx);
                dUdx_t = fl.dUdx(yIdx,xIdx,tIdx);
                dVdy_t = fl.dVdy(yIdx,xIdx,tIdx);
                dUdt_t = fl.dUdt(yIdx,xIdx,tIdx);
                dVdt_t = fl.dVdt(yIdx,xIdx,tIdx);
                
                % Flow velocity
                s_U(:,i) = interp3(X_t,Y_t,T_t,U_t,...
                    s_x(:,i),s_y(:,i),s_x(:,i).*0+t(i));
                s_V(:,i) = interp3(X_t,Y_t,T_t,V_t,...
                    s_x(:,i),s_y(:,i),s_x(:,i).*0+t(i));
                
                % Spatial gradient in velocity at segments
                s_dUdx(:,i) = interp3(X_t,Y_t,T_t,dUdx_t,...
                    s_x(:,i),s_y(:,i),s_x(:,i).*0+t(i));
                s_dVdy(:,i) = interp3(X_t,Y_t,T_t,dVdy_t,...
                    s_x(:,i),s_y(:,i),s_x(:,i).*0+t(i));
                
                % Flow acceleration at segments
                s_dUdt(:,i) = interp3(X_t,Y_t,T_t,dUdt_t,...
                    s_x(:,i),s_y(:,i),s_x(:,i).*0+t(i));
                s_dVdt(:,i) = interp3(X_t,Y_t,T_t,dVdt_t,...
                    s_x(:,i),s_y(:,i),s_x(:,i).*0+t(i));
                
                % Previous flow acceleration at segments
                t_prev = min(T_t(:));
                s_dUdt_prev(:,i) = interp3(X_t,Y_t,T_t,dUdt_t,...
                                   s_x(:,i),s_y(:,i),s_x(:,i).*0+t_prev);
                s_dVdt_prev(:,i) = interp3(X_t,Y_t,T_t,dVdt_t,...
                                   s_x(:,i),s_y(:,i),s_x(:,i).*0+t_prev);
                
            end
            
            % Clear for next iteration               
            clear tmp tIdx xOffset xMin xMax xIdx yOffset yMin yMax yIdx
            clear X_t Y_t T_t U_t V_t dUdx_t dVdy_t dUdt_t dVdt_t dt t_prev
        end
        
        % Relative velocity at segments
        s_relvel_x = s_U - s_dxdt;
        s_relvel_y = s_V - s_dydt;
        
        % Segment drag
        s_drag_x = 0.5 * prey.Cd .* sim.rho_water .* ...
            repmat(prey.wet_area,1,length(t)) .* ...
            s_relvel_x .* abs(s_relvel_x);
        
        s_drag_y = 0.5 * prey.Cd .* sim.rho_water .* ...
            repmat(prey.wet_area,1,length(t)) .* ...
            sim.rho_water .* s_relvel_y .* abs(s_relvel_y);
        
        % Total drag
        D = [sum(s_drag_x,1); sum(s_drag_y,1)];
        
        % Relative acceleration at segments
        s_relacc_x = s_dUdt - s_dUdt_prev;
        s_relacc_y = s_dVdt - s_dVdt_prev;
        
        % Pressure gradient at segments
        s_dPdx = -sim.rho_water * (s_dUdt + s_U.*s_dUdx);
        s_dPdy = -sim.rho_water * (s_dVdt + s_V.*s_dVdy);
        
        % Total pressure force
        PF(1,:) = sum(-s_dPdx.*repmat(prey.vol,1,length(t)),1);
        PF(2,:) = sum(-s_dPdy.*repmat(prey.vol,1,length(t)),1);
        
        % Acceleration reaction force
        AR(1,:) = sum(-prey.add_mass * sim.rho_water .* ...
            repmat(prey.vol,1,length(t)) .* s_relacc_x);
        
        AR(2,:) = sum(-prey.add_mass * sim.rho_water .* ...
            repmat(prey.vol,1,length(t)) .* s_relacc_y);
        
        % Body acceleration
        accelCOM = (D + PF + AR)./prey.mass;
        
        %t,D,PF,AR
        
        % Output: x speed
        dX(1,:) = velCOM_x;
        % Output: y speed
        dX(2,:) = velCOM_y;
        % Output: theta speed
        dX(3,:) = velCOM_ang;
        
        % Output: x acceleration
        dX(4,:) = accelCOM(1,:);
        % Output: y acceleration
        dX(5,:) = accelCOM(2,:);
        % Output: theta acceleration (torques not currently supported)
        dX(6,:) = 0.*accelCOM(2,:);
        
        %idx,D,PF,AR
        
        % Predator's position (used by status_check)
        pred_position(1,1) = interp1(fl.t,fl.pos(:,1),min(t));
        pred_position(1,2) = interp1(fl.t,fl.pos(:,2),min(t));
        
    end


end



%function [value,isterminal,direction] = evnts(t,X)
function status = status_check(t,y,flag,varargin)

global pred_position

if ~isempty(t)
    pos = y(1);
    if isnan(pos) || (pos <= pred_position(1))
        status = 1;
    else
        status = 0;
    end
else
    status = 1;
end
    
end


function esc = getThrust(t,p)
    esc = p.prey.esc;
end


