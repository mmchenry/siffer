function fl = flow_uniform_ramp(sim,pred,vis_flow)
% Creates a spatially uniform flow field that changes with time.
% All vectors set equal to gape speed.


%% Define time and spatial domain

% Time vector
t = linspace(0,sim.dur,sim.num_time)';

% Spatial vectors
xs = linspace(sim.flow_lim(1),sim.flow_lim(2),sim.num_x);
ys = sim.flow_lim(3):mean(diff(xs)):sim.flow_lim(4);

% Mesh position for flow in global FOR
[Xg,Yg] = meshgrid(xs,ys);

clear xs ys flow_lim


%% Predator gape and position

% Temporal variables for flow
pos      = [0.*t 0.*t];
gape_spd = getGapeSpeed_up(t,pred);


%% Interpolate feeding field in global FOR

% Store variables up to this point
fl.t = t;
fl.pos = pos;
fl.gape_spd = gape_spd;
fl.X = Xg;
fl.Y = Yg;

% Scale up the flow field for each iteration
for i = 1:sim.num_time

    % Store results
    %fl.T(:,:,i) = t(i) + 0.*Xg;
    fl.U(:,:,i) = Xg.*0 - gape_spd(i);
    fl.V(:,:,i) = Xg.*0;
    
end

% Find derivatives
[fl.dUdx,dUdy,fl.dUdt] = gradient(fl.U,fl.X(1,:),fl.Y(:,1),fl.t);
[dVdx,fl.dVdy,fl.dVdt] = gradient(fl.V,fl.X(1,:),fl.Y(:,1),fl.t);

clear t pos gape_spd gape Xg Yg f dUdy dVdx


%% Visualize & validate results


if vis_flow
    
% Make figure window    
f1 = figure;

% Line number to examine flow profile
l_num = round(size(fl.X,1)/2);

% Top speed to scale plots
spd_lim = 1.5*max([max(fl.V(:)) max(fl.U(:))]);

% Step through time
for i = 1:sim.num_time
    
    % Current speed
    spd_vals = sqrt((fl.U(:,:,i)).^2 + (fl.V(:,:,i)).^2);
    
    % Calculate a speed profile through center of field
    Uval = reshape(fl.U(:,l_num,i),size(fl.U,1),1,1);
    
    % Plot raw CFD data 
    subplot(3,1,1:2)
    h = pcolor(fl.X,fl.Y,spd_vals);
    set(h,'EdgeColor','none')
    caxis([0 spd_lim]);
    colorbar
    axis equal
    title(['t = ' num2str(fl.t(i))]);

    % Plot profiles 
    if 0
        % Verify derivatives
        subplot(3,1,3)
%         plot(fl.X(l_num,:),fl.dUdx(l_num,:,i),'k', ...
%              fl.X(l_num,:),[0 diff(fl.U(l_num,:,i))./diff(fl.X(l_num,:))],'r--')
          plot(fl.X(l_num,:),fl.dVdx(l_num,:,i),'k', ...
              fl.X(l_num,:),[0 diff(fl.V(l_num,:,i))./diff(fl.X(l_num,:))],'r--')
          
    % Speed profile
    else
        subplot(3,1,3)
        plot(fl.X(l_num,:),spd_vals(l_num,:),'k')
        ylim([0 spd_lim])
        xlabel('X');
        ylabel('U');
    end
     
    % Pause to display
    pause(.01)
    
    % Clear for next loop
    clear spd_vals

end

close

% Test gradident calculation over time
if 0
    figure;
    c_num = round(size(fl.X,2)/2);
    a = [0; diff(reshape(fl.U(l_num,c_num,:),size(fl.U,3),1,1))./mean(diff(fl.t))];
    a_test = reshape(fl.dUdt(l_num,c_num,:),size(fl.dUdt,3),1,1);
    plot(fl.t,a_test,'k',fl.t,a,'r--')
    ylabel('Acceleration')
end


end


function spd = getGapeSpeed_up(t,pred)
% Speed of flow (inertial FOR) at mouth with constant increase
% Keeps going beyond spd.max, if time permits

spd   = (pred.spd.max/pred.spd.t_max) * t;

end

function spd = getGapeSpeed_updown(t,pred)
% Speed of flow (inertial FOR) at mouth -- ramps up, then down at constant
% rate
idx        = t<=pred.spd.t_max;
spd        = 0.*t;
spd(idx)   = (pred.spd.max/pred.spd.t_max) * t(idx);
spd(~idx)  = -(pred.spd.max/pred.spd.t_max) * t(~idx) + 2*pred.spd.max;
spd(spd<0) = 0;

end


function gape = getGape(t,pred)
% Gape diameter
gape = pred.gape.max.*((t./pred.gape.t_max).*...
                  (exp(1-(t./pred.gape.t_max)))).^pred.gape.alpha;
end


function pos = getPredPos(t,pred)
    dist = pred.dist.init + pred.dist.max.*((t./pred.dist.t_max).*...
                         (exp(1-(t./pred.dist.t_max)))).^pred.dist.alpha;
    pos = [dist dist.*0];
end



end


