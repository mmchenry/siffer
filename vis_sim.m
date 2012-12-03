function vis_animate(fl,sim,prey,r)
% Creates animation of flow field or simulation data



% Load morphology data ('m')
load('adult morphology.mat');

% Frames per second for playing animation
fps = 60;

bod_clr = .7.*[1 1 1];


%% Visualize & validate results

    
% Make figure window    
figure;

% Line number to examine flow profile
l_num = round(size(fl.X,1)/2);

% Top speed to scale plots
spd_lim = 1.5*max([max(fl.V(:)) max(fl.U(:))]);

% Set frame advance to unity
frame_skip = 1;

% Start timer  
tstart = tic;

% Step through time
for i = 1:length(r.t)
    
    tmp = abs(r.t(i)-fl.t);
    idx = find(tmp==min(tmp),1,'first');
    
    % Current speed
    spd_vals = sqrt((fl.U(:,:,idx)).^2 + (fl.V(:,:,idx)).^2);
    
    % Plot raw CFD data 
    h = pcolor(fl.X,fl.Y,spd_vals);
    set(h,'EdgeColor','none')
    caxis([0 spd_lim]);
    colorbar
    axis equal
    xL = xlim;
    title(['t = ' num2str(r.t(i))]);
    
    % Get cooridnate of the prey body
    %[py_x,py_y] = prey_body(prey,r,i);
    
    % Add prey body
    hold on
    %h = fill(py_x,py_y,bod_clr);
    h = fill(r.marg_x(:,i),r.marg_y(:,i),bod_clr);
    set(h,'LineStyle','none')
    
    % Coordinates of predator body
    [pd_x,pd_y] = predator_body(m,r,i);
    h = fill(pd_x,pd_y,bod_clr);
    %set(h,'LineStyle','none')
    
    xlim(xL);
    hold off
    

    dur = toc(tstart);
    
    % Pause to display
    if dur<(1/fps)
        pause(1/fps - dur);
    else
        frame_skip = ceil((1/fps) / dur);
        pause(1e-3);
    end

    % Clear for next loop
    clear spd_vals
  
end

% Update title
if r.cap
    title('CAPTURED !')
else
    title('ESCAPED !')
end

pause(.5)


return

  
function [py_x,py_y] = prey_body(prey,r,i)
% Gives coordinate for prey body for a particular instant of a simulation
% in the inertial FOR

% Body coorindates in prey FOR
bod_x = [prey.s'  prey.s(end:-1:1)']-prey.sCOM;
bod_y = [prey.r' -prey.r(end:-1:1)'];

% Rotation matrix
M = [cos(r.pos(3,i)) -sin(r.pos(3,i)); ...
    sin(r.pos(3,i))  cos(r.pos(3,i))];

% Rotate
py_bod = M*[bod_x; bod_y];

% Translate position
py_x = py_bod(1,:) + r.pos(1,i);
py_y = py_bod(2,:) + r.pos(2,i);


function [pd_x,pd_y] = predator_body(m,r,i)

% Shape in pred FOR
bod_x = [m.s' m.s(end:-1:1)']-max(m.s);
bod_y = [m.w'./2 -m.w(end:-1:1)'./2];

% Translate in interial FOR
pd_x = bod_x + r.pred_pos(1,i);
pd_y = bod_y + r.pred_pos(2,i);

