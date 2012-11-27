function vis_field(fl,sim)
% Creates animation of flow field 

% Frames per second for playing animation
fps = 120;


%% Visualize & validate results

    
% Make figure window    
figure;

% Line number to examine flow profile
l_num = round(size(fl.X,1)/2);

% Top speed to scale plots
spd_lim = 1.5*max([max(fl.V(:)) max(fl.U(:))]);

% Set frame advance to unity
frame_skip = 1;

% Start frame index
i = 1;

% Step through time
while true
%for i = 1:sim.num_time
    
    % Start timer
    tic
    
    % Current speed
    spd_vals = sqrt((fl.U(:,:,i)).^2 + (fl.V(:,:,i)).^2);
    
    % Plot raw CFD data 
    h = pcolor(fl.X,fl.Y,spd_vals);
    set(h,'EdgeColor','none')
    caxis([0 spd_lim]);
    colorbar
    axis equal
    title(['t = ' num2str(fl.t(i))]);
     
    
    dur = toc;
    
    % Pause to display
    if dur<(1/fps)
        pause(1/fps - dur);
    else
        frame_skip = ceil((1/fps) / dur);
        pause(1e-3);
    end
    
    % Step to next frame
    i = i + frame_skip;
    
    % Check for end
    if i > sim.num_time
        break
    end

    % Clear for next loop
    clear spd_vals
  
end
    
    
    
end