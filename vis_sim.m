function vis_animate(fl,sim,prey,r)
% Creates animation of flow field or simulation data



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

% Define segment positions for all time values
s_x = repmat((prey.s-prey.sCOM) ,1,length(r.t)) .* ...
    repmat(cos(r.pos(3,:))      ,length(prey.s),1) + ...
    repmat(r.pos(1,:)           ,length(prey.s),1);

s_y = repmat((prey.s-prey.sCOM) ,1,length(r.t)) .* ...
    repmat(sin(r.pos(3,:))      ,length(prey.s),1) + ...
    repmat(r.pos(2,:)           ,length(prey.s),1);


% Step through time
for i = 1:length(r.t)
    
    
    

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




return




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
    pause(.5)
    
    % Clear for next loop
    clear spd_vals
    
    
    
    
end