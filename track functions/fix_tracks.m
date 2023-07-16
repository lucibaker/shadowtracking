function [tracks, tracklength] = fix_tracks(tracks, tracklength, searchrad, dt, skipframes)
% stitch together a track that is split up because the particle was not
% detected in one or more frames
%
% inputs:
% tracks array
% tracklengths arrays
% searchrad: search radius
% dt: time interval between frames (1/framerate) [s]
% skipframes: maximum number of skipped frames to allow

fprintf('\nrepairing broken tracks...')

% plot or not?
plot_on = true;

if plot_on
    figure;
    xlabel('x [m]'); ylabel('y [m]'); 
    c = jet(max(tracks(:,5)));
end

% fix broken tracks
for uid = 1:max(tracks(:,5))
    match_idx = 1;
    while numel(match_idx) > 0  % iterate until all matches are found
        uid_idx = tracks(:,5) == uid;
        if sum(uid_idx)
            track1 = tracks(uid_idx,:);        
        
            % endpoint of this track
            [~,idx_end] = max(track1(:,6));
            l_end = track1(idx_end,6);
            f_end = track1(idx_end,7);
            x_end = track1(idx_end,1);
            y_end = track1(idx_end,2);
            u_end = track1(idx_end,3);
            v_end = track1(idx_end,4);
            ax_end = track1(idx_end,8);
            ay_end = track1(idx_end,9);
        
            for df = 2:skipframes
                % predicted startpoint of next track
                f_start = f_end + df;
                x_start = x_end + u_end*dt*df;% + 0.5*ax_end*(dt*df)^2; % acceleration optional
                y_start = y_end + v_end*dt*df;% + 0.5*ay_end*(dt*df)^2;
            
                % find matching startpoints of other tracks
                match_idx = find( ~uid_idx & ... % not original track
                    tracks(:,6) == 0 & ... % check for start of track 
                    tracks(:,7) == f_start & ... % check for frame number
                    abs(tracks(:,1) - x_start) < searchrad & ... % check for x position
                    abs(tracks(:,2) - y_start) < searchrad ); % check for y position
                
                if numel(match_idx) > 0
        
                    % if more than one match, choose the one closest to the predicted startpoint
                    if numel(match_idx) > 1
                        [~,I] = min((tracks(match_idx,1) - x_start).^2 + (tracks(match_idx,2) - y_start).^2);
                        match_idx = match_idx(I);
                    end
                    match_uid = tracks(match_idx,5);
                    match_idx = tracks(:,5) == match_uid;
                    
                    % reassign uid of this track to the uid of the found match
                    tracks(match_idx,5) = uid;
                    
                    % update matched track's lifetime counter
                    tracks(match_idx,6) = ((l_end + df + 1):(l_end + df + sum(match_idx)))';
        
                    % update tracklengths
                    tracklength(uid) = tracklength(match_uid) + tracklength(uid);
                    tracklength(match_uid) = nan;
        
                    % interpolate missing values
                    fixed_idx = tracks(:,5) == uid; % fixed track indices
                    l_interp = ((l_end + 1):(l_end + df - 1))'; % frames to interpolate
                    
                    track_interp = zeros(length(l_interp),size(tracks,2));
                    track_interp(:,1:2) = interp1(tracks(fixed_idx,6),tracks(fixed_idx,1:2),l_interp); % interpolate positions
                    track_interp(:,3) = gradient(track_interp(:,1))/dt; % velocity components
                    track_interp(:,4) = gradient(track_interp(:,2))/dt; 
                    track_interp(:,5) = uid*ones(size(l_interp)); % uid
                    track_interp(:,6:7) = [l_interp, ((f_end+1):(f_start-1))']; % lifetime and frame counters
                    track_interp(:,8) = gradient(track_interp(:,3))/dt; % acceleration components
                    track_interp(:,9) = gradient(track_interp(:,4))/dt; % acceleration components
                    track_interp(:,10:end) = interp1(tracks(fixed_idx,6),tracks(fixed_idx,10:end),l_interp); % errchk and orientations
                 
                    tracks = vertcat(tracks,track_interp);
        
                    % plot
                    if plot_on 
                        fixed_idx = tracks(:,5) == uid;
                        plot(tracks(fixed_idx,1),tracks(fixed_idx,2),'.','color',c(randi(length(c)),:));
                        hold on; plot(tracks(find(match_idx,1),1),tracks(find(match_idx,1),2),'ko','markersize',8,'linewidth',1); %hold off
                        axis equal; axis([-.5 .5 -.45 .05]);
                        pause; 
                    end 
                    
                    break  % when match is found, move to the next iteration with new endpoint
                end  
            end
        else
            match_idx = []; % break the loop if no track is found at this uid
        end
    end
end

matches_found = sum(isnan(tracklength));

tracklength(isnan(tracklength)) = [];

% sort tracks by uid
counter = 1;
uid_new = 1;
tracks_unsorted = tracks;
tracks = zeros(size(tracks));
tracklength_unsorted = tracklength;
tracklength = zeros(size(tracklength));
for uid = 1:max(tracks_unsorted(:,5))
    idx = tracks_unsorted(:,5)==uid;
    if any(idx)
        tracks(counter:counter+sum(idx)-1,:) = sortrows(tracks_unsorted(idx,:),7,'ascend');
        tracks(counter:counter+sum(idx)-1, 5) = uid_new;
        tracklength(uid_new) = sum(idx);
        counter = counter + sum(idx);
        uid_new = uid_new + 1;
    end
end

fprintf('found %i matches\n',matches_found)

end